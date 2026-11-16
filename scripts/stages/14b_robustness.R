# Stage 14b：丰度指数（壁细胞、周细胞、平滑肌、内皮、Müller；节细胞/视杆作部位阳性对照）的稳健性
#   1) 主模型：指数 ~ 病程 + 部位 + 年龄 + 性别 + RIN + 死亡后间隔，供体内相关（duplicateCorrelation）
#   2) 供体层面置换：病程标签在供体间打乱（同一供体的样本一起移动），2000 次
#   3) 逐一剔除供体后的效应范围；4) 黄斑、周边分别拟合；
#   5) 糖尿病供体内的剂量关系：指数 ~ ETDRS 严重度等级（10–65）；病程（可解析为年数者）
#   6) Müller 状态：RLBP1、VIM、GLUL 相对 Müller 标记的比值
suppressPackageStartupMessages({library(limma); library(edgeR); library(data.table); library(parallel)})
source(file.path("scripts", "lib", "repro.R"))
cfg <- load_config(); set_global_seed(cfg$seed)
root <- p22_root(); setwd(root)
stage <- "14b_robustness"
out <- file.path("results", stage); dir.create(out, recursive = TRUE, showWarnings = FALSE)
n_perm <- 2000

meta_all <- fread(file.path("results", "06_retina_pilot", "sample_meta.tsv"))
cnt <- as.matrix(read.csv(gzfile(cfg$retina$counts), row.names = 1, check.names = FALSE))[, meta_all$sample]
gm <- fread(cfg$deconvolution$gene_map); setnames(gm, 1:2, c("id", "gene")); gm[, id := sub("\\..*$", "", id)]; gm <- gm[!duplicated(id)]
sym <- gm$gene[match(rownames(cnt), gm$id)]; cnt <- cnt[!is.na(sym), ]; cnt <- rowsum(cnt, sym[!is.na(sym)])
lc <- cpm(calcNormFactors(DGEList(cnt)), log = TRUE, prior.count = 2)

ms <- fread(file.path("results", "14a_marker_sets", "marker_sets.tsv"))
idx <- t(sapply(split(ms$gene, ms$set), function(g) colMeans(lc[intersect(g, rownames(lc)), , drop = FALSE])))
mul <- split(ms$gene, ms$set)[["muller_top8"]]
for (g in c("RLBP1", "VIM", "GLUL", "GFAP")) idx <- rbind(idx, setNames(list(lc[g, ] - colMeans(lc[mul, ])), paste0("ratio_", g, "_muller"))[[1]])
rownames(idx)[(nrow(idx) - 3):nrow(idx)] <- paste0("ratio_", c("RLBP1", "VIM", "GLUL", "GFAP"), "_muller")
fwrite(data.table(index = rownames(idx), idx), file.path(out, "indices_by_sample.tsv"), sep = "\t")

dur_years <- function(x) { v <- suppressWarnings(as.numeric(sub("^[^0-9]*([0-9]+).*$", "\\1", x))); v[!grepl("[0-9]", x)] <- NA; v }
meta_all[, `:=`(region = factor(region, levels = c("macula", "periphery")), dur = dur_years(duration))]
covs <- "region + age + sex + rin + pmi"

fit3 <- function(Y, d, rhs, cor = NULL) {
  X <- model.matrix(as.formula(paste("~", rhs)), data = d)
  if (is.null(cor)) { cor <- duplicateCorrelation(Y, X, block = d$donor)$consensus.correlation; if (!is.finite(cor)) cor <- 0 }
  list(f = eBayes(lmFit(Y, X, block = d$donor, correlation = cor)), cor = cor)
}
grab <- function(f, cf) { tt <- topTable(f, coef = cf, number = Inf, sort.by = "none", confint = TRUE)
  data.table(index = rownames(tt), coef = cf, diff = tt$logFC, ci_low = tt$CI.L, ci_high = tt$CI.R, t = tt$t, p = tt$P.Value) }

# ---- 1) 主模型 + 2) 置换 ----
m3 <- meta_all[group %in% c("Control", "Diabetic", "NPDR")][, group := factor(group, levels = c("Control", "Diabetic", "NPDR"))]
Y3 <- idx[, m3$sample]
main <- fit3(Y3, m3, paste("group +", covs))
res <- rbind(grab(main$f, "groupDiabetic"), grab(main$f, "groupNPDR"), grab(main$f, "regionperiphery"))
donors <- unique(m3[, .(donor, group)])
perm <- mclapply(seq_len(n_perm), function(b) {
  set.seed(cfg$seed * 100000 + b); sh <- copy(donors)[, group := sample(group)]
  dp <- copy(m3)[, group := factor(sh$group[match(donor, sh$donor)], levels = levels(m3$group))]
  f <- fit3(Y3, dp, paste("group +", covs), cor = main$cor)$f
  cbind(f$t[, "groupDiabetic"], f$t[, "groupNPDR"])
}, mc.cores = 20L)
pt <- simplify2array(perm)
res[coef == "groupDiabetic", p_perm := (1 + rowSums(abs(pt[, 1, ]) >= abs(t))) / (1 + n_perm)]
res[coef == "groupNPDR", p_perm := (1 + rowSums(abs(pt[, 2, ]) >= abs(t))) / (1 + n_perm)]
res[, donor_cor := main$cor]
fwrite(res, file.path(out, "index_stage_tests.tsv"), sep = "\t")

# ---- 3) 逐一剔除供体 ----
loo <- rbindlist(lapply(unique(m3$donor), function(dn) {
  d <- m3[donor != dn]; f <- fit3(Y3[, d$sample], d, paste("group +", covs), cor = main$cor)$f
  rbind(grab(f, "groupDiabetic"), grab(f, "groupNPDR"))[, dropped := dn]
}))
fwrite(loo[, .(min_diff = min(diff), max_diff = max(diff), max_p = max(p)), by = .(index, coef)], file.path(out, "leave_one_donor_out.tsv"), sep = "\t")

# ---- 4) 黄斑、周边分别 ----
reg <- rbindlist(lapply(c("macula", "periphery"), function(rg) {
  d <- m3[region == rg]; f <- fit3(Y3[, d$sample], d, "group + age + sex + rin + pmi")$f
  rbind(grab(f, "groupDiabetic"), grab(f, "groupNPDR"))[, region := rg]
}))
fwrite(reg, file.path(out, "index_by_region.tsv"), sep = "\t")

# ---- 5) 糖尿病供体内剂量关系（含 PDR 两位供体） ----
md <- meta_all[group != "Control" & !is.na(severity)]
fs <- fit3(idx[, md$sample], md, paste("severity +", covs))$f
dose <- grab(fs, "severity")[, `:=`(predictor = "ETDRS level", n_donors = uniqueN(md$donor))]
mdd <- md[!is.na(dur)]
if (uniqueN(mdd$donor) >= 10) {
  fd <- fit3(idx[, mdd$sample], mdd, paste("dur +", covs))$f
  dose <- rbind(dose, grab(fd, "dur")[, `:=`(predictor = "diabetes duration (years)", n_donors = uniqueN(mdd$donor))])
}
mnp <- md[group != "PDR"]                                   # 敏感性：去掉仅 2 位供体的 PDR
fn <- fit3(idx[, mnp$sample], mnp, paste("severity +", covs))$f
dose <- rbind(dose, grab(fn, "severity")[, `:=`(predictor = "ETDRS level, PDR excluded", n_donors = uniqueN(mnp$donor))])
fwrite(dose, file.path(out, "dose_response.tsv"), sep = "\t")

key <- c("mural_top8", "mural_top5", "mural_top12", "mural_top20", "pericyte_only_top8", "smc_only_top8",
         "endothelial_top8", "muller_top8", "ganglion_top8", "rod_top8", "ratio_RLBP1_muller", "ratio_VIM_muller", "ratio_GFAP_muller")
print(res[index %in% key, .(index, coef, diff = round(diff, 3), t = round(t, 2), p = signif(p, 2), p_perm = signif(p_perm, 2))][order(coef, index)])
print(fread(file.path(out, "leave_one_donor_out.tsv"))[index %in% key])
print(reg[index %in% key, .(index, coef, region, t = round(t, 2), p = signif(p, 2))][order(index, coef, region)])
print(dose[index %in% key, .(index, predictor, n_donors, diff = round(diff, 4), t = round(t, 2), p = signif(p, 2))])
write_provenance(stage, inputs = c(cfg$retina$counts, file.path("results", "14a_marker_sets", "marker_sets.tsv")),
                 outputs = list.files(out, "\\.tsv$", full.names = TRUE), seed = cfg$seed, extra = list(n_perm = n_perm))
