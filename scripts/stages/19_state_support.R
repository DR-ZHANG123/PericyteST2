# Stage 19：支撑主线的补充检验
#   1) 独立人源模型：糖尿病视网膜微血管芯片（GSE228900；内皮 + 周细胞 + 星形胶质共培养），糖尿病刺激 vs 甘露醇
#   2) 8 标记壁细胞指数逐个留出标记：DR 前的减少由哪些标记承担
#   3) IL-33（ST2 配体）相对 Müller 标记的组间差异
#   4) IL1RL1 相对壁细胞含量：加入糖尿病黄斑水肿（DME）为协变量
suppressPackageStartupMessages({library(limma); library(edgeR); library(data.table); library(parallel)})
source(file.path("scripts", "lib", "repro.R"))
cfg <- load_config(); set_global_seed(cfg$seed)
root <- p22_root(); setwd(root)
stage <- "19_state_support"; sc <- cfg$state_support; n_perm <- sc$n_perm
out <- file.path("results", stage); dir.create(out, recursive = TRUE, showWarnings = FALSE)
gm <- fread(cfg$deconvolution$gene_map); setnames(gm, 1:2, c("id", "gene")); gm[, id := sub("\\..*$", "", id)]; gm <- gm[!duplicated(id)]

# ---- 1) 芯片 ----
fs <- list.files(sc$chip_dir, "Exp2-Cond[23].*counts.txt.gz", full.names = TRUE)
M <- sapply(fs, function(f) { x <- read.delim(gzfile(f), header = FALSE); setNames(x$V2, x$V1) })
colnames(M) <- sub(".*_(Exp.*)\\.counts.txt.gz", "\\1", colnames(M)); M <- M[!grepl("^__", rownames(M)), ]
sym <- gm$gene[match(sub("\\..*$", "", rownames(M)), gm$id)]; M <- rowsum(M[!is.na(sym), ], sym[!is.na(sym)])
lcc <- cpm(calcNormFactors(DGEList(M)), log = TRUE, prior.count = 2)
info <- data.table(sample = colnames(lcc), condition = ifelse(grepl("Cond2", colnames(lcc)), "mannitol", "diabetic"),
                   day = ifelse(grepl("t2", colnames(lcc)), "d14", "d28"))
info[, condition := factor(condition, levels = c("mannitol", "diabetic"))]
fc <- eBayes(lmFit(lcc, model.matrix(~ day + condition, info)))
tc <- topTable(fc, coef = "conditiondiabetic", number = Inf, sort.by = "none", confint = TRUE)
chip <- data.table(gene = rownames(tc), lfc = tc$logFC, ci_low = tc$CI.L, ci_high = tc$CI.R, t = tc$t, p = tc$P.Value, ave_logcpm = tc$AveExpr)
chip[ave_logcpm >= sc$chip_min_logcpm, fdr := p.adjust(p, "BH")]
fwrite(chip[order(p)], file.path(out, "chip_all_genes.tsv"), sep = "\t")
fwrite(chip[gene %in% sc$chip_genes][order(match(gene, sc$chip_genes))], file.path(out, "chip_selected.tsv"), sep = "\t")
fwrite(info, file.path(out, "chip_samples.tsv"), sep = "\t")

# ---- 供体视网膜 ----
meta_all <- fread(file.path("results", "06_retina_pilot", "sample_meta.tsv"))
cnt <- as.matrix(read.csv(gzfile(cfg$retina$counts), row.names = 1, check.names = FALSE))[, meta_all$sample]
s2 <- gm$gene[match(rownames(cnt), gm$id)]; cnt <- cnt[!is.na(s2), ]; cnt <- rowsum(cnt, s2[!is.na(s2)])
d <- DGEList(cnt); d <- calcNormFactors(d[filterByExpr(d, group = meta_all$group), , keep.lib.sizes = FALSE])
lc <- cpm(d, log = TRUE, prior.count = 2)
d_all <- DGEList(cnt); d_all$samples$norm.factors <- d$samples$norm.factors      # 标记指数：不按表达过滤，归一化因子同上
lc_all <- cpm(d_all, log = TRUE, prior.count = 2)
ms <- fread(file.path("results", "14a_marker_sets", "marker_sets.tsv"))
share <- fread(file.path("results", "09c_ratio_validation", "bulk_share_by_celltype.tsv")); setnames(share, 1, "gene")
mur <- share[mural >= cfg$mural_state$share_min & gene %in% rownames(lc)]$gene
content <- apply(lc[mur, ], 2, median)
covs <- "region + age + sex + rin + pmi"
m3 <- meta_all[group %in% c("Control", "Diabetic", "NPDR")][, group := factor(group, levels = c("Control", "Diabetic", "NPDR"))]
fitm <- function(Y, dd, rhs, cor = NULL) {
  X <- model.matrix(as.formula(paste("~", rhs)), data = dd)
  if (is.null(cor)) { cor <- duplicateCorrelation(Y, X, block = dd$donor)$consensus.correlation; if (!is.finite(cor)) cor <- 0 }
  list(f = eBayes(lmFit(Y, X, block = dd$donor, correlation = cor)), cor = cor)
}
grab <- function(f, cf) { tt <- topTable(f, coef = cf, number = Inf, sort.by = "none", confint = TRUE)
  data.table(measure = rownames(tt), coef = cf, diff = tt$logFC, ci_low = tt$CI.L, ci_high = tt$CI.R, t = tt$t, p = tt$P.Value) }
donors <- unique(m3[, .(donor, group)])
perm_test <- function(Y, dd, rhs) {               # 供体层面置换：组别在供体间交换，样本跟随供体
  main <- fitm(Y, dd, rhs)
  pt <- simplify2array(mclapply(seq_len(n_perm), function(b) {
    set.seed(cfg$seed * 100000 + b); sh <- copy(donors)[, group := sample(group)]
    dp <- copy(dd)[, group := factor(sh$group[match(donor, sh$donor)], levels = levels(dd$group))]
    fitm(Y, dp, rhs, main$cor)$f$t[, c("groupDiabetic", "groupNPDR"), drop = FALSE]
  }, mc.cores = 16L))
  res <- rbind(grab(main$f, "groupDiabetic"), grab(main$f, "groupNPDR"))
  res[, p_perm := mapply(function(m, cf, tv) (1 + sum(abs(pt[m, cf, ]) >= abs(tv))) / (1 + n_perm), measure, coef, t)]
}

# ---- 2) 8 标记壁细胞指数逐个留出 ----
m8 <- intersect(ms[set == "mural_top8"]$gene, rownames(lc_all))
L <- rbind(colMeans(lc_all[m8, ]), t(sapply(m8, function(g) colMeans(lc_all[setdiff(m8, g), ]))),
           colMeans(lc_all[setdiff(m8, c("IL1RL1", "ADAMTS4")), ]), content)
rownames(L) <- c("all eight", paste("without", m8), "without IL1RL1 and ADAMTS4", "median of mural-enriched genes")
loo <- perm_test(L[, m3$sample], m3, paste("group +", covs))
per_marker <- perm_test(lc_all[m8, m3$sample], m3, paste("group +", covs))
fwrite(loo, file.path(out, "mural_index_leave_marker_out.tsv"), sep = "\t")
bys <- data.table(sample = colnames(lc), mural_content = content, mural_top8 = L["all eight", ],
                  mural_top8_without_il1rl1_adamts4 = L["without IL1RL1 and ADAMTS4", ], il1rl1_rel_content = lc["IL1RL1", ] - content)
fwrite(merge(meta_all[, .(sample, donor, group, region, severity)], bys, by = "sample"), file.path(out, "mural_measures_by_sample.tsv"), sep = "\t")
fwrite(per_marker, file.path(out, "mural_markers_single.tsv"), sep = "\t")

# 供体层面 IL1RL1 与胶原程序：去掉组别与协变量后的残差（与 stage 18a 的关联检验同一模型），供作图
dd <- fread(file.path("results", "18a_mural_state", "il1rl1_by_donor.tsv"))
dd[, il_resid := resid(lm(il ~ group + age + sex + rin + pmi, dd))][, col_resid := resid(lm(col ~ group + age + sex + rin + pmi, dd))]
fwrite(dd, file.path(out, "il1rl1_collagen_residuals.tsv"), sep = "\t")

# ---- 3) IL-33 相对 Müller 标记 ----
mul <- colMeans(lc[intersect(ms[set == "muller_top8"]$gene, rownames(lc)), ])
Y3 <- rbind(IL33_relative_to_muller = lc["IL33", ] - mul, IL33_bulk = lc["IL33", ])
il33 <- perm_test(Y3[, m3$sample], m3, paste("group +", covs))
il33 <- merge(il33, data.table(measure = rownames(Y3), share_muller = share[gene == "IL33"]$muller,
                               share_endothelial = share[gene == "IL33"]$endothelial), by = "measure")
fwrite(il33, file.path(out, "il33_tests.tsv"), sep = "\t")

# ---- 4) IL1RL1 相对含量 + DME ----
il <- lc["IL1RL1", ] - content
m3[, dme := as.logical(dme)]
Yd <- rbind(il[m3$sample], il[m3$sample]); rownames(Yd) <- c("IL1RL1_rel", "IL1RL1_rel_dup")
fd <- fitm(Yd, m3, paste("group + dme +", covs))$f
dme_fit <- rbind(grab(fd, "groupDiabetic"), grab(fd, "groupNPDR"))[measure == "IL1RL1_rel"]
dme_fit[, model := "IL1RL1 relative to mural content, with diabetic macular oedema in the model"]
fwrite(dme_fit, file.path(out, "il1rl1_dme.tsv"), sep = "\t")

print(chip[gene %in% sc$chip_genes, .(gene, lfc = round(lfc, 2), p = signif(p, 2), fdr = signif(fdr, 2))])
print(loo[, .(measure, coef, diff = round(diff, 2), p_perm = signif(p_perm, 2))])
print(il33[, .(measure, coef, diff = round(diff, 2), p_perm = signif(p_perm, 2))]); print(dme_fit)
write_provenance(stage, inputs = c(cfg$retina$counts, fs), outputs = list.files(out, "\\.tsv$", full.names = TRUE),
                 seed = cfg$seed, extra = list(n_perm = n_perm, chip_samples = nrow(info)))
