# Stage 18a：糖尿病供体中壁细胞的“状态”改变（区别于数量减少）
#   壁细胞含量 = 壁细胞富集基因（bulk 信号 ≥60% 来自壁细胞）的中位数——单个基因的状态改变对中位数影响小
#   1) 含量本身的组间差异；2) 每个壁细胞基因相对含量的组间差异（供体置换 + BH）；
#   3) IL1RL1 相对含量：校正内皮指数（IL1RL1 有 19% 内皮来源）、糖尿病供体内与 ETDRS 等级、与胶原程序的供体内关联（去掉组别后）；
#   4) 濒死期/应激即刻早期基因的组间差异（排除死亡过程假象）；5) 输出 Geneformer 用的状态签名
suppressPackageStartupMessages({library(limma); library(edgeR); library(data.table); library(parallel)})
source(file.path("scripts", "lib", "repro.R"))
cfg <- load_config(); set_global_seed(cfg$seed)
root <- p22_root(); setwd(root)
stage <- "18a_mural_state"; pc <- cfg$mural_state; n_perm <- pc$n_perm
out <- file.path("results", stage); dir.create(out, recursive = TRUE, showWarnings = FALSE)

meta_all <- fread(file.path("results", "06_retina_pilot", "sample_meta.tsv"))
cnt <- as.matrix(read.csv(gzfile(cfg$retina$counts), row.names = 1, check.names = FALSE))[, meta_all$sample]
gm <- fread(cfg$deconvolution$gene_map); setnames(gm, 1:2, c("id", "gene")); gm[, id := sub("\\..*$", "", id)]; gm <- gm[!duplicated(id)]
sym <- gm$gene[match(rownames(cnt), gm$id)]; cnt <- cnt[!is.na(sym), ]; cnt <- rowsum(cnt, sym[!is.na(sym)])
d <- DGEList(cnt); d <- calcNormFactors(d[filterByExpr(d, group = meta_all$group), , keep.lib.sizes = FALSE])
lc <- cpm(d, log = TRUE, prior.count = 2)

share <- fread(file.path("results", "09c_ratio_validation", "bulk_share_by_celltype.tsv")); setnames(share, 1, "gene")
mur <- share[mural >= pc$share_min & gene %in% rownames(lc)]$gene
ms <- fread(file.path("results", "14a_marker_sets", "marker_sets.tsv"))
content <- apply(lc[mur, ], 2, median)
endo <- colMeans(lc[intersect(ms[set == "endothelial_top8"]$gene, rownames(lc)), ])

covs <- "region + age + sex + rin + pmi"
m3 <- meta_all[group %in% c("Control", "Diabetic", "NPDR")][, group := factor(group, levels = c("Control", "Diabetic", "NPDR"))]
m3[, endo_idx := endo[sample]]
fitm <- function(Y, dd, rhs, cor = NULL) {
  X <- model.matrix(as.formula(paste("~", rhs)), data = dd)
  if (is.null(cor)) { cor <- duplicateCorrelation(Y, X, block = dd$donor)$consensus.correlation; if (!is.finite(cor)) cor <- 0 }
  list(f = eBayes(lmFit(Y, X, block = dd$donor, correlation = cor)), cor = cor)
}
grab <- function(f, cf) { tt <- topTable(f, coef = cf, number = Inf, sort.by = "none", confint = TRUE)
  data.table(measure = rownames(tt), coef = cf, diff = tt$logFC, ci_low = tt$CI.L, ci_high = tt$CI.R, t = tt$t, p = tt$P.Value) }
donors <- unique(m3[, .(donor, group)])
perm_test <- function(Y, rhs, label) {           # 供体层面置换：在供体间交换组别，样本跟随其供体
  main <- fitm(Y, m3, rhs)
  pt <- simplify2array(mclapply(seq_len(n_perm), function(b) {
    set.seed(cfg$seed * 100000 + b); sh <- copy(donors)[, group := sample(group)]
    dp <- copy(m3)[, group := factor(sh$group[match(donor, sh$donor)], levels = levels(m3$group))]
    fitm(Y, dp, rhs, main$cor)$f$t[, c("groupDiabetic", "groupNPDR"), drop = FALSE]
  }, mc.cores = 16L))
  res <- rbind(grab(main$f, "groupDiabetic"), grab(main$f, "groupNPDR"))
  res[, p_perm := mapply(function(m, cf, tv) (1 + sum(abs(pt[m, cf, ]) >= abs(tv))) / (1 + n_perm), measure, coef, t)]
  res[, model := label]
}
dup <- function(v, nm) { Y <- rbind(v, v); rownames(Y) <- c(nm, paste0(nm, "_dup")); Y }

# ---- 1) 壁细胞含量 ----
r_content <- perm_test(dup(content[m3$sample], "mural_content"), paste("group +", covs), "group + covariates")[measure == "mural_content"]

# ---- 2) 每个壁细胞基因相对含量 ----
R <- sweep(lc[mur, m3$sample], 2, content[m3$sample])
r_gene <- perm_test(R, paste("group +", covs), "relative to mural content")
r_gene[, fdr := p.adjust(p, "BH"), by = coef][, fdr_perm := p.adjust(p_perm, "BH"), by = coef]
r_bulk <- rbind(grab(fitm(lc[mur, m3$sample], m3, paste("group +", covs))$f, "groupDiabetic"),
                grab(fitm(lc[mur, m3$sample], m3, paste("group +", covs))$f, "groupNPDR"))[, .(measure, coef, bulk_diff = diff, bulk_p = p)]
r_gene <- merge(r_gene, r_bulk, by = c("measure", "coef"))
r_gene <- merge(r_gene, share[, .(measure = gene, share_mural = mural, share_pericyte = pericyte, share_smc = smooth_muscle, share_endothelial = endothelial)], by = "measure")
fwrite(r_gene[order(coef, p)], file.path(out, "mural_gene_state.tsv"), sep = "\t")

# ---- 3) IL1RL1 ----
il <- lc["IL1RL1", ] - content
r_il <- rbind(perm_test(dup(il[m3$sample], "IL1RL1_rel"), paste("group +", covs), "relative to mural content"),
              perm_test(dup(il[m3$sample], "IL1RL1_rel"), paste("group + endo_idx +", covs), "relative to mural content, adjusted for endothelial index"),
              perm_test(dup((lc["IL1RL1", ] - endo)[m3$sample], "IL1RL1_rel"), paste("group +", covs), "relative to endothelial markers"))[measure == "IL1RL1_rel"]
md <- meta_all[group != "Control" & !is.na(severity)]
sev <- rbind(grab(fitm(dup(il[md$sample], "IL1RL1_rel"), md, paste("severity +", covs))$f, "severity")[measure == "IL1RL1_rel"][, predictor := "ETDRS level"],
             { mn <- md[group != "PDR"]; grab(fitm(dup(il[mn$sample], "IL1RL1_rel"), mn, paste("severity +", covs))$f, "severity")[measure == "IL1RL1_rel"][, predictor := "ETDRS level, PDR excluded"] })
# 供体层面：每位供体取样本平均；与胶原程序的关联在去掉组别与协变量后计算（避免“两者都随组别变”造成的关联）
colp <- colMeans(lc[intersect(pc$program_up, rownames(lc)), ]) - content
dd <- meta_all[, .(il = mean(il[sample]), col = mean(colp[sample]), group = group[1], age = age[1], sex = sex[1],
                   rin = mean(rin), pmi = mean(pmi)), by = donor]
res_il <- resid(lm(il ~ group + age + sex + rin + pmi, dd)); res_col <- resid(lm(col ~ group + age + sex + rin + pmi, dd))
ct <- cor.test(res_il, res_col, method = "spearman", exact = FALSE)
assoc <- data.table(x = "IL1RL1 relative to mural content", y = "collagen genes relative to mural content",
                    adjusted_for = "group, age, sex, RIN, PMI (donor means)", rho = ct$estimate, p = ct$p.value, n_donors = nrow(dd))
by_group <- dd[, .(n = .N, il_mean = mean(il), il_sd = sd(il)), by = group]
fwrite(rbind(r_il, r_content, fill = TRUE), file.path(out, "il1rl1_tests.tsv"), sep = "\t")
fwrite(sev, file.path(out, "il1rl1_severity.tsv"), sep = "\t")
fwrite(assoc, file.path(out, "il1rl1_collagen_association.tsv"), sep = "\t")
fwrite(dd, file.path(out, "il1rl1_by_donor.tsv"), sep = "\t")

# ---- 4) 濒死期 / 即刻早期基因 ----
ieg <- intersect(c("FOS", "FOSB", "JUN", "JUNB", "EGR1", "ATF3", "IER2", "NR4A1", "HSPA1A", "DUSP1", "SOCS3", "HIF1A", "VEGFA", "ADM"), rownames(lc))
r_ieg <- rbind(grab(fitm(lc[ieg, m3$sample], m3, paste("group +", covs))$f, "groupDiabetic"),
               grab(fitm(lc[ieg, m3$sample], m3, paste("group +", covs))$f, "groupNPDR"))
fwrite(r_ieg, file.path(out, "agonal_genes.tsv"), sep = "\t")

# ---- 5) Geneformer 状态签名：相对含量的组间差异（所有壁细胞基因，不按显著性筛选）----
sig <- dcast(r_gene[, .(gene = measure, coef, diff)], gene ~ coef, value.var = "diff")
setnames(sig, c("groupDiabetic", "groupNPDR"), c("lfc_diabetic", "lfc_npdr"))
fwrite(sig, file.path(out, "state_signature.tsv"), sep = "\t")

print(r_content[, .(coef, diff = round(diff, 3), p = signif(p, 2), p_perm = signif(p_perm, 2))])
print(r_gene[fdr_perm < 0.1, .(measure, coef, diff = round(diff, 2), p_perm = signif(p_perm, 2), fdr_perm = signif(fdr_perm, 2))][order(coef, p_perm)])
print(r_il[, .(model, coef, diff = round(diff, 2), p_perm = signif(p_perm, 2))]); print(sev[, .(predictor, diff = round(diff, 3), p = signif(p, 2))])
print(assoc); print(by_group); print(r_ieg[p < 0.05])
write_provenance(stage, inputs = c(cfg$retina$counts, file.path("results", "09c_ratio_validation", "bulk_share_by_celltype.tsv")),
                 outputs = list.files(out, "\\.tsv$", full.names = TRUE), seed = cfg$seed, extra = list(n_perm = n_perm, n_mural_genes = length(mur)))
