# Stage 21：壁细胞状态基因在供体间如何共变
#   每个样本：基因 log-CPM − 壁细胞含量（52 个壁细胞富集基因中位数）；每位供体取黄斑、周边平均
#   1) 供体 × 基因矩阵（减去对照均值），供热图
#   2) 组内相关：每个基因对 组别 + 年龄 + 性别 + RIN + PMI 回归取残差，基因两两 Spearman 相关（43 位供体），BH 校正
suppressPackageStartupMessages({library(limma); library(edgeR); library(data.table)})
source(file.path("scripts", "lib", "repro.R"))
cfg <- load_config(); set_global_seed(cfg$seed)
root <- p22_root(); setwd(root)
stage <- "21_state_structure"; genes <- cfg$state_structure$genes
out <- file.path("results", stage); dir.create(out, recursive = TRUE, showWarnings = FALSE)

meta_all <- fread(file.path("results", "06_retina_pilot", "sample_meta.tsv"))
cnt <- as.matrix(read.csv(gzfile(cfg$retina$counts), row.names = 1, check.names = FALSE))[, meta_all$sample]
gm <- fread(cfg$deconvolution$gene_map); setnames(gm, 1:2, c("id", "gene")); gm[, id := sub("\\..*$", "", id)]; gm <- gm[!duplicated(id)]
s2 <- gm$gene[match(rownames(cnt), gm$id)]; cnt <- cnt[!is.na(s2), ]; cnt <- rowsum(cnt, s2[!is.na(s2)])
d <- DGEList(cnt); d <- calcNormFactors(d[filterByExpr(d, group = meta_all$group), , keep.lib.sizes = FALSE])
lc <- cpm(d, log = TRUE, prior.count = 2)
d_all <- DGEList(cnt); d_all$samples$norm.factors <- d$samples$norm.factors
la <- cpm(d_all, log = TRUE, prior.count = 2)
share <- fread(file.path("results", "09c_ratio_validation", "bulk_share_by_celltype.tsv")); setnames(share, 1, "gene")
mur <- share[mural >= cfg$mural_state$share_min & gene %in% rownames(lc)]$gene
content <- apply(lc[mur, ], 2, median)

rel <- t(la[genes, meta_all$sample] - matrix(content[meta_all$sample], nrow = length(genes), ncol = nrow(meta_all), byrow = TRUE))
dt <- cbind(meta_all[, .(sample, donor, group, region, severity, age, sex, rin, pmi)], as.data.table(rel))
dn <- dt[, c(list(group = group[1], severity = severity[1], age = age[1], sex = sex[1], rin = mean(rin), pmi = mean(pmi)),
             lapply(.SD, mean)), by = donor, .SDcols = genes]
ctrl <- dn[group == "Control", lapply(.SD, mean), .SDcols = genes]
cen <- copy(dn); for (g in genes) set(cen, j = g, value = cen[[g]] - ctrl[[g]])
fwrite(cen, file.path(out, "state_genes_by_donor.tsv"), sep = "\t")

dn[, group := factor(group, levels = c("Control", "Diabetic", "NPDR", "PDR"))]
res <- sapply(genes, function(g) residuals(lm(as.formula(paste(g, "~ group + age + sex + rin + pmi")), data = dn)))
pairs <- CJ(gene_a = genes, gene_b = genes)[match(gene_a, genes) < match(gene_b, genes)]
pairs[, c("rho", "p") := {
  ct <- mapply(function(a, b) unlist(cor.test(res[, a], res[, b], method = "spearman", exact = FALSE)[c("estimate", "p.value")]), gene_a, gene_b)
  list(ct[1, ], ct[2, ])
}]
pairs[, fdr := p.adjust(p, "BH")][, n_donors := nrow(dn)][, adjusted_for := "donor group, age, sex, RNA integrity, post-mortem interval"]
fwrite(pairs, file.path(out, "state_gene_correlations.tsv"), sep = "\t")
print(pairs[order(p)])

write_provenance(stage, inputs = c(cfg$retina$counts, file.path("results", "06_retina_pilot", "sample_meta.tsv"),
                                    file.path("results", "09c_ratio_validation", "bulk_share_by_celltype.tsv")),
                 outputs = list.files(out, pattern = "tsv$", full.names = TRUE))
