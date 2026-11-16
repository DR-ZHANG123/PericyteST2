# Stage 11b：人供体视网膜（GSE160306）中内皮、壁细胞、Müller 的“自身改变”与“数量改变”
# 用 stage 09c 确定的基因集：每个细胞类群的稳定特异标记（分母）与主要来自该类群的待检基因。
#   丰度指数 = 标记基因 log2CPM 的均值（反映该类群转录本在组织中的多少）
#   比值     = 待检基因 log2CPM − 丰度指数（抵消细胞数量变化，反映细胞自身表达；内皮、Müller 已在模拟中验证，壁细胞检出力不足）
# 模型：~ 病程 + 部位 + 年龄 + 性别 + RIN + 死亡后间隔，供体内相关用 duplicateCorrelation；PDR（2 人）不参与。
suppressPackageStartupMessages({library(limma); library(edgeR); library(data.table)})
source(file.path("scripts", "lib", "repro.R"))
cfg <- load_config(); set_global_seed(cfg$seed)
root <- p22_root(); setwd(root)
stage <- "11b_ratio_retina"
out <- file.path("results", stage); dir.create(out, recursive = TRUE, showWarnings = FALSE)

meta <- fread(file.path("results", "06_retina_pilot", "sample_meta.tsv"))[group %in% c("Control", "Diabetic", "NPDR")]
meta[, group := factor(group, levels = c("Control", "Diabetic", "NPDR"))]
meta[, region := factor(region, levels = c("macula", "periphery"))]
cnt <- as.matrix(read.csv(gzfile(cfg$retina$counts), row.names = 1, check.names = FALSE))[, meta$sample]
gm <- fread(cfg$deconvolution$gene_map); setnames(gm, 1:2, c("id", "gene")); gm[, id := sub("\\..*$", "", id)]; gm <- gm[!duplicated(id)]
sym <- gm$gene[match(rownames(cnt), gm$id)]; cnt <- cnt[!is.na(sym), ]; cnt <- rowsum(cnt, sym[!is.na(sym)])
dge <- calcNormFactors(DGEList(cnt)); lc <- cpm(dge, log = TRUE, prior.count = 2)

sets <- fread(file.path("results", "09c_ratio_validation", "gene_sets.tsv"))
X <- model.matrix(~ group + region + age + sex + rin + pmi, data = meta)
fit <- function(Y) {
  cor <- duplicateCorrelation(Y, X, block = meta$donor)$consensus.correlation; if (!is.finite(cor)) cor <- 0
  eBayes(lmFit(Y, X, block = meta$donor, correlation = cor))
}
res <- list(); abund <- list()
for (g in unique(sets$group)) {
  mk <- intersect(sets[group == g & role == "marker", gene], rownames(lc))
  tg <- intersect(sets[group == g & role == "tested", gene], rownames(lc))
  tg <- tg[rowMeans(lc[tg, , drop = FALSE]) > 0]                     # bulk 中有表达
  idx <- colMeans(lc[mk, , drop = FALSE])
  fa <- fit(rbind(abundance = idx, abundance_dup = idx))
  for (cf in c("groupDiabetic", "groupNPDR")) {
    tt <- topTable(fa, coef = cf, number = Inf, sort.by = "none")["abundance", ]
    abund[[length(abund) + 1]] <- data.table(group = g, contrast = cf, n_markers = length(mk), diff = tt$logFC, t = tt$t, p = tt$P.Value)
  }
  for (method in c("bulk", "ratio")) {
    Y <- lc[tg, , drop = FALSE]; if (method == "ratio") Y <- sweep(Y, 2, idx)
    f <- fit(Y)
    for (cf in c("groupDiabetic", "groupNPDR")) {
      tt <- topTable(f, coef = cf, number = Inf, sort.by = "none")
      res[[length(res) + 1]] <- data.table(group = g, method = method, contrast = cf, gene = rownames(tt),
                                           diff = tt$logFC, t = tt$t, p = tt$P.Value, fdr = tt$adj.P.Val)
    }
  }
}
res <- rbindlist(res); abund <- rbindlist(abund)
fwrite(res[order(group, method, contrast, p)], file.path(out, "ratio_tests.tsv"), sep = "\t")
fwrite(abund, file.path(out, "abundance_index_tests.tsv"), sep = "\t")
summ <- res[, .(n_genes = .N, n_fdr05 = sum(fdr < 0.05), n_p01 = sum(p < 0.01), pi1 = 1 - propTrueNull(p),
                frac_down = mean(t < 0)), by = .(group, method, contrast)]
fwrite(summ, file.path(out, "ratio_summary.tsv"), sep = "\t")
print(abund); print(summ)
goi <- cfg$deconvolution$ratio_exclude_markers
print(dcast(res[gene %in% goi], group + gene ~ method + contrast, value.var = "t")[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])
write_provenance(stage, inputs = c(cfg$retina$counts, file.path("results", "09c_ratio_validation", "gene_sets.tsv")),
                 outputs = list.files(out, "\\.tsv$", full.names = TRUE), seed = cfg$seed)
