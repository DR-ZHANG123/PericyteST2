# Stage 02：MAGNet 糖尿病效应的通路层面分析（先导）
# 单基因效应小时，糖尿病的影响可能分散在整条通路上。这里用 stage 01 的 t 统计量做预排序 GSEA，
# 分别在 NF、DCM 和合并模型中估计每条通路的富集分数，并看 NF 与 DCM（互不重叠的个体）是否一致。
# 基因集：MSigDB Hallmark 与 Reactome（msigdbr），所用基因集快照写到 data/interim/genesets.tsv。
suppressPackageStartupMessages({library(data.table); library(fgsea); library(msigdbr)})
source(file.path("scripts", "lib", "repro.R"))
cfg <- load_config(); set_global_seed(cfg$seed)
root <- p22_root(); setwd(root)
stage <- "02_magnet_genesets"; in_stage <- "01_magnet_diabetes"
out <- file.path("results", stage); dir.create(out, recursive = TRUE, showWarnings = FALSE)

gs_file <- file.path(cfg$paths$interim, "genesets.tsv")
if (!file.exists(gs_file)) {
  h <- msigdbr(species = "Homo sapiens", collection = "H")
  r <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
  gs <- rbind(as.data.table(h)[, .(collection = "H", gs_name, gene_symbol)],
              as.data.table(r)[, .(collection = "REACTOME", gs_name, gene_symbol)])
  dir.create(dirname(gs_file), recursive = TRUE, showWarnings = FALSE)
  fwrite(unique(gs), gs_file, sep = "\t")
}
gs <- fread(gs_file)
pathways <- split(gs$gene_symbol, gs$gs_name)
coll <- unique(gs[, .(gs_name, collection)])

rank_from <- function(model) {
  d <- fread(file.path("results", in_stage, sprintf("de_%s.tsv", model)))[!is.na(symbol) & symbol != ""]
  d <- d[order(-abs(t))][!duplicated(symbol)]
  setNames(d$t, d$symbol)
}
run <- function(model) {
  st <- rank_from(model)
  res <- fgseaMultilevel(pathways, st, minSize = 15, maxSize = 500, eps = 0, nproc = 1)
  data.table(model = model, pathway = res$pathway, size = res$size, NES = res$NES, p = res$pval, padj = res$padj,
             leading = vapply(res$leadingEdge, function(x) paste(head(x, 15), collapse = ","), ""))
}
res <- rbindlist(lapply(c("pool", "nf", "dcm"), run))
res <- merge(res, coll, by.x = "pathway", by.y = "gs_name")
res[, padj := p.adjust(p, "BH"), by = .(model, collection)]
fwrite(res[order(model, p)], file.path(out, "gsea_all.tsv"), sep = "\t")

w <- dcast(res, pathway + collection ~ model, value.var = c("NES", "padj"))
conc <- w[, .(n = .N, rho_nf_dcm = cor(NES_nf, NES_dcm, method = "spearman", use = "complete.obs")), by = collection]
both <- w[padj_nf < 0.05 & padj_dcm < 0.05]
both[, same_direction := sign(NES_nf) == sign(NES_dcm)]
fwrite(conc, file.path(out, "nes_concordance.tsv"), sep = "\t")
fwrite(both[order(-abs(NES_pool))], file.path(out, "pathways_sig_both.tsv"), sep = "\t")
fwrite(w[order(padj_pool)], file.path(out, "nes_wide.tsv"), sep = "\t")
print(conc)
print(w[collection == "H"][order(padj_pool)][1:25, .(pathway, NES_nf, padj_nf, NES_dcm, padj_dcm, NES_pool, padj_pool)])
print(both[order(-abs(NES_pool))][1:40, .(pathway, NES_nf, padj_nf, NES_dcm, padj_dcm)])

write_provenance(stage, inputs = c(list.files(file.path("results", in_stage), "^de_.*tsv$", full.names = TRUE), gs_file),
                 outputs = list.files(out, "\\.tsv$", full.names = TRUE), seed = cfg$seed,
                 extra = list(msigdbr = as.character(packageVersion("msigdbr"))))
