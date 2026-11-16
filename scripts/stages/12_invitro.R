# Stage 12：人视网膜细胞体外模型中，糖尿病相关刺激能否重现供体视网膜的改变
#   - GSE204880：人视网膜周细胞（HRPC）BSA 对照 vs 晚期糖基化终产物（AGEs），各 3 个（标准化表达值）
#   - GSE290024：人 iPSC 视网膜类器官，17 mM vs 25 mM 葡萄糖 28 天，各 3 个（计数）
# 输出目标基因的组间差异（limma），以及全基因组 pi1。
suppressPackageStartupMessages({library(limma); library(edgeR); library(data.table)})
source(file.path("scripts", "lib", "repro.R"))
cfg <- load_config(); set_global_seed(cfg$seed)
root <- p22_root(); setwd(root)
stage <- "12_invitro"
out <- file.path("results", stage); dir.create(out, recursive = TRUE, showWarnings = FALSE)
gm <- fread(cfg$deconvolution$gene_map); setnames(gm, 1:2, c("id", "gene")); gm[, id := sub("\\..*$", "", id)]
gm <- gm[!duplicated(id)]
to_sym <- function(M) { s <- gm$gene[match(rownames(M), gm$id)]; M <- M[!is.na(s), , drop = FALSE]; rowsum(M, s[!is.na(s)]) }
targets <- c("RGS5", "NOTCH3", "PDGFRB", "KCNJ8", "ABCC9", "DES", "TXNIP", "RLBP1", "VIM", "GFAP", "GLUL", "SLC1A3",
             "CRABP1", "DKK3", "VEGFA", "ARR3", "RHO", "RBPMS")

tab <- function(f, coef, model) { tt <- topTable(f, coef = coef, number = Inf, sort.by = "none")
  list(t = data.table(model = model, gene = rownames(tt), diff = tt$logFC, t = tt$t, p = tt$P.Value, fdr = tt$adj.P.Val),
       s = data.table(model = model, n_genes = nrow(tt), n_fdr05 = sum(tt$adj.P.Val < 0.05), pi1 = 1 - propTrueNull(tt$P.Value))) }

# 周细胞 + AGEs
fs <- sort(list.files(file.path(cfg$paths$raw, "gse204880"), "HRPC.*txt.gz$", full.names = TRUE))
M <- do.call(cbind, lapply(fs, function(f) { x <- read.delim(gzfile(f), row.names = 1); setNames(x[[1]], rownames(x)) }))
colnames(M) <- sub("^GSM[0-9]+_|\\.txt\\.gz$", "", basename(fs))
M <- to_sym(M); L <- log2(M + 1); L <- L[rowMeans(L) > 1, ]
grp <- factor(ifelse(grepl("AGEs", colnames(L)), "AGEs", "control"), levels = c("control", "AGEs"))
f1 <- eBayes(lmFit(L, model.matrix(~ grp)), trend = TRUE); r1 <- tab(f1, "grpAGEs", "pericyte_AGEs")

# 类器官 + 高糖
o <- read.delim(gzfile(file.path(cfg$paths$raw, "gse290024", "GSE290024_RAW_counts.txt.gz")), row.names = 1, check.names = FALSE)
o <- to_sym(as.matrix(o)); g2 <- factor(ifelse(grepl("glucose", colnames(o)), "glucose25", "control"), levels = c("control", "glucose25"))
dge <- DGEList(o); dge <- calcNormFactors(dge[filterByExpr(dge, group = g2), , keep.lib.sizes = FALSE])
X2 <- model.matrix(~ g2); f2 <- eBayes(lmFit(voom(dge, X2), X2)); r2 <- tab(f2, "g2glucose25", "organoid_glucose")

all_t <- rbind(r1$t, r2$t); fwrite(all_t[order(model, p)], file.path(out, "invitro_de.tsv"), sep = "\t")
fwrite(rbind(r1$s, r2$s), file.path(out, "invitro_summary.tsv"), sep = "\t")
fwrite(all_t[gene %in% targets][order(model, gene)], file.path(out, "invitro_targets.tsv"), sep = "\t")
print(rbind(r1$s, r2$s)); print(dcast(all_t[gene %in% targets], gene ~ model, value.var = "t"))
write_provenance(stage, inputs = c(fs, file.path(cfg$paths$raw, "gse290024", "GSE290024_RAW_counts.txt.gz")),
                 outputs = list.files(out, "\\.tsv$", full.names = TRUE), seed = cfg$seed)
