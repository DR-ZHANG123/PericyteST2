# Stage 06：人视网膜（GSE160306）各病程阶段的细胞程序与缺氧程序（先导）
# 问题：糖尿病在出现视网膜病变之前，视网膜哪些细胞程序已经改变？缺氧/HIF 程序何时出现？
# 设计：每位供体有黄斑与周边两份样本；病程按供体分组（对照 / 糖尿病无 DR / NPDR / PDR），
#       个体内相关用 limma duplicateCorrelation 处理（检验单位是供体）。
# 读出：细胞类型标记基因得分（组内 z 分数均值）、HIF 靶基因得分，以及基因层面的早期差异。
suppressPackageStartupMessages({library(limma); library(edgeR); library(data.table)})
source(file.path("scripts", "lib", "repro.R"))
cfg <- load_config(); set_global_seed(cfg$seed)
root <- p22_root(); setwd(root)
stage <- "06_retina_pilot"
out <- file.path("results", stage); dir.create(out, recursive = TRUE, showWarnings = FALSE)
raw <- file.path(cfg$paths$raw, "gse160306")

# ---------- 样本信息：series matrix 中每个样本的 key: value 不对齐，逐样本解析 ----------
sm <- readLines(gzfile(file.path(raw, "GSE160306_series_matrix.txt.gz")))
ch <- grep("^!Sample_characteristics_ch1", sm, value = TRUE)
cells <- lapply(strsplit(ch, "\t"), function(x) gsub('^"|"$', "", x[-1]))
n_s <- length(cells[[1]])
meta <- rbindlist(lapply(seq_len(n_s), function(j) {
  kv <- vapply(cells, `[`, "", j); kv <- kv[grepl(": ", kv)]
  as.list(setNames(sub("^[^:]+: ", "", kv), sub(": .*$", "", kv)))
}), fill = TRUE)
meta <- meta[, .(sample = sampleID, donor, tissue, sex = gender, age = as.numeric(age), group = disease_group,
                 dme = dme, severity = suppressWarnings(as.numeric(dr_severity_score)),
                 duration = duration_diabetis, rin = suppressWarnings(as.numeric(rin)),
                 pmi = suppressWarnings(as.numeric(post_mortem_interval)))]
meta[, region := ifelse(grepl("Macula", tissue), "macula", "periphery")]
meta[, group := factor(group, levels = c("Control", "Diabetic", "NPDR", "PDR"))]
meta[is.na(rin), rin := median(meta$rin, na.rm = TRUE)]
meta[is.na(pmi), pmi := median(meta$pmi, na.rm = TRUE)]
fwrite(meta, file.path(out, "sample_meta.tsv"), sep = "\t")
fwrite(meta[, .(n_samples = .N, n_donors = uniqueN(donor)), by = group][order(group)], file.path(out, "group_counts.tsv"), sep = "\t")

# ---------- 表达 ----------
cnt <- as.matrix(read.csv(gzfile(file.path(raw, "GSE160306_counts.txt.gz")), row.names = 1, check.names = FALSE))
meta <- meta[sample %in% colnames(cnt) & !is.na(group)]; cnt <- cnt[, meta$sample]
dge <- DGEList(cnt); keep <- filterByExpr(dge, group = meta$group); dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
lc <- cpm(dge, log = TRUE, prior.count = 2)
sym <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, rownames(lc), "SYMBOL", "ENSEMBL", multiVals = "first")
ok <- !is.na(sym) & !duplicated(sym); lcs <- lc[ok, ]; rownames(lcs) <- sym[ok]

markers <- list(
  rod = c("RHO", "GNAT1", "PDE6A", "PDE6B", "CNGA1", "SAG", "NRL", "NR2E3"),
  cone = c("ARR3", "GNAT2", "PDE6C", "PDE6H", "OPN1SW", "CNGB3", "GUCA1C"),
  bipolar = c("VSX2", "CABP5", "PRKCA", "GRM6", "TRPM1", "SCGN", "VSX1"),
  amacrine = c("GAD1", "GAD2", "SLC6A9", "TFAP2B", "TFAP2A", "C1QL2"),
  ganglion = c("RBPMS", "POU4F1", "POU4F2", "SNCG", "NEFL", "NEFM", "SLC17A6"),
  horizontal = c("ONECUT1", "ONECUT2", "LHX1", "CALB1"),
  muller = c("RLBP1", "GLUL", "SLC1A3", "WIF1", "DKK3", "CRABP1"),
  gliosis = c("GFAP", "CHI3L1", "LCN2", "SERPINA3", "CD44"),
  microglia = c("C1QA", "C1QB", "P2RY12", "CX3CR1", "TMEM119", "CSF1R", "AIF1"),
  endothelial = c("PECAM1", "CDH5", "CLDN5", "VWF", "FLT1", "KDR", "ESAM"),
  pericyte = c("PDGFRB", "RGS5", "KCNJ8", "ABCC9", "NOTCH3"),
  hif_targets = c("VEGFA", "ANGPTL4", "ADM", "SLC2A1", "BNIP3", "PDK1", "NDRG1", "EGLN3", "P4HA1", "ANKRD37", "CA9"),
  interferon = c("IFIT1", "IFIT3", "ISG15", "MX1", "OAS1", "RSAD2", "IFI44L"))
z <- t(scale(t(lcs)))
S <- t(vapply(markers, function(m) colMeans(z[intersect(m, rownames(z)), , drop = FALSE]), numeric(ncol(z))))
fwrite(data.table(score = rownames(S), S), file.path(out, "cell_scores.tsv"), sep = "\t")

fit_scores <- function(Y, form, coefs) {
  X <- model.matrix(form, data = meta)
  dc <- duplicateCorrelation(Y, X, block = meta$donor)$consensus.correlation
  f <- eBayes(lmFit(Y, X, block = meta$donor, correlation = dc))
  rbindlist(lapply(coefs, function(cf) { tt <- topTable(f, coef = cf, number = Inf, sort.by = "none")
    data.table(feature = rownames(tt), contrast = cf, estimate = tt$logFC, t = tt$t, p = tt$P.Value, fdr = tt$adj.P.Val) }))[, icc := dc][]
}
# 各阶段相对对照；校正部位、年龄、性别、RIN、死亡后间隔
res_grp <- fit_scores(S, ~ group + region + age + sex + rin + pmi, c("groupDiabetic", "groupNPDR", "groupPDR"))
meta[, stage_num := as.integer(group) - 1L]
res_trend <- fit_scores(S, ~ stage_num + region + age + sex + rin + pmi, "stage_num")
# 部位差异：糖尿病（任何阶段）效应在黄斑与周边是否不同
meta[, anydm := as.integer(group != "Control")]
res_reg <- fit_scores(S, ~ anydm * region + age + sex + rin + pmi, c("anydm", "anydm:regionperiphery"))
res <- rbind(res_grp, res_trend, res_reg)
fwrite(res, file.path(out, "cell_score_tests.tsv"), sep = "\t")
print(dcast(res, feature ~ contrast, value.var = "t")[order(-abs(stage_num))])
print(dcast(res, feature ~ contrast, value.var = "p")[order(stage_num)])

# 基因层面：糖尿病无 DR vs 对照（最早阶段），以及阶段趋势
X <- model.matrix(~ group + region + age + sex + rin + pmi, data = meta)
v <- voom(dge, X); dc <- duplicateCorrelation(v, X, block = meta$donor)$consensus.correlation
v <- voom(dge, X, block = meta$donor, correlation = dc)
f <- eBayes(lmFit(v, X, block = meta$donor, correlation = dc))
gl <- rbindlist(lapply(c("groupDiabetic", "groupNPDR", "groupPDR"), function(cf) {
  tt <- topTable(f, coef = cf, number = Inf, sort.by = "none")
  data.table(gene_id = rownames(tt), symbol = sym[rownames(tt)], contrast = cf, logFC = tt$logFC, t = tt$t, p = tt$P.Value, fdr = tt$adj.P.Val) }))
fwrite(gl[order(contrast, p)], file.path(out, "gene_tests.tsv"), sep = "\t")
print(gl[, .(n_fdr05 = sum(fdr < 0.05), pi1 = 1 - propTrueNull(p)), by = contrast])
write_provenance(stage, inputs = list.files(raw, full.names = TRUE), outputs = list.files(out, "\\.tsv$", full.names = TRUE),
                 seed = cfg$seed, extra = list(donor_icc_genes = dc))
