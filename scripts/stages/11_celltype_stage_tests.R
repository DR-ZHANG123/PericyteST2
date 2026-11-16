# Stage 11：在供体层面检验各病程阶段中，哪些视网膜细胞类型的比例或自身表达发生改变
# 输入：stage 10 的 TAPE 细胞比例与细胞类型特异表达、stage 06 的样本信息
# 1) 阳性对照：黄斑 vs 周边的细胞比例差异（视锥、节细胞应在黄斑更高，视杆在周边更高）
# 2) 细胞比例 ~ 病程 + 部位 + 年龄 + 性别 + RIN + 死亡后间隔（供体内相关用 duplicateCorrelation）
# 3) 各细胞类型自身表达的全基因组检验：FDR<0.05 基因数、pi1；病程标签在供体间置换（同一供体两份样本一起换）校准 pi1
suppressPackageStartupMessages({library(limma); library(data.table); library(parallel)})
source(file.path("scripts", "lib", "repro.R"))
cfg <- load_config(); set_global_seed(cfg$seed)
root <- p22_root(); setwd(root)
stage <- "11_celltype_stage_tests"
out <- file.path("results", stage); dir.create(out, recursive = TRUE, showWarnings = FALSE)
n_perm <- 200; n_cores <- 20L

meta <- fread(file.path("results", "06_retina_pilot", "sample_meta.tsv"))
meta <- meta[group %in% c("Control", "Diabetic", "NPDR")]           # PDR 仅 2 位供体，不参与检验
meta[, group := factor(group, levels = c("Control", "Diabetic", "NPDR"))]
meta[, region := factor(region, levels = c("macula", "periphery"))]
covs <- "region + age + sex + rin + pmi"

fit_block <- function(Y, d, rhs, coefs) {
  X <- model.matrix(as.formula(paste("~", rhs)), data = d)
  cor <- duplicateCorrelation(Y, X, block = d$donor)$consensus.correlation
  if (!is.finite(cor)) cor <- 0
  f <- eBayes(lmFit(Y, X, block = d$donor, correlation = cor))
  list(fit = f, cor = cor, tabs = lapply(setNames(coefs, coefs), function(cf) topTable(f, coef = cf, number = Inf, sort.by = "none")))
}

# ---------- 细胞比例 ----------
frac_tests <- rbindlist(lapply(c("tape", "nnls"), function(m) {
  fr <- fread(file.path("results", "10_deconv_retina", sprintf("fractions_%s.tsv", m)))
  setnames(fr, 1, "sample"); fr <- fr[sample %in% meta$sample]
  d <- meta[match(fr$sample, meta$sample)]
  P <- as.matrix(fr[, -1]); P <- pmin(pmax(P, 1e-4), 1 - 1e-4)
  Y <- t(qlogis(P)); colnames(Y) <- fr$sample
  r <- fit_block(Y, d, paste("group +", covs), c("groupDiabetic", "groupNPDR", "regionperiphery"))
  rbindlist(lapply(names(r$tabs), function(cf) data.table(method = m, contrast = cf, cell_class = rownames(r$tabs[[cf]]),
    logit_diff = r$tabs[[cf]]$logFC, t = r$tabs[[cf]]$t, p = r$tabs[[cf]]$P.Value, fdr = r$tabs[[cf]]$adj.P.Val,
    mean_fraction = rowMeans(P)[rownames(r$tabs[[cf]])])))
}))
fwrite(frac_tests, file.path(out, "fraction_tests.tsv"), sep = "\t")
print(dcast(frac_tests, cell_class ~ method + contrast, value.var = "t")[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])

# ---------- 细胞类型自身表达 ----------
ctdir <- cfg$deconvolution$ctexpr_dir
classes <- sub("\\.tsv\\.gz$", "", list.files(ctdir, "\\.tsv\\.gz$"))
targets <- c("RGS5", "NOTCH3", "PDGFRB", "KCNJ8", "DLL4", "TEK", "CDH5", "CLDN5", "PLVAP", "PECAM1", "JAG1", "PDGFB",
             "RLBP1", "VIM", "GFAP", "GLUL", "TXNIP", "COL4A1", "NID1", "APLN")
donors <- unique(meta[, .(donor, group)])
summ <- list(); targ <- list()
for (cl in classes) {
  E <- fread(cmd = paste("zcat", shQuote(file.path(ctdir, paste0(cl, ".tsv.gz"))))); setnames(E, 1, "sample"); E <- E[sample %in% meta$sample]
  d <- meta[match(E$sample, meta$sample)]
  Y <- t(as.matrix(E[, -1])); colnames(Y) <- E$sample
  Y <- Y[apply(Y, 1, sd) > 0, ]
  r <- fit_block(Y, d, paste("group +", covs), c("groupDiabetic", "groupNPDR"))
  obs_pi1 <- sapply(r$tabs, function(tt) 1 - propTrueNull(tt$P.Value))
  # 置换：在供体间打乱病程标签，同一供体的黄斑与周边一起移动
  perm <- do.call(rbind, mclapply(seq_len(n_perm), function(b) {
    set.seed(cfg$seed * 10000 + b)
    sh <- copy(donors); sh[, group := sample(group)]
    dp <- copy(d); dp[, group := sh$group[match(dp$donor, sh$donor)]]
    X <- model.matrix(as.formula(paste("~ group +", covs)), data = dp)
    f <- eBayes(lmFit(Y, X, block = dp$donor, correlation = r$cor))
    sapply(c("groupDiabetic", "groupNPDR"), function(cf) 1 - propTrueNull(topTable(f, coef = cf, number = Inf, sort.by = "none")$P.Value))
  }, mc.cores = n_cores))
  for (cf in names(r$tabs)) {
    tt <- r$tabs[[cf]]
    summ[[length(summ) + 1]] <- data.table(cell_class = cl, contrast = cf, n_genes = nrow(tt),
      n_fdr05 = sum(tt$adj.P.Val < 0.05), pi1 = obs_pi1[cf], null_pi1_mean = mean(perm[, cf]),
      null_pi1_q95 = quantile(perm[, cf], 0.95), p_perm = (1 + sum(perm[, cf] >= obs_pi1[cf])) / (1 + n_perm), donor_cor = r$cor)
    hit <- intersect(targets, rownames(tt))
    if (length(hit)) targ[[length(targ) + 1]] <- data.table(cell_class = cl, contrast = cf, gene = hit,
      diff = tt[hit, "logFC"], t = tt[hit, "t"], p = tt[hit, "P.Value"], fdr = tt[hit, "adj.P.Val"])
    fwrite(data.table(gene = rownames(tt), diff = tt$logFC, t = tt$t, p = tt$P.Value, fdr = tt$adj.P.Val)[order(p)],
           file.path(out, sprintf("ctexpr_%s_%s.tsv", cl, sub("group", "", cf))), sep = "\t")
  }
}
summ <- rbindlist(summ); targ <- rbindlist(targ)
fwrite(summ, file.path(out, "ctexpr_summary.tsv"), sep = "\t"); fwrite(targ, file.path(out, "ctexpr_target_genes.tsv"), sep = "\t")
print(summ[order(contrast, -pi1)])
print(dcast(targ, gene ~ cell_class + contrast, value.var = "t")[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 1) else x)])
write_provenance(stage, inputs = c(list.files(file.path("results", "10_deconv_retina"), "\\.tsv$", full.names = TRUE),
                                   list.files(ctdir, full.names = TRUE)),
                 outputs = list.files(out, "\\.tsv$", full.names = TRUE), seed = cfg$seed, extra = list(n_perm = n_perm))
