# Stage 04：左心室转录组对各类个体特征的敏感度（先导）
# 目的：判断糖尿病信号弱是“数据检测不到”还是“糖尿病本身影响小”。
#   - 阳性对照：性别、年龄、BMI、高血压在同一模型、同一样本中能检测到多少差异
#   - 糖尿病的三种模型：完全校正 / 不校正 BMI 与高血压 / 只校正技术变量
# 指标：FDR<0.05 基因数、pi1（1 - 真零假设比例）、以及效应量分布（|logFC| 的 95% 分位）
suppressPackageStartupMessages({library(limma); library(edgeR); library(data.table)})
source(file.path("scripts", "lib", "repro.R"))
cfg <- load_config(); set_global_seed(cfg$seed)
root <- p22_root(); setwd(root)
stage <- "04_effect_landscape"
out <- file.path("results", stage); dir.create(out, recursive = TRUE, showWarnings = FALSE)

cnt <- fread(cfg$magnet$counts); genes <- cnt[[1]]; cnt <- as.matrix(cnt[, -1]); rownames(cnt) <- genes
ph <- fread(cfg$magnet$pheno)
ph <- ph[etiology %in% cfg$magnet$etiologies & Diabetes %in% c("Yes", "No")]
ph[, bmi := weight / (height / 100)^2]
ph[, `:=`(dm = factor(Diabetes, levels = c("No", "Yes")), sex = factor(gender, levels = c("Male", "Female")),
          race = factor(race), htn = factor(Hypertension, levels = c("No", "Yes")), rin = RIN, tin = TIN.median.,
          pool = factor(Library.Pool))]
ph <- ph[!is.na(bmi) & !is.na(htn)]
ph[is.na(rin), rin := median(ph$rin, na.rm = TRUE)]
cnt <- cnt[, ph$sample_name]
dge <- DGEList(cnt); keep <- filterByExpr(dge, group = interaction(ph$etiology, ph$dm))
dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])

tech <- "rin + tin + pool"
models <- list(
  dm_full      = list(rhs = paste("dm + sex + age + race + bmi + htn +", tech), coef = "dmYes"),
  dm_no_bmihtn = list(rhs = paste("dm + sex + age + race +", tech), coef = "dmYes"),
  dm_tech_only = list(rhs = paste("dm +", tech), coef = "dmYes"),
  sex          = list(rhs = paste("sex + dm + age + race + bmi + htn +", tech), coef = "sexFemale"),
  age_10y      = list(rhs = paste("age10 + sex + dm + race + bmi + htn +", tech), coef = "age10"),
  bmi_5u       = list(rhs = paste("bmi5 + sex + age + dm + race + htn +", tech), coef = "bmi5"),
  htn          = list(rhs = paste("htn + sex + age + dm + race + bmi +", tech), coef = "htnYes"),
  race         = list(rhs = paste("race + sex + age + dm + bmi + htn +", tech), coef = "raceCaucasian")
)
ph[, `:=`(age10 = age / 10, bmi5 = bmi / 5)]

rows <- list(); keep_t <- list()
for (et in c("NF", "DCM")) {
  idx <- which(ph$etiology == et); d <- ph[idx]
  for (nm in names(models)) {
    X <- model.matrix(as.formula(paste("~", models[[nm]]$rhs)), data = d)
    f <- eBayes(lmFit(voom(dge[, idx], X), X)); tt <- topTable(f, coef = models[[nm]]$coef, number = Inf, sort.by = "none")
    rows[[length(rows) + 1]] <- data.table(stratum = et, model = nm, n = nrow(d),
      n_fdr05 = sum(tt$adj.P.Val < 0.05), pi1 = 1 - propTrueNull(tt$P.Value),
      q95_abs_logfc = quantile(abs(tt$logFC), 0.95), df_resid = f$df.residual[1])
    keep_t[[paste(et, nm, sep = ":")]] <- setNames(tt$t, rownames(tt))
  }
}
land <- rbindlist(rows)
fwrite(land, file.path(out, "effect_landscape.tsv"), sep = "\t")
tmat <- as.data.table(keep_t); tmat[, gene_id := names(keep_t[[1]])]
fwrite(tmat, file.path(out, "t_by_model.tsv"), sep = "\t")
# 同一特征在 NF 与 DCM 两组个体中的全基因组一致性（作为可复现性的阳性对照）
rep <- rbindlist(lapply(names(models), function(nm)
  data.table(model = nm, rho_nf_dcm = cor(keep_t[[paste("NF", nm, sep = ":")]], keep_t[[paste("DCM", nm, sep = ":")]], method = "spearman"))))
fwrite(rep, file.path(out, "nf_dcm_replication.tsv"), sep = "\t")
print(land); print(rep)
write_provenance(stage, inputs = c(cfg$magnet$counts, cfg$magnet$pheno),
                 outputs = list.files(out, "\\.tsv$", full.names = TRUE), seed = cfg$seed)
