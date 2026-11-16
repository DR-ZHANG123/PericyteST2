# Stage 01：MAGNet 左心室中糖尿病相关的转录差异（先导分析）
# 问题：糖尿病是否在未衰竭和衰竭的人类左心室中留下同一套转录程序？
#   - 合并模型：糖尿病 + 病因(NF/DCM) + 性别 + 年龄 + 种族 + BMI + 高血压 + RIN + TIN + 文库批次
#   - 分层模型：NF 与 DCM 各自单独拟合，两组个体互不重叠，用来互相复现
#   - 交互模型：糖尿病效应在 NF 与 DCM 之间是否不同
#   - 置换零分布：NF 内部打乱糖尿病标签，检验 NF→DCM 的方向一致性是否超出偶然
# 输出：results/01_magnet_diabetes/*.tsv + PROVENANCE.json
suppressPackageStartupMessages({
  library(limma); library(edgeR); library(data.table); library(AnnotationDbi); library(org.Hs.eg.db)
})
source(file.path("scripts", "lib", "repro.R"))
cfg <- load_config(); set_global_seed(cfg$seed)
root <- p22_root(); setwd(root)
stage <- "01_magnet_diabetes"
out <- file.path("results", stage); dir.create(out, recursive = TRUE, showWarnings = FALSE)

# ---------- 数据 ----------
cnt <- fread(cfg$magnet$counts); genes <- cnt[[1]]; cnt <- as.matrix(cnt[, -1]); rownames(cnt) <- genes
ph <- fread(cfg$magnet$pheno)
ph <- ph[etiology %in% cfg$magnet$etiologies & Diabetes %in% c("Yes", "No")]
ph[, bmi := weight / (height / 100)^2]
ph[, `:=`(dm = factor(Diabetes, levels = c("No", "Yes")), etiology = factor(etiology, levels = c("NF", "DCM")),
          sex = factor(gender), race = factor(race), htn = factor(Hypertension, levels = c("No", "Yes")),
          rin = RIN, tin = TIN.median., pool = factor(Library.Pool))]
ph <- ph[!is.na(bmi) & !is.na(htn)]
ph[is.na(rin), rin := median(ph$rin, na.rm = TRUE)]   # 6 例 RIN 缺失，用中位数填补并在补充方法说明
cnt <- cnt[, ph$sample_name]
sym <- mapIds(org.Hs.eg.db, keys = rownames(cnt), keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")

cat(sprintf("样本 %d（NF %d：DM %d；DCM %d：DM %d）\n", nrow(ph), sum(ph$etiology == "NF"),
            sum(ph$etiology == "NF" & ph$dm == "Yes"), sum(ph$etiology == "DCM"),
            sum(ph$etiology == "DCM" & ph$dm == "Yes")))
fwrite(ph[, .N, by = .(etiology, dm)][order(etiology, dm)], file.path(out, "sample_counts.tsv"), sep = "\t")

# 每组临床特征（供体层面描述）
desc <- ph[, .(n = .N, age = median(age), female = mean(sex == "Female"), bmi = median(bmi),
               htn = mean(htn == "Yes"), black = mean(race == "AA"), lvef = median(LVEF, na.rm = TRUE),
               rin = median(rin)), by = .(etiology, dm)][order(etiology, dm)]
fwrite(desc, file.path(out, "group_characteristics.tsv"), sep = "\t")

# ---------- 过滤与归一化 ----------
dge <- DGEList(cnt)
keep <- filterByExpr(dge, group = interaction(ph$etiology, ph$dm))
dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
cat(sprintf("保留基因 %d / %d\n", sum(keep), length(keep)))

covs <- "sex + age + race + bmi + htn + rin + tin + pool"
fit_dm <- function(idx, formula_rhs, coef = "dmYes") {
  d <- ph[idx]
  X <- model.matrix(as.formula(paste("~", formula_rhs)), data = d)
  v <- voom(dge[, idx], X)
  f <- eBayes(lmFit(v, X))
  tt <- topTable(f, coef = coef, number = Inf, sort.by = "none")
  data.table(gene_id = rownames(tt), symbol = sym[rownames(tt)], logFC = tt$logFC,
             t = tt$t, p = tt$P.Value, fdr = tt$adj.P.Val, ave = tt$AveExpr)
}

all_idx <- seq_len(nrow(ph)); nf_idx <- which(ph$etiology == "NF"); dcm_idx <- which(ph$etiology == "DCM")
res_pool <- fit_dm(all_idx, paste("dm + etiology +", covs))
res_nf   <- fit_dm(nf_idx,  paste("dm +", covs))
res_dcm  <- fit_dm(dcm_idx, paste("dm +", covs))
res_int  <- fit_dm(all_idx, paste("dm * etiology +", covs), coef = "dmYes:etiologyDCM")
for (nm in c("pool", "nf", "dcm", "int")) {
  r <- get(paste0("res_", nm)); setorder(r, p)
  fwrite(r, file.path(out, sprintf("de_%s.tsv", nm)), sep = "\t")
}

# ---------- NF 与 DCM 两组独立个体之间的复现 ----------
m <- merge(res_nf[, .(gene_id, symbol, t_nf = t, p_nf = p, fc_nf = logFC)],
           res_dcm[, .(gene_id, t_dcm = t, p_dcm = p, fc_dcm = logFC)], by = "gene_id")
rho_all <- cor(m$t_nf, m$t_dcm, method = "spearman")
concord <- function(mm, top) {
  s <- mm[order(p_nf)][seq_len(top)]
  c(sign_agree = mean(sign(s$t_nf) == sign(s$t_dcm)), rho = cor(s$t_nf, s$t_dcm, method = "spearman"))
}
tops <- c(50, 100, 200, 500)
obs <- t(sapply(tops, function(k) concord(m, k)))

# 置换：只打乱 NF 内的糖尿病标签（保留协变量），DCM 结果不动
n_perm <- 200
perm <- replicate(n_perm, {
  ph_p <- copy(ph); ph_p$dm[nf_idx] <- sample(ph_p$dm[nf_idx])
  d <- ph_p[nf_idx]; X <- model.matrix(as.formula(paste("~ dm +", covs)), data = d)
  f <- eBayes(lmFit(voom(dge[, nf_idx], X), X)); tt <- topTable(f, coef = "dmYes", number = Inf, sort.by = "none")
  mm <- merge(data.table(gene_id = rownames(tt), t_nf = tt$t, p_nf = tt$P.Value), m[, .(gene_id, t_dcm)], by = "gene_id")
  c(sapply(tops, function(k) concord(mm, k)["sign_agree"]), cor(mm$t_nf, mm$t_dcm, method = "spearman"))
})
conc <- data.table(top_n = c(tops, NA), statistic = c(rep("sign_agree", length(tops)), "spearman_all"),
                   observed = c(obs[, "sign_agree"], rho_all),
                   null_mean = rowMeans(perm),
                   null_p = (1 + rowSums(perm >= c(obs[, "sign_agree"], rho_all))) / (1 + n_perm))
fwrite(conc, file.path(out, "nf_dcm_concordance.tsv"), sep = "\t")
fwrite(m, file.path(out, "nf_vs_dcm_t.tsv"), sep = "\t")

# ---------- 候选基因 ----------
cand <- merge(m[symbol %in% cfg$magnet$candidate_genes],
              res_pool[, .(gene_id, fc_pool = logFC, p_pool = p, fdr_pool = fdr)], by = "gene_id")
fwrite(cand[order(p_pool)], file.path(out, "candidate_genes.tsv"), sep = "\t")

# ---------- 汇总 ----------
summ <- data.table(model = c("pool", "nf", "dcm", "interaction"),
                   n_fdr05 = sapply(list(res_pool, res_nf, res_dcm, res_int), function(r) sum(r$fdr < cfg$magnet$fdr)),
                   n_p001 = sapply(list(res_pool, res_nf, res_dcm, res_int), function(r) sum(r$p < 0.001)),
                   pi1 = sapply(list(res_pool, res_nf, res_dcm, res_int), function(r) 1 - limma::propTrueNull(r$p)))
fwrite(summ, file.path(out, "model_summary.tsv"), sep = "\t")
print(summ); print(conc); print(head(res_pool[order(p)], 30)); print(cand[order(p_pool)])

outs <- list.files(out, pattern = "\\.tsv$", full.names = TRUE)
write_provenance(stage, inputs = c(cfg$magnet$counts, cfg$magnet$pheno), outputs = outs, seed = cfg$seed,
                 extra = list(n_samples = nrow(ph), n_genes = sum(keep), n_perm = n_perm))
