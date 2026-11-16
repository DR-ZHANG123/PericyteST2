# Stage 05：去除隐变量（svaseq，保护糖尿病与已知协变量）后的糖尿病效应，并用标签置换校准（先导）
# 供体心脏的转录组被脑死亡/重症监护应激等未测量因素主导，可能掩盖糖尿病信号。
# MAGNet 原始分析用 svaseq 隐变量；这里同法估计隐变量，再检验糖尿病；
# 同时把糖尿病标签在层内打乱、整套流程（含 svaseq）重做，得到“流程本身”会产生多少假信号的零分布。
suppressPackageStartupMessages({library(limma); library(edgeR); library(sva); library(data.table)})
source(file.path("scripts", "lib", "repro.R"))
cfg <- load_config(); set_global_seed(cfg$seed)
root <- p22_root(); setwd(root)
stage <- "05_sva_calibrated"
out <- file.path("results", stage); dir.create(out, recursive = TRUE, showWarnings = FALSE)
n_perm <- 20

cnt <- fread(cfg$magnet$counts); genes <- cnt[[1]]; cnt <- as.matrix(cnt[, -1]); rownames(cnt) <- genes
ph <- fread(cfg$magnet$pheno)
ph <- ph[etiology %in% cfg$magnet$etiologies & Diabetes %in% c("Yes", "No")]
ph[, bmi := weight / (height / 100)^2]
ph[, `:=`(dm = factor(Diabetes, levels = c("No", "Yes")), sex = factor(gender), race = factor(race),
          htn = factor(Hypertension, levels = c("No", "Yes")), rin = RIN, tin = TIN.median., pool = factor(Library.Pool))]
ph <- ph[!is.na(bmi) & !is.na(htn)]
ph[is.na(rin), rin := median(ph$rin, na.rm = TRUE)]
cnt <- cnt[, ph$sample_name]
dge <- DGEList(cnt); keep <- filterByExpr(dge, group = interaction(ph$etiology, ph$dm))
dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
norm_cpm <- cpm(dge)

covs <- "sex + age + race + bmi + htn + rin + tin + pool"
run_one <- function(idx, dm_vec) {
  d <- ph[idx]; d$dm <- dm_vec
  mod <- model.matrix(as.formula(paste("~ dm +", covs)), data = d)
  mod0 <- model.matrix(as.formula(paste("~", covs)), data = d)
  sv <- svaseq(norm_cpm[, idx], mod, mod0)$sv
  X <- cbind(mod, sv); colnames(X)[(ncol(mod) + 1):ncol(X)] <- paste0("SV", seq_len(ncol(sv)))
  f <- eBayes(lmFit(voom(dge[, idx], X), X)); tt <- topTable(f, coef = "dmYes", number = Inf, sort.by = "none")
  list(tt = tt, n_sv = ncol(sv))
}
n_cores <- min(20L, parallel::detectCores())
summ <- list(); tabs <- list()
for (et in c("NF", "DCM")) {   # 先算真实标签，结果立即写出
  idx <- which(ph$etiology == et); r <- run_one(idx, ph$dm[idx]); tabs[[et]] <- r$tt
  fwrite(data.table(gene_id = rownames(r$tt), r$tt)[order(P.Value)], file.path(out, sprintf("de_sva_%s.tsv", tolower(et))), sep = "\t")
  cat(sprintf("%s: n_sv=%d FDR<0.05=%d pi1=%.3f\n", et, r$n_sv, sum(r$tt$adj.P.Val < 0.05), 1 - propTrueNull(r$tt$P.Value)))
}
for (et in c("NF", "DCM")) {   # 置换：每次重新抽签、重做 svaseq；各次置换的种子由主种子派生
  idx <- which(ph$etiology == et); tt <- tabs[[et]]
  obs <- c(n_fdr05 = sum(tt$adj.P.Val < 0.05), pi1 = 1 - propTrueNull(tt$P.Value))
  perm_seeds <- cfg$seed * 1000 + seq_len(n_perm) + if (et == "NF") 0 else 500
  nul <- do.call(rbind, parallel::mclapply(perm_seeds, function(s) { set.seed(s); rp <- run_one(idx, sample(ph$dm[idx]))
    c(n_fdr05 = sum(rp$tt$adj.P.Val < 0.05), pi1 = 1 - propTrueNull(rp$tt$P.Value)) }, mc.cores = n_cores))
  summ[[et]] <- data.table(stratum = et, n_fdr05 = obs["n_fdr05"], pi1 = obs["pi1"],
                           null_n_fdr05_mean = mean(nul[, "n_fdr05"]), null_n_fdr05_max = max(nul[, "n_fdr05"]),
                           null_pi1_mean = mean(nul[, "pi1"]), null_pi1_max = max(nul[, "pi1"]),
                           p_pi1 = (1 + sum(nul[, "pi1"] >= obs["pi1"])) / (1 + n_perm))
}
summ <- rbindlist(summ); fwrite(summ, file.path(out, "sva_summary.tsv"), sep = "\t"); print(summ)
m <- merge(data.table(gene_id = rownames(tabs$NF), t_nf = tabs$NF$t, p_nf = tabs$NF$P.Value),
           data.table(gene_id = rownames(tabs$DCM), t_dcm = tabs$DCM$t), by = "gene_id")
rep <- data.table(rho_all = cor(m$t_nf, m$t_dcm, method = "spearman"),
                  sign_top200 = m[order(p_nf)][1:200, mean(sign(t_nf) == sign(t_dcm))])
fwrite(rep, file.path(out, "sva_nf_dcm_replication.tsv"), sep = "\t"); print(rep)
write_provenance(stage, inputs = c(cfg$magnet$counts, cfg$magnet$pheno),
                 outputs = list.files(out, "\\.tsv$", full.names = TRUE), seed = cfg$seed, extra = list(n_perm = n_perm))
