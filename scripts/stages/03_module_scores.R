# Stage 03：以个体为单位检验通路程序与糖尿病的关系（先导）
# 预排序 GSEA 对高度共表达的基因集（核糖体、氧化磷酸化）给出偏小的 p 值，因为它把基因当独立单位。
# 这里改用每个人的通路得分（组内 z 分数后取均值），在 NF 与 DCM 中分别做协变量校正的线性模型，
# 显著性用 Freedman–Lane 置换（打乱简化模型残差），检验单位是个体。
suppressPackageStartupMessages({library(limma); library(edgeR); library(data.table)})
source(file.path("scripts", "lib", "repro.R"))
cfg <- load_config(); set_global_seed(cfg$seed)
root <- p22_root(); setwd(root)
stage <- "03_module_scores"
out <- file.path("results", stage); dir.create(out, recursive = TRUE, showWarnings = FALSE)
n_perm <- 2000

cnt <- fread(cfg$magnet$counts); genes <- cnt[[1]]; cnt <- as.matrix(cnt[, -1]); rownames(cnt) <- genes
ph <- fread(cfg$magnet$pheno)
ph <- ph[etiology %in% cfg$magnet$etiologies & Diabetes %in% c("Yes", "No")]
ph[, bmi := weight / (height / 100)^2]
ph[, `:=`(dm = as.numeric(Diabetes == "Yes"), sex = factor(gender), race = factor(race),
          htn = factor(Hypertension, levels = c("No", "Yes")), rin = RIN, tin = TIN.median., pool = factor(Library.Pool))]
ph <- ph[!is.na(bmi) & !is.na(htn)]
ph[is.na(rin), rin := median(ph$rin, na.rm = TRUE)]
cnt <- cnt[, ph$sample_name]

sym <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, rownames(cnt), "SYMBOL", "ENSEMBL", multiVals = "first")
dge <- DGEList(cnt); keep <- filterByExpr(dge, group = interaction(ph$etiology, ph$dm))
dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
lc <- cpm(dge, log = TRUE, prior.count = 2)
s <- sym[rownames(lc)]; ok <- !is.na(s) & !duplicated(s); lc <- lc[ok, ]; rownames(lc) <- s[ok]

gs <- fread(file.path(cfg$paths$interim, "genesets.tsv"))
sets <- split(gs$gene_symbol, gs$gs_name)
sets <- lapply(sets, intersect, rownames(lc)); sets <- sets[lengths(sets) >= 15 & lengths(sets) <= 500]

score_stratum <- function(et) {
  idx <- ph$etiology == et; d <- ph[idx]; z <- t(scale(t(lc[, idx])))
  S <- t(vapply(sets, function(g) colMeans(z[g, , drop = FALSE]), numeric(sum(idx))))   # 通路 × 个体
  X0 <- model.matrix(~ sex + age + race + bmi + htn + rin + tin + pool, data = d)
  X <- cbind(X0, dm = d$dm)
  tstat <- function(Y) { f <- eBayes(lmFit(Y, X), trend = FALSE); f$t[, "dm"] }   # 通路数 > 个体数不成问题
  t_obs <- tstat(S)
  fit0 <- lmFit(S, X0); fitted0 <- fitted(fit0); res0 <- S - fitted0
  null_max <- numeric(n_perm); ge <- numeric(length(t_obs))
  for (b in seq_len(n_perm)) {
    tb <- tstat(fitted0 + res0[, sample(ncol(res0))])
    ge <- ge + (abs(tb) >= abs(t_obs)); null_max[b] <- max(abs(tb))
  }
  data.table(stratum = et, pathway = names(t_obs), size = lengths(sets)[names(t_obs)], t = t_obs,
             p_perm = (ge + 1) / (n_perm + 1),
             p_fwer = vapply(abs(t_obs), function(x) (sum(null_max >= x) + 1) / (n_perm + 1), 0),
             delta_sd = coef(lmFit(S, X))[, "dm"])
}
res <- rbindlist(lapply(c("NF", "DCM"), score_stratum))
res[, q_perm := p.adjust(p_perm, "BH"), by = stratum]
w <- dcast(res, pathway + size ~ stratum, value.var = c("t", "p_perm", "q_perm", "p_fwer", "delta_sd"))
fwrite(res, file.path(out, "module_tests_long.tsv"), sep = "\t")
fwrite(w[order(pmax(p_perm_NF, p_perm_DCM))], file.path(out, "module_tests_wide.tsv"), sep = "\t")

# 两层之间的一致性：所有通路 t 的秩相关；两层都 p<0.05 且方向一致的通路数
conc <- data.table(rho_t = cor(w$t_NF, w$t_DCM, method = "spearman"),
                   n_both_p05_same = w[p_perm_NF < 0.05 & p_perm_DCM < 0.05 & sign(t_NF) == sign(t_DCM), .N],
                   n_both_p05_opp = w[p_perm_NF < 0.05 & p_perm_DCM < 0.05 & sign(t_NF) != sign(t_DCM), .N],
                   n_nf_p05 = w[p_perm_NF < 0.05, .N], n_dcm_p05 = w[p_perm_DCM < 0.05, .N], n_sets = nrow(w))
fwrite(conc, file.path(out, "stratum_concordance.tsv"), sep = "\t")
print(conc)
print(w[order(pmax(p_perm_NF, p_perm_DCM))][1:30, .(pathway = substr(pathway, 1, 60), t_NF, p_perm_NF, t_DCM, p_perm_DCM, p_fwer_NF, p_fwer_DCM)])

write_provenance(stage, inputs = c(cfg$magnet$counts, cfg$magnet$pheno, file.path(cfg$paths$interim, "genesets.tsv")),
                 outputs = list.files(out, "\\.tsv$", full.names = TRUE), seed = cfg$seed,
                 extra = list(n_perm = n_perm, n_sets = length(sets)))
