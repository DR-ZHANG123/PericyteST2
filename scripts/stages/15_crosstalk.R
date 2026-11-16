# Stage 15：周细胞转录本减少与哪些内皮 / Müller 来源信号同步（为体外细胞实验选择靶点）
#   1) 基因 bulk 信号 ≥ 一半来自内皮或 Müller 时，取该基因相对来源细胞标记的比值（抵消来源细胞数量变化）
#   2) 组间差异：比值 ~ 组 + 部位 + 年龄 + 性别 + RIN + 死亡后间隔，供体内相关；供体层面置换
#   3) 供体层面关联：每位供体两部位取均值，周细胞指数 ~ 比值 + 来源细胞标记指数 + 协变量（全部 41 位、仅糖尿病供体）；置换比值
#      来源标记指数必须进模型：内皮标记与周细胞指数正相关（血管含量），比值的分母会造成假负相关；Müller 反之
suppressPackageStartupMessages({library(limma); library(edgeR); library(data.table); library(parallel)})
source(file.path("scripts", "lib", "repro.R"))
cfg <- load_config(); set_global_seed(cfg$seed)
root <- p22_root(); setwd(root)
stage <- "15_crosstalk"; cc <- cfg$crosstalk
out <- file.path("results", stage); dir.create(out, recursive = TRUE, showWarnings = FALSE)

meta_all <- fread(file.path("results", "06_retina_pilot", "sample_meta.tsv"))
cnt <- as.matrix(read.csv(gzfile(cfg$retina$counts), row.names = 1, check.names = FALSE))[, meta_all$sample]
gm <- fread(cfg$deconvolution$gene_map); setnames(gm, 1:2, c("id", "gene")); gm[, id := sub("\\..*$", "", id)]; gm <- gm[!duplicated(id)]
sym <- gm$gene[match(rownames(cnt), gm$id)]; cnt <- cnt[!is.na(sym), ]; cnt <- rowsum(cnt, sym[!is.na(sym)])
lc <- cpm(calcNormFactors(DGEList(cnt)), log = TRUE, prior.count = 2)

ms <- split(fread(file.path("results", "14a_marker_sets", "marker_sets.tsv"))$gene,
            fread(file.path("results", "14a_marker_sets", "marker_sets.tsv"))$set)
mk <- list(endothelial = ms$endothelial_top8, muller = ms$muller_top8)
share <- fread(file.path("results", "09c_ratio_validation", "bulk_share_by_celltype.tsv")); setnames(share, 1, "gene")
expr_ok <- rownames(lc)[rowMeans(lc) >= cc$min_mean_logcpm]

# ---- 候选：来源份额达标的全部基因 + 文献配体 ----
cand <- rbindlist(lapply(names(mk), function(s) share[get(s) >= cc$source_share_min & gene %in% expr_ok & !gene %in% unlist(mk),
                                                     .(gene, source = s, source_share = get(s))]))
lig <- share[gene %in% intersect(cc$ligands, expr_ok)]
lig <- lig[, .(gene, source = ifelse(endothelial >= muller, "endothelial", "muller"), source_share = pmax(endothelial, muller))]
cand <- unique(rbind(cand[, literature := gene %in% cc$ligands], lig[!gene %in% cand$gene][, literature := TRUE]), by = "gene")
R <- t(sapply(seq_len(nrow(cand)), function(i) lc[cand$gene[i], ] - colMeans(lc[mk[[cand$source[i]]], ])))
rownames(R) <- cand$gene
peri <- colMeans(lc[ms$pericyte_only_top8, ]); mural <- colMeans(lc[ms$mural_top8, ])

# ---- 2) 组间差异 ----
covs <- "region + age + sex + rin + pmi"
m3 <- meta_all[group %in% c("Control", "Diabetic", "NPDR")][, group := factor(group, levels = c("Control", "Diabetic", "NPDR"))]
fitg <- function(Y, d, cor = NULL) {
  X <- model.matrix(as.formula(paste("~ group +", covs)), data = d)
  if (is.null(cor)) { cor <- duplicateCorrelation(Y, X, block = d$donor)$consensus.correlation; if (!is.finite(cor)) cor <- 0 }
  list(f = eBayes(lmFit(Y, X, block = d$donor, correlation = cor)), cor = cor)
}
Y3 <- R[, m3$sample]; main <- fitg(Y3, m3)
donors <- unique(m3[, .(donor, group)])
pt <- simplify2array(mclapply(seq_len(cc$n_perm), function(b) {
  set.seed(cfg$seed * 100000 + b); sh <- copy(donors)[, group := sample(group)]
  dp <- copy(m3)[, group := factor(sh$group[match(donor, sh$donor)], levels = levels(m3$group))]
  fitg(Y3, dp, main$cor)$f$t[, c("groupDiabetic", "groupNPDR")]
}, mc.cores = 20L))
grp <- rbindlist(lapply(c("groupDiabetic", "groupNPDR"), function(cf) {
  tt <- topTable(main$f, coef = cf, number = Inf, sort.by = "none")
  data.table(gene = rownames(tt), coef = cf, diff = tt$logFC, t = tt$t, p = tt$P.Value,
             p_perm = (1 + rowSums(abs(pt[, cf, ]) >= abs(tt$t))) / (1 + cc$n_perm))
}))
grp[, fdr := p.adjust(p, "BH"), by = coef]
grp <- merge(cand, grp, by = "gene")

# ---- 3) 供体层面关联 ----
dmeta <- m3[, .(group = group[1], age = age[1], sex = sex[1], rin = mean(rin), pmi = mean(pmi)), by = donor]
dmean <- function(v) tapply(v[m3$sample], m3$donor, mean)[dmeta$donor]
Rd <- t(apply(R, 1, dmean)); pd <- dmean(peri)
Sd <- sapply(mk, function(g) dmean(colMeans(lc[g, ])))
# 线性模型中 x 的 t 值（Frisch–Waugh：先把 x、y 对协变量取残差）；所有置换一次矩阵计算
assoc_t <- function(Xp, y, Z) {
  M <- diag(length(y)) - Z %*% solve(crossprod(Z), t(Z))
  rx <- M %*% Xp; ry <- drop(M %*% y); sxx <- colSums(rx^2)
  beta <- drop(crossprod(rx, ry)) / sxx; df <- length(y) - ncol(Z) - 1
  rbind(beta = beta, t = beta / sqrt((sum(ry^2) - beta^2 * sxx) / df / sxx))
}
assoc <- function(sel, rhs, label) {
  d <- droplevels(dmeta[sel]); y <- pd[sel]; Z <- model.matrix(as.formula(paste("~", rhs)), data = d)
  rbindlist(lapply(rownames(Rd), function(g) {
    x <- Rd[g, sel]; Zg <- cbind(Z, source_index = Sd[sel, cand$source[cand$gene == g]])
    obs <- assoc_t(matrix(x), y, Zg)
    set.seed(cfg$seed); null <- assoc_t(replicate(cc$n_perm, sample(x)), y, Zg)["t", ]
    data.table(gene = g, set = label, n_donors = sum(sel), beta = obs["beta", 1], t = obs["t", 1],
               r_partial = cor(resid(lm(x ~ Zg - 1)), resid(lm(y ~ Zg - 1))), p_perm = (1 + sum(abs(null) >= abs(obs["t", 1]))) / (1 + cc$n_perm))
  }))
}
as_all <- assoc(rep(TRUE, nrow(dmeta)), "group + age + sex + rin + pmi", "all donors, adjusted for group")
as_all2 <- assoc(rep(TRUE, nrow(dmeta)), "age + sex + rin + pmi", "all donors")
as_dm <- assoc(dmeta$group != "Control", "group + age + sex + rin + pmi", "diabetic donors")
as <- rbind(as_all2, as_all, as_dm)[, fdr := p.adjust(p_perm, "BH"), by = set]
as <- merge(cand, as, by = "gene")

fwrite(grp[order(coef, p)], file.path(out, "source_ratio_group_tests.tsv"), sep = "\t")
fwrite(as[order(set, p_perm)], file.path(out, "pericyte_association.tsv"), sep = "\t")
fwrite(data.table(donor = dmeta$donor, group = dmeta$group, pericyte_index = pd, t(Rd[intersect(cc$ligands, rownames(Rd)), ])),
       file.path(out, "donor_ligand_ratios.tsv"), sep = "\t")

# ---- 4) 复核：换内皮标记数（5/8/12/20）、bulk 本身、黄斑与周边分别 ----
fit1 <- function(y, d, rhs) {
  X <- model.matrix(as.formula(paste("~", rhs)), data = d); Y <- rbind(y, y)
  cr <- if (uniqueN(d$donor) < nrow(d)) duplicateCorrelation(Y, X, block = d$donor)$consensus.correlation else 0
  f <- if (cr != 0) eBayes(lmFit(Y, X, block = d$donor, correlation = cr)) else eBayes(lmFit(Y, X))
  rbindlist(lapply(c("groupDiabetic", "groupNPDR"), function(cf) data.table(coef = cf, diff = f$coef[1, cf], p = f$p.value[1, cf])))
}
rob <- rbindlist(lapply(cc$followup, function(g) {
  sets <- c(paste0("endothelial_top", c(5, 8, 12, 20)), "bulk")
  r1 <- rbindlist(lapply(sets, function(st) {
    y <- if (st == "bulk") lc[g, m3$sample] else lc[g, m3$sample] - colMeans(lc[ms[[st]], m3$sample])
    fit1(y, m3, paste("group +", covs))[, `:=`(measure = st, region = "both")]
  }))
  r2 <- rbindlist(lapply(c("macula", "periphery"), function(rg) { d <- m3[region == rg]
    fit1(lc[g, d$sample] - colMeans(lc[ms$endothelial_top8, d$sample]), d, "group + age + sex + rin + pmi")[, `:=`(measure = "endothelial_top8", region = rg)] }))
  rbind(r1, r2)[, gene := g]
}))
fwrite(rob, file.path(out, "followup_robustness.tsv"), sep = "\t")

lg <- cand[literature == TRUE, gene]
w <- dcast(grp[gene %in% lg], gene + source + source_share ~ coef, value.var = c("diff", "p_perm"))
w <- merge(w, dcast(as[gene %in% lg], gene ~ set, value.var = "p_perm"), by = "gene")
print(w[order(source, gene)], digits = 2)
cat("\n全基因组：组间 FDR<0.1 的来源基因数\n"); print(grp[fdr < 0.1, .N, by = .(source, coef)])
cat("\n与周细胞指数关联最强（全部供体，未校正组）\n"); print(as[set == "all donors"][order(p_perm)][1:20, .(gene, source, source_share, beta, r_partial, p_perm, fdr)], digits = 2)
cat("\n糖尿病供体内\n"); print(as[set == "diabetic donors"][order(p_perm)][1:15, .(gene, source, source_share, beta, r_partial, p_perm, fdr)], digits = 2)
write_provenance(stage, inputs = c(cfg$retina$counts, file.path("results", "14a_marker_sets", "marker_sets.tsv"),
                                   file.path("results", "09c_ratio_validation", "bulk_share_by_celltype.tsv")),
                 outputs = list.files(out, "\\.tsv$", full.names = TRUE), seed = cfg$seed, extra = list(n_perm = cc$n_perm))
