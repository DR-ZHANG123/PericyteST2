# Stage 15b：内皮 PDGF-B 上调是否真的来自内皮
#   单核参考中 PDGFB 的 bulk 信号 38% 来自内皮、33% 来自无长突细胞、12% 来自小胶质（09c bulk_share）。
#   因此在“PDGFB 相对内皮标记”的模型中再校正无长突细胞、小胶质（以及周细胞）标记指数；
#   并检验无长突细胞、小胶质指数本身是否随组别改变。供体层面置换（两部位样本一起移动）。
suppressPackageStartupMessages({library(limma); library(edgeR); library(data.table); library(parallel)})
source(file.path("scripts", "lib", "repro.R"))
cfg <- load_config(); set_global_seed(cfg$seed)
root <- p22_root(); setwd(root)
stage <- "15b_pdgfb_source"; cc <- cfg$crosstalk
out <- file.path("results", stage); dir.create(out, recursive = TRUE, showWarnings = FALSE)

meta_all <- fread(file.path("results", "06_retina_pilot", "sample_meta.tsv"))
cnt <- as.matrix(read.csv(gzfile(cfg$retina$counts), row.names = 1, check.names = FALSE))[, meta_all$sample]
gm <- fread(cfg$deconvolution$gene_map); setnames(gm, 1:2, c("id", "gene")); gm[, id := sub("\\..*$", "", id)]; gm <- gm[!duplicated(id)]
sym <- gm$gene[match(rownames(cnt), gm$id)]; cnt <- cnt[!is.na(sym), ]; cnt <- rowsum(cnt, sym[!is.na(sym)])
lc <- cpm(calcNormFactors(DGEList(cnt)), log = TRUE, prior.count = 2)
msd <- fread(file.path("results", "14a_marker_sets", "marker_sets.tsv")); ms <- split(msd$gene, msd$set)
idx <- function(s) colMeans(lc[intersect(ms[[s]], rownames(lc)), , drop = FALSE])

m3 <- meta_all[group %in% c("Control", "Diabetic", "NPDR")][, group := factor(group, levels = c("Control", "Diabetic", "NPDR"))]
s <- m3$sample
m3[, `:=`(endo = idx("endothelial_top8")[s], amacrine = idx("amacrine_top8")[s], microglia = idx("microglia_top8")[s],
          pericyte = idx("pericyte_only_top8")[s], pdgfb = lc["PDGFB", s])]
m3[, pdgfb_ratio := pdgfb - endo]
covs <- "region + age + sex + rin + pmi"

fit_one <- function(y, d, rhs, cor = NULL) {
  X <- model.matrix(as.formula(paste("~", rhs)), data = d); Y <- rbind(y, y)
  if (is.null(cor)) { cor <- duplicateCorrelation(Y, X, block = d$donor)$consensus.correlation; if (!is.finite(cor)) cor <- 0 }
  f <- eBayes(lmFit(Y, X, block = d$donor, correlation = cor))
  list(tab = rbindlist(lapply(c("groupDiabetic", "groupNPDR"), function(cf) {
    ci <- f$coef[1, cf] + c(-1, 1) * qt(0.975, f$df.total[1]) * f$stdev.unscaled[1, cf] * sqrt(f$s2.post[1])
    data.table(coef = cf, diff = f$coef[1, cf], ci_low = ci[1], ci_high = ci[2], t = f$t[1, cf], p = f$p.value[1, cf])
  })), cor = cor)
}
donors <- unique(m3[, .(donor, group)])
perm_p <- function(y, rhs, obs, cor) {
  tt <- simplify2array(mclapply(seq_len(cc$n_perm), function(b) {
    set.seed(cfg$seed * 100000 + b); sh <- copy(donors)[, group := sample(group)]
    dp <- copy(m3)[, group := factor(sh$group[match(donor, sh$donor)], levels = levels(m3$group))]
    fit_one(y, dp, rhs, cor)$tab$t
  }, mc.cores = 20L))
  sapply(1:2, function(k) (1 + sum(abs(tt[k, ]) >= abs(obs$t[k]))) / (1 + cc$n_perm))
}
models <- list(
  list(y = "pdgfb_ratio", rhs = paste("group +", covs), label = "PDGFB relative to endothelial markers"),
  list(y = "pdgfb_ratio", rhs = paste("group +", covs, "+ amacrine + microglia"), label = "relative to endothelial markers, adjusted for amacrine and microglial indices"),
  list(y = "pdgfb_ratio", rhs = paste("group +", covs, "+ amacrine + microglia + pericyte"), label = "relative to endothelial markers, adjusted for amacrine, microglial and pericyte indices"),
  list(y = "pdgfb", rhs = paste("group +", covs, "+ endo + amacrine + microglia"), label = "bulk PDGFB adjusted for endothelial, amacrine and microglial indices"),
  list(y = "amacrine", rhs = paste("group +", covs), label = "amacrine index"),
  list(y = "microglia", rhs = paste("group +", covs), label = "microglial index"),
  list(y = "endo", rhs = paste("group +", covs), label = "endothelial index"))
res <- rbindlist(lapply(models, function(m) {
  r <- fit_one(m3[[m$y]], m3, m$rhs)
  r$tab[, `:=`(measure = m$label, p_perm = perm_p(m3[[m$y]], m$rhs, r$tab, r$cor))]
}))
fwrite(res, file.path(out, "pdgfb_source_tests.tsv"), sep = "\t")

# 供体层面：每位供体两部位均值（作图用），与周细胞指数的相关（偏相关，校正组别与协变量）
dd <- m3[, .(group = group[1], pdgfb_ratio = mean(pdgfb_ratio), pericyte = mean(pericyte), endo = mean(endo),
             age = age[1], sex = sex[1], rin = mean(rin), pmi = mean(pmi)), by = donor]
fwrite(dd, file.path(out, "pdgfb_by_donor.tsv"), sep = "\t")
print(res[, .(measure, coef, diff = round(diff, 3), p = signif(p, 2), p_perm = signif(p_perm, 2))])
write_provenance(stage, inputs = c(cfg$retina$counts, file.path("results", "14a_marker_sets", "marker_sets.tsv")),
                 outputs = list.files(out, "\\.tsv$", full.names = TRUE), seed = cfg$seed, extra = list(n_perm = cc$n_perm))
