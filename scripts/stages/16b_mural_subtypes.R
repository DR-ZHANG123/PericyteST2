# Stage 16b：周细胞与平滑肌转录本在供体中是否分开改变
#   1) 各标记定义的指数：组间差异 + 供体层面置换
#   2) 周细胞 − 平滑肌 指数差：周细胞减少是否显著多于平滑肌
#   3) 周细胞指数校正平滑肌指数、校正渗漏细胞（节细胞、视杆、内皮、Müller）指数后的组间差异
#   4) 单基因效应；5) 糖尿病供体内与 ETDRS 等级的关系
suppressPackageStartupMessages({library(limma); library(edgeR); library(data.table); library(parallel)})
source(file.path("scripts", "lib", "repro.R"))
cfg <- load_config(); set_global_seed(cfg$seed)
root <- p22_root(); setwd(root)
stage <- "16b_mural_subtypes"; n_perm <- cfg$mural_subtypes$n_perm
out <- file.path("results", stage); dir.create(out, recursive = TRUE, showWarnings = FALSE)

meta_all <- fread(file.path("results", "06_retina_pilot", "sample_meta.tsv"))
cnt <- as.matrix(read.csv(gzfile(cfg$retina$counts), row.names = 1, check.names = FALSE))[, meta_all$sample]
gm <- fread(cfg$deconvolution$gene_map); setnames(gm, 1:2, c("id", "gene")); gm[, id := sub("\\..*$", "", id)]; gm <- gm[!duplicated(id)]
sym <- gm$gene[match(rownames(cnt), gm$id)]; cnt <- cnt[!is.na(sym), ]; cnt <- rowsum(cnt, sym[!is.na(sym)])
lc <- cpm(calcNormFactors(DGEList(cnt)), log = TRUE, prior.count = 2)

mm <- fread(file.path("results", "16a_mural_markers", "mural_marker_sets.tsv"))
base <- fread(file.path("results", "14a_marker_sets", "marker_sets.tsv"))
sets <- c(split(mm$gene, mm$set), split(base$gene, base$set)[c("pericyte_only_top8", "smc_only_top8", "ganglion_top8", "rod_top8", "endothelial_top8", "muller_top8")])
idx <- t(sapply(sets, function(g) colMeans(lc[intersect(g, rownames(lc)), , drop = FALSE])))
pairs <- list(strict = c("pericyte_strict", "smc_strict"), relaxed = c("pericyte_relaxed", "smc_relaxed"),
              classic = c("pericyte_classic", "smc_classic"), main = c("pericyte_only_top8", "smc_only_top8"))
for (p in names(pairs)) idx <- rbind(idx, setNames(list(idx[pairs[[p]][1], ] - idx[pairs[[p]][2], ]), paste0("pericyte_minus_smc_", p))[[1]])
rownames(idx)[(nrow(idx) - length(pairs) + 1):nrow(idx)] <- paste0("pericyte_minus_smc_", names(pairs))
fwrite(data.table(index = rownames(idx), idx), file.path(out, "indices_by_sample.tsv"), sep = "\t")

covs <- "region + age + sex + rin + pmi"
m3 <- meta_all[group %in% c("Control", "Diabetic", "NPDR")][, group := factor(group, levels = c("Control", "Diabetic", "NPDR"))]
fitm <- function(Y, d, rhs, cor = NULL) {
  X <- model.matrix(as.formula(paste("~", rhs)), data = d)
  if (is.null(cor)) { cor <- duplicateCorrelation(Y, X, block = d$donor)$consensus.correlation; if (!is.finite(cor)) cor <- 0 }
  list(f = eBayes(lmFit(Y, X, block = d$donor, correlation = cor)), cor = cor)
}
grab <- function(f, cf) { tt <- topTable(f, coef = cf, number = Inf, sort.by = "none", confint = TRUE)
  data.table(measure = rownames(tt), coef = cf, diff = tt$logFC, ci_low = tt$CI.L, ci_high = tt$CI.R, t = tt$t, p = tt$P.Value) }
donors <- unique(m3[, .(donor, group)])
perm_test <- function(Y, rhs, label) {
  main <- fitm(Y, m3, rhs)
  pt <- simplify2array(mclapply(seq_len(n_perm), function(b) {
    set.seed(cfg$seed * 100000 + b); sh <- copy(donors)[, group := sample(group)]
    dp <- copy(m3)[, group := factor(sh$group[match(donor, sh$donor)], levels = levels(m3$group))]
    fitm(Y, dp, rhs, main$cor)$f$t[, c("groupDiabetic", "groupNPDR"), drop = FALSE]
  }, mc.cores = 16L))
  res <- rbind(grab(main$f, "groupDiabetic"), grab(main$f, "groupNPDR"))
  res[, p_perm := mapply(function(m, cf, tv) (1 + sum(abs(pt[m, cf, ]) >= abs(tv))) / (1 + n_perm), measure, coef, t)]
  res[, model := label]
}

# ---- 1) + 2) ----
r12 <- perm_test(idx[, m3$sample], paste("group +", covs), "group + covariates")
# ---- 3) 校正平滑肌、渗漏细胞 ----
r3 <- rbindlist(lapply(c("pericyte_strict", "pericyte_relaxed", "pericyte_only_top8"), function(pc) {
  smc <- sub("pericyte_only_top8", "smc_only_top8", sub("pericyte", "smc", pc))
  m3[, `:=`(smc_idx = idx[smc, sample], gang = idx["ganglion_top8", sample], rod = idx["rod_top8", sample],
            endo = idx["endothelial_top8", sample], mul = idx["muller_top8", sample])]
  Y <- rbind(idx[pc, m3$sample], idx[pc, m3$sample]); rownames(Y) <- c(pc, paste0(pc, "_dup"))
  rbind(perm_test(Y, paste("group + smc_idx +", covs), "adjusted for smooth-muscle index"),
        perm_test(Y, paste("group + gang + rod + endo + mul +", covs), "adjusted for ganglion, rod, endothelial, Müller indices"))[measure == pc]
}))
fwrite(rbind(r12, r3), file.path(out, "group_tests.tsv"), sep = "\t")

# ---- 4) 单基因 ----
gs <- unique(mm[set %in% c("pericyte_strict", "smc_strict", "pericyte_classic", "smc_classic"), .(set, gene)])
gs <- gs[gene %in% rownames(lc)]; G <- lc[unique(gs$gene), m3$sample]; fg <- fitm(G, m3, paste("group +", covs))$f
per_gene <- merge(gs, rbind(grab(fg, "groupDiabetic"), grab(fg, "groupNPDR")), by.x = "gene", by.y = "measure", allow.cartesian = TRUE)
fwrite(per_gene[order(set, coef, gene)], file.path(out, "per_gene_tests.tsv"), sep = "\t")

# ---- 5) 糖尿病供体内与 ETDRS 等级 ----
md <- meta_all[group != "Control" & !is.na(severity)]
dose <- rbind(grab(fitm(idx[, md$sample], md, paste("severity +", covs))$f, "severity")[, predictor := "ETDRS level"],
              { mn <- md[group != "PDR"]; grab(fitm(idx[, mn$sample], mn, paste("severity +", covs))$f, "severity")[, predictor := "ETDRS level, PDR excluded"] })
fwrite(dose, file.path(out, "dose_response.tsv"), sep = "\t")

key <- c("pericyte_strict", "smc_strict", "pericyte_relaxed", "smc_relaxed", "pericyte_classic", "smc_classic",
         "pericyte_only_top8", "smc_only_top8", grep("minus", rownames(idx), value = TRUE))
print(r12[measure %in% key, .(measure, coef, diff = round(diff, 3), t = round(t, 2), p = signif(p, 2), p_perm = signif(p_perm, 2))][order(coef, measure)])
print(r3[, .(measure, model, coef, diff = round(diff, 3), p = signif(p, 2), p_perm = signif(p_perm, 2))])
print(per_gene[, .(gene, set, coef, diff = round(diff, 2), p = signif(p, 2))][order(set, coef)])
print(dose[measure %in% key, .(measure, predictor, diff = round(diff, 4), p = signif(p, 2))])
write_provenance(stage, inputs = c(cfg$retina$counts, file.path("results", "16a_mural_markers", "mural_marker_sets.tsv"),
                                   file.path("results", "14a_marker_sets", "marker_sets.tsv")),
                 outputs = list.files(out, "\\.tsv$", full.names = TRUE), seed = cfg$seed, extra = list(n_perm = n_perm))
