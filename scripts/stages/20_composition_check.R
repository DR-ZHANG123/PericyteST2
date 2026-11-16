# Stage 20：周细胞/平滑肌组成与“状态”比值的区分
#   1) 周细胞指数不含 IL1RL1/ADAMTS4（状态基因不参与定义“细胞减少”）；周细胞 − 平滑肌
#   2) IL1RL1 相对周细胞专属标记（不含 IL1RL1）；IL1RL1 相对壁细胞含量 + 周细胞−平滑肌组成指数
#   3) 胶原：相对平滑肌标记、相对周细胞标记；与 IL1RL1 的供体关联加入组成指数后
#   4) NPDR 与无 DR 之间壁细胞含量的差异（早晚两阶段效应是否可区分）
suppressPackageStartupMessages({library(limma); library(edgeR); library(data.table); library(parallel)})
source(file.path("scripts", "lib", "repro.R"))
cfg <- load_config(); set_global_seed(cfg$seed)
root <- p22_root(); setwd(root)
stage <- "20_composition_check"; cc <- cfg$composition_check; n_perm <- cc$n_perm
out <- file.path("results", stage); dir.create(out, recursive = TRUE, showWarnings = FALSE)

meta_all <- fread(file.path("results", "06_retina_pilot", "sample_meta.tsv"))
cnt <- as.matrix(read.csv(gzfile(cfg$retina$counts), row.names = 1, check.names = FALSE))[, meta_all$sample]
gm <- fread(cfg$deconvolution$gene_map); setnames(gm, 1:2, c("id", "gene")); gm[, id := sub("\\..*$", "", id)]; gm <- gm[!duplicated(id)]
s2 <- gm$gene[match(rownames(cnt), gm$id)]; cnt <- cnt[!is.na(s2), ]; cnt <- rowsum(cnt, s2[!is.na(s2)])
d <- DGEList(cnt); d <- calcNormFactors(d[filterByExpr(d, group = meta_all$group), , keep.lib.sizes = FALSE])
lc <- cpm(d, log = TRUE, prior.count = 2)
d_all <- DGEList(cnt); d_all$samples$norm.factors <- d$samples$norm.factors
la <- cpm(d_all, log = TRUE, prior.count = 2)
share <- fread(file.path("results", "09c_ratio_validation", "bulk_share_by_celltype.tsv")); setnames(share, 1, "gene")
mur <- share[mural >= cfg$mural_state$share_min & gene %in% rownames(lc)]$gene
content <- apply(lc[mur, ], 2, median)
ms <- rbind(fread(file.path("results", "14a_marker_sets", "marker_sets.tsv"))[, .(set, gene)],
            fread(file.path("results", "16a_mural_markers", "mural_marker_sets.tsv"))[, .(set, gene)])
ex <- cc$exclude_state_genes
idx <- function(s, drop = ex) colMeans(la[intersect(setdiff(ms$gene[ms$set == s], drop), rownames(la)), , drop = FALSE])
peri <- sapply(c("pericyte_only_top8", "pericyte_strict", "pericyte_relaxed"), idx)
smc <- sapply(c("smc_only_top8", "smc_strict", "smc_relaxed"), idx)
comp <- peri[, "pericyte_only_top8"] - smc[, "smc_only_top8"]          # 周细胞−平滑肌组成指数（不含状态基因）
col <- colMeans(la[cc$collagen, ])

covs <- "region + age + sex + rin + pmi"
m3 <- meta_all[group %in% c("Control", "Diabetic", "NPDR")][, group := factor(group, levels = c("Control", "Diabetic", "NPDR"))]
m3[, comp := comp[sample]][, content := content[sample]]
fitm <- function(Y, dd, rhs, cor = NULL) {
  X <- model.matrix(as.formula(paste("~", rhs)), data = dd)
  if (is.null(cor)) { cor <- duplicateCorrelation(Y, X, block = dd$donor)$consensus.correlation; if (!is.finite(cor)) cor <- 0 }
  list(f = eBayes(lmFit(Y, X, block = dd$donor, correlation = cor)), cor = cor)
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
  res[, model := label][]
}
S <- m3$sample
Y1 <- rbind(t(peri[S, ]), t(smc[S, ]),
            pericyte_minus_smc_top8 = comp[S], pericyte_minus_smc_strict = peri[S, "pericyte_strict"] - smc[S, "smc_strict"],
            pericyte_minus_smc_relaxed = peri[S, "pericyte_relaxed"] - smc[S, "smc_relaxed"],
            IL1RL1_rel_pericyte_top8 = la["IL1RL1", S] - peri[S, "pericyte_only_top8"],
            IL1RL1_rel_pericyte_strict = la["IL1RL1", S] - peri[S, "pericyte_strict"],
            collagen_rel_smc = col[S] - smc[S, "smc_only_top8"], collagen_rel_pericyte = col[S] - peri[S, "pericyte_only_top8"],
            collagen_rel_content = col[S] - content[S])
rownames(Y1)[1:6] <- paste0(rownames(Y1)[1:6], "_without_IL1RL1_ADAMTS4")
r1 <- perm_test(Y1, paste("group +", covs), "group + covariates")
Y2 <- rbind(IL1RL1_rel_content = la["IL1RL1", S] - content[S], collagen_rel_content = col[S] - content[S])
r2 <- perm_test(Y2, paste("group + comp +", covs), "with pericyte-minus-smooth-muscle composition index")
fwrite(rbind(r1, r2), file.path(out, "composition_tests.tsv"), sep = "\t")
bys <- data.table(sample = colnames(la), pericyte_top8 = peri[, "pericyte_only_top8"], smc_top8 = smc[, "smc_only_top8"],
                  il1rl1_rel_pericyte = la["IL1RL1", ] - peri[, "pericyte_only_top8"], comp = comp)
fwrite(merge(meta_all[, .(sample, donor, group, region, severity)], bys, by = "sample"), file.path(out, "measures_by_sample.tsv"), sep = "\t")

# NPDR 与无 DR 之间：壁细胞含量、IL1RL1 相对含量
m3b <- copy(m3)[, group := relevel(group, ref = "Diabetic")]
f3 <- fitm(rbind(content = content[S], IL1RL1_rel_content = la["IL1RL1", S] - content[S]), m3b, paste("group +", covs))$f
stage_diff <- grab(f3, "groupNPDR")[, comparison := "NPDR vs diabetes without DR"]
fwrite(stage_diff, file.path(out, "npdr_vs_diabetic.tsv"), sep = "\t")

md <- meta_all[group != "Control" & !is.na(severity)]
sv <- function(dd) { Y <- rbind(v = (la["IL1RL1", ] - peri[, "pericyte_only_top8"])[dd$sample], v2 = (la["IL1RL1", ] - content)[dd$sample])
  rownames(Y) <- c("IL1RL1_rel_pericyte", "IL1RL1_rel_content"); grab(fitm(Y, dd, paste("severity +", covs))$f, "severity") }
sev <- rbind(sv(md)[, predictor := "ETDRS level"], sv(md[group != "PDR"])[, predictor := "ETDRS level, PDR excluded"])
fwrite(sev, file.path(out, "il1rl1_severity.tsv"), sep = "\t")

# 供体层面：IL1RL1 与胶原的关联，加入组成指数
dd <- meta_all[, .(il = mean((la["IL1RL1", ] - content)[sample]), colc = mean((col - content)[sample]), comp = mean(comp[sample]),
                   group = group[1], age = age[1], sex = sex[1], rin = mean(rin), pmi = mean(pmi)), by = donor]
assoc <- rbindlist(lapply(list(c("group + age + sex + rin + pmi", "group, age, sex, RIN, PMI"),
                               c("group + comp + age + sex + rin + pmi", "group, composition index, age, sex, RIN, PMI")), function(m) {
  ri <- resid(lm(as.formula(paste("il ~", m[1])), dd)); rc <- resid(lm(as.formula(paste("colc ~", m[1])), dd))
  ct <- cor.test(ri, rc, method = "spearman", exact = FALSE)
  data.table(adjusted_for = m[2], rho = ct$estimate, p = ct$p.value, n_donors = nrow(dd))
}))
gene_assoc <- rbindlist(lapply(cc$collagen, function(g) {
  x <- meta_all[, .(cg = mean((la[g, ] - content)[sample])), by = donor]; dd2 <- merge(dd, x, by = "donor")
  ri <- resid(lm(il ~ group + comp + age + sex + rin + pmi, dd2)); rc <- resid(lm(cg ~ group + comp + age + sex + rin + pmi, dd2))
  ct <- cor.test(ri, rc, method = "spearman", exact = FALSE)
  data.table(adjusted_for = paste0(g, " only; group, composition index, age, sex, RIN, PMI"), rho = ct$estimate, p = ct$p.value, n_donors = nrow(dd2))
}))
dd[, il_resid := resid(lm(il ~ group + comp + age + sex + rin + pmi, dd))][, col_resid := resid(lm(colc ~ group + comp + age + sex + rin + pmi, dd))]
fwrite(dd, file.path(out, "il1rl1_collagen_residuals.tsv"), sep = "\t")
cs <- cor.test(dd$il, dd$comp, method = "spearman", exact = FALSE)
assoc <- rbind(assoc, gene_assoc, data.table(adjusted_for = "IL1RL1 vs composition index, unadjusted", rho = cs$estimate, p = cs$p.value, n_donors = nrow(dd)))
fwrite(assoc, file.path(out, "il1rl1_collagen_composition.tsv"), sep = "\t")

print(rbind(r1, r2)[, .(measure, model = substr(model, 1, 20), coef, diff = round(diff, 2), p_perm = signif(p_perm, 2))], nrows = 100)
print(stage_diff); print(assoc)
write_provenance(stage, inputs = c(cfg$retina$counts), outputs = list.files(out, "\\.tsv$", full.names = TRUE), seed = cfg$seed,
                 extra = list(n_perm = n_perm))
