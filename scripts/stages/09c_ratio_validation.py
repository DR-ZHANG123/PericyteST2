"""Stage 09c：区分“细胞变少”与“细胞自身改变”——同类稳定标记作分母的比值检验（留出供体模拟验证）。

基因 g 若在视网膜 bulk 中的信号主要来自细胞类群 G（份额 ≥ 0.5，G 内平均计数 ≥ 0.2），则 g 与 G 的稳定特异标记之比
抵消 G 的细胞数量变化。G：壁细胞（周细胞 + 平滑肌）、内皮、Müller。标记 = 份额最高、G 内表达 ≥ 0.5 的 8 个基因，
预先排除生物学关注基因（config 的 ratio_exclude_markers）。两种情景，各 40 vs 40 份测试供体视网膜样伪 bulk：
  intrinsic：在周细胞、内皮、Müller 中把各自 30 个“主要来自本类”的基因乘以 0.5（细胞比例不变）
  cell_loss：周细胞与内皮比例减半（表达不变）
比较 bulk 直接检验与比值检验的检出率、AUC 与假阳性率。
"""
from __future__ import annotations

import sys
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from deconv import VASC, make_pseudobulk, retina_props, split_donors  # noqa: E402
from repro import ROOT, load_config, set_global_seed, write_provenance  # noqa: E402

STAGE = "09c_ratio_validation"
GROUPS = {"mural": ["pericyte", "smooth_muscle"], "endothelial": ["endothelial"], "muller": ["muller"]}
SPIKE_IN = {"mural": "pericyte", "endothelial": "endothelial", "muller": "muller"}
N, N_SPIKE, N_MARKERS = 40, 30, 8


def share_tables(X, labels, classes, genes, comp):
    f = comp.mean().reindex(classes).fillna(0).values
    for v in VASC:
        f[classes.index(v)] = 0.01
    f /= f.sum()
    raw = pd.DataFrame(np.vstack([X[labels == c].mean(axis=0) for c in classes]).T, index=genes, columns=classes)
    w = raw * f; share = w.div(w.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    for g, m in GROUPS.items():
        share[g] = share[m].sum(axis=1); raw[g] = raw[m].max(axis=1)
    return share, raw


def define_sets(share, raw, exclude):
    sets = {}
    for g in GROUPS:
        dom = share.index[(share[g] >= 0.5) & (raw[g] >= 0.2)]
        cand = [x for x in share.loc[dom, g].sort_values(ascending=False).index if x not in exclude and raw.loc[x, g] >= 0.5]
        markers = cand[:N_MARKERS]
        sets[g] = (list(markers), [x for x in dom if x not in markers])
    return sets


def welch(a, b):
    t, p = stats.ttest_ind(b, a, axis=0, equal_var=False)
    t = np.nan_to_num(t); p = np.nan_to_num(p, nan=1.0)
    return t, p, multipletests(p, method="fdr_bh")[1]


def evaluate(lb, ga, sets, spiked):
    rows = []
    for g, (markers, tests) in sets.items():
        denom = lb[markers].mean(axis=1)
        for method, M in [("bulk", lb[tests]), ("ratio", lb[tests].sub(denom, axis=0))]:
            t, p, q = welch(M.values[ga], M.values[~ga])
            s = np.isin(tests, spiked.get(g, []))
            rows.append(dict(group=g, method=method, n_tested=len(tests), n_spiked=int(s.sum()),
                             auc=roc_auc_score(s, np.abs(t)) if 0 < s.sum() < len(s) else np.nan,
                             detected_fdr05=float((q[s] < 0.05).mean()) if s.any() else np.nan,
                             right_direction=float((t[s] < 0).mean()) if s.any() else np.nan,
                             false_pos_rate=float((q[~s] < 0.05).mean())))
    return rows


def main() -> None:
    cfg = load_config(); dc = cfg["deconvolution"]; set_global_seed(cfg["seed"])
    rng = np.random.default_rng(cfg["seed"] + 2)
    out = ROOT / "results" / STAGE; out.mkdir(parents=True, exist_ok=True)
    ref = anndata.read_h5ad(ROOT / cfg["reference"]["out_h5ad"])
    X = ref.X.toarray().astype(np.float32) if hasattr(ref.X, "toarray") else np.asarray(ref.X, dtype=np.float32)
    labels = ref.obs["cell_class"].astype(str).values; classes = sorted(set(labels)); genes = np.asarray(ref.var_names)
    train_d, test_d = split_donors(ref.obs.assign(cell_class=labels), dc["test_donor_fraction"], cfg["seed"])
    tr = ref.obs.donor_id.isin(train_d).values; te = ~tr
    comp = pd.read_csv(ROOT / cfg["reference"]["composition"], sep="\t", index_col=[0, 1])
    comp_test = comp[comp.index.get_level_values(0).isin(test_d)]
    share, raw = share_tables(X[tr], labels[tr], classes, genes, comp)
    sets = define_sets(share, raw, set(dc["ratio_exclude_markers"]))
    pd.DataFrame([(g, "marker", x) for g, (m, t) in sets.items() for x in m] + [(g, "tested", x) for g, (m, t) in sets.items() for x in t],
                 columns=["group", "role", "gene"]).to_csv(out / "gene_sets.tsv", sep="\t", index=False)
    share[list(GROUPS)].round(4).to_csv(out / "bulk_share.tsv", sep="\t")

    lt = labels[te]; gidx = pd.Index(genes); rows = []
    # 情景 1：细胞自身改变
    spiked = {g: list(rng.choice(sets[g][1], min(N_SPIKE, len(sets[g][1]) // 2), replace=False)) for g in GROUPS}
    Xb = X[te].copy()
    for g, gs in spiked.items():
        Xb[np.ix_(lt == SPIKE_IN[g], gidx.get_indexer(gs))] *= 0.5
    pa = retina_props(comp_test, classes, N, rng); pb = retina_props(comp_test, classes, N, rng)
    ba, _ = make_pseudobulk(X[te], lt, classes, N, dc["cells_per_pseudobulk"], rng, props=pa)
    bb, _ = make_pseudobulk(Xb, lt, classes, N, dc["cells_per_pseudobulk"], rng, props=pb)
    lb = np.log2(pd.DataFrame(np.vstack([ba, bb]), columns=genes).pipe(lambda d: d.div(d.sum(axis=1), axis=0) * 1e6) + 1)
    ga = np.r_[np.ones(N, bool), np.zeros(N, bool)]
    rows += [dict(r, scenario="intrinsic") for r in evaluate(lb, ga, sets, spiked)]
    # 情景 2：周细胞与内皮变少，表达不变
    pc = retina_props(comp_test, classes, N, rng)
    for v in ["pericyte", "endothelial"]:
        pc[:, classes.index(v)] *= 0.5
    pc /= pc.sum(axis=1, keepdims=True)
    bc, _ = make_pseudobulk(X[te], lt, classes, N, dc["cells_per_pseudobulk"], rng, props=pc)
    lb2 = np.log2(pd.DataFrame(np.vstack([ba, bc]), columns=genes).pipe(lambda d: d.div(d.sum(axis=1), axis=0) * 1e6) + 1)
    rows += [dict(r, scenario="cell_loss") for r in evaluate(lb2, ga, sets, {})]
    res = pd.DataFrame(rows); res.to_csv(out / "ratio_vs_bulk.tsv", sep="\t", index=False)
    print(res.round(3).to_string(index=False))
    for g, (m, t) in sets.items():
        print(g, "markers:", m, "| n tested:", len(t), "| genes of interest tested:", [x for x in dc["ratio_exclude_markers"] if x in t])
    write_provenance(STAGE, [ROOT / cfg["reference"]["out_h5ad"]], sorted(out.glob("*.tsv")), seed=cfg["seed"],
                     extra={"n_per_group": N, "n_spiked": N_SPIKE, "n_markers": N_MARKERS})


if __name__ == "__main__":
    main()
