"""Stage 09b：深度学习反卷积能否把“某一细胞类型自身的表达改变”找出来并归到正确的细胞类型（留出供体模拟）。

用测试供体的细胞构造两组视网膜样伪 bulk（各 40 份，组成独立抽取）。第二组中预先植入已知改变：
周细胞 RGS5、NOTCH3 + 48 个随机周细胞表达基因 ×0.5；内皮 DLL4 + 49 个随机内皮表达基因 ×0.5；Müller 50 个随机基因 ×2。
TAPE（mixed 训练）逐样本自适应得到各细胞类型表达，比较两组：
  - 植入基因在其所属细胞类型中的检出（AUC、FDR<0.05 检出率）与方向
  - 被错误归到其他细胞类型的比例；未植入基因的假阳性率
  - 对照：直接在 bulk 上比较（只能检出、不能归属）
"""
from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
import anndata
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from deconv import make_pseudobulk, retina_props, split_donors, to_tape_anndata, training_props  # noqa: E402
from repro import ROOT, load_config, set_global_seed, write_provenance  # noqa: E402

STAGE = "09b_spikein"
SPIKES = {"pericyte": (["RGS5", "NOTCH3"], 48, 0.5), "endothelial": (["DLL4"], 49, 0.5), "muller": ([], 50, 2.0)}
N_PER_GROUP = 40


def pick_spikes(X, labels, genes, rng):
    chosen, used = {}, set()
    for cls, (fixed, n_rand, fac) in SPIKES.items():
        m = X[labels == cls].mean(axis=0)
        pool = [g for g in genes[m >= 0.5] if g not in used and g not in fixed]
        sel = list(fixed) + list(rng.choice(pool, n_rand, replace=False))
        used |= set(sel); chosen[cls] = (sel, fac)
    return chosen


def welch(a: pd.DataFrame, b: pd.DataFrame) -> pd.DataFrame:
    t, p = stats.ttest_ind(b.values, a.values, axis=0, equal_var=False)
    t = np.nan_to_num(t); p = np.nan_to_num(p, nan=1.0)
    return pd.DataFrame({"t": t, "p": p, "fdr": multipletests(p, method="fdr_bh")[1]}, index=a.columns)


def main() -> None:
    from TAPE import Deconvolution
    cfg = load_config(); dc = cfg["deconvolution"]; set_global_seed(cfg["seed"])
    rng = np.random.default_rng(cfg["seed"] + 1)
    out = ROOT / "results" / STAGE; out.mkdir(parents=True, exist_ok=True)
    ref = anndata.read_h5ad(ROOT / cfg["reference"]["out_h5ad"])
    X = ref.X.toarray().astype(np.float32) if hasattr(ref.X, "toarray") else np.asarray(ref.X, dtype=np.float32)
    labels = ref.obs["cell_class"].astype(str).values; classes = sorted(set(labels)); genes = np.asarray(ref.var_names)
    train_d, test_d = split_donors(ref.obs.assign(cell_class=labels), dc["test_donor_fraction"], cfg["seed"])
    tr = ref.obs.donor_id.isin(train_d).values; te = ~tr
    comp = pd.read_csv(ROOT / cfg["reference"]["composition"], sep="\t", index_col=[0, 1])
    cdon = comp.index.get_level_values(0)
    comp_train, comp_test = comp[cdon.isin(train_d)], comp[cdon.isin(test_d)]

    spikes = pick_spikes(X[te], labels[te], genes, rng)
    pd.DataFrame([(c, g, f) for c, (gs, f) in spikes.items() for g in gs], columns=["cell_class", "gene", "factor"]
                 ).to_csv(out / "spiked_genes.tsv", sep="\t", index=False)
    Xb = X[te].copy(); lt = labels[te]; gidx = pd.Index(genes)
    for cls, (gs, fac) in spikes.items():
        Xb[np.ix_(lt == cls, gidx.get_indexer(gs))] *= fac
    pa = retina_props(comp_test, classes, N_PER_GROUP, rng); pb = retina_props(comp_test, classes, N_PER_GROUP, rng)
    ba, _ = make_pseudobulk(X[te], lt, classes, N_PER_GROUP, dc["cells_per_pseudobulk"], rng, props=pa)
    bb, _ = make_pseudobulk(Xb, lt, classes, N_PER_GROUP, dc["cells_per_pseudobulk"], rng, props=pb)
    bulk = pd.DataFrame(np.vstack([ba, bb]), columns=genes, index=[f"A{i}" for i in range(N_PER_GROUP)] + [f"B{i}" for i in range(N_PER_GROUP)])

    props = training_props(dc["apply_training"], comp_train, classes, dc["n_pseudobulk_train"], rng)
    sim_X, sim_p = make_pseudobulk(X[tr], labels[tr], classes, dc["n_pseudobulk_train"], dc["cells_per_pseudobulk"], rng, props=props, sparse_prob=0.5)
    sim = to_tape_anndata(sim_X, sim_p, genes, classes)
    sig_sum = None; cwd = os.getcwd()
    for seed in dc["seeds"]:
        with tempfile.TemporaryDirectory() as tmp:
            os.chdir(tmp)
            sigm, _ = Deconvolution(sim, bulk, datatype="counts", mode="high-resolution", adaptive=True,
                                    batch_size=dc["tape_batch_size"], epochs=dc["tape_epochs"], seed=seed)
            os.chdir(cwd)
        sig_sum = {c: sigm[c].copy() for c in classes} if sig_sum is None else {c: sig_sum[c] + sigm[c] for c in classes}

    ga = bulk.index.str.startswith("A")
    rows, attrib = [], []
    for cls in classes:
        S = sig_sum[cls] / len(dc["seeds"])
        w = welch(S[ga], S[~ga])
        for scls, (gs, fac) in spikes.items():
            g_in = w.index.intersection(gs); g_out = w.index.difference(sum((v[0] for v in spikes.values()), []))
            y = np.r_[np.ones(len(g_in)), np.zeros(len(g_out))]
            score = np.r_[w.loc[g_in, "t"].abs(), w.loc[g_out, "t"].abs()]
            exp_sign = np.sign(np.log(fac))
            attrib.append(dict(estimated_in=cls, spiked_in=scls, n_spiked=len(g_in),
                               auc=roc_auc_score(y, score) if len(g_in) else np.nan,
                               detected_fdr05=float((w.loc[g_in, "fdr"] < 0.05).mean()) if len(g_in) else np.nan,
                               right_direction=float((np.sign(w.loc[g_in, "t"]) == exp_sign).mean()) if len(g_in) else np.nan,
                               false_pos_rate_unspiked=float((w.loc[g_out, "fdr"] < 0.05).mean())))
        for g in ["RGS5", "NOTCH3", "DLL4"]:
            if g in w.index:
                rows.append(dict(estimated_in=cls, gene=g, t=w.loc[g, "t"], fdr=w.loc[g, "fdr"]))
    lb = np.log2(bulk.div(bulk.sum(axis=1), axis=0) * 1e6 + 1)
    wb = welch(lb[ga], lb[~ga])
    for scls, (gs, fac) in spikes.items():
        g_in = wb.index.intersection(gs); g_out = wb.index.difference(sum((v[0] for v in spikes.values()), []))
        attrib.append(dict(estimated_in="bulk", spiked_in=scls, n_spiked=len(g_in),
                           auc=roc_auc_score(np.r_[np.ones(len(g_in)), np.zeros(len(g_out))], np.r_[wb.loc[g_in, "t"].abs(), wb.loc[g_out, "t"].abs()]),
                           detected_fdr05=float((wb.loc[g_in, "fdr"] < 0.05).mean()),
                           right_direction=float((np.sign(wb.loc[g_in, "t"]) == np.sign(np.log(fac))).mean()),
                           false_pos_rate_unspiked=float((wb.loc[g_out, "fdr"] < 0.05).mean())))
    at = pd.DataFrame(attrib); at.to_csv(out / "spikein_attribution.tsv", sep="\t", index=False)
    pd.DataFrame(rows).to_csv(out / "spikein_key_genes.tsv", sep="\t", index=False)
    print(at.round(3).to_string(index=False)); print(pd.DataFrame(rows).round(3).to_string(index=False))
    write_provenance(STAGE, [ROOT / cfg["reference"]["out_h5ad"]], sorted(out.glob("*.tsv")), seed=cfg["seed"],
                     extra={"n_per_group": N_PER_GROUP, "training": dc["apply_training"]})


if __name__ == "__main__":
    main()
