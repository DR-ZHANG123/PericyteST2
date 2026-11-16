"""Stage 09：反卷积方法在留出供体上的准确度（深度学习 TAPE、Scaden vs 经典 NNLS）。

参考细胞按供体分为训练/测试；训练伪 bulk 只用训练供体细胞，测试伪 bulk 只用测试供体细胞，真实比例已知。
测试两种组成：随机 Dirichlet，以及“视网膜样”组成（神经/胶质比例取测试供体各部位的实际组成，
血管三类比例在 0.2%–3% 均匀抽取）。深度学习方法比较两种训练组成：uniform（随机 Dirichlet，软件默认）与
mixed（一半随机 + 一半训练供体的视网膜样组成），各重复 3 个种子。
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
import anndata
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from deconv import (evaluate, make_pseudobulk, marker_signature, nnls_fractions, retina_props, split_donors,  # noqa: E402
                    to_tape_anndata, training_props)
from repro import ROOT, load_config, set_global_seed, write_provenance  # noqa: E402
from retina_ref import hrca_class  # noqa: E402

STAGE = "09_deconv_benchmark"


def main() -> None:
    from TAPE import Deconvolution
    from TAPE.deconvolution import ScadenDeconvolution
    cfg = load_config(); dc = cfg["deconvolution"]; set_global_seed(cfg["seed"])
    rng = np.random.default_rng(cfg["seed"])
    out = ROOT / "results" / STAGE; out.mkdir(parents=True, exist_ok=True)
    ref = anndata.read_h5ad(ROOT / cfg["reference"]["out_h5ad"])
    X = ref.X.toarray().astype(np.float32) if hasattr(ref.X, "toarray") else np.asarray(ref.X, dtype=np.float32)
    labels = ref.obs["cell_class"].astype(str).values
    classes = sorted(set(labels)); genes = np.asarray(ref.var_names)
    train_d, test_d = split_donors(ref.obs.assign(cell_class=labels), dc["test_donor_fraction"], cfg["seed"])
    tr = ref.obs.donor_id.isin(train_d).values; te = ~tr
    pd.DataFrame({"donor_id": sorted(train_d | test_d), "set": ["test" if d in test_d else "train" for d in sorted(train_d | test_d)]}
                 ).to_csv(out / "donor_split.tsv", sep="\t", index=False)

    comp = pd.read_csv(ROOT / cfg["reference"]["composition"], sep="\t", index_col=[0, 1])
    cdon = comp.index.get_level_values("donor_id")
    comp_train, comp_test = comp[cdon.isin(train_d)], comp[cdon.isin(test_d)]
    tests = {}
    for name, props in [("random", None), ("retina_like", retina_props(comp_test, classes, dc["n_pseudobulk_test"], rng))]:
        bx, bp = make_pseudobulk(X[te], labels[te], classes, dc["n_pseudobulk_test"], dc["cells_per_pseudobulk"], rng, props=props)
        idx = [f"{name}_{i}" for i in range(len(bx))]
        tests[name] = (pd.DataFrame(bx, index=idx, columns=genes), pd.DataFrame(bp, index=idx, columns=classes))
        tests[name][1].assign(test=name).to_csv(out / f"truth_{name}.tsv", sep="\t")

    sig = marker_signature(X[tr], labels[tr], classes, genes)
    rows, preds = [], []
    for tname, (bulk, truth) in tests.items():
        p = nnls_fractions(bulk, sig)
        rows.append(evaluate(p, truth).assign(method="NNLS", training="none", test=tname, seed=np.nan))
        preds.append(p[classes].assign(method="NNLS", training="none", test=tname))
    import torch
    for scheme in dc["training_schemes"]:
        props = training_props(scheme, comp_train, classes, dc["n_pseudobulk_train"], rng)
        sim_X, sim_p = make_pseudobulk(X[tr], labels[tr], classes, dc["n_pseudobulk_train"], dc["cells_per_pseudobulk"], rng,
                                       props=props, sparse_prob=0.5)
        sim = to_tape_anndata(sim_X, sim_p, genes, classes)
        for tname, (bulk, truth) in tests.items():
            for seed in dc["seeds"]:
                _, pt = Deconvolution(sim, bulk, datatype="counts", mode="overall", adaptive=False,
                                      batch_size=dc["tape_batch_size"], epochs=dc["tape_epochs"], seed=seed)
                pt.index = bulk.index
                np.random.seed(seed); torch.manual_seed(seed)
                ps = ScadenDeconvolution(sim, bulk, batch_size=dc["tape_batch_size"], epochs=dc["tape_epochs"])
                ps.index = bulk.index
                for m, pr in [("TAPE", pt), ("Scaden", ps)]:
                    rows.append(evaluate(pr[classes], truth).assign(method=m, training=scheme, test=tname, seed=seed))
                    if seed == dc["seeds"][0]:
                        preds.append(pr[classes].assign(method=m, training=scheme, test=tname))
    res = pd.concat(rows); res.to_csv(out / "benchmark_metrics.tsv", sep="\t", index=False)
    pd.concat(preds).to_csv(out / "benchmark_predictions_seed1.tsv", sep="\t")
    summ = res.groupby(["test", "method", "training"])[["pearson", "ccc", "rmse"]].mean().round(3)
    summ.to_csv(out / "benchmark_summary.tsv", sep="\t"); print(summ)
    print(res[res.test == "retina_like"].groupby(["cell_class", "method", "training"]).ccc.mean().unstack([1, 2]).round(3))
    write_provenance(STAGE, [ROOT / cfg["reference"]["out_h5ad"]], sorted(out.glob("*.tsv")), seed=cfg["seed"],
                     extra={"n_train_donors": len(train_d), "n_test_donors": len(test_d)})


if __name__ == "__main__":
    main()
