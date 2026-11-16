"""Stage 10：用深度学习反卷积（TAPE）估计每份人视网膜样本的细胞比例与各细胞类型自身的表达。

训练伪 bulk 用全部参考细胞，训练组成按 config 的 apply_training（mixed：一半随机 + 一半视网膜样；留出供体上的比较见 stage 09/09b）；对 GSE160306 每份样本做 TAPE 的逐样本自适应，
得到 细胞比例 与 细胞类型×基因 的样本特异表达。重复 5 个种子取平均，并报告种子间一致性。
经典 NNLS 比例一并输出作对照。细胞类型特异表达矩阵较大，存 data/interim（不入库）。
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

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from deconv import make_pseudobulk, marker_signature, nnls_fractions, to_tape_anndata, training_props  # noqa: E402
from repro import ROOT, load_config, set_global_seed, write_provenance  # noqa: E402

STAGE = "10_deconv_retina"


def load_bulk(cfg) -> pd.DataFrame:
    """GSE160306 计数 → 样本×基因名（同名基因计数相加）。"""
    cnt = pd.read_csv(ROOT / cfg["retina"]["counts"], index_col=0)
    gm = pd.read_csv(ROOT / cfg["deconvolution"]["gene_map"], sep="\t", index_col=0)["gene"]
    gm.index = gm.index.str.split(".").str[0]; gm = gm[~gm.index.duplicated()]
    cnt.index = gm.reindex(cnt.index).values
    cnt = cnt[pd.notna(cnt.index)].groupby(level=0).sum()
    return cnt.T


def main() -> None:
    from TAPE import Deconvolution
    cfg = load_config(); dc = cfg["deconvolution"]; set_global_seed(cfg["seed"])
    rng = np.random.default_rng(cfg["seed"])
    out = ROOT / "results" / STAGE; out.mkdir(parents=True, exist_ok=True)
    ctdir = ROOT / dc["ctexpr_dir"]; ctdir.mkdir(parents=True, exist_ok=True)
    ref = anndata.read_h5ad(ROOT / cfg["reference"]["out_h5ad"])
    X = ref.X.toarray().astype(np.float32) if hasattr(ref.X, "toarray") else np.asarray(ref.X, dtype=np.float32)
    labels = ref.obs["cell_class"].astype(str).values
    classes = sorted(set(labels)); genes = np.asarray(ref.var_names)
    bulk = load_bulk(cfg)
    shared = np.array([g for g in genes if g in set(bulk.columns)])
    gi = pd.Index(genes).get_indexer(shared)
    X, genes, bulk = X[:, gi], shared, bulk[shared]

    comp = pd.read_csv(ROOT / cfg["reference"]["composition"], sep="\t", index_col=[0, 1])
    props = training_props(dc["apply_training"], comp, classes, dc["n_pseudobulk_train"], rng)
    sim_X, sim_p = make_pseudobulk(X, labels, classes, dc["n_pseudobulk_train"], dc["cells_per_pseudobulk"], rng,
                                   props=props, sparse_prob=0.5)
    sim = to_tape_anndata(sim_X, sim_p, genes, classes)
    fr, sig_sum = [], None
    cwd = os.getcwd()
    for seed in dc["apply_seeds"]:
        with tempfile.TemporaryDirectory() as tmp:
            os.chdir(tmp)                       # TAPE 把模型写成当前目录下的 model.pth，每个种子各用一个目录
            sigm, pred = Deconvolution(sim, bulk, datatype="counts", mode="high-resolution", adaptive=True,
                                       batch_size=dc["tape_batch_size"], epochs=dc["tape_epochs"], seed=seed)
            os.chdir(cwd)
        fr.append(pred[classes].assign(seed=seed))
        if sig_sum is None:
            sig_sum = {c: sigm[c].copy() for c in classes}
        else:
            for c in classes:
                sig_sum[c] += sigm[c]
    n = len(dc["apply_seeds"])
    for c in classes:
        (sig_sum[c] / n).to_csv(ctdir / f"{c}.tsv.gz", sep="\t", compression="gzip")
    fr = pd.concat(fr); fr.index.name = "sample"
    fr.to_csv(out / "fractions_by_seed.tsv", sep="\t")
    mean = fr.groupby(level=0)[classes].mean(); mean.to_csv(out / "fractions_tape.tsv", sep="\t")
    # 种子间一致性：每个细胞类型，各种子估计与其余种子平均的相关
    cons = []
    for c in classes:
        m = fr.reset_index().pivot(index="sample", columns="seed", values=c)
        cons.append(dict(cell_class=c, mean_pairwise_r=np.nanmean(np.corrcoef(m.T.values)[np.triu_indices(n, 1)]),
                         mean_fraction=m.values.mean()))
    pd.DataFrame(cons).to_csv(out / "seed_consistency.tsv", sep="\t", index=False)
    sig = marker_signature(X, labels, classes, genes)
    nnls_fractions(bulk, sig).to_csv(out / "fractions_nnls.tsv", sep="\t")
    print(mean.describe().T.round(4)); print(pd.DataFrame(cons).round(3))
    write_provenance(STAGE, [ROOT / cfg["reference"]["out_h5ad"], ROOT / cfg["retina"]["counts"]],
                     sorted(out.glob("*.tsv")), seed=cfg["seed"],
                     extra={"seeds": dc["apply_seeds"], "n_genes": int(len(genes)), "ctexpr_dir": dc["ctexpr_dir"]})


if __name__ == "__main__":
    main()
