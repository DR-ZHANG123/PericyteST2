"""Stage 22：单细胞参考（48,863 个细胞核，13 类）的 UMAP，以及关键基因在各细胞类型中的表达（点图）。

仅用于展示参考中的细胞类型与基因归属（IL1RL1、PDGFRB 在周细胞，IL33 在 Müller 胶质，PDGFB 在内皮/无长突等），
不参与任何检验。每个细胞核按总计数归一化到 1e4 后 log1p；高变基因 → PCA（30 维）→ UMAP（umap-learn，固定随机种子）。
输出：umap_cells.tsv（坐标 + 类别 + 所选基因表达）、gene_by_class.tsv（各类平均表达与表达细胞比例）。
"""
from __future__ import annotations

import sys
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import umap

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from repro import ROOT, load_config, set_global_seed, write_provenance  # noqa: E402

STAGE = "22_reference_umap"


def main() -> None:
    cfg = load_config(); seed = cfg["seed"]; set_global_seed(seed); pc = cfg["reference_umap"]
    out = ROOT / "results" / STAGE; out.mkdir(parents=True, exist_ok=True)
    src = ROOT / "data" / "interim" / "retina_reference.h5ad"
    a = anndata.read_h5ad(src)
    a.var_names = a.var.feature_name.astype(str).values; a.var_names_make_unique()
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    genes = [g for g in pc["genes"] if g in a.var_names]
    expr = pd.DataFrame(a[:, genes].X.toarray(), columns=genes, index=a.obs_names)
    rows = []
    for cl, idx in a.obs.groupby("cell_class", observed=True).indices.items():
        sub = expr.iloc[idx]
        for g in genes:
            rows.append(dict(cell_class=cl, gene=g, mean_log_expr=float(sub[g].mean()), frac_expressing=float((sub[g] > 0).mean()),
                             n_nuclei=len(idx)))
    pd.DataFrame(rows).to_csv(out / "gene_by_class.tsv", sep="\t", index=False)
    sc.pp.highly_variable_genes(a, n_top_genes=pc["n_hvg"], flavor="seurat")
    b = a[:, a.var.highly_variable].copy(); sc.pp.scale(b, max_value=10)
    sc.tl.pca(b, n_comps=pc["n_pcs"], random_state=seed)
    # 直接用 umap-learn（本环境 scanpy 的 kNN 与 scipy 稀疏数组接口不兼容）
    emb = umap.UMAP(n_neighbors=pc["n_neighbors"], min_dist=0.3, random_state=seed).fit_transform(b.obsm["X_pca"])
    u = pd.DataFrame(emb, columns=["umap1", "umap2"], index=b.obs_names)
    u["cell_class"] = b.obs.cell_class.astype(str).values; u["dataset"] = np.where(b.obs.dataset_id == cfg["reference"]["hrca_dataset"], "retina", "optic nerve")
    u = u.join(expr.round(3))
    u.to_csv(out / "umap_cells.tsv", sep="\t", index_label="cell")
    print(u.cell_class.value_counts())
    write_provenance(STAGE, inputs=[src], outputs=[out / "umap_cells.tsv", out / "gene_by_class.tsv"], seed=seed)


if __name__ == "__main__":
    main()
