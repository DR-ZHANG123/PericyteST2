"""Stage 13：小鼠视网膜血管细胞单细胞（GSE204880，野生型 vs db/db 21 周，各 1 个样本）中的内皮–周细胞 Notch 轴。

只作跨物种描述（每组 1 只动物，不做统计检验）。该样本几乎全为内皮细胞（Pdgfrb+ 非内皮细胞约 70 个），
只能读内皮侧：Dll4/Tek/Plvap/Cdh5 等（周细胞行若细胞不足则为空）的
平均表达与表达细胞比例。细胞类别用标记基因聚类注释；周细胞注释不使用 Rgs5/Notch3（它们是被检验的基因）。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from repro import ROOT, load_config, set_global_seed, write_provenance  # noqa: E402

STAGE = "13_mouse_vascular"
ANNOT = {"endothelial": ["Pecam1", "Cdh5", "Cldn5", "Flt1", "Ly6c1"],
         "pericyte": ["Pdgfrb", "Kcnj8", "Abcc9", "Higd1b", "Vtn"],
         "smooth_muscle": ["Acta2", "Myh11", "Tagln", "Cnn1"],
         "microglia": ["P2ry12", "Cx3cr1", "C1qa", "Tmem119"],
         "muller_astro": ["Rlbp1", "Glul", "Gfap", "Aqp4"]}
TESTS = {"endothelial": ["Dll4", "Tek", "Plvap", "Cdh5", "Cldn5", "Esm1", "Apln", "Angpt2", "Kdr", "Hey1"],
         "pericyte": ["Notch3", "Rgs5", "Heyl", "Hey1", "Nrarp", "Pdgfrb", "Kcnj8", "Txnip", "Angpt1"]}


def main() -> None:
    cfg = load_config(); set_global_seed(cfg["seed"])
    d = ROOT / cfg["paths"]["raw"] / "gse204880" / "mouse"
    ads = []
    for tag, cond in [("GSM6199212_WTR_1", "wild_type"), ("GSM6199213_DBR_1", "db_db")]:
        a = sc.read_10x_mtx(d, prefix=f"{tag}_", var_names="gene_symbols", make_unique=True)
        a.obs["condition"] = cond; ads.append(a)
    ad = sc.concat(ads, label=None, index_unique="-")
    ad.var["mt"] = ad.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(ad, qc_vars=["mt"], inplace=True, percent_top=None)
    ad = ad[(ad.obs.n_genes_by_counts >= 300) & (ad.obs.pct_counts_mt < 20)].copy()
    ad.layers["counts"] = ad.X.copy()
    sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad)
    sc.pp.highly_variable_genes(ad, n_top_genes=2000, batch_key="condition")
    sc.pp.pca(ad, n_comps=30, random_state=cfg["seed"])
    from sklearn.cluster import KMeans          # 聚类用 PCA 上的 KMeans（k=25），再按标记基因得分注释各簇
    ad.obs["leiden"] = KMeans(n_clusters=25, n_init=10, random_state=cfg["seed"]).fit_predict(ad.obsm["X_pca"]).astype(str)
    for k, g in ANNOT.items():
        sc.tl.score_genes(ad, [x for x in g if x in ad.var_names], score_name=f"s_{k}", random_state=cfg["seed"])
    cl = ad.obs.groupby("leiden")[[f"s_{k}" for k in ANNOT]].mean()
    lab = cl.idxmax(axis=1).str[2:].where(cl.max(axis=1) > 0.3, "other")
    ad.obs["cell_class"] = ad.obs.leiden.map(lab).astype(str)

    out = ROOT / "results" / STAGE; out.mkdir(parents=True, exist_ok=True)
    comp = ad.obs.groupby(["condition", "cell_class"]).size().unstack(fill_value=0)
    comp.to_csv(out / "cell_counts.tsv", sep="\t")
    rows = []
    for cls, genes in TESTS.items():
        for cond in ["wild_type", "db_db"]:
            s = ad[(ad.obs.cell_class == cls) & (ad.obs.condition == cond)]
            for g in genes:
                if g in s.var_names:
                    x = s[:, g].X.toarray().ravel()
                    rows.append(dict(cell_class=cls, condition=cond, gene=g, n_cells=s.n_obs,
                                     mean_log_expr=float(x.mean()), frac_expressing=float((x > 0).mean())))
    t = pd.DataFrame(rows); t.to_csv(out / "notch_axis_by_condition.tsv", sep="\t", index=False)
    w = t.pivot_table(index=["cell_class", "gene"], columns="condition", values=["mean_log_expr", "frac_expressing"])
    w.columns = [f"{a}_{b}" for a, b in w.columns]
    w["log_ratio_dbdb_vs_wt"] = np.log2((w["mean_log_expr_db_db"] + 0.01) / (w["mean_log_expr_wild_type"] + 0.01))
    w.round(4).to_csv(out / "notch_axis_wide.tsv", sep="\t")
    print(comp); print(w.round(3).to_string())
    write_provenance(STAGE, sorted(d.glob("*")), sorted(out.glob("*.tsv")), seed=cfg["seed"],
                     extra={"note": "one animal per condition; descriptive only"})


if __name__ == "__main__":
    main()
