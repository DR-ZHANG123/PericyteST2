"""Stage 08：构建人视网膜单细胞参考（神经元/胶质来自 HRCA snRNA，血管细胞来自同实验室同化学的视神经 snRNA）。

HRCA 神经视网膜 snRNA 不含血管细胞；视神经图谱的毛细血管/静脉内皮、周细胞、平滑肌属于中枢型血管，
与视网膜血管（血–视网膜屏障）最接近，且与 HRCA 同为 10x 3' v3 细胞核测序。有孔内皮（脉络膜型）不纳入。
每个大类按固定种子随机抽取至多 cells_per_class 个细胞核；抽到的 soma_joinid 写入结果表，Census 版本固定。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from repro import ROOT, load_config, set_global_seed, write_provenance  # noqa: E402
from retina_ref import hrca_class  # noqa: E402

STAGE = "08_reference"


def main() -> None:
    import cellxgene_census
    cfg = load_config(); rc = cfg["reference"]; set_global_seed(cfg["seed"])
    rng = np.random.default_rng(cfg["seed"])
    cols = ["soma_joinid", "dataset_id", "cell_type", "tissue", "donor_id", "sex", "development_stage"]
    with cellxgene_census.open_soma(census_version=rc["census_version"]) as census:
        hs = census["census_data"]["homo_sapiens"]
        obs = hs.obs.read(value_filter=f"dataset_id in ['{rc['hrca_dataset']}', '{rc['vascular_dataset']}']",
                          column_names=cols).concat().to_pandas()
        for c in ["cell_type", "tissue", "donor_id", "dataset_id"]:
            obs[c] = obs[c].astype(str)
        vmap = {ct: k for k, v in rc["vascular_classes"].items() for ct in v}
        obs["cell_class"] = np.where(obs.dataset_id == rc["hrca_dataset"], obs.cell_type.map(hrca_class),
                                     obs.cell_type.map(vmap))
        obs = obs[obs.cell_class.notna()]
        picked = []
        for cls, g in obs.groupby("cell_class"):
            n = min(rc["cells_per_class"], len(g))
            picked.append(g.iloc[np.sort(rng.choice(len(g), n, replace=False))])
        sel = pd.concat(picked).sort_values("soma_joinid")
        ad = cellxgene_census.get_anndata(census, organism="Homo sapiens", obs_coords=sel.soma_joinid.to_numpy(),
                                          column_names={"obs": cols, "var": ["feature_id", "feature_name"]})
    ad.obs = ad.obs.astype({c: str for c in ["cell_type", "tissue", "donor_id", "dataset_id"]})
    ad.obs = ad.obs.merge(sel[["soma_joinid", "cell_class"]], on="soma_joinid", how="left").set_index(ad.obs.index)
    ad.var_names = ad.var["feature_name"].astype(str).values
    ad.var_names_make_unique()
    ad = ad[:, np.asarray((ad.X > 0).sum(axis=0)).ravel() >= 10].copy()
    out_h5 = ROOT / rc["out_h5ad"]; out_h5.parent.mkdir(parents=True, exist_ok=True)
    ad.write_h5ad(out_h5, compression="gzip")

    out = ROOT / "results" / STAGE; out.mkdir(parents=True, exist_ok=True)
    ad.obs[["soma_joinid", "dataset_id", "cell_type", "cell_class", "tissue", "donor_id"]].to_csv(
        out / "reference_cells.tsv", sep="\t", index=False)
    summ = ad.obs.groupby("cell_class").agg(n_cells=("soma_joinid", "size"), n_donors=("donor_id", "nunique")).reset_index()
    summ.to_csv(out / "reference_summary.tsv", sep="\t", index=False)
    print(summ.to_string(index=False)); print(ad)
    write_provenance(STAGE, [], [out / "reference_cells.tsv", out / "reference_summary.tsv"], seed=cfg["seed"],
                     extra={"census_version": rc["census_version"], "n_genes": int(ad.n_vars), "h5ad": rc["out_h5ad"]})


if __name__ == "__main__":
    main()
