"""视网膜参考的细胞大类映射（stage 08 与 09 共用）。"""
from __future__ import annotations


def hrca_class(ct: str) -> str | None:
    """把 HRCA 细胞类型映射到大类；OFFx 在灵长类是弥散双极细胞亚型。"""
    if ct == "retinal rod cell":
        return "rod"
    if ct in ("retinal cone cell", "S cone cell"):
        return "cone"
    if "bipolar" in ct or ct == "OFFx cell":
        return "bipolar"
    if "amacrine" in ct:
        return "amacrine"
    if "ganglion" in ct:
        return "ganglion"
    if "horizontal" in ct:
        return "horizontal"
    return {"Mueller cell": "muller", "astrocyte": "astrocyte", "microglial cell": "microglia",
            "retinal pigment epithelial cell": "rpe"}.get(ct)


def hrca_composition(cfg):
    import cellxgene_census
    import pandas as pd
    rc = cfg["reference"]
    with cellxgene_census.open_soma(census_version=rc["census_version"]) as census:
        obs = census["census_data"]["homo_sapiens"].obs.read(
            value_filter=f"dataset_id == '{rc['hrca_dataset']}'", column_names=["cell_type", "donor_id", "tissue"]).concat().to_pandas()
    obs = obs.astype(str); obs["cell_class"] = obs.cell_type.map(hrca_class)
    comp = obs.dropna(subset=["cell_class"]).groupby(["donor_id", "tissue", "cell_class"]).size().unstack(fill_value=0)
    comp = comp[comp.sum(axis=1) >= 2000]
    return comp.div(comp.sum(axis=1), axis=0)
