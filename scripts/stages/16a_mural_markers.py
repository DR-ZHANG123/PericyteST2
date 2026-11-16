"""Stage 16a：周细胞 / 平滑肌标记的多种定义（严格、宽松、经典），并导出壁细胞相关基因的参考谱。

用于检验周细胞转录本减少是否独立于平滑肌（活体成像报道 DR 前先减少的是小动脉平滑肌）。
"""
from __future__ import annotations

import sys
from pathlib import Path

import anndata
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from repro import ROOT, load_config, set_global_seed, write_provenance  # noqa: E402

STAGE = "16a_mural_markers"
OTHER = {"pericyte": "smooth_muscle", "smooth_muscle": "pericyte"}


def main() -> None:
    cfg = load_config(); set_global_seed(cfg["seed"]); mc = cfg["mural_subtypes"]
    out = ROOT / "results" / STAGE; out.mkdir(parents=True, exist_ok=True)
    ref = anndata.read_h5ad(ROOT / cfg["reference"]["out_h5ad"])
    lab = ref.obs["cell_class"].astype(str).values
    means = pd.DataFrame({c: np.asarray(ref.X[lab == c].mean(axis=0)).ravel() for c in ("pericyte", "smooth_muscle")},
                         index=ref.var_names)
    share_path = ROOT / "results" / "09c_ratio_validation" / "bulk_share_by_celltype.tsv"
    share = pd.read_csv(share_path, sep="\t", index_col=0)
    prof = means.join(share.add_prefix("share_"), how="inner")
    coding = ~prof.index.str.contains(r"^(?:LINC|LOC|ENSG|MIR|SNOR|AC\d|AL\d|AP\d)|-(?:AS|DT|IT|OT)\d*$", regex=True)
    rows = []
    for name, d in mc["sets"].items():
        c, o = d["cell"], OTHER[d["cell"]]
        ok = (prof[f"share_{c}"] >= d["share"]) & (prof[c] >= mc["min_group_mean"]) \
            & (prof[c] >= d["fold"] * prof[o]) & coding
        for g in prof.index[ok]:
            rows.append(dict(set=name, gene=g))
    for name, genes in mc["classic"].items():
        rows += [dict(set=name, gene=g) for g in genes if g in prof.index]
    ms = pd.DataFrame(rows)
    ms = ms.merge(prof, left_on="gene", right_index=True)
    ms.to_csv(out / "mural_marker_sets.tsv", sep="\t", index=False)
    print(ms.groupby("set").gene.apply(lambda x: f"{len(x)}: " + ", ".join(x)).to_string())
    write_provenance(STAGE, [ROOT / cfg["reference"]["out_h5ad"], share_path], [out / "mural_marker_sets.tsv"], seed=cfg["seed"])


if __name__ == "__main__":
    main()
