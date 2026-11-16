"""Stage 14a：定义各细胞类群的稳定特异标记集合（用于丰度指数与稳健性检验）。

标记来自单细胞参考（与 bulk 无关）：在平均视网膜组成下 bulk 贡献份额 ≥ 0.5、该类细胞内平均计数 ≥ 0.5、
蛋白编码基因（排除 LINC/LOC/ENSG/-AS/-DT/MIR/SNOR 命名），并预先排除生物学关注基因。
集合：壁细胞、内皮、Müller 各取份额最高的前 5/8/12/20 个；另设只来自周细胞（周细胞计数 ≥ 2 倍平滑肌）、
只来自平滑肌（反之）的集合；节细胞与视杆作部位阳性对照；无长突细胞与小胶质是 PDGFB 的另两个来源（stage 15b 校正用）。
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import anndata
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from deconv import VASC  # noqa: E402
from repro import ROOT, load_config, set_global_seed, write_provenance  # noqa: E402

STAGE = "14a_marker_sets"
NONCODING = re.compile(r"^(LINC|LOC|ENSG|MIR|SNOR|RNU|RN7|Y_RNA)|-(AS|DT|IT|OT)\d*$|^AC\d|^AL\d|^AP\d")
GROUPS = {"mural": ["pericyte", "smooth_muscle"], "endothelial": ["endothelial"], "muller": ["muller"],
          "pericyte_only": ["pericyte"], "smc_only": ["smooth_muscle"], "ganglion": ["ganglion"], "rod": ["rod"],
          "amacrine": ["amacrine"], "microglia": ["microglia"]}   # 后两者：PDGFB 的其他来源细胞（stage 15b）


def main() -> None:
    cfg = load_config(); set_global_seed(cfg["seed"])
    out = ROOT / "results" / STAGE; out.mkdir(parents=True, exist_ok=True)
    ref = anndata.read_h5ad(ROOT / cfg["reference"]["out_h5ad"])
    labels = ref.obs["cell_class"].astype(str).values; classes = sorted(set(labels)); genes = np.asarray(ref.var_names)
    X = ref.X
    raw = pd.DataFrame(np.vstack([np.asarray(X[labels == c].mean(axis=0)).ravel() for c in classes]).T, index=genes, columns=classes)
    comp = pd.read_csv(ROOT / cfg["reference"]["composition"], sep="\t", index_col=[0, 1])
    f = comp.mean().reindex(classes).fillna(0)
    for v in VASC:
        f[v] = 0.01
    f /= f.sum()
    w = raw * f; share = w.div(w.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    excl = set(cfg["deconvolution"]["ratio_exclude_markers"])
    coding = pd.Series([not NONCODING.search(g) for g in genes], index=genes)
    rows = []
    for g, members in GROUPS.items():
        sh = share[members].sum(axis=1); rw = raw[members].max(axis=1)
        ok = (sh >= 0.5) & (rw >= 0.5) & coding & ~share.index.isin(excl)
        if g == "pericyte_only":
            ok &= raw["pericyte"] >= 2 * raw["smooth_muscle"]
        if g == "smc_only":
            ok &= raw["smooth_muscle"] >= 2 * raw["pericyte"]
        ranked = sh[ok].sort_values(ascending=False)
        sizes = [5, 8, 12, 20] if g in ("mural", "endothelial", "muller") else [8]
        for n in sizes:
            for gene in ranked.index[:n]:
                rows.append(dict(set=f"{g}_top{n}", group=g, gene=gene, share=round(float(sh[gene]), 4),
                                 mean_count_in_group=round(float(rw[gene]), 3)))
    ms = pd.DataFrame(rows); ms.to_csv(out / "marker_sets.tsv", sep="\t", index=False)
    print(ms.groupby("set").gene.apply(lambda x: ", ".join(x)).to_string())
    write_provenance(STAGE, [ROOT / cfg["reference"]["out_h5ad"], ROOT / cfg["reference"]["composition"]],
                     [out / "marker_sets.tsv"], seed=cfg["seed"])


if __name__ == "__main__":
    main()
