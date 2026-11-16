"""Figure 2A/B：单细胞参考的 UMAP（13 类细胞）与关键基因在各细胞类型中的表达点图。只读 results/22_reference_umap。"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from figstyle import CLASS_LABEL, ROOT

R = ROOT / "results" / "22_reference_umap"
ORDER = ["rod", "cone", "bipolar", "amacrine", "horizontal", "ganglion", "muller", "astrocyte", "microglia", "rpe",
         "endothelial", "pericyte", "smooth_muscle"]
PALETTE = {"rod": "#9ECAE1", "cone": "#4292C6", "bipolar": "#08519C", "amacrine": "#A1D99B", "horizontal": "#41AB5D",
           "ganglion": "#006D2C", "muller": "#FDAE6B", "astrocyte": "#E6550D", "microglia": "#BCBDDC", "rpe": "#756BB1",
           "endothelial": "#969696", "pericyte": "#C44E52", "smooth_muscle": "#7B3F00"}
GENES = ["PDGFRB", "RGS5", "IL1RL1", "ADAMTS4", "COL1A2", "COL3A1", "MYH11", "PECAM1", "PDGFB", "TXNIP", "RLBP1", "IL33", "GFAP"]


def panel_umap(ax):
    u = pd.read_csv(R / "umap_cells.tsv", sep="\t", usecols=["cell", "umap1", "umap2", "cell_class"])
    u = u.sample(frac=1, random_state=0)
    ax.scatter(u.umap1, u.umap2, s=0.4, c=u.cell_class.map(PALETTE), lw=0, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    ax.set_title(f"Single-nucleus reference, {len(u):,} nuclei", loc="left")
    return u[["cell", "umap1", "umap2", "cell_class"]]


def panel_dot(ax, cax):
    d = pd.read_csv(R / "gene_by_class.tsv", sep="\t")
    d["scaled"] = d.groupby("gene").mean_log_expr.transform(lambda v: v / v.max())
    for i, cl in enumerate(ORDER):
        for j, g in enumerate(GENES):
            r = d[(d.cell_class == cl) & (d.gene == g)].iloc[0]
            ax.scatter(j, i, s=2 + 55 * r.frac_expressing, c=[[r.scaled]], cmap="Reds", vmin=0, vmax=1, edgecolor="#555555", lw=0.25)
    ax.set_xticks(range(len(GENES))); ax.set_xticklabels(GENES, rotation=45, ha="right", rotation_mode="anchor", fontstyle="italic")
    ax.set_yticks(range(len(ORDER))); ax.set_yticklabels([CLASS_LABEL[c] for c in ORDER])
    ax.set_xlim(-0.6, len(GENES) - 0.4); ax.set_ylim(len(ORDER) - 0.5, -0.5)
    # 细胞类型颜色键（与 A 的 UMAP 颜色一致）：放在 y 轴刻度标签右侧、坐标轴左侧
    for i, cl in enumerate(ORDER):
        ax.scatter(-1.15, i, s=26, marker="s", color=PALETTE[cl], clip_on=False, lw=0)
    ax.tick_params(axis="y", pad=11)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    sm = plt.cm.ScalarMappable(cmap="Reds", norm=plt.Normalize(0, 1))
    cb = plt.colorbar(sm, cax=cax); cb.set_label("Mean expression\n(scaled to maximum)"); cb.outline.set_linewidth(0.5)
    for k, f in enumerate([0.1, 0.5, 1.0]):
        ax.scatter([], [], s=2 + 55 * f, color="#BBBBBB", edgecolor="#555555", lw=0.25, label=f"{int(f * 100)}%")
    ax.legend(title="Nuclei expressing", loc="upper left", bbox_to_anchor=(1.02, -0.02), frameon=False, borderaxespad=0,
              handletextpad=0.2, labelspacing=0.6, title_fontsize=8)
    ax.set_title("Genes of interest across cell types", loc="left")
    return d
