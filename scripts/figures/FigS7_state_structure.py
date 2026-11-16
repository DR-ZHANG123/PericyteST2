"""补充图 S7：壁细胞状态基因在供体间的共变。
A 供体 × 基因热图（相对壁细胞含量，减去对照均值）；B 组内相关矩阵（Spearman，校正组别与协变量）。
只读 results/21_state_structure。"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from figstyle import FULL_W, MM, ROOT, STAGE_COLORS, STAGE_LABEL, apply_style, num, panel_label, save_fig, save_source

R = ROOT / "results" / "21_state_structure"
GENES = ["IL1RL1", "ADAMTS4", "EDNRA", "GJA4", "COL1A1", "COL1A2", "COL3A1", "COL6A3"]
GROUPS = [("Control", "Control"), ("Diabetic", "Diabetes without DR"), ("NPDR", "NPDR"), ("PDR", "PDR")]


def panel_heat(ax, cax):
    d = pd.read_csv(R / "state_genes_by_donor.tsv", sep="\t")
    d["g"] = pd.Categorical(d.group, [g for g, _ in GROUPS], ordered=True)
    d = d.sort_values(["g", "severity", "donor"]).reset_index(drop=True)
    M = d[GENES].T.values
    im = ax.imshow(np.clip(M, -2, 2), cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto", interpolation="nearest")
    ax.set_yticks(range(len(GENES))); ax.set_yticklabels(GENES, fontstyle="italic"); ax.set_xticks([])
    edges = np.cumsum([0] + [int((d.group == g).sum()) for g, _ in GROUPS])
    for e in edges[1:-1]:
        ax.axvline(e - 0.5, color="black", lw=0.8)
    ax.axhline(3.5, color="white", lw=1.5)
    for (g, lab), a, b in zip(GROUPS, edges[:-1], edges[1:]):
        ax.add_patch(plt.Rectangle((a - 0.5, -1.45), b - a, 0.8, color=STAGE_COLORS[STAGE_LABEL[g]], clip_on=False, lw=0))
        ax.text((a + b - 1) / 2, -1.05, lab.replace("\n", " "), ha="center", va="center", color="white", fontsize=7)
    ax.set_xlim(-0.5, len(d) - 0.5); ax.set_ylim(len(GENES) - 0.5, -0.5)
    ax.set_xlabel("Donors, ordered by group and ETDRS level")
    cb = plt.colorbar(im, cax=cax); cb.set_label("vs control (log$_2$)"); cb.outline.set_linewidth(0.5)
    ax.set_title("State genes relative to mural content, each donor", loc="left", pad=30)
    return d[["donor", "group", "severity"] + GENES]


def panel_corr(ax, cax):
    c = pd.read_csv(R / "state_gene_correlations.tsv", sep="\t")
    n = len(GENES); M = np.full((n, n), np.nan); F = np.full((n, n), np.nan)
    for _, r in c.iterrows():
        i, j = GENES.index(r.gene_a), GENES.index(r.gene_b)
        a, b = max(i, j), min(i, j); M[a, b] = r.rho; F[a, b] = r.fdr
    im = ax.imshow(np.ma.masked_invalid(M), cmap="RdBu_r", vmin=-0.8, vmax=0.8)
    for a in range(n):
        for b in range(a):
            ax.text(b, a, num(M[a, b], 2) + ("*" if F[a, b] < 0.05 else ""), ha="center", va="center", fontsize=7,
                    color="white" if abs(M[a, b]) > 0.55 else "black")
    ax.set_xticks(range(n - 1)); ax.set_xticklabels(GENES[:-1], rotation=45, ha="right", fontstyle="italic")
    ax.set_yticks(range(1, n)); ax.set_yticklabels(GENES[1:], fontstyle="italic")
    ax.set_xlim(-0.5, n - 1.5); ax.set_ylim(n - 0.5, 0.5)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = plt.colorbar(im, cax=cax); cb.set_label("Spearman ρ"); cb.outline.set_linewidth(0.5)
    ax.set_title("Co-variation across donors", loc="left")
    return c


def main() -> None:
    apply_style()
    fig = plt.figure(figsize=(FULL_W, 150 * MM))
    g1 = fig.add_gridspec(1, 2, width_ratios=[1, 0.025], wspace=0.04, left=0.1, right=0.9, top=0.89, bottom=0.6)
    g2 = fig.add_gridspec(1, 2, width_ratios=[1, 0.04], wspace=0.06, left=0.1, right=0.56, top=0.48, bottom=0.09)
    ax = fig.add_subplot(g1[0, 0]); save_source(panel_heat(ax, fig.add_subplot(g1[0, 1])), "FigS7", "A"); panel_label(ax, "A")
    ax = fig.add_subplot(g2[0, 0]); save_source(panel_corr(ax, fig.add_subplot(g2[0, 1])), "FigS7", "B"); panel_label(ax, "B")
    save_fig(fig, "FigureS7", supplementary=True)


if __name__ == "__main__":
    main()
