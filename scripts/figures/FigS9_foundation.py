"""Supplementary Figure S9：Geneformer in silico 删除——周细胞 vs 平滑肌依赖、PDGFRB 分位数、版本间一致性。只读 results/17_foundation_model。"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from figstyle import FULL_W, MM, ROOT, apply_style, panel_label, save_fig, save_source
from Fig4_pdgf_foundation import panel_c as deletion_scatter, panel_d as pdgfrb_percentiles

R = ROOT / "results" / "17_foundation_model"
LAB = {"V1-10M": "V1, 10M", "V2-104M": "V2, 104M", "V2-316M": "V2, 316M"}


def main() -> None:
    apply_style()
    fig = plt.figure(figsize=(FULL_W, 150 * MM))
    top = fig.add_gridspec(1, 2, width_ratios=[1.2, 1], wspace=0.5, top=0.95, bottom=0.6)
    ax = fig.add_subplot(top[0, 0]); save_source(deletion_scatter(ax), "FigS9", "A"); panel_label(ax, "A")
    ax = fig.add_subplot(top[0, 1]); save_source(pdgfrb_percentiles(ax), "FigS9", "B"); panel_label(ax, "B")
    ax = fig.add_axes([0.3, 0.13, 0.55, 0.32])
    c = pd.read_csv(R / "model_concordance.tsv", sep="\t")
    c["pair"] = c.model_a.map(LAB) + " vs " + c.model_b.map(LAB)
    pairs = list(dict.fromkeys(c.pair))
    cols = {("pericyte", "shift_z"): "#F2D0D1", ("pericyte", "identity_z"): "#C44E52",
            ("smooth_muscle", "shift_z"): "#D6E2F2", ("smooth_muscle", "identity_z"): "#4C72B0"}
    labels = {("pericyte", "shift_z"): "Pericytes, embedding shift", ("pericyte", "identity_z"): "Pericytes, identity",
              ("smooth_muscle", "shift_z"): "Smooth muscle, embedding shift", ("smooth_muscle", "identity_z"): "Smooth muscle, identity"}
    for k, key in enumerate(cols):
        d = c[(c.cell_class == key[0]) & (c.measure == key[1])].set_index("pair").reindex(pairs)
        ax.barh(np.arange(len(pairs)) + (k - 1.5) * 0.2, d.spearman, height=0.19, color=cols[key], label=labels[key])
    ax.set_yticks(range(len(pairs))); ax.set_yticklabels(pairs); ax.invert_yaxis(); ax.axvline(0, color="black", lw=0.6)
    ax.set_xlabel("Spearman correlation of gene scores"); ax.set_title("Agreement between model versions", loc="left")
    ax.legend(loc="upper center", bbox_to_anchor=(0.4, -0.3), ncol=2, borderpad=0.3, handletextpad=0.3, columnspacing=0.8)
    save_source(c, "FigS9", "C"); panel_label(ax, "C")
    save_fig(fig, "FigureS9", supplementary=True)


if __name__ == "__main__":
    main()
