"""Figure 5：人视网膜周细胞（HRPC）中高糖削弱 PDGF-BB 信号与 ST2 诱导；敲低 IL1RL1 重现供体中的胶原高、ADAMTS4 低状态。
只读 results/23_wetlab/（measurements.tsv、contrasts.tsv、blots/）；作图输入写到 figures/source_data/Fig5_*.csv。
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

from figstyle import FULL_W, MM, ROOT, apply_style, panel_label, save_fig, save_source

R = ROOT / "results" / "23_wetlab"
GROUPS = ["NG", "NG+PDGF-BB", "HG", "HG+PDGF-BB", "MAN", "MAN+PDGF-BB"]
XPOS = [0, 1, 2.6, 3.6, 5.2, 6.2]
COND_COL = {"NG": "#8C8C8C", "HG": "#C44E52", "MAN": "#4C72B0"}
SI_GROUPS = ["siNC", "siIL1RL1-1", "siIL1RL1-2"]
SI_COL = {"siNC": "#8C8C8C", "siIL1RL1-1": "#DD8452", "siIL1RL1-2": "#937860"}

M = pd.read_csv(R / "measurements.tsv", sep="\t")
C = pd.read_csv(R / "contrasts.tsv", sep="\t")


def stars(p: float) -> str:
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."


def pval(exp, rd, contrast):
    return C[(C.experiment == exp) & (C.readout == rd) & (C.contrast == contrast)].p.iloc[0]


def bracket(ax, x0, x1, y, label, h):
    ax.plot([x0, x0, x1, x1], [y, y + h, y + h, y], color="black", lw=0.6)
    ax.text((x0 + x1) / 2, y + h, label, ha="center", va="bottom", fontsize=7 if label == "n.s." else 8)


def factorial(ax, exp, rd, ylabel, title, brackets):
    d = M[(M.experiment == exp) & (M.readout == rd)]
    rng = np.random.default_rng(23)
    tops = []
    for x, g in zip(XPOS, GROUPS):
        v = d[d.group == g].value.values; cond = g.split("+")[0]; pdgf = "PDGF" in g
        ax.bar(x, v.mean(), width=0.78, color=COND_COL[cond] if pdgf else "white", edgecolor=COND_COL[cond], lw=0.9, zorder=2)
        ax.errorbar(x, v.mean(), yerr=v.std(ddof=1), color="black", lw=0.6, capsize=1.5, zorder=3)
        ax.scatter(x + rng.uniform(-0.16, 0.16, len(v)), v, s=6, color="black", lw=0, zorder=4)
        tops.append(max(v.max(), v.mean() + v.std(ddof=1)))
    ymax = max(tops); h = 0.035 * ymax
    lvl = ymax * 1.08
    for (i0, i1, contrast) in brackets:
        lo = max(tops[min(i0, i1):max(i0, i1) + 1]) + 0.06 * ymax
        y = max(lo, lvl if contrast != brackets[0][2] else lo)
        bracket(ax, XPOS[i0], XPOS[i1], y, stars(pval(exp, rd, contrast)), h)
        lvl = y + 0.2 * ymax
    ax.set_ylim(0, lvl + 0.02 * ymax)
    ax.set_xticks(XPOS); ax.set_xticklabels(["−", "+"] * 3)
    ax.set_xticks([0.5, 3.1, 5.7], minor=True); ax.set_xticklabels(["NG", "HG", "MAN"], minor=True)
    ax.tick_params(axis="x", which="minor", length=0, pad=10); ax.tick_params(axis="x", which="major", length=0, pad=1.5)
    ax.set_xlim(-0.7, 6.9); ax.set_ylabel(ylabel); ax.set_title(title, loc="left")
    return d.assign(panel_readout=rd)


def sirna(ax, exp, genes, ylabel, title):
    d = M[(M.experiment == exp) & (M.readout.isin(genes))]
    rng = np.random.default_rng(23); w = 0.26; ymax = d.value.max()
    for i, gene in enumerate(genes):
        for j, g in enumerate(SI_GROUPS):
            v = d[(d.readout == gene) & (d.group == g)].value.values; x = i + (j - 1) * w
            ax.bar(x, v.mean(), width=w * 0.9, color=SI_COL[g], lw=0, zorder=2)
            ax.errorbar(x, v.mean(), yerr=v.std(ddof=1), color="black", lw=0.6, capsize=1.2, zorder=3)
            ax.scatter(x + rng.uniform(-0.05, 0.05, len(v)), v, s=5, color="black", lw=0, zorder=4)
            if g != "siNC":
                ax.text(x, max(v.max(), v.mean() + v.std(ddof=1)) + 0.04 * ymax, stars(pval(exp, gene, f"{g} vs siNC")),
                        ha="center", va="bottom", fontsize=7, rotation=90)
    ax.axhline(1, color="#BBBBBB", lw=0.6, ls="--", zorder=0)
    ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes, fontstyle="italic"); ax.tick_params(axis="x", length=0)
    ax.set_ylim(0, ymax * 1.32); ax.set_ylabel(ylabel); ax.set_title(title, loc="left")
    return d


def blots(ax):
    names = [("p_PDGFRB", "p-PDGFRβ", "170 kDa"), ("PDGFRB", "PDGFRβ", "170 kDa"), ("GAPDH", "GAPDH", "37 kDa")]
    W = 600; y = 0; rows = []
    for key, lab, kda in names:
        im = Image.open(R / "blots" / f"strip_{key}.png").convert("L")
        hh = round(im.height * W / im.width * 0.9)
        arr = np.asarray(im.resize((W, hh)))
        ax.imshow(arr, cmap="gray", vmin=0, vmax=255, extent=(0, W, y + hh, y), aspect="auto")
        ax.text(-12, y + hh / 2, lab, ha="right", va="center"); ax.text(W + 12, y + hh / 2, kda, ha="left", va="center")
        rows.append(dict(blot=lab, size=kda, source=f"results/23_wetlab/blots/strip_{key}.png")); y += hh + 14
    for k in range(6):
        ax.text((k + 0.5) * W / 6, -8, "−+"[k % 2], ha="center", va="bottom")
    for k, cond in enumerate(["NG", "HG", "MAN"]):
        ax.text((2 * k + 1) * W / 6, -48, cond, ha="center", va="bottom")
        ax.plot([(2 * k + 0.15) * W / 6, (2 * k + 1.85) * W / 6], [-42, -42], color="black", lw=0.6)
    ax.text(-12, -8, "PDGF-BB", ha="right", va="bottom")
    ax.set_xlim(0, W); ax.set_ylim(y, -2); ax.set_xticks([]); ax.set_yticks([]); ax.axis("off")
    ax.set_title("Western blot, PDGF-BB 10 min", loc="left", pad=40)
    return pd.DataFrame(rows)


def main() -> None:
    apply_style()
    plt.rcParams.update({"mathtext.fontset": "custom", "mathtext.it": "Liberation Sans:italic", "mathtext.rm": "Liberation Sans"})
    fig = plt.figure(figsize=(FULL_W, 182 * MM))
    g1 = fig.add_gridspec(1, 3, width_ratios=[1.3, 0.9, 0.9], wspace=0.75, left=0.11, right=0.99, top=0.9, bottom=0.71)
    g2 = fig.add_gridspec(1, 5, wspace=0.62, left=0.07, right=0.99, top=0.6, bottom=0.42)
    g3 = fig.add_gridspec(1, 3, width_ratios=[0.8, 1.1, 1.1], wspace=0.42, left=0.07, right=0.99, top=0.3, bottom=0.08)
    src = {}
    ax = fig.add_subplot(g1[0, 0]); src["A"] = blots(ax); panel_label(ax, "A")
    ax = fig.add_subplot(g1[0, 1])
    src["B"] = factorial(ax, "western_blot", "p-PDGFRβ/PDGFRβ", "p-PDGFRβ / PDGFRβ\n(relative to NG)", "Receptor phosphorylation",
                         [(1, 3, "Interaction HG (PDGF-BB effect in HG minus NG)"), (3, 5, "Interaction HG vs MAN (PDGF-BB effect in HG minus MAN)")])
    panel_label(ax, "B")
    ax = fig.add_subplot(g1[0, 2])
    src["C"] = factorial(ax, "qpcr_1h", "TXNIP", "Relative mRNA", r"$\it{TXNIP}$, 1 h",
                         [(0, 2, "HG vs NG (no PDGF-BB)"), (0, 4, "MAN vs NG (no PDGF-BB)")])
    panel_label(ax, "C")
    br = [(1, 3, "Interaction HG (PDGF-BB effect in HG minus NG)"), (3, 5, "Interaction HG vs MAN (PDGF-BB effect in HG minus MAN)")]
    parts = []
    for k, (exp, rd) in enumerate([("qpcr_1h", "EGR1"), ("qpcr_1h", "FOS"), ("qpcr_24h", "sST2"), ("qpcr_24h", "ST2L"), ("qpcr_24h", "IL1RL1")]):
        ax = fig.add_subplot(g2[0, k])
        title = {"IL1RL1": r"$\it{IL1RL1}$ total, 24 h"}.get(rd, rf"$\it{{{rd}}}$, {'1' if exp == 'qpcr_1h' else '24'} h")
        parts.append(factorial(ax, exp, rd, "Relative mRNA" if k in (0, 2) else "", title, br))
        if k == 0:
            panel_label(ax, "D")
        if k == 2:
            panel_label(ax, "E")
    src["D"] = pd.concat(parts[:2]); src["E"] = pd.concat(parts[2:])
    ax = fig.add_subplot(g3[0, 0])
    src["F"] = factorial(ax, "elisa_24h", "sST2 (pg/mL)", "sST2 in medium (pg/mL)", "Secreted sST2, 24 h", br); panel_label(ax, "F")
    ax = fig.add_subplot(g3[0, 1])
    src["G"] = sirna(ax, "sirna_knockdown", ["IL1RL1", "ST2L", "sST2"], "Relative mRNA", "Knock-down, 72 h"); panel_label(ax, "G")
    ax = fig.add_subplot(g3[0, 2])
    src["H"] = sirna(ax, "sirna_targets", ["COL1A2", "COL3A1", "ADAMTS4"], "Relative mRNA", "Matrix genes after knock-down, 72 h"); panel_label(ax, "H")
    hs = [plt.Rectangle((0, 0), 1, 1, color=SI_COL[g]) for g in SI_GROUPS]
    ax.legend(hs, SI_GROUPS, loc="upper right", handlelength=0.9, borderpad=0.3, fontsize=7)
    for k, d in src.items():
        save_source(d, "Fig5", k)
    save_fig(fig, "Figure5")


if __name__ == "__main__":
    main()
