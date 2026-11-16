"""Figure 4：糖尿病供体视网膜中相对内皮标记的 PDGFB 升高；其他内皮来源、与壁细胞相关的信号基因不变（TXNIP 在 NPDR 升高）。
panel_c/panel_d（Geneformer）由补充图 S7 调用。
只读 results/15_crosstalk、results/15b_pdgfb_source、results/15_crosstalk、results/17_foundation_model。"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from figstyle import CLASS_LABEL, FULL_W, MM, ROOT, STAGE_COLORS, STAGE_LABEL, apply_style, panel_label, save_fig, save_source, stage_order

R = ROOT / "results"
CONTRASTS = [("groupDiabetic", "Diabetes, no DR"), ("groupNPDR", "NPDR")]
C_COL = [STAGE_COLORS["Diabetes without retinopathy"], STAGE_COLORS["NPDR"]]
MODELS = [("V1-10M", "V1, 10M"), ("V2-104M", "V2, 104M"), ("V2-316M", "V2, 316M")]
CELL_COL = {"pericyte": "#C44E52", "smooth_muscle": "#4C72B0"}


def fmt_p(p):
    return "P < 0.001" if p < 0.001 else f"P = {p:.3f}" if p < 0.1 else f"P = {p:.2f}"


def panel_a(ax):
    d = pd.read_csv(R / "15b_pdgfb_source" / "pdgfb_by_donor.tsv", sep="\t")
    d["centred"] = d.pdgfb_ratio - d.loc[d.group == "Control", "pdgfb_ratio"].mean()
    d["stage"] = d.group.map(STAGE_LABEL)
    t = pd.read_csv(R / "15b_pdgfb_source" / "pdgfb_source_tests.tsv", sep="\t")
    t = t[t.measure == "PDGFB relative to endothelial markers"]
    for k, s in enumerate(stage_order()):
        v = d[d.stage == s].centred.values
        ax.boxplot(v, positions=[k], widths=0.5, showfliers=False, medianprops=dict(color="black", lw=1),
                   boxprops=dict(lw=0.6), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
        ax.scatter(k + np.random.default_rng(k).uniform(-0.15, 0.15, len(v)), v, s=11, color=STAGE_COLORS[s], edgecolor="white", lw=0.3, zorder=3)
    top = d.centred.max()
    for k, (cf, _) in enumerate(CONTRASTS, start=1):
        ax.text(k, top + 0.12 + (k - 1) * 0.36, fmt_p(t[t.coef == cf].p_perm.iloc[0]), ha="center", va="bottom")
    ax.set_ylim(d.centred.min() - 0.2, top + 1.0); ax.axhline(0, color="#CCCCCC", lw=0.6, zorder=0)
    ax.set_xticks(range(3)); ax.set_xticklabels(["Ctrl", "No\nDR", "NPDR"])
    ax.set_ylabel("PDGFB / endothelial markers\nvs control (log$_2$)"); ax.set_title("PDGFB", loc="left")
    return d


def panel_b(ax):
    t = pd.read_csv(R / "15b_pdgfb_source" / "pdgfb_source_tests.tsv", sep="\t")
    f = pd.read_csv(R / "15_crosstalk" / "followup_robustness.tsv", sep="\t"); f = f[f.gene == "PDGFB"]
    items = [("PDGFB relative to endothelial markers", "Ratio, 8 markers"),
             ("relative to endothelial markers, adjusted for amacrine and microglial indices", "+ amacrine, microglia"),
             ("relative to endothelial markers, adjusted for amacrine, microglial and pericyte indices", "+ also pericytes"),
             ("bulk PDGFB adjusted for endothelial, amacrine and microglial indices", "Bulk, cell-adjusted")]
    rows = []
    for meas, lab in items:
        for cf, cl in CONTRASTS:
            r = t[(t.measure == meas) & (t.coef == cf)].iloc[0]
            rows.append(dict(item=lab, contrast=cl, diff=r["diff"], lo=r.ci_low, hi=r.ci_high, p=r.p_perm))
    for st, lab in [("endothelial_top5", "Ratio, 5 markers"), ("endothelial_top12", "Ratio, 12 markers"), ("endothelial_top20", "Ratio, 20 markers")]:
        for cf, cl in CONTRASTS:
            r = f[(f.measure == st) & (f.region == "both") & (f.coef == cf)].iloc[0]
            rows.append(dict(item=lab, contrast=cl, diff=r["diff"], lo=np.nan, hi=np.nan, p=r.p))
    for rg in ["macula", "periphery"]:
        for cf, cl in CONTRASTS:
            r = f[(f.measure == "endothelial_top8") & (f.region == rg) & (f.coef == cf)].iloc[0]
            rows.append(dict(item=f"{rg.capitalize()} only", contrast=cl, diff=r["diff"], lo=np.nan, hi=np.nan, p=r.p))
    df = pd.DataFrame(rows)
    labs = list(dict.fromkeys(df.item))
    for i, lab in enumerate(labs):
        for j, (cf, cl) in enumerate(CONTRASTS):
            r = df[(df.item == lab) & (df.contrast == cl)].iloc[0]; y = i + (j - 0.5) * 0.32
            if np.isfinite(r.lo):
                ax.errorbar(r["diff"], y, xerr=[[r["diff"] - r.lo], [r.hi - r["diff"]]], fmt="o", ms=3.4, color=C_COL[j], lw=0.8, capsize=1.2)
            else:
                ax.scatter(r["diff"], y, s=14, color=C_COL[j], marker="s", zorder=3)
    ax.axvline(0, color="black", lw=0.6); ax.set_yticks(range(len(labs))); ax.set_yticklabels(labs); ax.invert_yaxis()
    ax.set_xlabel("Difference vs control (log$_2$)"); ax.set_title("PDGFB across measures", loc="left")
    return df


def panel_c(ax):
    g = pd.read_csv(R / "17_foundation_model" / "gene_dependency.tsv", sep="\t")
    g = g[g.model == "V2-104M"]
    p = g[g.cell_class == "pericyte"].set_index("gene_id"); s = g[g.cell_class == "smooth_muscle"].set_index("gene_id")
    both = p.index.intersection(s.index)
    d = pd.DataFrame({"symbol": p.loc[both, "symbol"], "receptor": p.loc[both, "is_receptor"],
                      "pericyte": p.loc[both, "identity_z"], "smooth_muscle": s.loc[both, "identity_z"]})
    oth = d[~d.receptor]; rec = d[d.receptor]
    ax.scatter(oth.smooth_muscle, oth.pericyte, s=4, color="#D0D0D0", lw=0, zorder=1)
    ax.scatter(rec.smooth_muscle, rec.pericyte, s=10, color="#4C72B0", lw=0, zorder=2)
    lim = [min(d.min(numeric_only=True)) - 0.1, max(d.max(numeric_only=True)) + 0.1]
    ax.plot(lim, lim, color="#999999", lw=0.6, ls="--", zorder=0)
    pd_row = d[d.symbol == "PDGFRB"].iloc[0]
    ax.scatter(pd_row.smooth_muscle, pd_row.pericyte, s=26, color="#C44E52", edgecolor="black", lw=0.5, zorder=4)
    ax.annotate("PDGFRB", (pd_row.smooth_muscle, pd_row.pericyte), xytext=(-38, 22), textcoords="offset points",
                fontstyle="italic", arrowprops=dict(arrowstyle="-", lw=0.5, color="black", shrinkA=0, shrinkB=3))
    ax.set_xlabel("Smooth-muscle identity dependency (z)"); ax.set_ylabel("Pericyte identity dependency (z)")
    ax.set_title("In silico deletion, V2-104M", loc="left")
    return d


def panel_d(ax):
    g = pd.read_csv(R / "17_foundation_model" / "gene_dependency.tsv", sep="\t")
    g = g[g.symbol == "PDGFRB"]
    rows = []
    for k, (m, lab) in enumerate(MODELS):
        for c, off in [("smooth_muscle", -0.12), ("pericyte", 0.12)]:
            r = g[(g.model == m) & (g.cell_class == c)]
            if len(r):
                v = 100 * r.identity_pct.iloc[0]
                ax.scatter(k + off, v, s=26, color=CELL_COL[c], zorder=3)
                rows.append(dict(model=lab, cell_class=c, percentile=v))
        sub = [x for x in rows if x["model"] == lab]
        if len(sub) == 2:
            ax.plot([k - 0.12, k + 0.12], [sub[0]["percentile"], sub[1]["percentile"]], color="#999999", lw=0.7, zorder=2)
    ax.axhline(50, color="#CCCCCC", lw=0.6, ls="--", zorder=0)
    ax.set_xticks(range(len(MODELS))); ax.set_xticklabels([m[1].replace(", ", "\n") for m in MODELS]); ax.set_xlim(-0.5, len(MODELS) - 0.5)
    ax.set_ylim(0, 135); ax.set_yticks([0, 25, 50, 75, 100]); ax.set_ylabel("PDGFRB percentile among genes")
    ax.set_title("PDGFRB dependency by model", loc="left")
    ax.scatter([], [], s=26, color=CELL_COL["pericyte"], label="Pericytes"); ax.scatter([], [], s=26, color=CELL_COL["smooth_muscle"], label="Smooth muscle")
    ax.legend(loc="upper left", ncol=1, borderpad=0.3, handletextpad=0.2)
    return pd.DataFrame(rows)


def panel_ligands(ax):
    """文献中作用于壁细胞或血管应激的内皮来源基因（参考中内皮占比最高者），相对内皮标记；PDGFB 取 A 的同一模型。"""
    s = pd.read_csv(R / "15_crosstalk" / "source_ratio_group_tests.tsv", sep="\t")
    s = s[(s.source == "endothelial") & s.literature & (s.gene != "PDGFB")][["gene", "coef", "diff", "p_perm"]]
    t = pd.read_csv(R / "15b_pdgfb_source" / "pdgfb_source_tests.tsv", sep="\t")
    t = t[t.measure == "PDGFB relative to endothelial markers"].assign(gene="PDGFB")[["gene", "coef", "diff", "p_perm"]]
    d = pd.concat([t, s], ignore_index=True)
    order = ["PDGFB"] + list(s[s.coef == "groupNPDR"].sort_values("diff", ascending=False).gene)
    for i, g in enumerate(order):
        for j, (cf, _) in enumerate(CONTRASTS):
            r = d[(d.gene == g) & (d.coef == cf)].iloc[0]; y = i + (j - 0.5) * 0.34
            filled = r.p_perm < 0.05
            ax.scatter(r["diff"], y, s=22, color=C_COL[j] if filled else "white", edgecolor=C_COL[j], lw=0.9, zorder=3)
    ax.axvline(0, color="black", lw=0.6); ax.set_yticks(range(len(order))); ax.set_yticklabels(order, fontstyle="italic")
    ax.invert_yaxis(); ax.set_xlabel("Difference vs control (log$_2$)")
    ax.set_title("Endothelial genes, own markers", loc="left")
    ax.scatter([], [], s=22, color="#777777", label="P < 0.05"); ax.scatter([], [], s=22, color="white", edgecolor="#777777", lw=0.9, label="P ≥ 0.05")
    ax.legend(loc="lower right", handletextpad=0.1, borderpad=0.3)
    return d


def main() -> None:
    apply_style()
    fig = plt.figure(figsize=(FULL_W, 80 * MM))
    r1 = fig.add_gridspec(1, 3, width_ratios=[0.72, 1.0, 0.8], wspace=1.05, left=0.08, right=0.99, top=0.88, bottom=0.2)
    for spec, fn, letter in [(r1[0, 0], panel_a, "A"), (r1[0, 1], panel_b, "B"), (r1[0, 2], panel_ligands, "C")]:
        ax = fig.add_subplot(spec)
        save_source(fn(ax), "Fig4", letter)
        panel_label(ax, letter)
    hs = [plt.Line2D([], [], marker="o", ls="", color=c, ms=4.5) for c in C_COL]
    fig.legend(hs, ["Diabetes without DR", "NPDR"], loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.0))
    save_fig(fig, "Figure4")


if __name__ == "__main__":
    main()
