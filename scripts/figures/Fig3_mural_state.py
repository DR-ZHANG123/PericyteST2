"""Figure 3：NPDR 中壁细胞含量下降（周细胞与平滑肌相近）；IL1RL1 相对周细胞自身标记下降、胶原上升，不由周细胞/平滑肌组成解释。
只读 results/18a_mural_state、results/19_state_support、results/20_composition_check。"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from figstyle import FULL_W, MM, ROOT, STAGE_COLORS, STAGE_LABEL, apply_style, panel_label, save_fig, save_source, stage_order

R = ROOT / "results"
CONTRASTS = [("groupDiabetic", "Diabetes, no DR"), ("groupNPDR", "NPDR")]
C_COL = [STAGE_COLORS["Diabetes without retinopathy"], STAGE_COLORS["NPDR"]]
LEVELS = ["Control", "10", "20", "35", "43", "53", "65"]


def fmt_p(p):
    return "P < 0.001" if p < 0.001 else f"P = {p:.3f}" if p < 0.1 else f"P = {p:.2f}"


def by_donor(d, col, groups=("Control", "Diabetic", "NPDR")):
    d = d[d.group.isin(groups)].copy()
    ref = d[d.group == "Control"].groupby("region")[col].mean()
    d["centred"] = d[col] - d.region.map(ref)
    out = d.groupby(["donor", "group"], as_index=False).agg(centred=("centred", "mean"), severity=("severity", "first"))
    out["stage"] = out.group.map(STAGE_LABEL)
    return out


def strip(ax, dv, ylabel, title, ptext, step):
    for k, s in enumerate(stage_order()):
        v = dv[dv.stage == s].centred.values
        ax.boxplot(v, positions=[k], widths=0.5, showfliers=False, medianprops=dict(color="black", lw=1),
                   boxprops=dict(lw=0.6), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
        ax.scatter(k + np.random.default_rng(k).uniform(-0.15, 0.15, len(v)), v, s=11, color=STAGE_COLORS[s],
                   edgecolor="white", lw=0.3, zorder=3)
    top, lo = dv.centred.max(), dv.centred.min()
    for k, txt in enumerate(ptext, start=1):
        ax.text(k, top + step * (0.3 + (k - 1)), txt, ha="center", va="bottom")
    ax.set_ylim(lo - 0.3 * step, top + 2.6 * step); ax.axhline(0, color="#CCCCCC", lw=0.6, zorder=0)
    ax.set_xticks(range(3)); ax.set_xticklabels(["Ctrl", "No\nDR", "NPDR"])
    ax.set_ylabel(ylabel); ax.set_title(title, loc="left")


def forest(ax, df, labels, title, xlabel):
    for i, lab in enumerate(labels):
        for j, (cf, cl) in enumerate(CONTRASTS):
            r = df[(df.label == lab) & (df.coef == cf)].iloc[0]; y = i + (j - 0.5) * 0.34
            ax.errorbar(r["diff"], y, xerr=[[r["diff"] - r.ci_low], [r.ci_high - r["diff"]]], fmt="o", ms=3.2,
                        color=C_COL[j], lw=0.8, capsize=1.2)
    ax.axvline(0, color="black", lw=0.6)
    ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels); ax.invert_yaxis()
    ax.set_xlabel(xlabel); ax.set_title(title, loc="left")


def panel_a(ax):
    s = pd.read_csv(R / "19_state_support" / "mural_measures_by_sample.tsv", sep="\t")
    t = pd.read_csv(R / "18a_mural_state" / "il1rl1_tests.tsv", sep="\t"); t = t[t.measure == "mural_content"]
    dv = by_donor(s, "mural_content")
    strip(ax, dv, "Mural content vs control (log$_2$)", "Mural-cell content", [fmt_p(t[t.coef == cf].p_perm.iloc[0]) for cf, _ in CONTRASTS], 0.22)
    return dv


def panel_b(ax):
    t = pd.read_csv(R / "19_state_support" / "mural_index_leave_marker_out.tsv", sep="\t")
    lab = {"all eight": "All eight markers", "without IL1RL1": "Without IL1RL1", "without ADAMTS4": "Without ADAMTS4",
           "without IL1RL1 and ADAMTS4": "Without both", "median of mural-enriched genes": "Median of 52 genes"}
    others = [m for m in t.measure.unique() if m.startswith("without ") and m not in lab]
    t["label"] = t.measure.map(lambda m: lab.get(m, "Without one other marker"))
    keep = t[t.measure.isin(lab)].copy()
    forest(ax, keep, list(lab.values()), "Mural index, markers left out", "Difference vs control (log$_2$, 95% CI)")
    ot = t[t.measure.isin(others)]
    for j, (cf, _) in enumerate(CONTRASTS):                      # 其余 6 个标记逐个留出：范围
        v = ot[ot.coef == cf]["diff"]
        ax.plot([v.min(), v.max()], [0 + (j - 0.5) * 0.34] * 2, color=C_COL[j], lw=3.2, alpha=0.35, solid_capstyle="butt", zorder=1)
    ax.set_xlim(-1.3, 0.35)
    return t[["measure", "coef", "diff", "ci_low", "ci_high", "p", "p_perm"]]


def panel_c(ax):
    g = pd.read_csv(R / "18a_mural_state" / "mural_gene_state.tsv", sep="\t")
    genes = ["IL1RL1", "ADAMTS4", "EDNRA", "GJA4", "COL1A1", "COL1A2", "COL3A1", "COL6A3"]
    d = g[g.measure.isin(genes)].copy(); d["label"] = d.measure
    forest(ax, d, genes, "Genes relative to mural content", "Difference vs control (log$_2$, 95% CI)")
    ax.set_yticklabels(genes, fontstyle="italic"); ax.axhline(3.5, color="#DDDDDD", lw=0.6)
    return d[["measure", "coef", "diff", "ci_low", "ci_high", "p", "p_perm", "fdr_perm", "share_pericyte", "share_smc"]]


def panel_comp(ax):
    t = pd.read_csv(R / "20_composition_check" / "composition_tests.tsv", sep="\t")
    t = t[t.model == "group + covariates"]
    lab = {"pericyte_only_top8_without_IL1RL1_ADAMTS4": "Pericyte, top 8", "pericyte_strict_without_IL1RL1_ADAMTS4": "Pericyte, strict",
           "pericyte_relaxed_without_IL1RL1_ADAMTS4": "Pericyte, relaxed", "smc_only_top8_without_IL1RL1_ADAMTS4": "Smooth muscle, top 8",
           "smc_strict_without_IL1RL1_ADAMTS4": "Smooth muscle, strict", "smc_relaxed_without_IL1RL1_ADAMTS4": "Smooth muscle, relaxed"}
    d = t[t.measure.isin(lab)].copy(); d["label"] = d.measure.map(lab)
    forest(ax, d, list(lab.values()), "Marker indices without IL1RL1", "Difference vs control (log$_2$, 95% CI)")
    ax.axhline(2.5, color="#DDDDDD", lw=0.6)
    return d


def panel_il_peri(ax):
    s = pd.read_csv(R / "20_composition_check" / "measures_by_sample.tsv", sep="\t")
    t = pd.read_csv(R / "20_composition_check" / "composition_tests.tsv", sep="\t")
    t = t[t.measure == "IL1RL1_rel_pericyte_top8"]
    dv = by_donor(s, "il1rl1_rel_pericyte")
    strip(ax, dv, "IL1RL1 / pericyte markers\nvs control (log$_2$)", "IL1RL1 within pericytes", [fmt_p(t[t.coef == cf].p_perm.iloc[0]) for cf, _ in CONTRASTS], 0.7)
    return dv


def panel_e(ax):
    d = pd.read_csv(R / "20_composition_check" / "il1rl1_collagen_residuals.tsv", sep="\t")
    a = pd.read_csv(R / "20_composition_check" / "il1rl1_collagen_composition.tsv", sep="\t")
    a = a[a.adjusted_for == "group, composition index, age, sex, RIN, PMI"].iloc[0]
    for g, lab in STAGE_LABEL.items():
        x = d[d.group == g]
        ax.scatter(x.il_resid, x.col_resid, s=12, color=STAGE_COLORS[lab], edgecolor="white", lw=0.3, zorder=3)
    b = np.polyfit(d.il_resid, d.col_resid, 1); xs = np.linspace(d.il_resid.min(), d.il_resid.max(), 10)
    ax.plot(xs, np.polyval(b, xs), color="black", lw=0.8, zorder=2)
    lo, hi = ax.get_ylim(); ax.set_ylim(lo, hi + 0.3 * (hi - lo))
    rho = f"ρ = −{abs(a.rho):.2f}" if a.rho < 0 else f"ρ = {a.rho:.2f}"
    ax.text(0.98, 0.97, rho + f"\nP < 0.001, {int(a.n_donors)} donors", transform=ax.transAxes, ha="right", va="top")
    ax.set_xlabel("IL1RL1, donor residual (log$_2$)"); ax.set_ylabel("Collagen genes, donor residual (log$_2$)")
    ax.set_title("IL1RL1 and collagen within groups", loc="left")
    return d


def main() -> None:
    apply_style()
    fig = plt.figure(figsize=(FULL_W, 152 * MM))
    r1 = fig.add_gridspec(1, 3, width_ratios=[0.8, 1.15, 1.15], wspace=1.15, top=0.95, bottom=0.6)
    r2 = fig.add_gridspec(1, 3, width_ratios=[1.15, 0.85, 1.1], wspace=0.95, top=0.43, bottom=0.11)
    for spec, fn, letter in [(r1[0, 0], panel_a, "A"), (r1[0, 1], panel_b, "B"), (r1[0, 2], panel_comp, "C"),
                             (r2[0, 0], panel_c, "D"), (r2[0, 1], panel_il_peri, "E"), (r2[0, 2], panel_e, "F")]:
        ax = fig.add_subplot(spec)
        save_source(fn(ax), "Fig3", letter); panel_label(ax, letter)
    labs = ["Control", "Diabetes without DR", "NPDR", "PDR"]
    hs = [plt.Line2D([], [], marker="o", ls="", color=STAGE_COLORS[k if k != "Diabetes without DR" else "Diabetes without retinopathy"], ms=4.5) for k in labs]
    fig.legend(hs, labs, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.0))
    save_fig(fig, "Figure3")


if __name__ == "__main__":
    main()
