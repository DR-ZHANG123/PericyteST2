"""Supplementary Figure S8：Müller 胶质早期状态改变、随等级加深的胶质增生与视杆减少，以及体外单一刺激不能重现供体改变。"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from figstyle import FULL_W, MM, ROOT, STAGE_COLORS, apply_style, panel_label, save_fig, save_source, stage_order, num
from FigS_mural_robustness import donor_centred, load
from FigS4_composition import donor_values, strip

LEVELS = ["Control", "10", "20", "35", "43", "53", "65"]


def level_color(lv):
    return (STAGE_COLORS["Control"] if lv == "Control" else STAGE_COLORS["Diabetes without retinopathy"] if lv == "10"
            else STAGE_COLORS["PDR"] if lv == "65" else STAGE_COLORS["NPDR"])


def by_level(ax, d, col, ylabel, title, p):
    dv = donor_centred(d, col, groups=("Control", "Diabetic", "NPDR", "PDR"))
    dv["level"] = np.where(dv.group == "Control", "Control", dv.severity.fillna(-1).astype(int).astype(str))
    for k, lv in enumerate(LEVELS):
        v = dv[dv.level == lv].centred.values
        ax.scatter(k + np.random.default_rng(k + 7).uniform(-0.12, 0.12, len(v)), v, s=14, color=level_color(lv),
                   edgecolor="white", lw=0.3, zorder=3)
        if len(v):
            ax.hlines(np.mean(v), k - 0.25, k + 0.25, color="black", lw=1)
    ax.set_xticks(range(len(LEVELS))); ax.set_xticklabels(["Ctrl"] + LEVELS[1:]); ax.set_xlabel("ETDRS severity level")
    ax.axhline(0, color="#CCCCCC", lw=0.6, zorder=0); ax.set_ylabel(ylabel); ax.set_title(title, loc="left")
    ax.text(0.98, 0.97, f"Trend: P = {p:.3f}", transform=ax.transAxes, ha="right", va="top")
    lo, hi = ax.get_ylim(); ax.set_ylim(lo, hi + 0.18 * (hi - lo))
    return dv


def main() -> None:
    apply_style()
    d, T = load()
    st, dr = T["index_stage_tests"], T["dose_response"]
    fig = plt.figure(figsize=(FULL_W, 160 * MM))
    top = fig.add_gridspec(1, 2, width_ratios=[1, 1.5], wspace=0.45, top=0.95, bottom=0.6)
    bot = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.15], wspace=0.6, top=0.43, bottom=0.07)

    # A：TAPE 估计的 Müller 比例（供体层面）
    ax = fig.add_subplot(top[0, 0])
    res = ROOT / "results"
    fr = pd.read_csv(res / "10_deconv_retina" / "fractions_tape.tsv", sep="\t")
    fr = fr.rename(columns={fr.columns[0]: "sample"}) if "sample" not in fr.columns else fr
    meta = pd.read_csv(res / "06_retina_pilot" / "sample_meta.tsv", sep="\t")
    dv = donor_values(fr, meta, "muller")
    strip(ax, dv, "Logit fraction vs control", "Müller glial fraction (TAPE)")
    ax.set_xticklabels(["Ctrl", "No DR", "NPDR"])
    save_source(dv, "FigS8", "A"); panel_label(ax, "A")

    # B：Müller 状态（相对 Müller 稳定标记的比值）
    ax = fig.add_subplot(top[0, 1])
    rows = []
    top = max(donor_centred(d, c).centred.max() for c in ["ratio_RLBP1_muller", "ratio_VIM_muller"])
    for g_i, (col, lab) in enumerate([("ratio_RLBP1_muller", "RLBP1 (CRALBP)"), ("ratio_VIM_muller", "VIM (vimentin)")]):
        dv = donor_centred(d, col)
        for k, s in enumerate(stage_order()):
            v = dv[dv.stage == s].centred.values; x = g_i * 3.6 + k
            ax.boxplot(v, positions=[x], widths=0.55, showfliers=False, medianprops=dict(color="black", lw=1),
                       boxprops=dict(lw=0.6), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
            ax.scatter(x + np.random.default_rng(k + 3 * g_i).uniform(-0.15, 0.15, len(v)), v, s=10, color=STAGE_COLORS[s],
                       edgecolor="white", lw=0.3, zorder=3)
        r = st[(st["index"] == col) & (st.coef == "groupDiabetic")].iloc[0]
        ax.text(g_i * 3.6 + 1, top + 0.15, f"{lab.split()[0]}\nP = {r.p_perm:.3f}", ha="center", va="bottom", fontstyle="normal")
        rows.append(dv.assign(gene=lab))
    ax.set_xticks([0, 1, 2, 3.6, 4.6, 5.6]); ax.set_xticklabels(["Ctrl", "No DR", "NPDR"] * 2)
    ax.set_ylim(-1.3, top + 0.85); ax.axhline(0, color="#CCCCCC", lw=0.6, zorder=0)
    ax.set_ylabel("Ratio to Müller markers vs control (log$_2$)"); ax.set_title("Müller glial state, before retinopathy", loc="left")
    save_source(pd.concat(rows), "FigS8", "B")
    panel_label(ax, "B")

    # C：GFAP（反应性胶质增生）随等级上升
    ax = fig.add_subplot(bot[0, 0])
    p = dr[(dr["index"] == "ratio_GFAP_muller") & (dr.predictor == "ETDRS level")].p.iloc[0]
    save_source(by_level(ax, d, "ratio_GFAP_muller", "GFAP ratio vs control (log$_2$)", "Gliosis rises with grade", p), "FigS8", "C")
    panel_label(ax, "C")

    # D：GFAP 比值与糖尿病病程（可解析为年数的糖尿病供体）
    ax = fig.add_subplot(bot[0, 1])
    import re
    dd = d[d.group.isin(["Diabetic", "NPDR", "PDR"])].copy()
    ref = d[d.group == "Control"].groupby("region")["ratio_GFAP_muller"].mean()
    dd["centred"] = dd["ratio_GFAP_muller"] - dd.region.map(ref)
    dd["years"] = dd.duration.astype(str).map(lambda x: float(re.search(r"(\d+)", x).group(1)) if re.search(r"\d", x) else np.nan)
    dv = dd.dropna(subset=["years"]).groupby(["donor", "group"], as_index=False).agg(years=("years", "first"), centred=("centred", "mean"))
    for g, col in [("Diabetic", STAGE_COLORS["Diabetes without retinopathy"]), ("NPDR", STAGE_COLORS["NPDR"]), ("PDR", STAGE_COLORS["PDR"])]:
        v = dv[dv.group == g]
        ax.scatter(v.years, v.centred, s=16, color=col, edgecolor="white", lw=0.3, zorder=3,
                   label={"Diabetic": "Diabetes, no DR"}.get(g, g))
    b = np.polyfit(dv.years, dv.centred, 1); xs = np.linspace(dv.years.min(), dv.years.max(), 50)
    ax.plot(xs, np.polyval(b, xs), color="black", lw=0.8)
    ax.axhline(0, color="#CCCCCC", lw=0.6, zorder=0)
    p = dr[(dr["index"] == "ratio_GFAP_muller") & (dr.predictor == "diabetes duration (years)")].p.iloc[0]
    ax.set_xlabel("Diabetes duration (years)"); ax.set_ylabel("GFAP ratio vs control (log$_2$)")
    ax.set_title("Gliosis and duration", loc="left")
    lo, hi = ax.get_ylim(); ax.set_ylim(lo, hi + 0.75 * (hi - lo))
    ax.text(0.02, 0.97, f"P = {p:.3f} ({dv.donor.nunique()} donors)", transform=ax.transAxes, ha="left", va="top")
    ax.legend(loc="upper right", bbox_to_anchor=(1, 0.9), borderpad=0.3, handletextpad=0.2)
    save_source(dv, "FigS8", "D")
    panel_label(ax, "D")

    # E：供体视网膜与体外单一刺激的方向比较（t 值）
    ax = fig.add_subplot(bot[0, 2])
    g6 = pd.read_csv(ROOT / "results" / "06_retina_pilot" / "gene_tests.tsv", sep="\t")
    iv = pd.read_csv(ROOT / "results" / "12_invitro" / "invitro_targets.tsv", sep="\t")
    genes = ["RGS5", "NOTCH3", "PDGFRB", "KCNJ8", "RLBP1", "VIM", "GLUL", "GFAP", "TXNIP", "VEGFA"]
    cols = [("donor", "groupDiabetic", "Donors, no DR"), ("donor", "groupNPDR", "Donors, NPDR"),
            ("iv", "pericyte_AGEs", "Pericytes + AGEs"), ("iv", "organoid_glucose", "Organoids + glucose")]
    M = np.full((len(genes), len(cols)), np.nan)
    for j, (src, key, _) in enumerate(cols):
        for i, g in enumerate(genes):
            s = g6[(g6.symbol == g) & (g6.contrast == key)] if src == "donor" else iv[(iv.gene == g) & (iv.model == key)]
            if len(s):
                M[i, j] = s.t.iloc[0]
    Mc = np.clip(M, -6, 6)
    im = ax.imshow(np.ma.masked_invalid(Mc), cmap="RdBu_r", vmin=-6, vmax=6, aspect="auto")
    ax.set_facecolor("#EEEEEE")
    for i in range(len(genes)):
        for j in range(len(cols)):
            txt = "n.d." if np.isnan(M[i, j]) else num(M[i, j])
            t = ax.text(j, i, txt, ha="center", va="center", fontsize=8, color="white" if abs(Mc[i, j]) > 4 else "black")
            t._exempt = True
    ax.set_xticks(range(len(cols))); ax.set_xticklabels([c[2] for c in cols], rotation=40, ha="right", rotation_mode="anchor"); ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontstyle="italic"); ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.set_label("t (clipped at ±6)"); cb.outline.set_linewidth(0.4)
    ax.set_title("Single exposures in vitro", loc="left")
    save_source(pd.DataFrame(M, index=genes, columns=[c[2].replace("\n", " ") for c in cols]).reset_index(names="gene"), "FigS8", "E")
    panel_label(ax, "E")
    save_fig(fig, "FigureS8", supplementary=True)


if __name__ == "__main__":
    main()
