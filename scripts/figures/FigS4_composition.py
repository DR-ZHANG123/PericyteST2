"""Supplementary Figure S4：人供体视网膜的细胞组成（TAPE 与 NNLS）。Müller 比例在 Figure 5A。只读 results/ 表。"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from figstyle import (CLASS_LABEL, FULL_W, MM, ROOT, STAGE_COLORS, STAGE_LABEL, apply_style, panel_label, save_fig,
                      save_source, stage_order, num)

STABLE_R = 0.7      # 5 个种子间平均相关 ≥ 0.7 的细胞类型才在正文图中展示


def donor_values(fr: pd.DataFrame, meta: pd.DataFrame, cls: str, region: str | None = None) -> pd.DataFrame:
    """每份样本 logit(比例) 减去同部位对照均值，再按供体平均。"""
    d = meta.merge(fr[["sample", cls]], on="sample")
    d = d[d.group.isin(["Control", "Diabetic", "NPDR"])]
    if region:
        d = d[d.region == region]
    p = d[cls].clip(1e-4, 1 - 1e-4); d["logit"] = np.log(p / (1 - p))
    d["centred"] = d["logit"] - d.groupby("region")["logit"].transform(lambda x: x[d.loc[x.index, "group"] == "Control"].mean())
    out = d.groupby(["donor", "group"], as_index=False)["centred"].mean()
    out["stage"] = out.group.map(STAGE_LABEL)
    return out


def strip(ax, df, ylabel, title):
    order = stage_order()
    for k, s in enumerate(order):
        v = df[df.stage == s].centred.values
        ax.boxplot(v, positions=[k], widths=0.5, showfliers=False, medianprops=dict(color="black", lw=1),
                   boxprops=dict(lw=0.6), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
        jit = np.random.default_rng(k).uniform(-0.15, 0.15, len(v))
        ax.scatter(k + jit, v, s=12, color=STAGE_COLORS[s], edgecolor="white", lw=0.3, zorder=3)
    ax.set_xticks(range(3)); ax.set_xticklabels(["Control", "Diabetes,\nno DR", "NPDR"])
    ax.set_ylabel(ylabel); ax.set_title(title, loc="left")
    ax.axhline(0, color="#CCCCCC", lw=0.6, zorder=0)


def main() -> None:
    apply_style()
    res = ROOT / "results"
    ft = pd.read_csv(res / "11_celltype_stage_tests" / "fraction_tests.tsv", sep="\t")
    sc = pd.read_csv(res / "10_deconv_retina" / "seed_consistency.tsv", sep="\t")
    fr = pd.read_csv(res / "10_deconv_retina" / "fractions_tape.tsv", sep="\t").rename(columns={"Unnamed: 0": "sample"})
    if "sample" not in fr.columns:
        fr = fr.rename(columns={fr.columns[0]: "sample"})
    meta = pd.read_csv(res / "06_retina_pilot" / "sample_meta.tsv", sep="\t")
    stable = sc[sc.mean_pairwise_r >= STABLE_R].cell_class.tolist()
    order = [c for c in ["rod", "bipolar", "amacrine", "ganglion", "muller", "astrocyte", "microglia", "rpe"] if c in stable]

    fig = plt.figure(figsize=(FULL_W, 150 * MM))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], width_ratios=[1.25, 1.25, 1], hspace=0.62, wspace=0.7)

    # A：部位阳性对照（周边 vs 黄斑的 t 值）
    ax = fig.add_subplot(gs[0, 0])
    cls_a = ["ganglion", "cone", "amacrine", "rod", "muller"]
    rows = []
    for k, c in enumerate(cls_a):
        for j, (m, col) in enumerate([("tape", "#4C72B0"), ("nnls", "#8C8C8C")]):
            t = ft[(ft.method == m) & (ft.contrast == "regionperiphery") & (ft.cell_class == c)].t.iloc[0]
            ax.barh(k + (j - 0.5) * 0.36, t, height=0.34, color=col, edgecolor="white", lw=0.3)
            rows.append(dict(cell_class=c, method=m.upper(), t_periphery_vs_macula=t))
    ax.set_yticks(range(len(cls_a))); ax.set_yticklabels([CLASS_LABEL[c] for c in cls_a]); ax.invert_yaxis()
    ax.axvline(0, color="black", lw=0.6); ax.set_xlabel("t, periphery vs macula")
    ax.set_title("Anatomy recovered", loc="left")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color="#4C72B0", label="TAPE"), Patch(color="#8C8C8C", label="NNLS")], loc="lower left", borderpad=0.3)
    save_source(pd.DataFrame(rows), "FigS4", "A")
    panel_label(ax, "A")

    # B：各阶段相对对照的 t 值（稳定可估的细胞类型）
    ax = fig.add_subplot(gs[0, 1:])
    ax.set_anchor("E")
    cols = [("tape", "groupDiabetic"), ("tape", "groupNPDR"), ("nnls", "groupDiabetic"), ("nnls", "groupNPDR")]
    M = np.array([[ft[(ft.method == m) & (ft.contrast == cf) & (ft.cell_class == c)].t.iloc[0] for (m, cf) in cols] for c in order])
    im = ax.imshow(M.T, cmap="RdBu_r", vmin=-3.5, vmax=3.5, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            t = ax.text(i, j, num(M[i, j]), ha="center", va="center", fontsize=8, color="white" if abs(M[i, j]) > 2.4 else "black")
            t._exempt = True
    ax.set_xticks(range(len(order))); ax.set_xticklabels([CLASS_LABEL[c] for c in order], rotation=45, ha="right", rotation_mode="anchor")
    ax.set_yticks(range(4)); ax.set_yticklabels(["TAPE\ndiabetes, no DR", "TAPE\nNPDR", "NNLS\ndiabetes, no DR", "NNLS\nNPDR"])
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02); cb.set_label("t vs control"); cb.outline.set_linewidth(0.4)
    ax.set_title("Change in cell fraction with diabetes", loc="left")
    save_source(pd.DataFrame(M, index=order, columns=[f"{m}_{cf}" for m, cf in cols]).reset_index(names="cell_class"), "FigS4", "B")
    panel_label(ax, "B")

    # C：节细胞比例（黄斑，供体层面）
    ax = fig.add_subplot(gs[1, 0])
    dv = donor_values(fr, meta, "ganglion", "macula")
    strip(ax, dv, "Logit fraction vs control", "Ganglion cells (macula)")
    save_source(dv, "FigS4", "C")
    panel_label(ax, "C")

    # E：TAPE 5 个种子间一致性（决定哪些细胞类型可解读）
    ax = fig.add_subplot(gs[1, 1:])
    s2 = sc.sort_values("mean_pairwise_r")
    ax.barh(range(len(s2)), s2.mean_pairwise_r, color=["#4C72B0" if r >= STABLE_R else "#CCCCCC" for r in s2.mean_pairwise_r], height=0.7)
    ax.set_yticks(range(len(s2))); ax.set_yticklabels([CLASS_LABEL[c] for c in s2.cell_class])
    ax.axvline(STABLE_R, color="#C44E52", lw=0.6, ls="--"); ax.set_xlim(-0.1, 1.0)
    ax.set_xlabel("Agreement across 5 training runs (r)"); ax.set_title("Which fractions are stable", loc="left")
    save_source(s2, "FigS4", "D")
    panel_label(ax, "D")
    save_fig(fig, "FigureS4", supplementary=True)


if __name__ == "__main__":
    main()
