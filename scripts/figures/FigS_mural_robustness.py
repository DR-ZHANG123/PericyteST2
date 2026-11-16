"""补充图：壁细胞指数减少的稳健性；donor_centred/load 也供 Figure 5 使用。只读 results/14b_robustness 与样本信息。"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from figstyle import (FULL_W, MM, ROOT, STAGE_COLORS, STAGE_LABEL, apply_style, panel_label, save_fig, save_source,
                      stage_order)

CONTRASTS = [("groupDiabetic", "Diabetes, no DR"), ("groupNPDR", "NPDR")]


def load():
    r = ROOT / "results" / "14b_robustness"
    idx = pd.read_csv(r / "indices_by_sample.tsv", sep="\t").set_index("index").T
    meta = pd.read_csv(ROOT / "results" / "06_retina_pilot" / "sample_meta.tsv", sep="\t").set_index("sample")
    return idx.join(meta, how="inner"), {k: pd.read_csv(r / f"{k}.tsv", sep="\t") for k in
                                          ["index_stage_tests", "leave_one_donor_out", "index_by_region", "dose_response"]}


def donor_centred(d: pd.DataFrame, col: str, groups=("Control", "Diabetic", "NPDR")) -> pd.DataFrame:
    """样本指数减去同部位对照均值，再按供体平均（展示用；检验见主模型）。"""
    d = d[d.group.isin(groups)].copy()
    ref = d[d.group == "Control"].groupby("region")[col].mean()
    d["centred"] = d[col] - d.region.map(ref)
    out = d.groupby(["donor", "group"], as_index=False).agg(centred=("centred", "mean"), severity=("severity", "first"))
    out["stage"] = out.group.map(STAGE_LABEL)
    return out


def main() -> None:
    """补充图：壁细胞指数的稳健性（标记数、逐一剔除供体、黄斑/周边分别拟合）。"""
    apply_style()
    d, T = load()
    st = T["index_stage_tests"]
    fig = plt.figure(figsize=(110 * MM, 95 * MM))
    gs = fig.add_gridspec(1, 1)
    # C：稳健性（壁细胞指数）：标记数、逐一剔除供体后的范围、两个部位分别拟合
    ax = fig.add_subplot(gs[0, 0])
    loo, reg = T["leave_one_donor_out"], T["index_by_region"]
    items = [("mural_top5", "5 markers"), ("mural_top8", "8 markers"), ("mural_top12", "12 markers"), ("mural_top20", "20 markers")]
    rows, ylab = [], []
    y = 0
    for cf_i, (cf, clab) in enumerate(CONTRASTS):
        col = list(STAGE_COLORS.values())[cf_i + 1]
        for k, lab in items:
            r = st[(st["index"] == k) & (st.coef == cf)].iloc[0]
            ax.errorbar(r["diff"], y, xerr=[[r["diff"] - r.ci_low], [r.ci_high - r["diff"]]], fmt="o", ms=3.5, color=col, lw=0.8, capsize=1.2)
            rows.append(dict(contrast=clab, check=lab, diff=r["diff"], low=r.ci_low, high=r.ci_high)); ylab.append(lab); y += 1
        l = loo[(loo["index"] == "mural_top8") & (loo.coef == cf)].iloc[0]
        ax.plot([l.min_diff, l.max_diff], [y, y], color=col, lw=2.2, solid_capstyle="butt")
        rows.append(dict(contrast=clab, check="Leave one donor out (range)", diff=np.nan, low=l.min_diff, high=l.max_diff)); ylab.append("Leave one donor out"); y += 1
        for rg in ["macula", "periphery"]:
            r = reg[(reg["index"] == "mural_top8") & (reg.coef == cf) & (reg.region == rg)].iloc[0]
            ax.errorbar(r["diff"], y, xerr=[[r["diff"] - r.ci_low], [r.ci_high - r["diff"]]], fmt="s", ms=3.5, color=col, lw=0.8, capsize=1.2)
            rows.append(dict(contrast=clab, check=f"{rg.capitalize()} only", diff=r["diff"], low=r.ci_low, high=r.ci_high)); ylab.append(f"{rg.capitalize()} only"); y += 1
        y += 0.6
    ys = [i for i, _ in enumerate(ylab)]
    pos = []; yy = 0
    for cf_i in range(2):
        n = len(items) + 3
        pos += [yy + i for i in range(n)]; yy += n + 0.6
    ax.set_yticks(pos); ax.set_yticklabels(ylab); ax.invert_yaxis(); ax.axvline(0, color="black", lw=0.6)
    ax.set_xlabel("Difference vs control (log$_2$)"); ax.set_title("Robustness of mural loss", loc="left")
    ax.text(1.02, (len(items) + 2) / 2, "Diabetes,\nno DR", transform=ax.get_yaxis_transform(), va="center", ha="left",
            color=STAGE_COLORS["Diabetes without retinopathy"])
    ax.text(1.02, len(items) + 3.6 + (len(items) + 2) / 2, "NPDR", transform=ax.get_yaxis_transform(), va="center", ha="left",
            color=STAGE_COLORS["NPDR"])
    save_source(pd.DataFrame(rows), "FigS5", "A")

    save_fig(fig, "FigureS5", supplementary=True)


if __name__ == "__main__":
    main()
