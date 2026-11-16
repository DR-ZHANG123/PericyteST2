"""Figure 2：深度学习反卷积在留出供体上的表现、对植入改变的归属，以及同类标记比值的验证。只读 results/ 表。"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Fig2_reference_panels import panel_dot, panel_umap
from figstyle import (CLASS_LABEL, FULL_W, METHOD_COLORS, MM, ROOT, apply_style, panel_label, save_fig,
                      save_source)

ORDER = ["rod", "cone", "bipolar", "amacrine", "horizontal", "ganglion", "muller", "astrocyte", "microglia",
         "endothelial", "pericyte", "smooth_muscle", "rpe"]
METHODS = [("NNLS", "none", "NNLS"), ("Scaden", "uniform", "Scaden\nuniform"), ("Scaden", "mixed", "Scaden\nretina-informed"),
           ("TAPE", "uniform", "TAPE\nuniform"), ("TAPE", "mixed", "TAPE\nretina-informed")]


def heat(ax, M, xl, yl, vmin, vmax, cmap, fmt="{:.2f}", cbar_label=None):
    im = ax.imshow(M, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(xl))); ax.set_xticklabels(xl); ax.set_yticks(range(len(yl))); ax.set_yticklabels(yl)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            t = ax.text(j, i, fmt.format(v), ha="center", va="center", fontsize=8,
                        color="white" if (v - vmin) / (vmax - vmin) > 0.6 else "black")
            t._exempt = True
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    if cbar_label:
        cb = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.set_label(cbar_label); cb.outline.set_linewidth(0.4)
    return im


def main() -> None:
    apply_style()
    res = ROOT / "results"
    bm = pd.read_csv(res / "09_deconv_benchmark" / "benchmark_metrics.tsv", sep="\t")
    sp = pd.read_csv(res / "09b_spikein" / "spikein_attribution.tsv", sep="\t")
    rv = pd.read_csv(res / "09c_ratio_validation" / "ratio_vs_bulk.tsv", sep="\t")

    fig = plt.figure(figsize=(FULL_W, 232 * MM))
    top = fig.add_gridspec(1, 3, width_ratios=[1, 1.05, 0.035], wspace=0.42, left=0.04, right=0.83, top=0.965, bottom=0.73)
    ax = fig.add_subplot(top[0, 0]); save_source(panel_umap(ax), "Fig2", "A"); panel_label(ax, "A")
    ax = fig.add_subplot(top[0, 1]); save_source(panel_dot(ax, fig.add_subplot(top[0, 2])), "Fig2", "B"); panel_label(ax, "B")
    gs = fig.add_gridspec(2, 2, width_ratios=[1.35, 1], height_ratios=[1.15, 1], hspace=0.55, wspace=0.55, top=0.62, bottom=0.05)

    # A：视网膜样组成下各方法、各细胞类型的一致性相关系数（3 个种子平均）
    ax = fig.add_subplot(gs[0, 0])
    rl = bm[bm.test == "retina_like"]
    M = np.array([[rl[(rl.method == m) & (rl.training == t) & (rl.cell_class == c)].ccc.mean() for (m, t, _) in METHODS] for c in ORDER])
    heat(ax, M, [lab for *_, lab in METHODS], [CLASS_LABEL[c] for c in ORDER], 0, 1, "Blues", cbar_label="Concordance (CCC)")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    ax.set_title("Held-out donors, retina-like composition", loc="left")
    save_source(pd.DataFrame(M, index=ORDER, columns=[lab.replace("\n", " ") for *_, lab in METHODS]).reset_index(names="cell_class"), "Fig2", "C")
    panel_label(ax, "C")

    # B：两种测试组成下的平均一致性（每个种子一个点）
    ax = fig.add_subplot(gs[0, 1])
    rows = []
    for k, (m, t, lab) in enumerate(METHODS):
        for j, test in enumerate(["random", "retina_like"]):
            d = bm[(bm.method == m) & (bm.training == t) & (bm.test == test)]
            per_seed = d.groupby(d.seed.fillna(-1)).ccc.mean()
            x = k + (j - 0.5) * 0.32
            ax.scatter(np.full(len(per_seed), x), per_seed.values, s=14, color=["#BBBBBB", "#4C72B0"][j], edgecolor="none", zorder=3)
            ax.hlines(per_seed.mean(), x - 0.12, x + 0.12, color="black", lw=1)
            rows += [dict(method=lab.replace("\n", " "), test=test, seed=s, mean_ccc=v) for s, v in per_seed.items()]
    ax.set_xticks(range(len(METHODS))); ax.set_xticklabels([lab for *_, lab in METHODS], rotation=45, ha="right", rotation_mode="anchor")
    ax.set_ylabel("Mean CCC across cell types"); ax.set_ylim(0.3, 1.0)
    ax.scatter([], [], s=14, color="#BBBBBB", label="Random composition"); ax.scatter([], [], s=14, color="#4C72B0", label="Retina-like composition")
    ax.legend(loc="lower left", handletextpad=0.2, borderpad=0.3)
    save_source(pd.DataFrame(rows), "Fig2", "D")
    panel_label(ax, "D")

    # C：植入改变的检出（AUC）：行 = TAPE 估计的各细胞类型表达及 bulk，列 = 改变植入的细胞类型
    ax = fig.add_subplot(gs[1, 0])
    est = ["pericyte", "endothelial", "muller", "bipolar", "rod", "ganglion", "bulk"]
    cols = ["pericyte", "endothelial", "muller"]
    A = np.array([[sp[(sp.estimated_in == e) & (sp.spiked_in == c)].auc.iloc[0] for c in cols] for e in est])
    heat(ax, A, ["Pericyte", "Endothelial", "Müller glia"], [CLASS_LABEL.get(e, "Bulk tissue") + (" (TAPE)" if e != "bulk" else "") for e in est],
         0.4, 1.0, "Purples", cbar_label="AUC for planted genes")
    for k, c in enumerate(cols):
        ax.add_patch(plt.Rectangle((k - 0.5, est.index(c) - 0.5), 1, 1, fill=False, ec="#C44E52", lw=1.4))
    ax.set_xlabel("Cell type carrying the planted change")
    save_source(pd.DataFrame(A, index=est, columns=cols).reset_index(names="estimated_in"), "Fig2", "E")
    panel_label(ax, "E")

    # D：细胞变少（无自身改变）时的假阳性与细胞自身改变时的检出：bulk 直接比较 vs 同类标记比值
    ax = fig.add_subplot(gs[1, 1])
    lab_g = {"mural": "Mural", "endothelial": "Endothelial", "muller": "Müller"}
    rows = []
    for k, g in enumerate(["endothelial", "muller", "mural"]):
        for j, (meth, colr) in enumerate([("bulk", "#BBBBBB"), ("ratio", "#4C72B0")]):
            fp = rv[(rv.group == g) & (rv.method == meth) & (rv.scenario == "cell_loss")].false_pos_rate.iloc[0]
            det = rv[(rv.group == g) & (rv.method == meth) & (rv.scenario == "intrinsic")].detected_fdr05.iloc[0]
            ax.bar(k + (j - 0.5) * 0.36, fp, width=0.34, color=colr, edgecolor="white", lw=0.3)
            ax.bar(k + 3.6 + (j - 0.5) * 0.36, det, width=0.34, color=colr, edgecolor="white", lw=0.3)
            rows.append(dict(group=g, method=meth, false_positive_cell_loss=fp, detected_intrinsic=det))
    ax.set_xticks([0, 1, 2, 3.6, 4.6, 5.6]); ax.set_xticklabels([lab_g[g] for g in ["endothelial", "muller", "mural"]] * 2, rotation=45, ha="right", rotation_mode="anchor")
    ax.axvline(2.8, color="#CCCCCC", lw=0.6)
    ax.text(1, 0.93, "Cell loss:\nfalse positives", ha="center", va="top")
    ax.text(4.6, 0.93, "Intrinsic change:\ndetected", ha="center", va="top")
    ax.set_ylim(0, 1.0); ax.set_ylabel("Fraction of genes (FDR < 0.05)")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color="#BBBBBB", label="Bulk expression"), Patch(color="#4C72B0", label="Ratio to own markers")],
              loc="center right", bbox_to_anchor=(1.0, 0.58), borderpad=0.3)
    save_source(pd.DataFrame(rows), "Fig2", "F")
    panel_label(ax, "F")
    save_fig(fig, "Figure2")


if __name__ == "__main__":
    main()
