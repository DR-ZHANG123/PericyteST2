"""Supplementary Figures S1–S3 与 S6（S4 组成、S5 壁细胞稳健性、S7 基础模型见各自脚本）。只读 results/ 表。"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from figstyle import CLASS_LABEL, FULL_W, MM, ROOT, STAGE_COLORS, apply_style, panel_label, save_fig, save_source
from Fig2_deconvolution import METHODS, ORDER, heat

R = ROOT / "results"


def s1():
    """S1：随机组成测试下的一致性（留出供体）。"""
    bm = pd.read_csv(R / "09_deconv_benchmark" / "benchmark_metrics.tsv", sep="\t")
    d = bm[bm.test == "random"]
    M = np.array([[d[(d.method == m) & (d.training == t) & (d.cell_class == c)].ccc.mean() for (m, t, _) in METHODS] for c in ORDER])
    fig, ax = plt.subplots(figsize=(110 * MM, 105 * MM))
    heat(ax, M, [lab for *_, lab in METHODS], [CLASS_LABEL[c] for c in ORDER], 0, 1, "Blues", cbar_label="Concordance (CCC)")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    save_source(pd.DataFrame(M, index=ORDER, columns=[lab.replace("\n", " ") for *_, lab in METHODS]).reset_index(names="cell_class"), "FigS1", "A")
    save_fig(fig, "FigureS1", supplementary=True)


def s2():
    """S3：关键基因与标记基因的 bulk 信号来源（平均视网膜组成下各细胞类群的份额）。"""
    sh = pd.read_csv(R / "09c_ratio_validation" / "bulk_share.tsv", sep="\t", index_col=0)
    ms = pd.read_csv(R / "14a_marker_sets" / "marker_sets.tsv", sep="\t")
    key = ["RGS5", "NOTCH3", "KCNJ8", "PDGFRB", "DLL4", "TEK", "CDH5", "PLVAP", "RLBP1", "VIM", "GLUL", "GFAP"]
    mk = []
    for s in ["pericyte_only_top8", "smc_only_top8", "endothelial_top8", "muller_top8"]:
        mk += [g for g in ms[ms.set == s].gene if g not in mk]
    genes = [g for g in key + mk if g in sh.index]
    M = sh.loc[genes, ["mural", "endothelial", "muller"]].values
    fig, ax = plt.subplots(figsize=(85 * MM, 190 * MM))
    heat(ax, M, ["Mural", "Endothelial", "Müller"], genes, 0, 1, "Greens", cbar_label="Share of bulk signal")
    ax.set_yticklabels(genes, fontstyle="italic")
    ax.axhline(len(key) - 0.5, color="black", lw=0.8)
    save_source(sh.loc[genes, ["mural", "endothelial", "muller"]].reset_index(names="gene"), "FigS3", "A")
    save_fig(fig, "FigureS3", supplementary=True)


def s3():
    """S6：各细胞类群指数在黄斑、周边分别拟合的效应（95% CI）。"""
    reg = pd.read_csv(R / "14b_robustness" / "index_by_region.tsv", sep="\t")
    sets = [("mural_top8", "Mural"), ("pericyte_only_top8", "Pericytes"), ("smc_only_top8", "Smooth muscle"),
            ("endothelial_top8", "Endothelium"), ("muller_top8", "Müller glia"), ("ratio_RLBP1_muller", "RLBP1 ratio"),
            ("ratio_VIM_muller", "VIM ratio"), ("ratio_GFAP_muller", "GFAP ratio")]
    fig, axes = plt.subplots(1, 2, figsize=(FULL_W, 90 * MM), sharey=True)
    rows = []
    for a_i, (cf, title) in enumerate([("groupDiabetic", "Diabetes without DR vs control"), ("groupNPDR", "NPDR vs control")]):
        ax = axes[a_i]
        for i, (k, lab) in enumerate(sets):
            for j, (rg, mk) in enumerate([("macula", "o"), ("periphery", "s")]):
                r = reg[(reg["index"] == k) & (reg.coef == cf) & (reg.region == rg)].iloc[0]
                y = i + (j - 0.5) * 0.3
                ax.errorbar(r["diff"], y, xerr=[[r["diff"] - r.ci_low], [r.ci_high - r["diff"]]], fmt=mk, ms=3.5,
                            color=["#4C72B0", "#DD8452"][j], lw=0.8, capsize=1.2)
                rows.append(dict(contrast=title, index=lab, region=rg, diff=r["diff"], ci_low=r.ci_low, ci_high=r.ci_high, p=r.p))
        ax.axvline(0, color="black", lw=0.6); ax.set_title(title, loc="left"); ax.set_xlabel("Difference vs control (log$_2$, 95% CI)")
        panel_label(ax, "AB"[a_i])
    axes[0].set_yticks(range(len(sets))); axes[0].set_yticklabels([s[1] for s in sets]); axes[0].invert_yaxis()
    axes[1].errorbar([], [], fmt="o", color="#4C72B0", label="Macula"); axes[1].errorbar([], [], fmt="s", color="#DD8452", label="Periphery")
    fig.legend(*axes[1].get_legend_handles_labels(), loc="upper center", bbox_to_anchor=(0.5, -0.01), ncol=2, borderpad=0.3)
    save_source(pd.DataFrame(rows), "FigS6", "AB")
    save_fig(fig, "FigureS6", supplementary=True)


def s4():
    """S2：TAPE 逐样本细胞类型表达的病程检验——观察到的非零比例与供体间置换零分布。"""
    cs = pd.read_csv(R / "11_celltype_stage_tests" / "ctexpr_summary.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(FULL_W, 80 * MM), sharey=True)
    for a_i, (cf, title) in enumerate([("groupDiabetic", "Diabetes without DR vs control"), ("groupNPDR", "NPDR vs control")]):
        ax = axes[a_i]; d = cs[cs.contrast == cf].sort_values("cell_class")
        y = np.arange(len(d))
        ax.barh(y, d.null_pi1_q95, color="#E3E3E3", height=0.7, label="95th percentile, permuted labels")
        ax.scatter(d.pi1, y, color="#C44E52", s=16, zorder=3, label="Observed")
        ax.set_yticks(y); ax.set_yticklabels([CLASS_LABEL[c] for c in d.cell_class]); ax.invert_yaxis()
        ax.set_xlabel("Estimated proportion of non-null genes"); ax.set_title(title, loc="left"); ax.set_xlim(0, 1.0)
        panel_label(ax, "AB"[a_i])
    fig.legend(*axes[1].get_legend_handles_labels(), loc="upper center", bbox_to_anchor=(0.5, -0.01), ncol=2, borderpad=0.3)
    save_source(cs, "FigS2", "AB")
    save_fig(fig, "FigureS2", supplementary=True)


if __name__ == "__main__":
    apply_style()
    for f in (s1, s2, s3, s4):
        f()
