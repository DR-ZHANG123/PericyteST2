"""Shared figure style: 8 pt minimum type, upper-case panel letters, 180 mm or 85 mm width."""

每个作图脚本只读 results/ 中的表，把每个面板的作图输入写到 figures/source_data/，并保存 PDF + PNG。
面板字母按各面板实际绘制范围统一放置：同一行字母在同一水平线（高于该行最高内容），同一列字母左缘对齐。
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
MM = 1 / 25.4
FULL_W = 180 * MM
STAGE_COLORS = {"Control": "#8C8C8C", "Diabetes without retinopathy": "#4C72B0", "NPDR": "#C44E52", "PDR": "#8172B3"}
STAGE_LABEL = {"Control": "Control", "Diabetic": "Diabetes without retinopathy", "NPDR": "NPDR", "PDR": "PDR"}
METHOD_COLORS = {"NNLS": "#8C8C8C", "Scaden": "#55A868", "TAPE": "#4C72B0"}
CLASS_LABEL = {"amacrine": "Amacrine", "astrocyte": "Astrocyte", "bipolar": "Bipolar", "cone": "Cone",
               "endothelial": "Endothelial", "ganglion": "Ganglion", "horizontal": "Horizontal", "microglia": "Microglia",
               "muller": "Müller glia", "pericyte": "Pericyte", "rod": "Rod", "rpe": "RPE", "smooth_muscle": "Smooth muscle"}


def apply_style(base: float = 8.0) -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
        "font.size": base, "axes.titlesize": base, "axes.labelsize": base, "xtick.labelsize": base,
        "ytick.labelsize": base, "legend.fontsize": base, "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5,
        "ytick.major.size": 2.5, "lines.linewidth": 1.0, "legend.frameon": True, "legend.fancybox": False,
        "legend.edgecolor": "#CCCCCC", "savefig.bbox": "tight", "savefig.pad_inches": 0.05, "savefig.dpi": 300,
        "pdf.fonttype": 42, "ps.fonttype": 42, "axes.titleweight": "normal"})


def panel_label(ax, letter: str) -> None:
    fig = ax.figure
    if not hasattr(fig, "_panels"):
        fig._panels = []
    fig._panels.append((ax, letter))


def _overlap(a0, a1, b0, b1) -> bool:
    return min(a1, b1) - max(a0, b0) > 0.5 * min(a1 - a0, b1 - b0)


def align_panel_labels(fig, pad_pt: float = 3.0, col_tol: float = 0.06) -> None:
    items = getattr(fig, "_panels", [])
    if not items:
        return
    fig.canvas.draw(); r = fig.canvas.get_renderer(); inv = fig.transFigure.inverted()
    rec = []
    for ax, letter in items:
        tb = inv.transform_bbox(ax.get_tightbbox(r)); pos = ax.get_position()
        rec.append(dict(letter=letter, x0=tb.x0, top=tb.y1, y0=pos.y0, y1=pos.y1, px0=pos.x0, px1=pos.x1))
    rows = []
    for it in sorted(rec, key=lambda d: -d["y1"]):
        for row in rows:
            if any(_overlap(it["y0"], it["y1"], o["y0"], o["y1"]) for o in row):
                row.append(it); break
        else:
            rows.append([it])
    pad = pad_pt / 72 / fig.get_figheight()
    for row in rows:
        y = max(o["top"] for o in row) + pad
        for o in row:
            o["y"] = y
    parent = list(range(len(rec)))

    def find(k):
        while parent[k] != k:
            parent[k] = parent[parent[k]]; k = parent[k]
        return k
    for i in range(len(rec)):
        for j in range(i + 1, len(rec)):
            a, b = rec[i], rec[j]
            inter = min(a["px1"], b["px1"]) - max(a["px0"], b["px0"])
            wide = max(a["px1"] - a["px0"], b["px1"] - b["px0"])
            if abs(a["x0"] - b["x0"]) < col_tol or (inter > 0.5 * wide and abs(a["x0"] - b["x0"]) < 1.5 * col_tol):
                parent[find(i)] = find(j)
    groups = {}
    for k, o in enumerate(rec):
        groups.setdefault(find(k), []).append(o)
    for col in groups.values():
        x = min(o["x0"] for o in col)
        for o in col:
            o["x"] = x
    for o in rec:
        fig.text(o["x"], o["y"], o["letter"], fontsize=10, fontweight="bold", ha="left", va="bottom")
    fig._panels = []


def overlap_report(fig) -> list[str]:
    """文字与文字相交、文字压住散点，逐条报告（旋转文字与热图格内数值除外）。"""
    from matplotlib.collections import PathCollection
    fig.canvas.draw(); r = fig.canvas.get_renderer()
    texts = [t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text().strip()
             and not getattr(t, "_exempt", False)]
    drawn = set()
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            for tk in axis._update_ticks():
                drawn |= {id(tk.label1), id(tk.label2)}
    ticks = {id(tk.label1) for ax in fig.axes for axis in (ax.xaxis, ax.yaxis) for tk in axis.get_major_ticks() + axis.get_minor_ticks()}
    texts = [t for t in texts if id(t) not in ticks or id(t) in drawn]
    boxes = [(t, mpl.text.Text.get_window_extent(t, r)) for t in texts]
    out = []
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            (ti, bi), (tj, bj) = boxes[i], boxes[j]
            if ti.get_rotation() or tj.get_rotation():
                continue
            if bi.width > 0 and bj.width > 0 and bi.overlaps(bj):
                out.append(f"'{ti.get_text()[:25]}' overlaps '{tj.get_text()[:25]}'")
    for col in fig.findobj(PathCollection):
        if col.axes is None or not col.get_visible() or col.get_offsets().size == 0:
            continue
        pts = col.get_offset_transform().transform(np.asarray(col.get_offsets()))
        for t, b in boxes:
            if t.axes is not col.axes:
                continue
            inside = (pts[:, 0] > b.x0 + 1) & (pts[:, 0] < b.x1 - 1) & (pts[:, 1] > b.y0 + 1) & (pts[:, 1] < b.y1 - 1)
            if inside.any():
                out.append(f"'{t.get_text()[:25]}' covers {int(inside.sum())} point(s)")
    return sorted(set(out))


def save_fig(fig, name: str, supplementary: bool = False) -> list[Path]:
    align_panel_labels(fig)
    issues = overlap_report(fig)
    for msg in issues:
        print(f"OVERLAP {name}: {msg}")
    outdir = ROOT / "figures" / ("supplementary" if supplementary else "main"); outdir.mkdir(parents=True, exist_ok=True)
    paths = []
    for ext in ("pdf", "png"):
        p = outdir / f"{name}.{ext}"; fig.savefig(p); paths.append(p)
    plt.close(fig)
    return paths


def save_source(df: pd.DataFrame, figure: str, panel: str) -> Path:
    d = ROOT / "figures" / "source_data"; d.mkdir(parents=True, exist_ok=True)
    p = d / f"{figure}_{panel}.csv"; df.to_csv(p, index=False)
    return p


def num(v: float, d: int = 1) -> str:
    """热图格内数值：用真正的减号，避免出现 “-0.0”。"""
    s = f"{v:.{d}f}"
    if float(s) == 0:
        s = f"{0:.{d}f}"
    return s.replace("-", "\u2212")


def stage_order():
    return ["Control", "Diabetes without retinopathy", "NPDR"]
