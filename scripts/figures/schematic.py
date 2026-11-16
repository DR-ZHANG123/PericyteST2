"""示意图工具：读取矢量化图标（vector/*.svg，只含 M/L/C/Z 绝对坐标路径 + 填充色 + 平移）并画成 matplotlib 矢量路径；
带文字的方框；检查文字是否超出所在方框。"""
from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path

import matplotlib as mpl
import numpy as np
from matplotlib.patches import FancyBboxPatch, PathPatch
from matplotlib.path import Path as MPath

ICON_DIR = Path(__file__).resolve().parent / "assets" / "ai_icons" / "vector"
_TOK = re.compile(r"[MLCZ]|-?\d*\.?\d+(?:e-?\d+)?")


@lru_cache(maxsize=None)
def load_icon(name: str):
    svg = (ICON_DIR / f"{name}.svg").read_text()
    W = float(re.search(r'width="([\d.]+)"', svg).group(1)); H = float(re.search(r'height="([\d.]+)"', svg).group(1))
    shapes = []
    for m in re.finditer(r'<path d="([^"]+)" fill="(#[0-9A-Fa-f]{6})"(?: transform="translate\(([-\d.]+),([-\d.]+)\)")?', svg):
        d, col = m.group(1), m.group(2)
        tx, ty = float(m.group(3) or 0), float(m.group(4) or 0)
        verts, codes, cmd, buf = [], [], None, []
        start = None
        for tok in _TOK.findall(d):
            if tok in "MLCZ":
                cmd = tok
                if tok == "Z" and start is not None:
                    verts.append(start); codes.append(MPath.CLOSEPOLY)
                continue
            buf.append(float(tok))
            need = {"M": 2, "L": 2, "C": 6}[cmd]
            if len(buf) == need:
                pts = [(buf[i] + tx, buf[i + 1] + ty) for i in range(0, need, 2)]
                if cmd == "M":
                    verts.append(pts[0]); codes.append(MPath.MOVETO); start = pts[0]
                elif cmd == "L":
                    verts.append(pts[0]); codes.append(MPath.LINETO)
                else:
                    verts += pts; codes += [MPath.CURVE4] * 3
                buf = []
        if verts:
            shapes.append((np.array(verts), codes, col))
    return W, H, shapes


def draw_icon(ax, name: str, cx: float, cy: float, height: float, aspect: float):
    """以 (cx, cy) 为中心、高度 height（坐标轴单位）画图标；aspect = 坐标轴单位下 x/y 的物理比例修正。"""
    W, H, shapes = load_icon(name)
    s = height / H; w = W * s / aspect
    for v, codes, col in shapes:
        x = cx - w / 2 + v[:, 0] * s / aspect
        y = cy + height / 2 - v[:, 1] * s
        ax.add_patch(PathPatch(MPath(np.c_[x, y], codes), facecolor=col, edgecolor="none", lw=0))
    return w


def box(ax, x, y, w, h, fc="#F4F6F9", ec="#B8C2CF", lw=0.7, r=0.012):
    p = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}", fc=fc, ec=ec, lw=lw)
    ax.add_patch(p); p._is_box = True
    return p


def arrow(ax, x0, y0, x1, y1, color="#5A6878", lw=1.2):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, mutation_scale=9))


def text_outside_boxes(fig, ax, pairs) -> list[str]:
    """pairs：[(text, box_patch)]；文字边框必须完全落在方框内。"""
    fig.canvas.draw(); r = fig.canvas.get_renderer(); out = []
    for t, b in pairs:
        tb = mpl.text.Text.get_window_extent(t, r); bb = b.get_window_extent(r)
        if tb.x0 < bb.x0 + 1 or tb.x1 > bb.x1 - 1 or tb.y0 < bb.y0 + 1 or tb.y1 > bb.y1 - 1:
            out.append(f"text '{t.get_text()[:30]}' leaves its box")
    return out
