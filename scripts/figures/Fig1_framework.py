"""Figure 1：研究框架示意图（矢量）。

输入 figures/source_data/Figure1_vector_original.svg（四个插图为矢量路径，见 Figure1_icons/），
统一成与其他图一致的样式后输出 figures/main/Figure1.pdf/.png：
面板字母 10 pt 黑色粗体、位于该行最高内容之上并与标题左缘对齐；全图字体 Arial；英式拼写；图内脚注移入图注。
画布 1920 px 对应 180 mm（1 px = 0.2657 pt）。
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "figures" / "source_data" / "Figure1_vector_original.svg"
OUT = ROOT / "figures" / "main"
PT_PER_PX = 180 / 25.4 * 72 / 1920
LETTER_PX = round(10 / PT_PER_PX, 1)          # 10 pt，与 figstyle.align_panel_labels 一致
TITLE_PX = round(8.5 / PT_PER_PX, 1)
FONT = "Arial, Liberation Sans, Helvetica, sans-serif"
PANEL_X = {"A": 60, "B": 650, "C": 1370}
HEIGHT = 1012                                  # 删去脚注后裁去底部空白


def _t(id_, x, y, size, text, weight=400, fill="#334B5B", anchor="start", style="normal"):
    return (f'<text id="{id_}" x="{x}" y="{y}" font-family="{FONT}" font-size="{size}" font-weight="{weight}" '
            f'font-style="{style}" text-anchor="{anchor}" fill="{fill}">{text}</text>\n')


def _arrow(id_, x0, x1, y, col="#334B5B"):
    return (f'<path id="{id_}-shaft" d="M{x0} {y} L{x1 - 6} {y}" fill="none" stroke="{col}" stroke-width="3" stroke-linecap="round"/>\n'
            f'<path id="{id_}-head" d="M{x1} {y} L{x1 - 13} {y + 6.8} L{x1 - 13} {y - 6.8} Z" fill="{col}"/>\n')


def panel_c(svg: str) -> str:
    """面板 C：由“待做实验”改为人视网膜周细胞实验的结论（Figure 5）；虚线改实线。"""
    red, teal, grey = "#A8574C", "#387E83", "#647987"
    rep = {"Experimental hypothesis": "Experiments in pericytes", "Proposed experiments,": "25 mM glucose for 6 days,",
           "not yet performed": "vs 5.5 mM glucose or mannitol", "Does PDGF-BB regulate": "High glucose blunts",
           "ST2 isoforms under": "PDGF-BB signalling", "diabetic stress?": "and sST2 induction",
           "Control and high-glucose conditions": "Normal or high glucose, mannitol"}
    for a, b in rep.items():
        svg = svg.replace(f">{a}<", f">{b}<")
    svg = re.sub(r'(id="(hypothesis-boundary|findings-to-hypothesis-shaft)"[^>]*?) stroke-dasharray="[^"]*"', r"\1", svg)
    def varrow(id_, x, y0, y1, col):     # 竖直箭头，y0 → y1
        d = 1 if y1 > y0 else -1
        return (f'<path id="{id_}-shaft" d="M{x} {y0} L{x} {y1 - d * 8}" fill="none" stroke="{col}" stroke-width="3" stroke-linecap="round"/>\n'
                f'<path id="{id_}-head" d="M{x} {y1} L{x - 6.8} {y1 - d * 13} L{x + 6.8} {y1 - d * 13} Z" fill="{col}"/>\n')
    block = (_t("pdgf-bb-ligand", 1406, 842, 26, "PDGF-BB", 700, teal) + _arrow("ligand-receptor", 1530, 1572, 833)
             + _t("pdgfr-beta", 1584, 842, 26, "PDGFRβ", 700) + _arrow("receptor-sst2", 1696, 1738, 833, red)
             + _t("sst2-induced", 1750, 842, 26, "sST2", 700, red)
             + _t("hg-blunts", 1615, 878, 24, "both steps blunted by high glucose", 400, grey, "middle")
             + _t("knockdown-gene", 1592, 930, 26, "IL1RL1", 700, "#334B5B", "end", "italic")
             + _t("knockdown-label", 1600, 930, 26, "knock-down", 700)
             + _t("collagen-effect", 1572, 972, 27, "Collagen", 700, red, "end") + varrow("collagen-up", 1590, 976, 948, red)
             + _t("adamts4-effect", 1632, 972, 27, "ADAMTS4", 700, red) + varrow("adamts4-down", 1784, 946, 974, red))
    svg = re.sub(r'<text id="pdgf-bb-ligand".*?<text id="st2-isoforms"[^>]*>[^<]*</text>\s*', block, svg, flags=re.S)
    return svg


def restyle(svg: str) -> str:
    svg = panel_c(svg)
    svg = re.sub(r'font-family="[^"]*"', f'font-family="{FONT}"', svg)
    svg = svg.replace(">signaling<", ">signalling<")
    # 图内脚注与其分隔线删去（内容写入图注，避免重复）
    svg = re.sub(r'<text id="figure-note".*?</text>\s*', "", svg, flags=re.S)
    svg = re.sub(r'<path id="footer-rule"[^>]*/>\s*', "", svg)
    for k, x in PANEL_X.items():
        svg = re.sub(rf'<text id="panel-{k}-letter"[^>]*>{k}</text>',
                     f'<text id="panel-{k}-letter" x="{x}" y="44" font-family="{FONT}" font-size="{LETTER_PX}" '
                     f'font-weight="700" text-anchor="start" fill="#000000">{k}</text>', svg)
        svg = re.sub(rf'(<text id="panel-{k}-title") x="[\d.]+" y="[\d.]+"([^>]*?)font-size="[\d.]+"',
                     rf'\1 x="{x}" y="96"\2font-size="{TITLE_PX}"', svg)
    svg = svg.replace('height="1080" viewBox="0 0 1920 1080"', f'height="{HEIGHT}" viewBox="0 0 1920 {HEIGHT}"')
    svg = svg.replace('<rect width="1920" height="1080"', f'<rect width="1920" height="{HEIGHT}"')
    return svg


def main() -> None:
    import cairosvg
    svg = restyle(SRC.read_text())
    (ROOT / "figures" / "source_data" / "Figure1_vector.svg").write_text(svg)
    OUT.mkdir(parents=True, exist_ok=True)
    w_pt = 180 / 25.4 * 72
    cairosvg.svg2pdf(bytestring=svg.encode(), write_to=str(OUT / "Figure1.pdf"), output_width=w_pt, output_height=w_pt * HEIGHT / 1920)
    cairosvg.svg2png(bytestring=svg.encode(), write_to=str(OUT / "Figure1.png"), output_width=round(180 / 25.4 * 300))


if __name__ == "__main__":
    main()
