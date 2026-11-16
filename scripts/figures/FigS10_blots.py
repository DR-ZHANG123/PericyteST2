"""补充图 S10：Figure 5A 的未裁剪全膜（p-PDGFRβ、PDGFRβ、GAPDH），实验室原始图像，只读 results/23_wetlab/blots/。"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

from figstyle import FULL_W, MM, ROOT, apply_style, panel_label, save_fig, save_source

R = ROOT / "results" / "23_wetlab" / "blots"
BLOTS = [("p_PDGFRB", "p-PDGFRβ"), ("PDGFRB", "PDGFRβ"), ("GAPDH", "GAPDH")]


def main() -> None:
    apply_style()
    fig = plt.figure(figsize=(FULL_W, 165 * MM))
    gs = fig.add_gridspec(3, 1, hspace=0.12, left=0.02, right=0.98, top=0.95, bottom=0.02)
    rows = []
    for k, (key, lab) in enumerate(BLOTS):
        ax = fig.add_subplot(gs[k, 0]); p = R / f"full_{key}.png"
        ax.imshow(np.asarray(Image.open(p))); ax.set_xticks([]); ax.set_yticks([]); ax.axis("off")
        panel_label(ax, "ABC"[k]); rows.append(dict(panel="ABC"[k], blot=lab, image=str(p.relative_to(ROOT))))
    save_source(pd.DataFrame(rows), "FigS10", "ABC")
    save_fig(fig, "FigureS10", supplementary=True)


if __name__ == "__main__":
    main()
