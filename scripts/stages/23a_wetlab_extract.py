"""Stage 23a：从实验室交付的 pptx 中解析 HRPC 实验原始数值与 WB 图像（不手抄，避免录入错误）。

输入 data/raw/wetlab/wetlab_data_20261116.pptx（slide 1 WB；2 刺激 1 h qPCR；3 刺激 24 h qPCR；4 sST2 ELISA；5–6 siRNA qPCR）。
输出 results/23_wetlab/：measurements.tsv（长表：experiment, readout, group, replicate, value）与 blots/*.png。
qPCR 数值为 2^-ΔΔCt 相对表达（GAPDH 内参，以 NG 第 1 份为 1）；WB 为 p-PDGFRβ/PDGFRβ 灰度比值；ELISA 为 pg/mL。
"""
from __future__ import annotations

import io
import re
import sys
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from repro import ROOT, load_config, write_provenance  # noqa: E402

STAGE = "23_wetlab"
NS = {"a": "http://schemas.openxmlformats.org/drawingml/2006/main", "p": "http://schemas.openxmlformats.org/presentationml/2006/main"}


def tables(z: zipfile.ZipFile, slide: int) -> list[list[list[str]]]:
    root = ET.fromstring(z.read(f"ppt/slides/slide{slide}.xml"))
    out = []
    for tbl in root.iter(f"{{{NS['a']}}}tbl"):
        rows = []
        for tr in tbl.findall("a:tr", NS):
            rows.append(["".join(t.text or "" for t in tc.iter(f"{{{NS['a']}}}t")).strip() for tc in tr.findall("a:tc", NS)])
        out.append(rows)
    return out


def qpcr_rows(tab, experiment: str, groups: list[str]) -> list[dict]:
    """表头：分组 | <基因>相对表达 ...；每格为逗号分隔的重复值。"""
    head = [re.sub(r"\s*相对表达\s*", "", h).strip() for h in tab[0]]
    rec = []
    for row in tab[1:]:
        g = row[0].strip()
        if g not in groups:
            raise ValueError(f"unexpected group {g!r} in {experiment}")
        for gene, cell in zip(head[1:], row[1:]):
            for k, v in enumerate(re.split(r"[,，]\s*", cell), start=1):
                rec.append(dict(experiment=experiment, readout=gene, group=g, replicate=k, value=float(v)))
    return rec


def wb_rows(tab, groups):
    rec = []
    for row in tab[1:]:
        g = row[0].strip()
        if g not in groups:
            continue
        for k, v in enumerate(row[1:4], start=1):
            rec.append(dict(experiment="western_blot", readout="p-PDGFRβ/PDGFRβ", group=g, replicate=k, value=float(v)))
    return rec


def elisa_rows(tab, groups):
    rec = []
    for row in tab[1:]:
        m = re.match(r"(.+)-(\d+)$", row[1].strip())
        g, k = m.group(1), int(m.group(2))
        if g not in groups:
            raise ValueError(g)
        rec.append(dict(experiment="elisa_24h", readout="sST2 (pg/mL)", group=g, replicate=k, value=float(row[3])))
    return rec


def gapdh_crop(z) -> Image.Image:
    """正文 GAPDH 条带 = image3 按 slide 1 中 srcRect（千分之一百分比）裁剪。"""
    xml = z.read("ppt/slides/slide1.xml").decode()
    l, t, r, b = (int(x) / 1e5 for x in re.search(r'srcRect l="(\d+)" t="(\d+)" r="(\d+)" b="(\d+)"', xml).groups())
    im = Image.open(io.BytesIO(z.read("ppt/media/image3.png"))).convert("L")
    w, h = im.size
    return im.crop((round(l * w), round(t * h), round(w * (1 - r)), round(h * (1 - b))))


def main() -> None:
    cfg = load_config()["wetlab"]
    src = ROOT / cfg["source"]; out = ROOT / "results" / STAGE; (out / "blots").mkdir(parents=True, exist_ok=True)
    groups, sig = cfg["groups"], cfg["sirna_groups"]
    with zipfile.ZipFile(src) as z:
        rec = wb_rows(tables(z, 1)[0], groups)
        rec += qpcr_rows(tables(z, 2)[-1], "qpcr_1h", groups)
        rec += qpcr_rows(tables(z, 3)[-1], "qpcr_24h", groups)
        rec += elisa_rows(tables(z, 4)[0], groups)
        rec += qpcr_rows(tables(z, 5)[-1], "sirna_knockdown", sig)
        rec += qpcr_rows(tables(z, 6)[-1], "sirna_targets", sig)
        outs = []
        for name, f in cfg["blots"].items():
            im = gapdh_crop(z) if name == "GAPDH" else Image.open(io.BytesIO(z.read(f"ppt/media/{f}"))).convert("L")
            p = out / "blots" / f"strip_{name}.png"; im.save(p); outs.append(p)
        for name, f in cfg["blots_full"].items():
            p = out / "blots" / f"full_{name}.png"; Image.open(io.BytesIO(z.read(f"ppt/media/{f}"))).convert("RGB").save(p); outs.append(p)
    d = pd.DataFrame(rec)
    d.to_csv(out / "measurements.tsv", sep="\t", index=False)
    print(d.groupby(["experiment", "readout"]).size().to_string())
    write_provenance(STAGE + "/extract", [src], [out / "measurements.tsv", *outs])


if __name__ == "__main__":
    main()
