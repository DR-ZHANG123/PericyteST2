"""Stage 18b：Geneformer 预测——删除哪个受体最能把健康周细胞推向糖尿病供体中的壁细胞状态。

1) 状态方向：把每个健康周细胞中壁细胞基因的表达按供体中“相对壁细胞含量”的改变（stage 18a）缩放，重新排名后嵌入，
   方向 D = 嵌入(改变后) − 嵌入(原始)。评价某个签名内基因 g 时，用去掉 g 的签名重算 D（避免“删除 g ≈ 签名里 g 下降”的循环）。
2) 每个被删除基因的嵌入变化（stage 17a 已存）与 D 的余弦，细胞内按排名位置分箱标准化后，先按供体再在供体间平均。
3) 特异性：把签名的改变值在壁细胞基因间随机打乱，重复得到零分布；候选受体的对齐须超过打乱签名。
用法：python scripts/stages/18b_fm_state_match.py --model V2-104M --gpu 0
"""
from __future__ import annotations

import argparse, json, sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.lib.repro import ROOT, load_config, set_global_seed  # noqa: E402
from scripts.lib import fm  # noqa: E402

SCREEN_SET = {"V2-104M": "screen", "V1-10M": "replicate", "V2-316M": "replicate"}


def directions(model, ref_sub, dicts, version, input_size, layer, device, batch, scale_sets: list[dict]) -> np.ndarray:
    """对每个缩放方案返回 细胞 × 维度 的方向（改变后 − 原始）。"""
    special = version == "V2"
    base = fm.embed(model, fm.tokenize_cells(ref_sub, dicts, version, input_size), layer, special, device, batch)
    out = []
    for sc in scale_sets:
        e = fm.embed(model, fm.tokenize_cells(ref_sub, dicts, version, input_size, scale=sc), layer, special, device, batch)
        out.append(e - base)
    return np.stack(out)                                   # 方案 × 细胞 × 维度


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True); ap.add_argument("--gpu", default="0"); ap.add_argument("--batch", type=int, default=None)
    a = ap.parse_args()
    cfg = load_config(); fc = cfg["foundation_model"]; sm = cfg["fm_state_match"]; set_global_seed(cfg["seed"])
    sub, version, input_size = fc["models"][a.model]
    batch = a.batch or fc["batch_size"]; device = f"cuda:{a.gpu}"
    scr = ROOT / "data/interim/fm" / f"{a.model}_{SCREEN_SET[a.model]}"
    cells = json.loads((scr / "DONE").read_text())["cells"]["pericyte"]
    out_dir = ROOT / "data/interim/fm" / f"{a.model}_state"; out_dir.mkdir(parents=True, exist_ok=True)

    ref = ad.read_h5ad(ROOT / cfg["reference"]["out_h5ad"])
    sub_ad = ref[cells].copy()
    sym2ens = dict(zip(ref.var["feature_name"].astype(str), ref.var["feature_id"].astype(str).str.split(".").str[0]))
    sig = pd.read_csv(ROOT / "results/18a_mural_state/state_signature.tsv", sep="\t")
    sig["ens"] = sig["gene"].map(sym2ens); sig = sig.dropna(subset=["ens"])
    rng = np.random.default_rng(cfg["seed"])

    schemes, labels = [], []
    for col in ["lfc_npdr", "lfc_diabetic"]:
        full = dict(zip(sig["ens"], 2.0 ** sig[col]))
        schemes.append(full); labels.append((col, "full", ""))
        for g in sm["leave_out"]:                          # 签名内的候选基因：去掉自身后的签名
            e = sym2ens.get(g)
            if e in full:
                schemes.append({k: v for k, v in full.items() if k != e}); labels.append((col, "minus", g))
        for b in range(sm["n_shuffle"]):                   # 打乱签名：改变值在壁细胞基因间随机重排
            v = rng.permutation(sig[col].values)
            schemes.append(dict(zip(sig["ens"], 2.0 ** v))); labels.append((col, "shuffle", str(b)))

    dicts = fm.dictionaries(ROOT / fc["model_dir"], version)
    model = fm.load_model(ROOT / fc["model_dir"] / sub, device)
    D = directions(model, sub_ad, dicts, version, input_size, fc["emb_layer"], device, batch, schemes)
    np.save(out_dir / "directions.npy", D.astype(np.float32))
    pd.DataFrame(labels, columns=["signature", "kind", "tag"]).to_csv(out_dir / "directions.tsv", sep="\t", index=False)
    (out_dir / "DONE").write_text(json.dumps({"model": a.model, "cells": cells, "n_schemes": len(schemes)}))
    print(f"[{a.model}] directions {D.shape}", flush=True)


if __name__ == "__main__":
    main()
