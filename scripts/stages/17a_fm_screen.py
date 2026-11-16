"""Stage 17a：Geneformer 全基因组 in silico 删除筛查（健康人周细胞与平滑肌细胞）。

问题：周细胞的转录组状态最依赖哪些基因（尤其是接收旁分泌信号的受体）？PDGFRB 的依赖性在周细胞与平滑肌之间是否不同？
每个细胞逐个删除其检出的每个基因，记录细胞嵌入的位移与对各类细胞中心相似度的变化。
用法：python scripts/stages/17a_fm_screen.py --model V2-104M --set screen --gpu 0
"""
from __future__ import annotations

import argparse, json, sys, time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.lib.repro import ROOT, load_config, set_global_seed  # noqa: E402
from scripts.lib import fm  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True); ap.add_argument("--set", default="screen", choices=["screen", "replicate"])
    ap.add_argument("--gpu", default="0"); ap.add_argument("--batch", type=int, default=None)   # 只影响速度/显存，不影响结果
    ap.add_argument("--shard", default="0/1")            # k/n：多块 GPU 分摊细胞，结果按细胞写文件，与分片方式无关
    a = ap.parse_args()
    cfg = load_config(); fc = cfg["foundation_model"]; seed = cfg["seed"]; set_global_seed(seed)
    if a.batch:
        fc["batch_size"] = a.batch
    sub, version, input_size = fc["models"][a.model]
    gf_root = ROOT / fc["model_dir"]
    device = f"cuda:{a.gpu}"
    out_dir = ROOT / "data/interim/fm" / f"{a.model}_{a.set}"; out_dir.mkdir(parents=True, exist_ok=True)

    ref = ad.read_h5ad(ROOT / cfg["reference"]["out_h5ad"])
    rng = np.random.default_rng(seed)                         # 同一种子：三个模型用同一批细胞（replicate 是 screen 的子集）
    classes = sorted(ref.obs["cell_class"].unique())
    cent_idx = {c: fm.pick_cells(ref, c, fc["cells_centroid"], fc["min_genes"] // 2, rng) for c in classes}
    used = np.concatenate(list(cent_idx.values()))
    n_cells = fc["cells_screen"]
    scr_idx = {c: fm.pick_cells(ref, c, n, fc["min_genes"], rng, exclude=used) for c, n in n_cells.items()}
    if a.set == "replicate":
        scr_idx = {c: scr_idx[c][rng.permutation(len(scr_idx[c]))[:n]] for c, n in fc["cells_screen_replicate"].items()}

    dicts = fm.dictionaries(gf_root, version)
    special = version == "V2"
    model = fm.load_model(gf_root / sub, device)
    cents = []
    for c in classes:
        seqs = fm.tokenize_cells(ref[cent_idx[c]], dicts, version, input_size)
        e = fm.embed(model, seqs, fc["emb_layer"], special, device, fc["batch_size"])
        cents.append((e / np.linalg.norm(e, axis=1, keepdims=True)).mean(0))
    cents = np.stack(cents)
    np.save(out_dir / "centroids.npy", cents)
    (out_dir / "classes.json").write_text(json.dumps(classes))

    inv = {v: k for k, v in dicts["token"].items()}
    for c, idx in scr_idx.items():
        seqs = fm.tokenize_cells(ref[idx], dicts, version, input_size)
        k_sh, n_sh = map(int, a.shard.split("/"))
        for j, (cell, seq) in enumerate(zip(ref.obs_names[idx], seqs)):
            f = out_dir / f"{c}__{cell}.parquet"
            if j % n_sh != k_sh or f.exists():
                continue
            t0 = time.time()
            tok, rr, shift, dcent, delta = fm.delete_screen_cell(model, seq, fc["emb_layer"], special, device, fc["batch_size"], cents)
            df = pd.DataFrame({"cell": cell, "cell_class": c, "donor": ref.obs.loc[cell, "donor_id"],
                               "gene_id": [inv[t] for t in tok], "rel_rank": rr, "n_genes": len(tok), "shift": shift})
            for k, cl in enumerate(classes):
                df[f"d_{cl}"] = dcent[:, k]
            df.to_parquet(f)
            np.save(out_dir / f"{c}__{cell}__delta.npy", delta.astype(np.float16))
            print(f"[{a.model}] {c} {j + 1}/{len(idx)} genes={len(tok)} {time.time() - t0:.0f}s", flush=True)
    if not all((out_dir / f"{c}__{cell}.parquet").exists() for c, i in scr_idx.items() for cell in ref.obs_names[i]):
        return
    (out_dir / "DONE").write_text(json.dumps({"model": a.model, "set": a.set, "seed": seed,
                                              "cells": {c: list(map(str, ref.obs_names[i])) for c, i in scr_idx.items()}}))


if __name__ == "__main__":
    main()
