"""Stage 17b：高糖相关信号分子的 in silico 上调能否模拟周细胞失去 PDGFRB。

在 17a 主模型筛查用的同一批周细胞上：
  Δ_del  = 删除 PDGFRB 后的嵌入变化（只在检出 PDGFRB 的细胞上）；
  Δ_X    = 把候选基因 X 移到排名最前（过表达）后的嵌入变化；
  指标   = 每个细胞 cos(Δ_X, Δ_del) 的均值，以及对周细胞类中心相似度的变化。
零分布：同一检出频率区间内随机抽取的基因做同样的过表达。
用法：python scripts/stages/17b_fm_overexpress.py --model V2-104M --gpu 0
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

FREQ_BINS = [-0.01, 0.05, 0.2, 0.5, 1.0]


def main() -> None:
    ap = argparse.ArgumentParser(); ap.add_argument("--model", default=None); ap.add_argument("--gpu", default="0")
    ap.add_argument("--n_cells", type=int, default=100); ap.add_argument("--batch", type=int, default=None)  # 只影响显存
    a = ap.parse_args()
    cfg = load_config(); fc = cfg["foundation_model"]; seed = cfg["seed"]; set_global_seed(seed)
    if a.batch:
        fc["batch_size"] = a.batch
    name = a.model or fc["primary_model"]
    sub, version, input_size = fc["models"][name]
    gf_root = ROOT / fc["model_dir"]; device = f"cuda:{a.gpu}"; special = version == "V2"
    scr = ROOT / "data/interim/fm" / f"{name}_{'screen' if name == fc['primary_model'] else 'replicate'}"
    cells = json.loads((scr / "DONE").read_text())["cells"]["pericyte"][:a.n_cells]   # 与该模型 17a 同一批周细胞
    classes = json.loads((scr / "classes.json").read_text())
    out_dir = ROOT / "data/interim/fm" / f"{name}_overexpress"; out_dir.mkdir(parents=True, exist_ok=True)

    ref = ad.read_h5ad(ROOT / cfg["reference"]["out_h5ad"])
    dicts = fm.dictionaries(gf_root, version); tok = dicts["token"]
    sym2ens = dict(zip(ref.var["feature_name"].astype(str), ref.var["feature_id"].astype(str).str.split(".").str[0]))
    peri = ref[ref.obs["cell_class"] == "pericyte"]
    freq = pd.Series(np.asarray((peri.X > 0).mean(0)).ravel(), index=ref.var["feature_id"].astype(str).str.split(".").str[0].values)
    freq = freq[[g in tok for g in freq.index]]
    freq = freq[~freq.index.duplicated()]
    fbin = pd.cut(freq, FREQ_BINS, labels=False)

    rng = np.random.default_rng(seed)
    cand = {s: sym2ens[s] for s in fc["overexpress"] if s in sym2ens and sym2ens[s] in tok}
    pdgfrb = sym2ens[fc["delete_reference"]]
    null = {}
    for b in sorted(set(fbin[list(cand.values())])):
        pool = [g for g in fbin.index[fbin == b] if g not in cand.values() and g != pdgfrb]
        null[int(b)] = list(rng.choice(pool, size=min(fc["n_null_overexpress"], len(pool)), replace=False))
    genes = {**{f"cand:{s}": e for s, e in cand.items()},
             **{f"null{b}:{e}": e for b, lst in null.items() for e in lst}}

    model = fm.load_model(gf_root / sub, device)
    cents = np.load(scr / "centroids.npy")
    c_peri = cents[classes.index("pericyte")]
    seqs = fm.tokenize_cells(ref[cells], dicts, version, input_size)
    rows = []
    special_off = 1 if special else 0
    for cell, seq in zip(cells, seqs):
        H = fm.gene_states(model, seq, fc["emb_layer"], special, device)
        body = seq[special_off:special_off + H.shape[0]]
        has = tok[pdgfrb] in body
        d_del = None
        if has:                                             # 删除 PDGFRB：参照 = 原始细胞去掉 PDGFRB 的平均
            k = int(np.where(body == tok[pdgfrb])[0][0])
            ref_del = (H.sum(0) - H[k]) / (H.shape[0] - 1)
            d_del = fm.embed(model, [seq[seq != tok[pdgfrb]]], fc["emb_layer"], special, device)[0] - ref_del
        labs = list(genes)
        pert = fm.embed(model, [fm.overexpress(seq, tok[genes[l]], special) for l in labs], fc["emb_layer"], special, device,
                        fc["batch_size"], skip_first_gene=True)          # 过表达：被移动的基因不计入
        refs = []
        for l in labs:                                      # 参照 = 原始细胞去掉该基因（若检出）的平均
            hit = np.where(body == tok[genes[l]])[0]
            refs.append((H.sum(0) - H[hit[0]]) / (H.shape[0] - 1) if len(hit) else H.mean(0))
        refs = np.stack(refs)
        delta = pert - refs
        d_peri = fm._cos(pert, c_peri[None, :]) - fm._cos(refs, c_peri[None, :])
        sim = fm._cos(delta, d_del[None, :]) if has else np.full(len(labs), np.nan)
        for l, s_, dp in zip(labs, sim, d_peri):
            rows.append({"cell": cell, "label": l, "gene_id": genes[l], "has_pdgfrb": has, "cos_to_pdgfrb_deletion": s_,
                         "d_pericyte": dp, "freq_bin": int(fbin.get(genes[l], -1))})
    df = pd.DataFrame(rows)
    df.to_parquet(out_dir / "per_cell.parquet")
    print(df.groupby(df["label"].str.split(":").str[0])[["cos_to_pdgfrb_deletion", "d_pericyte"]].mean())


if __name__ == "__main__":
    main()
