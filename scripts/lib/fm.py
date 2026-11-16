"""Geneformer 工具：取细胞、分词、嵌入、逐基因删除 / 过表达的 in silico 扰动。

嵌入定义与 Geneformer 官方 InSilicoPerturber 一致：倒数第二层隐藏状态，对基因 token 取平均（不含 <cls>/<eos>）；
删除 = 从排名编码中去掉该基因；过表达 = 把该基因移到排名最前（未检出则插到最前）。
比较时与官方一致，把被扰动基因本身从两侧嵌入中剔除：原始细胞嵌入 = 除该基因外其余基因 token 的平均，
过表达后的嵌入 = 除被移动基因外的平均。这样测到的是扰动对其余基因语境的影响，而不是该基因自身向量被拿掉/加入。
"""
from __future__ import annotations

import pickle
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import torch
from transformers import BertForMaskedLM


def dictionaries(gf_root: Path, version: str) -> dict:
    pkg = gf_root / "geneformer"
    sub = pkg / "gene_dictionaries_30m" if version == "V1" else pkg
    tag = "gc30M" if version == "V1" else "gc104M"
    out = {}
    for key, name in [("token", "token_dictionary"), ("median", "gene_median_dictionary"),
                      ("name_id", "gene_name_id_dict"), ("mapping", "ensembl_mapping_dict")]:
        with open(sub / f"{name}_{tag}.pkl", "rb") as fh:
            out[key] = pickle.load(fh)
    return out


def tokenize_cells(adata: ad.AnnData, dicts: dict, version: str, input_size: int, scale: dict | None = None) -> list[np.ndarray]:
    """排名编码：每个细胞的计数按总数归一到 1e4，再除以该基因在预训练语料中的非零中位数，按降序排列。
    scale：{Ensembl ID: 倍数}，在排名前把这些基因的归一化表达乘以对应倍数（用于构造状态改变后的细胞）。"""
    tok, med = dicts["token"], dicts["median"]
    ens = adata.var["feature_id"].astype(str).str.split(".").str[0].values
    keep = np.array([(e in tok) and (e in med) for e in ens])
    X = adata.X[:, keep]
    ens_k = ens[keep]
    tokens = np.array([tok[e] for e in ens_k], dtype=np.int64)
    norm = np.array([med[e] for e in ens_k], dtype=np.float64)
    if scale:
        norm = norm / np.array([scale.get(e, 1.0) for e in ens_k], dtype=np.float64)
    n_counts = np.asarray(adata.X.sum(axis=1)).ravel()
    special = version == "V2"
    cls_id, eos_id = tok.get("<cls>"), tok.get("<eos>")
    room = input_size - 2 if special else input_size
    out = []
    for i in range(adata.n_obs):
        row = X[i]
        idx, val = row.indices, row.data.astype(np.float64)
        score = val / n_counts[i] * 1e4 / norm[idx]
        order = idx[np.argsort(-score, kind="stable")][:room]
        seq = tokens[order]
        if special:
            seq = np.concatenate([[cls_id], seq, [eos_id]])
        out.append(seq.astype(np.int64))
    return out


def load_model(path: Path, device: str) -> BertForMaskedLM:
    m = BertForMaskedLM.from_pretrained(str(path), output_hidden_states=True, attn_implementation="sdpa")
    return m.to(device).half().eval()


@torch.no_grad()
def gene_states(model, seq: np.ndarray, layer: int, special: bool, device: str) -> np.ndarray:
    """单个细胞各基因 token 的隐藏状态（不含特殊 token），顺序与 seq 中基因位置一致。"""
    ids = torch.from_numpy(seq)[None].to(device)
    hs = model.bert(input_ids=ids, attention_mask=torch.ones_like(ids), output_hidden_states=True).hidden_states[layer][0].float()
    return (hs[1:-1] if special else hs).cpu().numpy()


@torch.no_grad()
def _pool(model, seqs: list[np.ndarray], layer: int, special: bool, device: str, skip_first_gene: bool = False) -> np.ndarray:
    L = max(len(s) for s in seqs)
    ids = torch.zeros((len(seqs), L), dtype=torch.long)
    mask = torch.zeros((len(seqs), L), dtype=torch.long)
    for j, s in enumerate(seqs):
        ids[j, :len(s)] = torch.from_numpy(s); mask[j, :len(s)] = 1
    ids, mask = ids.to(device), mask.to(device)
    hs = model.bert(input_ids=ids, attention_mask=mask, output_hidden_states=True).hidden_states[layer].float()
    gm = mask.clone().float()
    if special:                                         # 去掉 <cls> 与 <eos>
        gm[:, 0] = 0
        lens = mask.sum(1)
        gm[torch.arange(len(seqs)), lens - 1] = 0
    if skip_first_gene:                                 # 过表达：被移到最前的基因本身不计入
        gm[:, 1 if special else 0] = 0
    emb = (hs * gm.unsqueeze(-1)).sum(1) / gm.sum(1, keepdim=True)
    return emb.cpu().numpy()


def embed(model, seqs, layer, special, device, batch=64, skip_first_gene=False) -> np.ndarray:
    order = np.argsort([len(s) for s in seqs])
    out = np.zeros((len(seqs), model.config.hidden_size), dtype=np.float32)
    for k in range(0, len(seqs), batch):
        sel = order[k:k + batch]
        out[sel] = _pool(model, [seqs[i] for i in sel], layer, special, device, skip_first_gene)
    return out


def _cos(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    a = a / np.linalg.norm(a, axis=-1, keepdims=True); b = b / np.linalg.norm(b, axis=-1, keepdims=True)
    return (a * b).sum(-1)


def delete_screen_cell(model, seq, layer, special, device, batch, centroids: np.ndarray):
    """对一个细胞逐个删除每个基因。返回：基因 token、相对排名、嵌入余弦位移、对各类中心余弦相似度的变化、嵌入变化向量。
    参照嵌入 = 原始细胞中除被删基因外其余基因 token 的平均（与官方定义一致）。"""
    start = 1 if special else 0
    stop = len(seq) - 1 if special else len(seq)
    pos = np.arange(start, stop)
    H = gene_states(model, seq, layer, special, device)            # n × d
    n = H.shape[0]
    ref = (H.sum(0)[None, :] - H) / (n - 1)                        # 第 k 行：去掉第 k 个基因后的原始平均
    pert = np.zeros((len(pos), H.shape[1]), dtype=np.float32)
    for k in range(0, len(pos), batch):
        chunk = [np.delete(seq, p) for p in pos[k:k + batch]]
        pert[k:k + len(chunk)] = _pool(model, chunk, layer, special, device)
    shift = 1 - _cos(pert, ref)
    d_cent = np.stack([_cos(pert, c[None, :]) - _cos(ref, c[None, :]) for c in centroids], axis=1)
    rel_rank = (pos - start) / max(stop - start - 1, 1)
    return seq[pos], rel_rank, shift, d_cent, pert - ref


def overexpress(seq: np.ndarray, token: int, special: bool) -> np.ndarray:
    body = seq[1:-1] if special else seq
    body = np.concatenate([[token], body[body != token]])[:len(body)]
    return np.concatenate([[seq[0]], body, [seq[-1]]]) if special else body


def pick_cells(adata: ad.AnnData, cell_class: str, n: int, min_genes: int, rng, exclude=None) -> np.ndarray:
    """按供体轮流抽取细胞，避免单一供体占多数。"""
    obs = adata.obs
    ng = np.diff(adata.X.indptr) if hasattr(adata.X, "indptr") else (adata.X > 0).sum(1)
    ok = (obs["cell_class"].values == cell_class) & (np.asarray(ng).ravel() >= min_genes)
    if exclude is not None:
        ok &= ~np.isin(np.arange(adata.n_obs), exclude)
    idx = np.where(ok)[0]
    donors = obs["donor_id"].values[idx]
    pools = {d: list(rng.permutation(idx[donors == d])) for d in pd.unique(donors)}
    chosen = []
    while len(chosen) < min(n, len(idx)):
        for d in list(pools):
            if pools[d]:
                chosen.append(pools[d].pop())
                if len(chosen) >= n:
                    break
    return np.array(sorted(chosen))
