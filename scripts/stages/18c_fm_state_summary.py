"""Stage 18c：汇总——每个基因的 in silico 删除与“糖尿病壁细胞状态方向”的对齐程度（CPU）。

对齐 = 删除引起的嵌入变化（stage 17a）与状态方向 D（stage 18b）的余弦；细胞内按排名位置分箱做中位数/MAD 标准化
（与 stage 17c 相同的位置校正），先在供体内平均，再在供体间平均。
签名内的基因只在用“去掉自身的签名”评价时进入排名（stage 18b leave_out），其余签名基因不参与排名。
特异性 P：同一基因对 100 个打乱签名方向的对齐构成零分布。
"""
from __future__ import annotations

import json, sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.lib.repro import ROOT, load_config, set_global_seed, write_provenance  # noqa: E402

OUT = ROOT / "results/18_fm_state"
SCREEN_SET = {"V2-104M": "screen", "V1-10M": "replicate", "V2-316M": "replicate"}


def receptors(cfg) -> set[str]:
    fc = cfg["foundation_model"]
    lr = pd.read_csv(ROOT / fc["receptor_db"], sep="\t")
    lr = lr[lr["sources"].astype(str).str.split(";").str.len() >= fc["receptor_min_sources"]]
    return {g for t in lr["target_genesymbol"].astype(str) for g in t.split("_")}


def bin_z(A: np.ndarray, rel_rank: np.ndarray, bins: int) -> np.ndarray:
    b = np.minimum((rel_rank * bins).astype(int), bins - 1)
    Z = np.empty_like(A)
    for k in np.unique(b):
        m = b == k
        med = np.median(A[m], axis=0)
        mad = np.median(np.abs(A[m] - med), axis=0) + 1e-12
        Z[m] = (A[m] - med) / (1.4826 * mad)
    return Z


def model_scores(model: str, cfg, sym: dict, sig_genes: set) -> pd.DataFrame:
    fc, sm = cfg["foundation_model"], cfg["fm_state_match"]
    scr = ROOT / "data/interim/fm" / f"{model}_{SCREEN_SET[model]}"
    st = ROOT / "data/interim/fm" / f"{model}_state"
    cells = json.loads((st / "DONE").read_text())["cells"]
    D = np.load(st / "directions.npy")                                  # S × C × d
    lab = pd.read_csv(st / "directions.tsv", sep="\t", keep_default_na=False)
    D = D / (np.linalg.norm(D, axis=-1, keepdims=True) + 1e-12)
    rows = []
    for i, cell in enumerate(cells):
        df = pd.read_parquet(scr / f"pericyte__{cell}.parquet")
        dl = np.load(scr / f"pericyte__{cell}__delta.npy").astype(np.float32)
        dl = dl / (np.linalg.norm(dl, axis=1, keepdims=True) + 1e-12)
        Z = bin_z(dl @ D[:, i, :].T, df["rel_rank"].values, fc["rank_bins"])   # 基因 × 方案
        df = df[["cell", "donor", "gene_id"]].copy()
        df["symbol"] = df["gene_id"].map(sym)
        for sg in ["lfc_npdr", "lfc_diabetic"]:
            full = np.where((lab.signature == sg) & (lab.kind == "full"))[0][0]
            shuf = np.where((lab.signature == sg) & (lab.kind == "shuffle"))[0]
            obs = Z[:, full].copy()
            for g in sm["leave_out"]:
                j = np.where((lab.signature == sg) & (lab.kind == "minus") & (lab.tag == g))[0]
                if len(j):
                    obs[df["symbol"].values == g] = Z[df["symbol"].values == g, j[0]]
            d = df.copy(); d["signature"] = sg; d["align_z"] = obs
            d = pd.concat([d, pd.DataFrame(Z[:, shuf], columns=[f"null_{k}" for k in range(len(shuf))], index=d.index)], axis=1)
            rows.append(d)
    x = pd.concat(rows, ignore_index=True)
    excl = sig_genes - set(sm["leave_out"])
    x = x[~x["symbol"].isin(excl)]
    vals = ["align_z"] + [c for c in x.columns if c.startswith("null_")]
    n_cells = len(cells)
    det = x.groupby(["signature", "gene_id"])["cell"].nunique().rename("n_cells")
    by_donor = x.groupby(["signature", "gene_id", "donor"])[vals].mean()
    g = by_donor.groupby(["signature", "gene_id"])[vals].mean().join(det)
    sym_idx = g.index.get_level_values("gene_id").map(sym)
    keep = (g["n_cells"] >= fc["min_cell_frac"] * n_cells) | (sym_idx.isin(sm["focus"]) & (g["n_cells"] >= sm["focus_min_cells"]))
    g = g[keep].copy()
    null = g[[c for c in vals if c != "align_z"]].values
    g["null_mean"] = null.mean(1); g["null_sd"] = null.std(1)
    g["p_shuffle"] = (1 + (null >= g["align_z"].values[:, None]).sum(1)) / (1 + null.shape[1])
    g = g.drop(columns=[c for c in vals if c != "align_z"]).reset_index()
    g["symbol"] = g["gene_id"].map(sym); g["model"] = model; g["frac_cells"] = g["n_cells"] / n_cells
    g["pct_all"] = g.groupby("signature")["align_z"].rank(pct=True)
    g["below_min_frac"] = g["frac_cells"] < fc["min_cell_frac"]
    return g


def main() -> None:
    cfg = load_config(); set_global_seed(cfg["seed"]); sm = cfg["fm_state_match"]
    OUT.mkdir(parents=True, exist_ok=True)
    ref = ad.read_h5ad(ROOT / cfg["reference"]["out_h5ad"], backed="r")
    sym = dict(zip(ref.var["feature_id"].astype(str).str.split(".").str[0], ref.var["feature_name"].astype(str)))
    sig = pd.read_csv(ROOT / "results/18a_mural_state/state_signature.tsv", sep="\t")
    rec = receptors(cfg)
    models = [m for m in cfg["foundation_model"]["models"] if (ROOT / "data/interim/fm" / f"{m}_state" / "DONE").exists()]
    allg = pd.concat([model_scores(m, cfg, sym, set(sig["gene"])) for m in models], ignore_index=True)
    allg["is_receptor"] = allg["symbol"].isin(rec)
    allg["pct_receptor"] = np.nan
    r = allg["is_receptor"]
    allg.loc[r, "pct_receptor"] = allg[r].groupby(["model", "signature"])["align_z"].rank(pct=True)
    allg.to_csv(OUT / "state_alignment_by_model.tsv", sep="\t", index=False)

    # 跨模型：三个版本都检出的受体，取各版本分位的最小值
    rr = allg[r].pivot_table(index=["signature", "symbol"], columns="model", values="pct_receptor")
    rr = rr.dropna()
    rr["min_pct"] = rr.min(1)
    ps = allg[r].pivot_table(index=["signature", "symbol"], columns="model", values="p_shuffle").add_prefix("p_shuffle_")
    cons = rr.join(ps).reset_index().sort_values(["signature", "min_pct"], ascending=[True, False])
    cons.to_csv(OUT / "receptor_consistency.tsv", sep="\t", index=False)
    conc = []
    for sg in ["lfc_npdr", "lfc_diabetic"]:
        p = allg[allg.signature == sg].pivot_table(index="gene_id", columns="model", values="align_z").dropna()
        for i, a in enumerate(models):
            for b in models[i + 1:]:
                conc.append(dict(signature=sg, model_a=a, model_b=b, n_genes=len(p), spearman=p[a].corr(p[b], method="spearman")))
    pd.DataFrame(conc).to_csv(OUT / "model_concordance.tsv", sep="\t", index=False)
    foc = allg[allg.symbol.isin(sm["focus"])].sort_values(["signature", "symbol", "model"])
    foc.to_csv(OUT / "focus_genes.tsv", sep="\t", index=False)
    pd.set_option("display.width", 220)
    print(foc[["signature", "symbol", "model", "align_z", "p_shuffle", "pct_all", "pct_receptor"]].round(3).to_string(index=False))
    print(cons.groupby("signature").head(12).round(3).to_string(index=False))
    print(pd.DataFrame(conc).round(3))
    write_provenance("18_fm_state", inputs=[ROOT / "results/18a_mural_state/state_signature.tsv"],
                     outputs=sorted(OUT.glob("*.tsv")), seed=cfg["seed"], extra={"models": models})


if __name__ == "__main__":
    main()
