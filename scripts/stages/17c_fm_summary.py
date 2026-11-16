"""Stage 17c：汇总 Geneformer 扰动结果（CPU）。

1) 基因层面依赖性：每个细胞内按排名位置分箱，删除效应减去同箱中位数（排名靠前的基因天然位移大，需校正），
   再在细胞间平均；对周细胞身份的影响用“与周细胞类中心相似度的下降”。
2) 受体排名：OmniPath 受体中，周细胞身份最依赖哪些；PDGFRB 的分位。
3) 周细胞 vs 平滑肌：同一基因删除对两类细胞身份的影响差异（按供体自助法给置信区间）。
4) 三个模型间的一致性。
5) 过表达：候选基因上调与删除 PDGFRB 的方向相似度，相对同检出频率随机基因的经验 P。
"""
from __future__ import annotations

import json, sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.lib.repro import ROOT, load_config, set_global_seed, write_provenance  # noqa: E402

OUT = ROOT / "results/17_foundation_model"


def receptors(cfg) -> set[str]:
    fc = cfg["foundation_model"]
    lr = pd.read_csv(ROOT / fc["receptor_db"], sep="\t")
    lr["n_src"] = lr["sources"].astype(str).str.split(";").str.len()
    lr = lr[lr["n_src"] >= fc["receptor_min_sources"]]
    return {g for t in lr["target_genesymbol"].astype(str) for g in t.split("_")}


def fc_focus(cfg) -> list[str]:
    return ["PDGFRB", "NOTCH3", "EDNRA", "MCAM", "OSMR", "ITGB1", "TGFBR2", "LRP1"]


def load_screen(d: Path) -> pd.DataFrame:
    """只读 DONE 中登记的细胞（目录里可能有更早、更大细胞集合留下的文件）。"""
    cells = json.loads((d / "DONE").read_text())["cells"]
    return pd.concat([pd.read_parquet(d / f"{c}__{x}.parquet") for c, xs in cells.items() for x in xs], ignore_index=True)


def position_adjust(df: pd.DataFrame, bins: int) -> pd.DataFrame:
    df = df.copy()
    df["bin"] = np.minimum((df["rel_rank"] * bins).astype(int), bins - 1)
    df["self"] = np.where(df["cell_class"] == "pericyte", df["d_pericyte"], df["d_smooth_muscle"])
    df["identity_loss"] = -df["self"]                       # 正值 = 离本类中心更远
    for col in ["shift", "identity_loss"]:
        med = df.groupby(["cell", "bin"])[col].transform("median")
        mad = df.groupby(["cell", "bin"])[col].transform(lambda x: (x - x.median()).abs().median() + 1e-12)
        df[f"{col}_z"] = (df[col] - med) / (1.4826 * mad)
    return df


def gene_summary(df: pd.DataFrame, n_cells: int, min_frac: float, sym: dict, rec: set) -> pd.DataFrame:
    g = df.groupby("gene_id").agg(n=("cell", "nunique"), mean_rel_rank=("rel_rank", "mean"),
                                  shift_z=("shift_z", "mean"), identity_z=("identity_loss_z", "mean"),
                                  identity_raw=("identity_loss", "mean"))
    g["frac_cells"] = g["n"] / n_cells
    g = g[g["frac_cells"] >= min_frac].copy()
    g["symbol"] = [sym.get(i, i) for i in g.index]
    g["is_receptor"] = g["symbol"].isin(rec)
    g["identity_pct"] = g["identity_z"].rank(pct=True)
    g["identity_pct_receptors"] = np.where(g["is_receptor"], g.loc[g["is_receptor"], "identity_z"].rank(pct=True)
                                           .reindex(g.index), np.nan)
    return g.sort_values("identity_z", ascending=False)


def donor_boot_diff(df: pd.DataFrame, gene: str, rng, n_boot=2000):
    sub = df[df["gene_id"] == gene]
    per = sub.groupby(["cell_class", "donor"])["identity_loss_z"].mean().reset_index()
    a = per.loc[per["cell_class"] == "pericyte", "identity_loss_z"].values
    b = per.loc[per["cell_class"] == "smooth_muscle", "identity_loss_z"].values
    if len(a) < 3 or len(b) < 3:
        return dict(diff=np.nan, lo=np.nan, hi=np.nan, p_mwu=np.nan, n_peri_donors=len(a), n_smc_donors=len(b))
    boots = [rng.choice(a, len(a)).mean() - rng.choice(b, len(b)).mean() for _ in range(n_boot)]
    return dict(diff=a.mean() - b.mean(), lo=np.percentile(boots, 2.5), hi=np.percentile(boots, 97.5),
                p_mwu=stats.mannwhitneyu(a, b).pvalue, n_peri_donors=len(a), n_smc_donors=len(b))


def main() -> None:
    cfg = load_config(); fc = cfg["foundation_model"]; seed = cfg["seed"]; set_global_seed(seed)
    rng = np.random.default_rng(seed); OUT.mkdir(parents=True, exist_ok=True)
    var = ad.read_h5ad(ROOT / cfg["reference"]["out_h5ad"], backed="r").var
    sym = dict(zip(var["feature_id"].astype(str).str.split(".").str[0], var["feature_name"].astype(str)))
    rec = receptors(cfg)
    fmdir = ROOT / "data/interim/fm"
    runs = {m: fmdir / f"{m}_{'screen' if m == fc['primary_model'] else 'replicate'}" for m in fc["models"]}

    summaries, contrasts = [], []
    for m, d in runs.items():
        df = position_adjust(load_screen(d), fc["rank_bins"])
        for c in ["pericyte", "smooth_muscle"]:
            sub = df[df["cell_class"] == c]
            g = gene_summary(sub, sub["cell"].nunique(), fc["min_cell_frac"], sym, rec)
            g.insert(0, "cell_class", c); g.insert(0, "model", m)
            summaries.append(g.reset_index())
        if m == fc["primary_model"]:
            # 删除后细胞向哪类细胞靠拢：对各类中心相似度变化，按细胞内同排名箱校正（与身份得分同法）
            dcols = [c for c in df.columns if c.startswith("d_")]
            for col in dcols:
                df[f"{col}_adj"] = df[col] - df.groupby(["cell", "bin"])[col].transform("median")
            dirn = df.groupby(["cell_class", "gene_id"])[[f"{c}_adj" for c in dcols]].mean()
            dirn.columns = [c[2:-4] for c in dirn.columns]
            dirn = dirn.reset_index(); dirn["symbol"] = dirn["gene_id"].map(sym)
            keep = summaries[-2]["gene_id"].tolist() + summaries[-1]["gene_id"].tolist()
            dirn[dirn["gene_id"].isin(keep)].to_csv(OUT / "deletion_direction.tsv", sep="\t", index=False)
            pg = summaries[-2].set_index("gene_id"); sg = summaries[-1].set_index("gene_id")
            both = pg.index.intersection(sg.index)
            top = set(pg.loc[both].sort_values("identity_z", ascending=False).index[:200]) | \
                  {i for i in both if pg.loc[i, "is_receptor"]}
            for gid in sorted(top):
                r = donor_boot_diff(df, gid, rng); r.update(gene_id=gid, symbol=sym.get(gid, gid),
                                                           is_receptor=bool(pg.loc[gid, "is_receptor"]))
                contrasts.append(r)
    summ = pd.concat(summaries, ignore_index=True)
    summ.to_csv(OUT / "gene_dependency.tsv", sep="\t", index=False)
    con = pd.DataFrame(contrasts)
    con["fdr"] = stats.false_discovery_control(con["p_mwu"].fillna(1))
    con.sort_values("diff", ascending=False).to_csv(OUT / "pericyte_vs_smc_dependency.tsv", sep="\t", index=False)

    # 模型间一致性：周细胞身份依赖得分的 Spearman 相关；各模型中 PDGFRB 的分位
    cons = []
    for cl in ["pericyte", "smooth_muscle"]:
        for val in ["identity_z", "shift_z"]:
            piv = summ[summ["cell_class"] == cl].pivot_table(index="gene_id", columns="model", values=val)
            ms = list(piv.columns)
            for i, a in enumerate(ms):
                for b in ms[i + 1:]:
                    ok = piv[[a, b]].dropna(); r, p = stats.spearmanr(ok[a], ok[b])
                    cons.append({"cell_class": cl, "measure": val, "model_a": a, "model_b": b, "n_genes": len(ok), "spearman": r, "p": p})
    pd.DataFrame(cons).to_csv(OUT / "model_concordance.tsv", sep="\t", index=False)

    # 跨模型一致性：每个基因在三个模型中周细胞身份依赖分位数的最小值；PDGFRB 的最小分位数在全部基因/受体中的经验位置
    pct = summ.pivot_table(index=["gene_id", "cell_class"], columns="model", values="identity_pct").dropna()
    pct["min_pct"] = pct.min(axis=1)
    pct = pct.reset_index(); pct["symbol"] = pct["gene_id"].map(sym); pct["is_receptor"] = pct["symbol"].isin(rec)
    rows = []
    for cl in ["pericyte", "smooth_muscle"]:
        sub = pct[pct.cell_class == cl]
        for gname in fc_focus(cfg):
            r = sub[sub.symbol == gname]
            if not len(r):
                continue
            v = r.min_pct.iloc[0]
            rows.append({"cell_class": cl, "symbol": gname, "min_percentile_across_models": v,
                         "genes": len(sub), "genes_at_or_above": int((sub.min_pct >= v).sum()),
                         "fraction_at_or_above": float((sub.min_pct >= v).mean()),
                         "receptors": int(sub.is_receptor.sum()), "receptors_at_or_above": int((sub[sub.is_receptor].min_pct >= v).sum())})
    pct.to_csv(OUT / "cross_model_percentiles.tsv", sep="\t", index=False)
    pd.DataFrame(rows).to_csv(OUT / "cross_model_consistency.tsv", sep="\t", index=False)
    focus = summ[summ["symbol"].isin(["PDGFRB", "NOTCH3", "TGFBR2", "EDNRA", "ITGB1", "LRP1", "TEK", "S1PR3",
                                      "IGF1R", "INSR", "EPHB4", "ACVRL1", "ITGA1", "AXL"])]
    focus.to_csv(OUT / "receptor_focus.tsv", sep="\t", index=False)

    # 过表达
    ov = []
    for m in fc["models"]:
        f = fmdir / f"{m}_overexpress" / "per_cell.parquet"
        if not f.exists():
            continue
        pc = pd.read_parquet(f)
        pc = pc[pc["has_pdgfrb"]]
        agg = pc.groupby(["label", "gene_id", "freq_bin"])[["cos_to_pdgfrb_deletion", "d_pericyte"]].mean().reset_index()
        agg["kind"] = agg["label"].str.split(":").str[0].str.replace(r"\d", "", regex=True)
        for _, r in agg[agg["kind"] == "cand"].iterrows():
            nul = agg[(agg["kind"] == "null") & (agg["freq_bin"] == r["freq_bin"])]
            ov.append({"model": m, "symbol": r["label"].split(":")[1], "freq_bin": r["freq_bin"],
                       "cos_to_pdgfrb_deletion": r["cos_to_pdgfrb_deletion"],
                       "null_mean": nul["cos_to_pdgfrb_deletion"].mean(), "null_sd": nul["cos_to_pdgfrb_deletion"].std(),
                       "p_emp": (1 + (nul["cos_to_pdgfrb_deletion"] >= r["cos_to_pdgfrb_deletion"]).sum()) / (1 + len(nul)),
                       "d_pericyte": r["d_pericyte"], "d_pericyte_null_mean": nul["d_pericyte"].mean(),
                       "p_emp_identity": (1 + (nul["d_pericyte"] <= r["d_pericyte"]).sum()) / (1 + len(nul)),
                       "n_null": len(nul), "n_cells": pc["cell"].nunique()})
        agg.to_csv(OUT / f"overexpress_all_{m}.tsv", sep="\t", index=False)
    if ov:
        pd.DataFrame(ov).to_csv(OUT / "overexpress_vs_pdgfrb.tsv", sep="\t", index=False)

    outs = sorted(OUT.glob("*.tsv"))
    write_provenance("17_foundation_model", [ROOT / fc["receptor_db"]], outs, seed,
                     extra={"models": fc["models"], "geneformer_revision": (ROOT / fc["model_dir"] / "SOURCE_SHA.txt").read_text().strip()})


if __name__ == "__main__":
    main()
