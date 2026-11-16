"""Stage 07：胰腺导管腺癌（TCGA-PAAD）中糖尿病与胰岛内分泌转录本（先导）。

问题：胰腺癌相关糖尿病是否伴随 β 细胞选择性减少？
读出：肿瘤样本中 β / α / δ 细胞与腺泡标记基因的 log-CPM；以腺泡含量校正取样带入的正常胰腺多少，
再比较有无糖尿病史（BCR Biotab history_of_diabetes）。只纳入导管腺癌（排除神经内分泌癌等）。
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from repro import ROOT, load_config, set_global_seed, write_provenance  # noqa: E402

STAGE = "07_paad_islet_pilot"
MARKERS = {
    "beta": ["INS", "IAPP", "MAFA", "NKX6-1", "SLC30A8", "G6PC2", "PCSK1"],
    "alpha": ["GCG", "IRX2", "ARX", "TTR"],
    "delta": ["SST", "HHEX", "LEPR"],
    "acinar": ["PRSS1", "CPA1", "CELA3A", "CTRB1", "PNLIP", "CPB1", "CTRC"],
}


def load_clinical(path: Path) -> pd.DataFrame:
    rows = list(csv.reader(open(path), delimiter="\t"))
    head = rows[0]
    cols = ["bcr_patient_barcode", "histologic_diagnosis", "diabetes_diagnosis_indicator",
            "age_at_initial_pathologic_diagnosis", "gender"]
    idx = [head.index(c) for c in cols]
    cl = pd.DataFrame([[r[i] for i in idx] for r in rows[3:]], columns=["case", "histology", "dm", "age", "sex"])
    cl = cl[cl["dm"].isin(["YES", "NO"]) & (cl["histology"] == "Pancreas-Adenocarcinoma Ductal Type")].copy()
    cl["age"] = pd.to_numeric(cl["age"], errors="coerce")
    cl["DM"] = (cl["dm"] == "YES").astype(int)
    return cl


def load_scores(counts: Path, probemap: Path) -> pd.DataFrame:
    ex = pd.read_csv(counts, sep="\t", index_col=0)
    ex = ex[[c for c in ex.columns if c[13:15] == "01"]]          # 原发肿瘤
    cnt = 2 ** ex - 1                                               # Xena 存的是 log2(count+1)
    lcpm = np.log2(cnt / cnt.sum() * 1e6 + 1)
    lcpm.index = pd.read_csv(probemap, sep="\t", index_col=0)["gene"].reindex(lcpm.index).values
    lcpm = lcpm.groupby(level=0).mean()
    sc = pd.DataFrame({k: lcpm.loc[[g for g in v if g in lcpm.index]].mean() for k, v in MARKERS.items()})
    for g in ["INS", "GCG", "SST"]:
        sc[g] = lcpm.loc[g]
    sc["case"] = [c[:12] for c in sc.index]
    return sc.drop_duplicates("case")


def main() -> None:
    cfg = load_config(); set_global_seed(cfg["seed"])
    raw = ROOT / cfg["paths"]["raw"] / "tcga_paad"
    files = [raw / "clinical_patient_paad.txt", raw / "TCGA-PAAD.star_counts.tsv.gz", raw / "gencode.v36.annotation.gtf.gene.probemap"]
    d = load_scores(files[1], files[2]).merge(load_clinical(files[0]), on="case")
    d["ins_minus_gcg"] = d["INS"] - d["GCG"]
    out = ROOT / "results" / STAGE; out.mkdir(parents=True, exist_ok=True)
    d.to_csv(out / "sample_scores.tsv", sep="\t", index=False)
    rows = []
    for v in ["beta", "alpha", "delta", "INS", "GCG", "SST", "ins_minus_gcg"]:
        m = smf.ols(f"{v} ~ DM + acinar + age + C(sex)", data=d).fit()
        rows.append(dict(readout=v, n=int(m.nobs), n_dm=int(d["DM"].sum()), coef_dm=m.params["DM"],
                         ci_low=m.conf_int().loc["DM", 0], ci_high=m.conf_int().loc["DM", 1], p=m.pvalues["DM"],
                         median_dm=d.loc[d.DM == 1, v].median(), median_nodm=d.loc[d.DM == 0, v].median(),
                         p_mannwhitney=mannwhitneyu(d.loc[d.DM == 1, v], d.loc[d.DM == 0, v]).pvalue))
    res = pd.DataFrame(rows); res.to_csv(out / "dm_islet_tests.tsv", sep="\t", index=False)
    print(res.round(3).to_string(index=False))
    write_provenance(STAGE, files, [out / "sample_scores.tsv", out / "dm_islet_tests.tsv"], seed=cfg["seed"])


if __name__ == "__main__":
    main()
