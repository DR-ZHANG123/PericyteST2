"""Stage 23b：HRPC 实验的统计（以独立培养孔/重复为单位；log2 尺度）。

六组设计（NG/HG/MAN × PDGF-BB）：每个读出拟合 log2(值) ~ 组（各组均值模型，残差方差合并），给出预先指定的对比：
  PDGF-BB 在 NG、HG、MAN 中各自的效应；HG、MAN 相对 NG 的基础差异；
  交互作用（HG 中的 PDGF-BB 效应 − NG 中的效应；MAN 同理）——即高糖是否削弱 PDGF-BB 应答，甘露醇是否同样削弱；
  以及 HG 与 MAN 直接比较（高糖效应是否超出等渗透压对照）。
siRNA：log2(值) ~ 组，两条 siIL1RL1 各自对 siNC。
输出 results/23_wetlab/contrasts.tsv、group_summary.tsv。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from repro import ROOT, load_config, write_provenance  # noqa: E402

STAGE = "23_wetlab"
FACTORIAL = {
    "PDGF-BB in NG": {"NG+PDGF-BB": 1, "NG": -1},
    "PDGF-BB in HG": {"HG+PDGF-BB": 1, "HG": -1},
    "PDGF-BB in MAN": {"MAN+PDGF-BB": 1, "MAN": -1},
    "HG vs NG (no PDGF-BB)": {"HG": 1, "NG": -1},
    "MAN vs NG (no PDGF-BB)": {"MAN": 1, "NG": -1},
    "HG+PDGF-BB vs NG+PDGF-BB": {"HG+PDGF-BB": 1, "NG+PDGF-BB": -1},
    "MAN+PDGF-BB vs NG+PDGF-BB": {"MAN+PDGF-BB": 1, "NG+PDGF-BB": -1},
    "Interaction HG (PDGF-BB effect in HG minus NG)": {"HG+PDGF-BB": 1, "HG": -1, "NG+PDGF-BB": -1, "NG": 1},
    "Interaction MAN (PDGF-BB effect in MAN minus NG)": {"MAN+PDGF-BB": 1, "MAN": -1, "NG+PDGF-BB": -1, "NG": 1},
    "HG+PDGF-BB vs MAN+PDGF-BB": {"HG+PDGF-BB": 1, "MAN+PDGF-BB": -1},
    "Interaction HG vs MAN (PDGF-BB effect in HG minus MAN)": {"HG+PDGF-BB": 1, "HG": -1, "MAN+PDGF-BB": -1, "MAN": 1},
}
SIRNA = {"siIL1RL1-1 vs siNC": {"siIL1RL1-1": 1, "siNC": -1}, "siIL1RL1-2 vs siNC": {"siIL1RL1-2": 1, "siNC": -1}}


def contrasts(d: pd.DataFrame, spec: dict) -> list[dict]:
    y = np.log2(d.value)
    g = d.group
    means, ns = y.groupby(g).mean(), y.groupby(g).size()
    resid = y - g.map(means)
    df = len(y) - len(means)
    s2 = float((resid ** 2).sum() / df)
    rows = []
    for name, w in spec.items():
        est = sum(c * means[k] for k, c in w.items())
        se = np.sqrt(s2 * sum(c ** 2 / ns[k] for k, c in w.items()))
        t = est / se; p = 2 * stats.t.sf(abs(t), df); q = stats.t.ppf(0.975, df)
        rows.append(dict(contrast=name, log2_diff=est, ci_low=est - q * se, ci_high=est + q * se,
                         fold=2 ** est, t=t, df=df, p=p))
    return rows


def main() -> None:
    out = ROOT / "results" / STAGE
    m = pd.read_csv(out / "measurements.tsv", sep="\t")
    rows = []
    for (exp, rd), d in m.groupby(["experiment", "readout"], sort=False):
        spec = SIRNA if exp.startswith("sirna") else FACTORIAL
        for r in contrasts(d, spec):
            rows.append(dict(experiment=exp, readout=rd, n_per_group=int(d.groupby("group").size().min()), **r))
    c = pd.DataFrame(rows)
    c.to_csv(out / "contrasts.tsv", sep="\t", index=False, float_format="%.6g")
    s = (m.assign(log2=np.log2(m.value)).groupby(["experiment", "readout", "group"], sort=False)
         .agg(n=("value", "size"), mean=("value", "mean"), sd=("value", "std"), geo_mean=("log2", lambda v: 2 ** v.mean()))
         .reset_index())
    s.to_csv(out / "group_summary.tsv", sep="\t", index=False, float_format="%.6g")
    pd.set_option("display.width", 220)
    print(c[["experiment", "readout", "contrast", "fold", "ci_low", "ci_high", "p"]].to_string())
    write_provenance(STAGE, [out / "measurements.tsv"], [out / "contrasts.tsv", out / "group_summary.tsv"])


if __name__ == "__main__":
    main()
