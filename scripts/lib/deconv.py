"""反卷积共用函数：参考按供体划分、构造伪 bulk、NNLS 基线、评价指标。"""
from __future__ import annotations

import anndata
import numpy as np
import pandas as pd
from scipy.optimize import nnls


def split_donors(obs: pd.DataFrame, frac_test: float, seed: int) -> tuple[set, set]:
    """按供体划分；每个细胞大类在训练与测试供体中都至少保留细胞（按大类分别抽供体）。"""
    rng = np.random.default_rng(seed)
    donors = np.array(sorted(obs.donor_id.unique()))
    test = set(rng.choice(donors, int(round(len(donors) * frac_test)), replace=False))
    for cls, g in obs.groupby("cell_class"):
        d = set(g.donor_id)
        if not d - test:          # 该类全部落在测试供体：移一个回训练
            test.discard(sorted(d)[0])
        if not d & test:          # 该类没有测试供体：补一个
            test.add(sorted(d)[0])
    return set(donors) - test, test


def make_pseudobulk(X: np.ndarray, labels: np.ndarray, classes: list[str], n_samples: int, n_cells: int,
                    rng: np.random.Generator, props: np.ndarray | None = None,
                    sparse_prob: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    """按给定（或随机 Dirichlet）细胞比例从细胞中有放回抽样，计数相加成伪 bulk；返回 (样本×基因, 真实比例)。
    sparse_prob：随机比例时，以此概率把一部分细胞类型的比例置零（TAPE/Scaden 的稀疏训练技巧）。"""
    idx = {c: np.where(labels == c)[0] for c in classes}
    if props is None:
        props = rng.dirichlet(np.ones(len(classes)), n_samples)
        for i in range(n_samples):
            if rng.random() < sparse_prob:
                drop = rng.choice(len(classes), rng.integers(1, len(classes) - 1), replace=False)
                props[i, drop] = 0
                props[i] /= props[i].sum()
    # 抽样结果写成稀疏“抽取计数”矩阵 W（样本×细胞），伪 bulk = W @ X，避免逐样本循环求和
    from scipy import sparse
    rows, cols, vals = [], [], []
    for i in range(n_samples):
        counts = rng.multinomial(n_cells, props[i])
        for k, c in enumerate(classes):
            if counts[k]:
                pick, cnt = np.unique(rng.choice(idx[c], counts[k], replace=True), return_counts=True)
                rows.append(np.full(len(pick), i)); cols.append(pick); vals.append(cnt)
    W = sparse.csr_matrix((np.concatenate(vals).astype(np.float32), (np.concatenate(rows), np.concatenate(cols))),
                          shape=(n_samples, X.shape[0]))
    out = np.asarray(W @ X, dtype=np.float32)
    return out, props


def to_tape_anndata(X: np.ndarray, props: np.ndarray, genes, classes) -> anndata.AnnData:
    """TAPE 训练数据格式：X 为样本×基因计数，obs 为各细胞类型比例。"""
    return anndata.AnnData(X=X, obs=pd.DataFrame(props, columns=classes), var=pd.DataFrame(index=list(genes)))


def marker_signature(X: np.ndarray, labels: np.ndarray, classes: list[str], genes: np.ndarray,
                     n_markers: int = 50) -> pd.DataFrame:
    """NNLS 基线的签名矩阵：每类平均原始计数；每类取相对其余类倍数最高的 n_markers 个基因。"""
    means = np.vstack([X[labels == c].mean(axis=0) for c in classes])          # 类 × 基因
    cpm = means / means.sum(axis=1, keepdims=True) * 1e6
    keep = set()
    for k in range(len(classes)):
        other = np.delete(cpm, k, axis=0).max(axis=0)
        fc = np.log2((cpm[k] + 1) / (other + 1))
        keep |= set(np.argsort(-fc)[:n_markers])
    keep = np.array(sorted(keep))
    return pd.DataFrame(means[:, keep].T, index=genes[keep], columns=classes)


def nnls_fractions(bulk: pd.DataFrame, sig: pd.DataFrame) -> pd.DataFrame:
    """bulk：样本×基因原始计数。在签名基因上做 NNLS，系数归一化为细胞比例。"""
    g = sig.index.intersection(bulk.columns)
    S = sig.loc[g].values
    out = []
    for _, row in bulk[g].iterrows():
        coef, _ = nnls(S, row.values.astype(float))
        out.append(coef / coef.sum() if coef.sum() > 0 else coef)
    return pd.DataFrame(out, index=bulk.index, columns=sig.columns)


def evaluate(pred: pd.DataFrame, truth: pd.DataFrame) -> pd.DataFrame:
    """逐细胞类型的 Pearson r、一致性相关系数（CCC）与 RMSE。"""
    rows = []
    for c in truth.columns:
        p, t = pred[c].values, truth[c].values
        r = np.corrcoef(p, t)[0, 1] if p.std() > 0 else np.nan
        ccc = 2 * np.cov(p, t, bias=True)[0, 1] / (p.var() + t.var() + (p.mean() - t.mean()) ** 2)
        rows.append(dict(cell_class=c, pearson=r, ccc=ccc, rmse=float(np.sqrt(np.mean((p - t) ** 2)))))
    return pd.DataFrame(rows)


VASC = ("endothelial", "pericyte", "smooth_muscle")


def retina_props(comp: pd.DataFrame, classes: list[str], n: int, rng: np.random.Generator,
                 vasc_range=(0.002, 0.03), concentration: float | None = None) -> np.ndarray:
    """视网膜样组成：神经/胶质比例取自 comp（供体×部位的实际组成），可用 Dirichlet(concentration×组成) 加扰动；
    血管三类比例在 vasc_range 内均匀抽取，其余按比例缩放。"""
    rows = comp.iloc[rng.integers(0, len(comp), n)]
    out = np.zeros((n, len(classes)))
    for k, c in enumerate(classes):
        if c in rows.columns:
            out[:, k] = rows[c].values
    if concentration is not None:
        out = np.vstack([rng.dirichlet(concentration * r + 0.05) for r in out])
        for k, c in enumerate(classes):
            if c in VASC:
                out[:, k] = 0
    vk = [classes.index(v) for v in VASC if v in classes]
    vf = rng.uniform(*vasc_range, size=(n, len(vk)))
    out *= (1 - vf.sum(axis=1, keepdims=True)) / out.sum(axis=1, keepdims=True)
    out[:, vk] = vf
    return out


def training_props(scheme: str, comp_train: pd.DataFrame, classes: list[str], n: int,
                   rng: np.random.Generator) -> np.ndarray | None:
    """训练组成：uniform = 随机 Dirichlet（返回 None，由 make_pseudobulk 生成并做稀疏化）；
    mixed = 一半随机稀疏 Dirichlet + 一半来自训练供体的视网膜样组成（加扰动，血管 0.1%–5%）。"""
    if scheme == "uniform":
        return None
    half = n // 2
    rnd = rng.dirichlet(np.ones(len(classes)), n - half)
    for i in range(len(rnd)):
        if rng.random() < 0.5:
            drop = rng.choice(len(classes), rng.integers(1, len(classes) - 1), replace=False)
            rnd[i, drop] = 0; rnd[i] /= rnd[i].sum()
    real = retina_props(comp_train, classes, half, rng, vasc_range=(0.001, 0.05), concentration=50.0)
    return np.vstack([rnd, real])[rng.permutation(n)]
