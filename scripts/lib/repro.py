"""可复现性工具（Python 版）：读配置、固定种子、写 results/<stage>/PROVENANCE.json。"""
from __future__ import annotations

import hashlib, json, os, platform, random, subprocess, sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config" / "run.yaml"


def load_config() -> dict:
    return yaml.safe_load(CONFIG.read_text())


def set_global_seed(seed: int) -> None:
    random.seed(seed); np.random.seed(seed); os.environ["PYTHONHASHSEED"] = str(seed)


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _rel(p: Path) -> str:
    p = Path(p).resolve()
    return str(p.relative_to(ROOT.resolve())) if p.is_relative_to(ROOT.resolve()) else str(p)


def write_provenance(stage: str, inputs, outputs, seed=None, extra: dict | None = None) -> Path:
    out_dir = ROOT / "results" / stage; out_dir.mkdir(parents=True, exist_ok=True)
    try:
        commit = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        commit = "unknown"
    rec = {"stage": stage, "finished_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"), "seed": seed,
           "config_sha": hashlib.sha256(CONFIG.read_bytes()).hexdigest()[:16], "git_commit": commit,
           "python": sys.version.split()[0], "platform": platform.platform(),
           "conda_env": Path(sys.prefix).name,
           "inputs": {_rel(p): md5(Path(p)) for p in inputs if Path(p).is_file()},
           "outputs": {_rel(p): md5(Path(p)) for p in outputs if Path(p).is_file()}}
    rec.update(extra or {})
    path = out_dir / "PROVENANCE.json"; path.write_text(json.dumps(rec, indent=2, ensure_ascii=False))
    return path
