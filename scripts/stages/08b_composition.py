"""Stage 08b：HRCA 各供体×部位的神经/胶质大类组成（用于构造视网膜样伪 bulk 组成）。"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from repro import ROOT, load_config, set_global_seed, write_provenance  # noqa: E402
from retina_ref import hrca_composition  # noqa: E402


def main() -> None:
    cfg = load_config(); set_global_seed(cfg["seed"])
    comp = hrca_composition(cfg)
    path = ROOT / cfg["reference"]["composition"]; path.parent.mkdir(parents=True, exist_ok=True)
    comp.to_csv(path, sep="\t")
    print(comp.shape); print(comp.mean().round(4))
    write_provenance("08b_composition", [], [path], seed=cfg["seed"], extra={"census_version": cfg["reference"]["census_version"]})


if __name__ == "__main__":
    main()
