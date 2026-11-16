#!/usr/bin/env bash
# Stage 17a/17b：Geneformer 扰动（pericytest2_fm 环境，需 GPU）。GPU 编号用 FM_GPU 指定（默认 0）。
# 共享 GPU 时可把 17a 拆成多个分片并行：--shard k/n（结果按细胞写文件，与分片方式无关）。
set -euo pipefail
G=${FM_GPU:-0}
for spec in "V2-104M screen" "V2-316M replicate" "V1-10M replicate"; do
  set -- $spec
  conda run -n pericytest2_fm python scripts/stages/17a_fm_screen.py --model "$1" --set "$2" --gpu "$G"
done
for m in V2-104M V2-316M V1-10M; do
  conda run -n pericytest2_fm python scripts/stages/17b_fm_overexpress.py --model "$m" --gpu "$G"
done
