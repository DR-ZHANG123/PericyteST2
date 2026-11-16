#!/usr/bin/env bash
# pericytest2 总控：按顺序运行分析 stage。
# 用法：bash run_all.sh [--resume|--force] [--skip=stageA,stageB]
#   --resume（默认）：results/<stage>/PROVENANCE.json 已存在则跳过
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"; cd "$ROOT"
MODE=resume; SKIP=""
for a in "$@"; do
  case "$a" in --force) MODE=force ;; --resume) MODE=resume ;; --skip=*) SKIP="${a#--skip=}" ;; esac
done
CONDA_BASE="${CONDA_BASE:-$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")}"
source "$CONDA_BASE/etc/profile.d/conda.sh"; conda activate pericytest2
mkdir -p logs
export CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-0}
export MPLBACKEND=Agg

STAGES=(
  "00_download|data/raw/MD5SUMS|bash scripts/stages/00_download.sh"
  "01_magnet_diabetes|results/01_magnet_diabetes/PROVENANCE.json|Rscript scripts/stages/01_magnet_diabetes.R"
  "02_magnet_genesets|results/02_magnet_genesets/PROVENANCE.json|Rscript scripts/stages/02_magnet_genesets.R"
  "03_module_scores|results/03_module_scores/PROVENANCE.json|Rscript scripts/stages/03_module_scores.R"
  "04_effect_landscape|results/04_effect_landscape/PROVENANCE.json|Rscript scripts/stages/04_effect_landscape.R"
  "05_sva_calibrated|results/05_sva_calibrated/PROVENANCE.json|Rscript scripts/stages/05_sva_calibrated.R"
  "06_retina_pilot|results/06_retina_pilot/PROVENANCE.json|Rscript scripts/stages/06_retina_pilot.R"
  "07_paad_islet_pilot|results/07_paad_islet_pilot/PROVENANCE.json|python scripts/stages/07_paad_islet_pilot.py"
  "08_reference|results/08_reference/PROVENANCE.json|python scripts/stages/08_reference.py"
  "08b_composition|results/08b_composition/PROVENANCE.json|python scripts/stages/08b_composition.py"
  "09_deconv_benchmark|results/09_deconv_benchmark/PROVENANCE.json|python scripts/stages/09_deconv_benchmark.py"
  "09b_spikein|results/09b_spikein/PROVENANCE.json|python scripts/stages/09b_spikein.py"
  "09c_ratio_validation|results/09c_ratio_validation/PROVENANCE.json|python scripts/stages/09c_ratio_validation.py"
  "10_deconv_retina|results/10_deconv_retina/PROVENANCE.json|python scripts/stages/10_deconv_retina.py"
  "11_celltype_stage_tests|results/11_celltype_stage_tests/PROVENANCE.json|Rscript scripts/stages/11_celltype_stage_tests.R"
  "11b_ratio_retina|results/11b_ratio_retina/PROVENANCE.json|Rscript scripts/stages/11b_ratio_retina.R"
  "14a_marker_sets|results/14a_marker_sets/PROVENANCE.json|python scripts/stages/14a_marker_sets.py"
  "14b_robustness|results/14b_robustness/PROVENANCE.json|Rscript scripts/stages/14b_robustness.R"
  "15_crosstalk|results/15_crosstalk/PROVENANCE.json|Rscript scripts/stages/15_crosstalk.R"
  "16a_mural_markers|results/16a_mural_markers/PROVENANCE.json|python scripts/stages/16a_mural_markers.py"
  "16b_mural_subtypes|results/16b_mural_subtypes/PROVENANCE.json|Rscript scripts/stages/16b_mural_subtypes.R"
  "15b_pdgfb_source|results/15b_pdgfb_source/PROVENANCE.json|Rscript scripts/stages/15b_pdgfb_source.R"
  "17a_fm_screen|data/interim/fm/V2-104M_screen/DONE|bash scripts/stages/17_fm_all.sh"
  "17c_fm_summary|results/17_foundation_model/PROVENANCE.json|conda run -n pericytest2_fm python scripts/stages/17c_fm_summary.py"
  "12_invitro|results/12_invitro/PROVENANCE.json|Rscript scripts/stages/12_invitro.R"
  "13_mouse_vascular|results/13_mouse_vascular/PROVENANCE.json|python scripts/stages/13_mouse_vascular.py"
  "18a_mural_state|results/18a_mural_state/PROVENANCE.json|Rscript scripts/stages/18a_mural_state.R"
  "19_state_support|results/19_state_support/PROVENANCE.json|Rscript scripts/stages/19_state_support.R"
  "20_composition_check|results/20_composition_check/PROVENANCE.json|Rscript scripts/stages/20_composition_check.R"
  "21_state_structure|results/21_state_structure/PROVENANCE.json|Rscript scripts/stages/21_state_structure.R"
  "22_reference_umap|results/22_reference_umap/PROVENANCE.json|python scripts/stages/22_reference_umap.py"
  "23a_wetlab_extract|results/23_wetlab/measurements.tsv|python scripts/stages/23a_wetlab_extract.py"
  "23b_wetlab_stats|results/23_wetlab/PROVENANCE.json|python scripts/stages/23b_wetlab_stats.py"
)
for s in "${STAGES[@]}"; do
  IFS='|' read -r name done_file cmd <<< "$s"
  if [[ ",$SKIP," == *",$name,"* ]]; then echo "[skip] $name"; continue; fi
  if [[ "$MODE" == resume && -e "$done_file" ]]; then echo "[done] $name"; continue; fi
  echo "[run ] $name"; bash -c "$cmd" 2>&1 | tee "logs/$name.log"
done
