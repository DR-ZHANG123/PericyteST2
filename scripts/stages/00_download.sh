#!/usr/bin/env bash
# 下载 pericytest2 用到的公开数据，并写 md5 校验表。只下载、不改动内容。
# 用法：bash scripts/stages/00_download.sh
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
RAW="$ROOT/data/raw"
ROOT_DATA="$ROOT/data"
mkdir -p "$RAW/magnet" "$RAW/gse26887"
UA="Mozilla/5.0"

fetch() {  # fetch <url> <dest>：已存在则跳过
  local url="$1" dest="$2"
  if [[ -s "$dest" ]]; then echo "[skip] $dest"; return; fi
  echo "[get ] $url"
  curl -fsSL -A "$UA" --retry 3 -o "$dest.part" "$url"
  mv "$dest.part" "$dest"
}

# MAGNet 左心室 RNA-seq（GEO GSE141910）：原始计数与临床表型（含糖尿病、高血压、体重身高）
# 表型文件固定到 MAGNet 仓库 commit 10c96c7（2020-09-30）
MAGNET_COMMIT=10c96c727d5e637366c5dcb653393fe614f6efb2
fetch "https://www.dropbox.com/s/i5dthgl5c5ij5gd/Counts.csv?dl=1" "$RAW/magnet/Counts.csv"
fetch "https://raw.githubusercontent.com/mpmorley/MAGNet/$MAGNET_COMMIT/phenoData.csv" "$RAW/magnet/phenoData.csv"

# GSE26887：缺血性心衰左心室活检，糖尿病 7 / 非糖尿病 12 / 对照 5（Affymetrix HuGene 1.0 ST）
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE26nnn/GSE26887/matrix/GSE26887_series_matrix.txt.gz" \
      "$RAW/gse26887/GSE26887_series_matrix.txt.gz"

# GSE160306：人尸检视网膜（黄斑 + 周边），对照 / 糖尿病无 DR / NPDR / PDR（含 DME），总 RNA-seq 计数
mkdir -p "$RAW/gse160306"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE160nnn/GSE160306/suppl/GSE160306_human_retina_DR_totalRNA_counts.txt.gz" \
      "$RAW/gse160306/GSE160306_counts.txt.gz"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE160nnn/GSE160306/matrix/GSE160306_series_matrix.txt.gz" \
      "$RAW/gse160306/GSE160306_series_matrix.txt.gz"

# TCGA-PAAD：STAR 计数（UCSC Xena GDC hub）与 BCR Biotab 病人临床表（含 history_of_diabetes）
mkdir -p "$RAW/tcga_paad"
XENA=https://gdc-hub.s3.us-east-1.amazonaws.com/download
fetch "$XENA/TCGA-PAAD.star_counts.tsv.gz" "$RAW/tcga_paad/TCGA-PAAD.star_counts.tsv.gz"
fetch "$XENA/gencode.v36.annotation.gtf.gene.probemap" "$RAW/tcga_paad/gencode.v36.annotation.gtf.gene.probemap"
fetch "https://api.gdc.cancer.gov/data/92be86dd-9f8c-4dfe-ba38-22cef2f6a9c7" "$RAW/tcga_paad/clinical_patient_paad.txt"

# GSE290024：人 iPSC 视网膜类器官，对照（17 mM 葡萄糖）vs 25 mM 葡萄糖 28 天，各 3 个
mkdir -p "$RAW/gse290024"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE290nnn/GSE290024/suppl/GSE290024_RAW_counts.txt.gz" "$RAW/gse290024/GSE290024_RAW_counts.txt.gz"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE290nnn/GSE290024/matrix/GSE290024_series_matrix.txt.gz" "$RAW/gse290024/GSE290024_series_matrix.txt.gz"

# GSE189712：原代人脑周细胞 PDGF-BB 10 ng/mL 1 h / 24 h vs 溶剂（3 位供体，配对）
mkdir -p "$RAW/gse189712"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189712/suppl/GSE189712_full_table_PDGFRBB.txt.gz" "$RAW/gse189712/GSE189712_full_table_PDGFRBB.txt.gz"
# GSE311133：人冠状动脉平滑肌细胞 PDGF-BB 20 ng/mL 48 h vs 溶剂（各 3）
mkdir -p "$RAW/gse311133"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE311nnn/GSE311133/suppl/GSE311133_pdgf_repsox_gene_expression.xls.gz" "$RAW/gse311133/GSE311133_pdgf_repsox_gene_expression.xls.gz"
# GSE228900：人视网膜内皮、周细胞、星形胶质细胞微血管芯片，糖尿病刺激 vs 甘露醇（第 14、28 天，各 3）
mkdir -p "$RAW/gse228900"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE228nnn/GSE228900/suppl/GSE228900_RAW.tar" "$RAW/gse228900/GSE228900_RAW.tar"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE228nnn/GSE228900/matrix/GSE228900_series_matrix.txt.gz" "$RAW/gse228900/GSE228900_series_matrix.txt.gz"
tar -xf "$RAW/gse228900/GSE228900_RAW.tar" -C "$RAW/gse228900"

# GSE204880：人视网膜周细胞（HRPC）BSA 对照 vs 晚期糖基化终产物（AGEs），各 3 个（只取 HRPC 样本）
mkdir -p "$RAW/gse204880"
for g in GSM6199220_HRPC_control_1 GSM6199221_HRPC_control_2 GSM6199222_HRPC_control_3 \
         GSM6199223_HRPC_AGEs_1 GSM6199224_HRPC_AGEs_2 GSM6199225_HRPC_AGEs_3; do
  fetch "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6199nnn/${g%%_*}/suppl/${g}.txt.gz" "$RAW/gse204880/${g}.txt.gz"
done

# 基因注释：GENCODE v36 基因 ID → 基因名（UCSC Xena 提供的 probemap）
mkdir -p "$RAW/annotation"
fetch "https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v36.annotation.gtf.gene.probemap" "$RAW/annotation/gencode.v36.gene.probemap"

# GSE204880 小鼠视网膜血管细胞单细胞（野生型 vs db/db，21 周，PECAM1/CD140b 富集内皮与周细胞），各 1 个样本，仅作跨物种描述
mkdir -p "$RAW/gse204880/mouse"
for g in GSM6199212_WTR_1 GSM6199213_DBR_1; do
  for f in barcodes.tsv.gz features.tsv.gz matrix.mtx.gz; do
    fetch "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6199nnn/${g%%_*}/suppl/${g}_${f}" "$RAW/gse204880/mouse/${g}_${f}"
  done
done

# OmniPath 配体–受体（ligrecextra，人）：受体名单，用于 stage 17 基础模型筛查中“受体”一类的界定
mkdir -p "$RAW/omnipath"
fetch "https://omnipathdb.org/interactions?datasets=ligrecextra&genesymbols=1&organisms=9606&fields=sources,references" "$RAW/omnipath/ligrec_interactions.tsv"

# Geneformer 预训练模型（Hugging Face ctheodoris/Geneformer，固定 revision；V1-10M、V2-104M、V2-316M + 词表），装入 pericytest2_fm 环境
GF_REV=1f7fbae4e469a5f4f1af8c111a529cfe1b3829f5
if [[ ! -f "$ROOT_DATA/models/geneformer/SOURCE_SHA.txt" ]]; then
  conda run -n pericytest2_fm python -c "from huggingface_hub import snapshot_download; snapshot_download('ctheodoris/Geneformer', revision='$GF_REV', local_dir='$ROOT_DATA/models/geneformer', allow_patterns=['geneformer/*','geneformer/**','setup.py','requirements.txt','MANIFEST.in','README.md','Geneformer-V1-10M/*','Geneformer-V2-104M/*','Geneformer-V2-316M/*'])"
  echo "$GF_REV" > "$ROOT_DATA/models/geneformer/SOURCE_SHA.txt"
  conda run -n pericytest2_fm pip install -q -e "$ROOT_DATA/models/geneformer" --no-deps
fi

cd "$RAW"
find . -type f ! -name 'MD5SUMS' ! -name '*.part' -print0 | sort -z | xargs -0 md5sum > MD5SUMS
echo "done; checksums in $RAW/MD5SUMS"
