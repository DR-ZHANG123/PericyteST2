# PericyteST2

Code for the study *High glucose uncouples PDGF-B from ST2 in retinal pericytes, matching an
IL1RL1-low, collagen-high mural-cell state in human diabetic retinopathy*.

The analysis has two halves. The first reads human donor retina (GEO GSE160306) against a
single-nucleus reference of 13 retinal cell classes, establishes by deconvolution which cells bulk
tissue can resolve, and measures rare vascular cells against stable markers of the same cell type.
The second analyses experiments in primary human retinal pericytes: PDGFRβ phosphorylation,
immediate-early genes, ST2 isoform transcripts, secreted soluble ST2 and knock-down of IL1RL1.

## Layout

```
config/run.yaml        every tunable parameter; scripts read it, nothing is hard-coded
run_all.sh             runs the stages in order, skipping those whose output exists
scripts/stages/        numbered analysis stages (Python, R, shell)
scripts/lib/repro.py   seeds, configuration and the PROVENANCE.json written by each stage
scripts/figures/       plotting only: reads result tables, writes panels and their input data
metadata/              public data sources with accessions
wetlab/                primers used for quantitative PCR
```

## Running it

```bash
conda env create -f environment.yml && conda activate pericytest2
bash scripts/stages/00_download.sh        # public data, with checksums
bash run_all.sh                           # analysis stages; a GPU is used by 09, 09b, 10 and 17
python scripts/figures/Fig3_mural_state.py   # and the other figure scripts
```

In silico perturbation with Geneformer runs in a second environment
(`environment_fm.yml`, `scripts/stages/17_fm_all.sh`) and needs a GPU.

Each stage writes `results/<stage>/PROVENANCE.json` with input and output checksums, the seed, a
hash of the configuration and an environment snapshot. Figures are drawn from the result tables
alone, and each panel's input is written to `figures/source_data/`, so a figure can be redrawn
without rerunning any model.

## Data

All donor and reference data are public and are listed with accessions in
`metadata/data_sources.tsv`: human donor retina (GEO GSE160306), the human retina and optic nerve
single-nucleus atlases (CZ CELLxGENE Discover), retinal pericytes exposed to advanced glycation end
products (GSE204880), retinal organoids in high glucose (GSE290024), Geneformer pretrained models
(Hugging Face) and ligand–receptor pairs (OmniPath). Measurements from the pericyte experiments are
parsed from the delivered laboratory file by `scripts/stages/23a_wetlab_extract.py`.

## Licence

MIT, see `LICENSE`.
