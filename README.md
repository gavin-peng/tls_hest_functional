# TLS Functional State Prediction from H&E (HEST-1k)

Characterises TLS (tertiary lymphoid structure) functional states at the instance
level using spatial transcriptomics (ST) as supervision. Pipeline: HookNet-TLS
detection -> ST neighbourhood extraction -> per-instance feature matrix ->
Leiden clustering -> morphological transfer (Random Forest + Graph Attention Network).

Manuscript: `paper/manuscript.md` / `paper/manuscript.pdf`


---

## Prerequisites

- WSL2 with NVIDIA GPU (GTX 1660 Ti or similar) for HookNet inference
- Python 3.11 venv at `.venv_hooknet/` for Stage 1 native inference
- Python 3.12 venv at `.venv/` for Stage 2 onwards
- `data/HEST_v1_3_0.csv` -- HEST metadata (1,276 samples)
- `data/ensembl_to_symbol.tsv` -- Ensembl -> gene symbol mapping
- H&E WSIs and ST h5ad files under `data/`

---

## Pipeline Overview

```
Stage 1: HookNet-TLS detection (H&E -> TLS polygons)
  scripts: src/stage1_tls_detect/
  generates:
  outputs/tls_results_cc354.tsv          raw polygon counts, cc≥354 filter
  outputs/tls_spot_coverage.tsv          per-sample instance counts after watershed split
  outputs/tls_masks/<sample_id>/         per-sample masks, heatmaps, coverage JSONs
  outputs/tls_detection_summary.tsv      all 747 instances × metadata (audit trail)

Stage 2: TLS characterisation (ST -> functional state clusters)
  scripts: src/stage2/
  generates:
  outputs/stage2_clustering/tls_neighborhood_spots.h5ad   ~64k spots × 40,270 genes
  outputs/stage2_clustering/instance_features_725.tsv     725 × 12 features (pre-BC23377)
  outputs/stage2_clustering/instance_features_670.tsv     670 × 12 features (canonical)
  outputs/stage2_clustering/leiden_670_r0.3.tsv           6 clusters (canonical)
  outputs/stage2_clustering/sample_patient_map.tsv

Stage 3: Morphological transfer (H&E -> C1/C4 classification)
  scripts: src/stage3/
  generates:
  outputs/stage3/  Random Forest AUROC=0.744; GAT AUROC=0.857

Validation: TCGA-BRCA / METABRIC / TCGA-CESC / Shiao TNBC
  scripts: src/validation_*.py
  generates:
  outputs/validation_tcga/
  outputs/validation_geo/
```

---

## Stage 1 -- TLS Detection

Scripts in `src/stage1_tls_detect/`. Run in order.

### 00_make_manifest.sh
Filter HEST samples by cancer type, TLS likelihood, and data availability.
Produces priority-tiered run manifest.

```bash
bash src/stage1_tls_detect/00_make_manifest.sh
```

**Output:** `outputs/run_manifest.tsv`

---

### 00_download_data.sh
Download HEST-1k H&E WSIs and ST h5ad files from HuggingFace. Skips existing files.

```bash
bash src/stage1_tls_detect/00_download_data.sh
```

**Output:** `data/<sample_id>/` directories

---

### 01_asap_preprocess.sh
Preprocess a WSI for HookNet using the ASAP Docker container. Produces
pyramidal TIFF at 0.5 µm/px + tissue mask.

```bash
bash src/stage1_tls_detect/01_asap_preprocess.sh <sample_id>
```

**Output:** `data/<sample_id>/asap_preprocessed/`
 
(They are cleaned up by  04_batch_pipeline.sh after hooknet consumed them)

---

### 02_build_hooknet_docker.sh / 03b_setup_hooknet_native.sh
Build the HookNet-TLS Docker image, or set up `.venv_hooknet/` for native GPU
inference (recommended on WSL2 where Docker GPU passthrough is unreliable).

```bash
bash src/stage1_tls_detect/02_build_hooknet_docker.sh
# or
bash src/stage1_tls_detect/03b_setup_hooknet_native.sh
```

---

### 03c_run_hooknet_native.sh
Run HookNet-TLS inference on a single sample. Produces per-sample
`_with_confidence.json` polygon file.

```bash
bash src/stage1_tls_detect/03c_run_hooknet_native.sh <sample_id>
```

**Output:** `outputs/tls_masks/<sample_id>/images/<sample_id>_asap_hooknettls_with_confidence.json`

---

### 04_batch_pipeline.sh
Orchestrator: loop over `outputs/run_manifest.tsv`, run ASAP preprocessing then
HookNet inference for all samples. Resumes from where it left off.

```bash
bash src/stage1_tls_detect/04_batch_pipeline.sh
```

---

### 05_refilter.sh
Re-run polygon filtering on existing `*_with_confidence.json` files without
re-running inference. Canonical threshold: `--tls-confidence-count 354`
(≈ 31,329 µm² minimum area, ≈ 200 µm diameter, grounded in published minimum
TLS size). Flags duplicate WSI slides and Xenium samples.

```bash
bash src/stage1_tls_detect/05_refilter.sh --tls-confidence-count 354
```

**Output:** `outputs/tls_results_cc354.tsv`

---

### 06_tls_spot_coverage.sh
Spatial join: count ST spots inside each TLS polygon and its 300 µm neighbourhood
ring. Splits large polygons (> 1 mm²) into sub-instances via watershed segmentation
(Gaussian σ = 100 µm, minimum inter-centre distance = 300 µm).

Excludes:
- Lymph node samples (n = 52) -- organised secondary lymphoid organs, not ectopic TLS
- Duplicate WSI slides (n = 13)
- Xenium samples (stored separately; incomparable spot density)

Instances with ≥ 1 co-registered ST spot are retained -> **747 instances across 119 samples**.
See `outputs/tls_detection_summary.tsv` for the full per-instance audit trail.

```bash
bash src/stage1_tls_detect/06_tls_spot_coverage.sh
```

**Output:** `outputs/tls_spot_coverage.tsv`,
`outputs/tls_masks/<sample_id>/tls_spot_coverage.json` (per-sample),
`outputs/tls_detection_summary.tsv` (all 747 instances × metadata)

---

### viz_tls_heatmap.py
Visualise TLS detection overlaid on H&E for QC.

```bash
source .venv/bin/activate
python3 src/stage1_tls_detect/viz_tls_heatmap.py <sample_id>
```

**Output:** `outputs/viz/<sample_id>_tls_viz.png`

---

## Stage 2 -- TLS Characterisation

Scripts in `src/stage2/`. Run in order.

### 00_make_patient_map.py
Build sample->patient map from HEST metadata. Disambiguates cross-study patient
ID collisions. SPA BC##### codes are globally unique and used directly.

```bash
source .venv/bin/activate
python3 src/stage2/00_make_patient_map.py
```

**Output:** `outputs/stage2_clustering/sample_patient_map.tsv`

---

### 03_extract_tls_neighborhood.py
Extract ST spots inside each TLS instance and its 300 µm neighbourhood ring.
Builds a unified AnnData with union gene space (40,270 genes) and zero-padding.
Uses ≥ 1 spot threshold; applies same lymph node / duplicate / MISC17 exclusions
as Stage 1.

```bash
python3 src/stage2/03_extract_tls_neighborhood.py
```

**Input:** `outputs/tls_spot_coverage.tsv`, `outputs/tls_masks/`,
per-sample ST h5ad files in `data/`  
**Output:** `outputs/stage2_clustering/tls_neighborhood_spots.h5ad`

---

### 03b_fix_oncotree.py
Fill missing `oncotree_code` values for samples not covered by HEST metadata.

```bash
python3 src/stage2/03b_fix_oncotree.py
```

**Input/Output:** `outputs/stage2_clustering/tls_neighborhood_spots.h5ad` (in-place)

---

### 04_instance_features.py
Build per-TLS feature matrix: 7 in-TLS gene module scores + 5 neighbourhood
scores (300 µm ring), computed on log(1+CPM)-normalised counts lazily per-instance.

```bash
python3 src/stage2/04_instance_features.py
```

**Input:** `outputs/stage2_clustering/tls_neighborhood_spots.h5ad`  
**Output:** `outputs/stage2_clustering/instance_features.tsv` (747 × 17)

---

### 05_batch_correct.py
Z-score scale features. Excludes SPA51/52/53 (patient BC23803 -- 22 instances,
single-patient dominance of a C5-like cluster). Canonical output: 725-instance
z-scaled matrix.

```bash
python3 src/stage2/05_batch_correct.py
```

**Input:** `outputs/stage2_clustering/instance_features.tsv`  
**Output:** `outputs/stage2_clustering/instance_features_725.tsv` (canonical pre-clustering)

---

### 04b_spatial_features.py
Extend feature matrix with polygon morphology (log area, circularity, convexity)
and radial expression gradients. Excludes BC23803 (SPA51–53) and BC23377
(SPA54–56). **Not used for primary clustering** (platform-confounded); available
for post-hoc analysis.

```bash
python3 src/stage2/04b_spatial_features.py
```

**Output:** `outputs/stage2_clustering/instance_features_raw_spatial.tsv` (670 × 22, unscaled)

---

### 06_cluster.py
Leiden clustering + UMAP on the 670-instance feature matrix (after additional
exclusion of SPA54/55/56, patient BC23377). Canonical resolution 0.3 -> 6 clusters
(ARI = 0.83 ± 0.14, 20 seeds).

```bash
python3 src/stage2/06_cluster.py \
    --features-file outputs/stage2_clustering/instance_features_670.tsv \
    --variant-name 670 --resolutions 0.1 0.3 0.5
```

**Input:** `outputs/stage2_clustering/instance_features_670.tsv`  
**Output:** `outputs/stage2_clustering/leiden_670_r0.3.tsv` **(canonical)**,
`outputs/stage2_clustering/cluster_summary.tsv`

---

## Stage 3 -- Morphological Transfer

Scripts in `src/stage3/`.

### 01_phase1_features.py / 02_phase1_classifier.py
Random Forest baseline: 9 aggregate H&E features, C1 vs C4 LOPO-CV.
Result: AUROC = 0.744.

### 03_phase2_embeddings.py / 04_phase2_gnn.py
Graph Attention Network on two-region spatial graph (in-TLS + neighbourhood).
Result: AUROC = 0.857 (+0.113 over RF baseline).

```bash
python3 src/stage3/01_phase1_features.py
python3 src/stage3/02_phase1_classifier.py
python3 src/stage3/03_phase2_embeddings.py
python3 src/stage3/04_phase2_gnn.py
```

**Output:** `outputs/stage3/`

---

## Validation

Scripts in `src/`. All read pre-computed expression and clinical data.

### validation_tcga.py
TCGA-BRCA progression-free survival (Cox PH, n=927).
Result: C1 HR = 0.915 [0.844–0.992], p = 0.031; C4 null (HR ≈ 1.00, p = 0.95).

```bash
python3 src/validation_tcga.py
```

**Input:** `outputs/validation_tcga/tcga_RSEM_gene_tpm.gz` (UCSC Xena TOIL RSEM,
707 MB -- download manually from UCSC Xena browser),
`outputs/validation_tcga/BRCA_clinical.tsv`,
`outputs/validation_tcga/BRCA_survival.txt`  
**Output:** `outputs/validation_tcga/tcga_results.txt`,
`outputs/validation_tcga/TCGA_BRCA_TLS_survival.png`

---

### validation_metabric.py
METABRIC overall survival (Cox PH, n=1,980).
Result: C1 HR = 0.855, p = 0.003; C4 null.

```bash
python3 src/validation_metabric.py
```

**Input:** `outputs/validation_tcga/METABRIC_clinical.tsv`,
`outputs/validation_tcga/METABRIC_expr_zscores.tsv` (fetched via cBioPortal API if absent)  
**Output:** `outputs/validation_tcga/metabric_results.txt`,
`outputs/validation_tcga/METABRIC_TLS_survival.png`

---

### validation_tcga_cesc.py
TCGA-CESC overall survival (Cox PH, n=304).
Result: C1 HR = 0.562 [0.345–0.918], p = 0.021.

```bash
python3 src/validation_tcga_cesc.py
```

**Output:** `outputs/validation_tcga/cesc_results.txt`

---

### validation_geo_shiao.py
Shiao et al. 2024 (GSE246613) scRNA-seq ICB response in TNBC (n=35,
pembrolizumab). Pre-treatment C1 scores predict pCR and RCB.
Result: C1 MWU p = 0.030; C1 ~ RCB ρ = −0.384, p = 0.040.

```bash
python3 src/validation_geo_shiao.py
```

**Input:** `outputs/validation_geo/GSE246613_immune.h5ad`  
**Output:** `outputs/validation_geo/shiao_C1C4_pseudobulk.tsv`,
`outputs/validation_geo/shiao_TLS_ICB.png`

---

## Manuscript

```bash
cd paper
source ../.venv/bin/activate
python3 assemble_figures.py        # writes paper/figures/fig1_overview.png, fig2_clustering.png
bash build_manuscript.sh           # pandoc + tectonic -> manuscript.pdf
```

Supplementary figures are generated by:
- `/tmp/make_s1_fig.py` -- Fig S1 (725-instance exclusion UMAP)
- `/tmp/make_platform_confound_fig.py` -- Fig S2 (platform distribution per cluster)

---

## Key Canonical Outputs

| File | Rows | Description |
|------|------|-------------|
| `outputs/tls_results_cc354.tsv` | 497 samples | Raw HookNet polygon counts, cc≥354 filter |
| `outputs/tls_spot_coverage.tsv` | 191 samples | Instances + spot counts after watershed split |
| `outputs/tls_detection_summary.tsv` | **747** | All instances × metadata; `exclusion` column marks BC23803/BC23377 |
| `outputs/stage2_clustering/tls_neighborhood_spots.h5ad` | ~64k spots | ST spots for 119 samples, 40,270 genes |
| `outputs/stage2_clustering/instance_features_725.tsv` | 725 | Z-scored 12-feature matrix (post-BC23803) |
| `outputs/stage2_clustering/instance_features_670.tsv` | **670** | Z-scored 12-feature matrix (canonical) |
| `outputs/stage2_clustering/instance_features_raw_spatial.tsv` | 670 | 22-feature matrix with morphology + gradients (unscaled) |
| `outputs/stage2_clustering/leiden_670_r0.3.tsv` | 670 | Cluster labels, 6 clusters, resolution 0.3 |
| `outputs/stage2_clustering/sample_patient_map.tsv` | 119 | Sample -> patient map |
| `outputs/validation_tcga/tcga_results.txt` | -- | TCGA-BRCA Cox results |
| `outputs/validation_tcga/metabric_results.txt` | -- | METABRIC Cox results |
| `outputs/validation_tcga/cesc_results.txt` | -- | TCGA-CESC Cox results |
| `outputs/validation_geo/shiao_C1C4_pseudobulk.tsv` | 35 | Shiao pseudobulk C1/C4 scores |

---

## Instance Exclusion Audit

| Stage | Count | Reason |
|-------|-------|--------|
| HookNet detections (cc≥354) | 860 polygons / 185 samples | -- |
| After lymph node exclusion | − 52 samples | Organised secondary lymphoid organs |
| After duplicate WSI exclusion | − 13 samples | Same physical slide fingerprint |
| After watershed split + ≥1 spot | **747 instances / 119 samples** | Stage 1 final output |
| After BC23803 exclusion | 725 instances / 116 samples | SPA51/52/53: single-patient C5 dominance |
| After BC23377 exclusion | **670 instances / 113 samples** | SPA54/55/56: single-patient C5 dominance |
