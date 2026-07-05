#!/bin/bash
# 00_download_data.sh
# Download HEST-1k data from HuggingFace. Always skips files that already exist.
#
# Usage:
#   ./00_download_data.sh [stage] [options]
#
# Stage (positional, sets default sample list and file types):
#   desktop    3 representative breast+melanoma samples, wsis+metadata (~0.66 GB)
#   pilot      3-sample pilot, all file types
#   1          full breast+melanoma cohort, wsis+metadata only
#   2          full cohort, wsis+metadata+st+cellvit_seg
#   all        full cohort, all file types
#
# Options:
#   --ids SID [SID ...]      Explicit sample IDs (overrides stage sample list)
#   --from-tsv FILE          Read sample_id column from a TSV file
#   --tls-positive-only      With --from-tsv: restrict to clean TLS-positive rows
#                            (status=ok, flag empty, n_tls>0)
#   --types TYPE[,TYPE ...]  File types to download. Comma-separated from:
#                              wsis, st, patches, metadata, thumbnails,
#                              patches_vis, spatial_plots, tissue_seg, cellvit_seg
#                            Default: wsis,metadata (for --ids/--from-tsv)
#                                     or determined by stage
#   --dir PATH               Output directory (default: ../../data)
#   --dry-run                Print what would be downloaded without downloading
#
# Examples:
#   # Desktop TLS pilot : 3 samples, wsi+metadata
#   ./00_download_data.sh desktop
#
#   # Download ST files for all clean TLS-positive samples
#   ./00_download_data.sh --from-tsv ../../outputs/tls_results_cc354.tsv \
#       --tls-positive-only --types st
#
#   # Download wsi+metadata for specific samples (used by 04_batch_pipeline.sh)
#   ./00_download_data.sh --ids SPA1 SPA2 SPA3 --types wsis,metadata --dir /path/to/data
#
#   # Full cohort ST for cluster Stage 2
#   ./00_download_data.sh 2 --dir /cluster/path/to/hest/data

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT="$(cd "$SCRIPT_DIR/../.." && pwd)"
VENV_PY="$PROJECT/.venv/bin/python3"

STAGE=""
EXPLICIT_IDS=()
FROM_TSV=""
TLS_POSITIVE_ONLY=0
TYPES=""
LOCAL="$PROJECT/data"
DRY_RUN=0

# First positional arg (if not a flag) is the stage
if [[ $# -gt 0 && "$1" != --* ]]; then
    STAGE="$1"
    shift
fi

while [[ $# -gt 0 ]]; do
    case $1 in
        --ids)
            shift
            while [[ $# -gt 0 && "$1" != --* ]]; do
                EXPLICIT_IDS+=("$1"); shift
            done ;;
        --from-tsv)         FROM_TSV="$2";  shift 2 ;;
        --tls-positive-only) TLS_POSITIVE_ONLY=1; shift ;;
        --types)            TYPES="$2";    shift 2 ;;
        --dir)              LOCAL="$2";    shift 2 ;;
        --dry-run)          DRY_RUN=1;     shift ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

STAGE="${STAGE}" \
FROM_TSV="${FROM_TSV}" \
TLS_POSITIVE_ONLY="${TLS_POSITIVE_ONLY}" \
TYPES="${TYPES}" \
LOCAL="${LOCAL}" \
DRY_RUN="${DRY_RUN}" \
"$VENV_PY" - "${EXPLICIT_IDS[@]}" <<'PYEOF'
import sys, os
from pathlib import Path
from huggingface_hub import hf_hub_download
import pandas as pd

REPO  = "MahmoodLab/hest"

stage             = os.environ.get('STAGE', '')
from_tsv          = os.environ.get('FROM_TSV', '')
tls_positive_only = os.environ.get('TLS_POSITIVE_ONLY', '0') == '1'
types_arg         = os.environ.get('TYPES', '')
local             = os.environ.get('LOCAL')
dry_run           = os.environ.get('DRY_RUN', '0') == '1'
explicit_ids      = sys.argv[1:]

ALL_TYPES = ["wsis", "st", "patches", "metadata", "thumbnails",
             "patches_vis", "spatial_plots", "tissue_seg", "cellvit_seg"]

PATTERNS = {
    "wsis":          lambda sid: [f"wsis/{sid}.tif"],
    "st":            lambda sid: [f"st/{sid}.h5ad"],
    "patches":       lambda sid: [f"patches/{sid}.h5"],
    "metadata":      lambda sid: [f"metadata/{sid}.json"],
    "thumbnails":    lambda sid: [f"thumbnails/{sid}_downscaled_fullres.jpeg"],
    "patches_vis":   lambda sid: [f"patches_vis/{sid}_patch_vis.jpg"],
    "spatial_plots": lambda sid: [f"spatial_plots/{sid}_spatial_plots.png"],
    "tissue_seg":    lambda sid: [f"tissue_seg/{sid}_contours.geojson",
                                   f"tissue_seg/{sid}_vis.jpg"],
    "cellvit_seg":   lambda sid: [f"cellvit_seg/{sid}_cellvit_seg.parquet"],
}

# Always download the metadata index
hf_hub_download(repo_id=REPO, repo_type="dataset",
                filename="HEST_v1_3_0.csv", local_dir=local)
csv = pd.read_csv(f"{local}/HEST_v1_3_0.csv", low_memory=False)

cohort_mask = (
    csv['tissue'].str.lower().str.contains('breast', na=False) |
    csv['oncotree_code'].isin(['SKCM', 'MEL'])
) & (csv['disease_state'].str.lower() == 'cancer')

# --- resolve sample list ---
if explicit_ids:
    ids = explicit_ids
    default_types = "wsis,metadata"
elif from_tsv:
    df = pd.read_csv(from_tsv, sep='\t')
    if tls_positive_only:
        df = df[(df.status == 'ok') &
                (df.flag.fillna('') == '') &
                (df.n_tls > 0)]
    ids = df['sample_id'].tolist()
    default_types = "st"
elif stage == 'desktop':
    ids = ["NCBI688", "NCBI776", "TENX14"]
    default_types = "wsis,metadata"
elif stage == 'pilot':
    ids = ["TENX68", "MEND90", "TENX71"]
    default_types = ",".join(ALL_TYPES)
elif stage == '1':
    ids = csv[cohort_mask]['id'].tolist()
    default_types = "wsis,metadata"
elif stage == '2':
    ids = csv[cohort_mask]['id'].tolist()
    default_types = "wsis,metadata,st,cellvit_seg"
elif stage == 'all':
    ids = csv[cohort_mask]['id'].tolist()
    default_types = ",".join(ALL_TYPES)
elif stage:
    print(f"Unknown stage: {stage}. Use desktop, pilot, 1, 2, or all.")
    sys.exit(1)
else:
    print("Specify a stage or --ids / --from-tsv.")
    sys.exit(1)

# --- resolve file types ---
types = [t.strip() for t in (types_arg or default_types).split(',') if t.strip()]
unknown = [t for t in types if t not in PATTERNS]
if unknown:
    print(f"Unknown file types: {unknown}. Valid: {list(PATTERNS)}")
    sys.exit(1)

# --- build file list ---
files = []
for sid in ids:
    for ftype in types:
        files.extend(PATTERNS[ftype](sid))

label = (f"--ids ({len(ids)} samples)" if explicit_ids else
         f"--from-tsv ({len(ids)} samples)" if from_tsv else
         f"stage={stage} ({len(ids)} samples)")
print(f"{label} | types={','.join(types)} | files={len(files)} | dir={local}")

if dry_run:
    for sid in ids[:5]:
        for ftype in types:
            for f in PATTERNS[ftype](sid):
                dest = Path(local) / f
                status = "EXISTS" if dest.exists() else "DOWNLOAD"
                print(f"  {status}  {f}")
    if len(ids) > 5:
        print(f"  ... and {len(ids)-5} more samples")
    sys.exit(0)

# --- download ---
ok = skip = err = 0
for fname in files:
    dest = Path(local) / fname
    if dest.exists():
        print(f"SKIP {fname} (exists)")
        skip += 1
        continue
    try:
        hf_hub_download(repo_id=REPO, repo_type="dataset",
                        filename=fname, local_dir=local)
        print(f"OK   {fname}")
        ok += 1
    except Exception as e:
        print(f"ERR  {fname}: {str(e)[:120]}")
        err += 1

print(f"Download: {ok} OK, {skip} skipped, {err} errors")
if err:
    sys.exit(1)
PYEOF
