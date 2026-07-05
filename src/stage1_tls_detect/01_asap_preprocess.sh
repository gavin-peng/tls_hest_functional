#!/bin/bash
# Preprocess a WSI for HookNet-TLS using the ASAP Docker.
#
# Produces ASAP-native pyramidal TIFF + tissue mask, readable by the
# WholeSlideData ASAP backend inside the HookNet inference container.
#
# What it does:
#   1. Extracts the WSI level closest to 0.5 µm/px using ASAP and saves as
#      a new pyramidal TIFF with spacing metadata embedded.
#   2. Creates a tissue-background mask via Otsu thresholding.
#
# Requires:
#   - hooknet-preprocess Docker image (build with 02_build_hooknet_docker.sh)
#   - HEST_v1_3_0.csv for pixel spacing lookup (or pass --spacing manually)
#
# Usage:
#   ./01_asap_preprocess.sh TENX68              # reads spacing from HEST CSV
#   ./01_asap_preprocess.sh TENX68 0.2738       # explicit spacing
#   ./01_asap_preprocess.sh MEND90 TENX71       # multiple samples

set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
WSI_DIR="$SCRIPT_DIR/../../data/wsis"
OUT_BASE="$SCRIPT_DIR/../../outputs/tls_masks"
CSV="$SCRIPT_DIR/../../data/HEST_v1_3_0.csv"
VENV="$SCRIPT_DIR/../../.venv/bin/python3"

lookup_spacing() {
    local sid=$1
    "$VENV" -c "
import pandas as pd
df = pd.read_csv('$CSV')
row = df[df['id']=='$sid']
if row.empty: exit(1)
print(row.iloc[0]['pixel_size_um_estimated'])
" 2>/dev/null
}

run_sample() {
    local SID=$1
    local SPACING=${2:-}

    if [ -z "$SPACING" ]; then
        SPACING=$(lookup_spacing "$SID") || {
            echo "ERROR: Could not find $SID in HEST CSV. Pass spacing manually."
            return 1
        }
    fi

    OUT_DIR="$OUT_BASE/$SID"
    mkdir -p "$OUT_DIR"

    echo "=== $SID (spacing: ${SPACING} µm/px) ==="
    docker run --rm \
        -v "$WSI_DIR:/input" \
        -v "$OUT_DIR:/output" \
        hooknet-preprocess \
        "/input/${SID}.tif" \
        "/output/${SID}_asap.tif" \
        "/output/${SID}_asap_mask.tif" \
        0.5 0.02 "$SPACING"
    echo ""
}

# Parse args: if second arg looks like a float, treat it as spacing for first sample
if [ $# -eq 2 ] && echo "$2" | grep -qE '^[0-9]+\.[0-9]+$'; then
    run_sample "$1" "$2"
else
    for SID in "${@:-TENX68}"; do
        run_sample "$SID"
    done
fi
