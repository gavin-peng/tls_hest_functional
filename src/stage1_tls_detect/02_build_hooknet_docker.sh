#!/bin/bash
# Build the hooknet-preprocess Docker image for WSI preprocessing.
#
# Only the preprocessing step uses Docker. Inference runs natively via
# 03b_setup_hooknet_native.sh + 03c_run_hooknet_native.sh.
#
# Images built:
#   hooknet-preprocess  -- ASAP-based WSI preprocessing (CPU, Ubuntu 20.04)
#                          Converts any pyramidal TIFF to ASAP-native format at 0.5 µm/px.
#                          Modified from upstream to support --force-spacing for TIFFs
#                          that have no embedded pixel spacing metadata (e.g. HEST TIFFs).
#
# Upstream repo: https://github.com/DIAGNijmegen/pathology-hooknet-tls
# Weights:       https://zenodo.org/records/10614942 (weights.h5, 49 MB)

set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../../models/pathology-hooknet-tls" && pwd)"
PATCH="$SCRIPT_DIR/hooknet_hest_compat.patch"

# Clone repo if missing
if [ ! -d "$REPO_DIR" ]; then
    echo "=== Cloning pathology-hooknet-tls ==="
    git clone https://github.com/DIAGNijmegen/pathology-hooknet-tls.git "$REPO_DIR"
fi

# Apply HEST compatibility patch if not already applied
cd "$REPO_DIR"
if git diff --quiet; then
    echo "=== Applying HEST compatibility patch ==="
    git apply "$PATCH"
else
    echo "=== Patch already applied, skipping ==="
fi

# Download weights if missing
WEIGHTS="$REPO_DIR/weights.h5"
if [ ! -f "$WEIGHTS" ]; then
    echo "Downloading HookNet weights from Zenodo (~49 MB)..."
    wget -q --show-progress \
        "https://zenodo.org/records/10614942/files/weights.h5" \
        -O "$WEIGHTS"
fi

echo "=== Building hooknet-preprocess ==="
docker build -t hooknet-preprocess "$REPO_DIR/preprocessing/"

echo ""
echo "Built images:"
docker images | grep hooknet
