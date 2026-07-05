#!/bin/bash
# Run HookNet-TLS inference natively (no Docker) using the .venv_hooknet environment.
#
# Prerequisite: run 03b_setup_hooknet_native.sh once first.
#
# Usage:
#   bash 03c_run_hooknet_native.sh NCBI688
#   bash 03c_run_hooknet_native.sh NCBI776
#   bash 03c_run_hooknet_native.sh TENX14

set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../../models/pathology-hooknet-tls" && pwd)"
OUT_BASE="$SCRIPT_DIR/../../outputs/tls_masks"
VENV="$SCRIPT_DIR/../../.venv_hooknet"
PYTHON="$VENV/bin/python3"
SID="${1:-TENX68}"
INPUT_DIR="$OUT_BASE/$SID"

if [ ! -f "$INPUT_DIR/${SID}_asap.tif" ] || [ ! -f "$INPUT_DIR/${SID}_asap_mask.tif" ]; then
    echo "ERROR: Preprocessed TIFFs not found in $INPUT_DIR"
    echo "Run 01_asap_preprocess.sh $SID first."
    exit 1
fi

# Remove stale lock file if present (left by a previously killed run)
LOCK="$INPUT_DIR/images/${SID}_asaphooknettls.lock"
if [ -f "$LOCK" ] && [ ! -f "$INPUT_DIR/images/${SID}_asap_hooknettls.tif" ]; then
    echo "Removing stale lock file: $LOCK"
    rm -f "$LOCK"
fi

mkdir -p "$INPUT_DIR/images"

echo "=== HookNet-TLS inference on $SID (native GPU) ==="

# Run hooknet in its own process group (setsid) so that after a segfault we
# can kill every descendant : resource_tracker, spawn_main workers, etc. :
# with a single kill -- -PGID, regardless of their command-line patterns.
# Without this, orphaned workers inherit the tee pipe fd and deadlock the
# calling batch pipeline indefinitely.
setsid env PYTHONPATH="$REPO_DIR" TF_FORCE_GPU_ALLOW_GROWTH=true "$PYTHON" -m hooknettls \
    "hooknettls.default.image_path=$INPUT_DIR/${SID}_asap.tif" \
    "hooknettls.default.mask_path=$INPUT_DIR/${SID}_asap_mask.tif" \
    "hooknettls.default.output_folder=$INPUT_DIR/images/" \
    "hooknettls.default.execute_inference.tmp_folder=/tmp/hooknet_${SID}/" \
    "hooknettls.default.iterator.backend=asap" \
    "hooknettls.default.model.model_weights=$REPO_DIR/weights.h5" \
    "hooknettls.default.model_weights=$REPO_DIR/weights.h5" \
    "hooknettls.default.iterator.cpus=1" &
HOOKNET_PID=$!
wait $HOOKNET_PID || true
# Kill the entire process group (PGID = HOOKNET_PID when launched via setsid).
# This cleans up all orphaned workers that would otherwise hold the pipe open.
kill -- -"$HOOKNET_PID" 2>/dev/null || true
sleep 1

# Re-check whether output was actually written; exit non-zero if not.
if [ ! -f "$INPUT_DIR/images/${SID}_asap_hooknettls.tif" ]; then
    echo "ERROR: inference output missing for $SID (likely segfault in ASAP writer)"
    exit 1
fi

echo ""
echo "=== Results ==="
ls -lh "$INPUT_DIR/images/" 2>/dev/null || echo "No output yet."
