#!/bin/bash
# 04_batch_pipeline.sh
# Orchestrator for Stage 1 TLS detection pipeline.
#
# Steps (run in order, each can be examined before the next):
#   1. Inference  : download WSIs -> preprocess (Docker) -> HookNet-TLS (GPU) -> cleanup
#   2. Refilter   : apply literature-grounded thresholds; flag duplicates + Xenium
#   3. Coverage   : download ST h5ad; count ST spots inside each TLS polygon (watershed split)
#
# By default only inference runs -- examine results before chaining further steps.
# Use flags to chain steps or run individual steps in isolation.
#
# Usage:
#   bash 04_batch_pipeline.sh [options]
#
# Inference options:
#   --manifest FILE     Run manifest from 00_make_manifest.sh (recommended for full dataset)
#   --tier TIERS        Comma-separated priority tiers to run (default: high)
#                       e.g. --tier high,medium  or  --tier xenium
#   --batch-size N      Samples per download chunk (default: 10)
#   --keep-wsi          Keep WSI + ASAP TIFF after inference (no cleanup)
#   --offset N          Skip first N pending samples (for manual resuming)
#   --dry-run           Print plan, change nothing
#
# Step-control options (default: inference only):
#   --refilter          Also run refilter after inference
#   --coverage          Also run refilter + ST download + spot coverage after inference
#   --end-to-end        Alias for --coverage (all steps, no stops)
#   --refilter-only     Skip inference; run refilter only
#   --coverage-only     Skip inference + refilter; run ST download + spot coverage only
#
# Refilter options:
#   --refilter-cc N     confidence_count for TLS/GC filtering (default: 354 ≈ 200 µm min)
#   --results FILE      Raw inference results TSV (default: outputs/tls_results.tsv)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT="$(cd "$SCRIPT_DIR/../.." && pwd)"
DATA_DIR="$PROJECT/data"
OUT_BASE="$PROJECT/outputs/tls_masks"
VENV_PY="$PROJECT/.venv/bin/python3"
RESULTS_TSV="$PROJECT/outputs/tls_results.tsv"
MANIFEST_TSV=""       # if set, read sample list from manifest instead of HEST CSV
MANIFEST_TIERS="high" # comma-separated tiers to run from manifest
BATCH_SIZE=10
KEEP_WSI=0
DRY_RUN=0
OFFSET=0
SKIP_INFER=0
RUN_REFILTER=0
RUN_COVERAGE=0
REFILTER_CC=354

# --- arg parsing ---
while [[ $# -gt 0 ]]; do
    case $1 in
        --batch-size)    BATCH_SIZE="$2";       shift 2 ;;
        --keep-wsi)      KEEP_WSI=1;            shift ;;
        --dry-run)       DRY_RUN=1;             shift ;;
        --results)       RESULTS_TSV="$2";      shift 2 ;;
        --offset)        OFFSET="$2";           shift 2 ;;
        --refilter-cc)   REFILTER_CC="$2";      shift 2 ;;
        --manifest)      MANIFEST_TSV="$2";     shift 2 ;;
        --tier)          MANIFEST_TIERS="$2";   shift 2 ;;
        --refilter)      RUN_REFILTER=1;        shift ;;
        --coverage)      RUN_REFILTER=1; RUN_COVERAGE=1; shift ;;
        --end-to-end)    RUN_REFILTER=1; RUN_COVERAGE=1; shift ;;
        --refilter-only) SKIP_INFER=1; RUN_REFILTER=1;  shift ;;
        --coverage-only) SKIP_INFER=1; RUN_COVERAGE=1;  shift ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

REFILTER_TSV="$PROJECT/outputs/tls_results_cc${REFILTER_CC}.tsv"

# --- helpers ---
log() { echo "[$(date '+%H:%M:%S')] $*"; }

count_polygons() {
    local json_file="$1"
    if [ ! -f "$json_file" ]; then echo 0; return; fi
    "$VENV_PY" -c "
import json, sys
try:
    d = json.load(open('$json_file'))
    feats = d.get('features', d) if isinstance(d, dict) else d
    print(len(feats))
except: print(0)
"
}

download_batch() {
    # Download WSI + metadata for a list of sample IDs.
    # Args: sid1 sid2 ...
    local ids=("$@")
    log "Downloading ${#ids[@]} samples: ${ids[*]}"
    bash "$SCRIPT_DIR/00_download_data.sh" \
        --ids "${ids[@]}" --types wsis,metadata --dir "$DATA_DIR"
}

append_result() {
    local sid="$1" code="$2" tissue="$3" n_tls="$4" n_gc="$5" status="$6"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$sid" "$code" "$tissue" "$n_tls" "$n_gc" "$status" >> "$RESULTS_TSV"
}

# ============================================================
# STEP 1 -- INFERENCE
# ============================================================
if [ "$SKIP_INFER" -eq 1 ]; then
    log "Skipping inference (--refilter-only or --coverage-only)"
else

# --- init results TSV ---
mkdir -p "$(dirname "$RESULTS_TSV")"
if [ ! -f "$RESULTS_TSV" ]; then
    printf 'sample_id\toncotree_code\ttissue\tn_tls\tn_gc\tstatus\n' > "$RESULTS_TSV"
    log "Created results file: $RESULTS_TSV"
fi

# --- get pending sample list ---
if [ -n "$MANIFEST_TSV" ]; then
    log "Reading cohort from manifest: $MANIFEST_TSV (tiers: $MANIFEST_TIERS)"
else
    log "Reading cohort from HEST CSV (no manifest -- using breast+melanoma default)"
fi
PENDING_JSON=$("$VENV_PY" - "$MANIFEST_TSV" "$MANIFEST_TIERS" <<PYEOF
import pandas as pd, os, json, sys

manifest_path = sys.argv[1]
tiers_arg     = sys.argv[2]

if manifest_path:
    mf = pd.read_csv(manifest_path, sep='\t', dtype=str)
    tiers = [t.strip() for t in tiers_arg.split(',')]
    # xenium samples run H&E detection -- include when tier explicitly requested
    sub = mf[mf['priority_tier'].isin(tiers)][['sample_id','oncotree_code','tissue']].copy()
    sub = sub.rename(columns={'sample_id': 'id'})
else:
    # Legacy fallback: breast + melanoma only
    df = pd.read_csv('$DATA_DIR/HEST_v1_3_0.csv', low_memory=False)
    mask = (
        df['tissue'].str.lower().str.contains('breast', na=False) |
        df['oncotree_code'].isin(['SKCM', 'MEL'])
    ) & (df['disease_state'].str.lower() == 'cancer')
    sub = df[mask][['id','oncotree_code','tissue']].copy()

out_base = '$OUT_BASE'
pending = []
for _, row in sub.iterrows():
    sid = row['id']
    # Support both old (filtered/) and new (images/filtered/) output layout
    tls_json_new = f"{out_base}/{sid}/images/filtered/{sid}_asap_hooknettls_tls_filtered.json"
    tls_json_old = f"{out_base}/{sid}/filtered/{sid}_asap_hooknettls_tls_filtered.json"
    if os.path.exists(tls_json_new) or os.path.exists(tls_json_old):
        continue
    pending.append({
            'id': sid,
            'oncotree_code': str(row.get('oncotree_code', '')),
            'tissue': str(row.get('tissue', '')),
        })

print(json.dumps(pending))
PYEOF
)

# Parse into parallel arrays
mapfile -t PENDING_IDS < <(echo "$PENDING_JSON" | "$VENV_PY" -c "
import json, sys
data = json.load(sys.stdin)
for d in data: print(d['id'])
")
mapfile -t PENDING_CODES < <(echo "$PENDING_JSON" | "$VENV_PY" -c "
import json, sys
data = json.load(sys.stdin)
for d in data: print(d['oncotree_code'])
")
mapfile -t PENDING_TISSUES < <(echo "$PENDING_JSON" | "$VENV_PY" -c "
import json, sys
data = json.load(sys.stdin)
for d in data: print(d['tissue'])
")

TOTAL=${#PENDING_IDS[@]}
log "Pending: $TOTAL samples (already done: skipped)"

if [ "$DRY_RUN" -eq 1 ]; then
    log "[DRY RUN] Would process ${PENDING_IDS[*]}"
    log "[DRY RUN] Batch size: $BATCH_SIZE | Keep WSI: $KEEP_WSI"
    exit 0
fi

if [ "$TOTAL" -eq 0 ]; then
    log "All samples already inferred. Nothing to do."
    exit 0
fi

# Apply offset
if [ "$OFFSET" -gt 0 ]; then
    PENDING_IDS=("${PENDING_IDS[@]:$OFFSET}")
    PENDING_CODES=("${PENDING_CODES[@]:$OFFSET}")
    PENDING_TISSUES=("${PENDING_TISSUES[@]:$OFFSET}")
    log "Offset $OFFSET applied; ${#PENDING_IDS[@]} samples remaining"
fi

# --- main loop ---
BATCH_NUM=0
IDX=0
TOTAL_REMAINING=${#PENDING_IDS[@]}
START_TIME=$(date +%s)
N_DONE=0

while [ "$IDX" -lt "${#PENDING_IDS[@]}" ]; do
    BATCH_NUM=$((BATCH_NUM + 1))
    END=$((IDX + BATCH_SIZE))
    if [ "$END" -gt "${#PENDING_IDS[@]}" ]; then END="${#PENDING_IDS[@]}"; fi

    BATCH_IDS=("${PENDING_IDS[@]:$IDX:$((END - IDX))}")
    BATCH_CODES=("${PENDING_CODES[@]:$IDX:$((END - IDX))}")
    BATCH_TISSUES=("${PENDING_TISSUES[@]:$IDX:$((END - IDX))}")

    log "=== Batch $BATCH_NUM: samples $((IDX + OFFSET + 1))–$((END + OFFSET)) of $((TOTAL + OFFSET)) ==="

    # Step 1: download batch
    DATA_DIR="$DATA_DIR" download_batch "${BATCH_IDS[@]}"

    # Step 2: preprocess all in batch (Docker, CPU)
    PREPROCESS_ARGS=()
    for SID in "${BATCH_IDS[@]}"; do
        WSI="$DATA_DIR/wsis/${SID}.tif"
        ASAP_OUT="$OUT_BASE/$SID/${SID}_asap.tif"
        if [ -f "$ASAP_OUT" ]; then
            log "  $SID: ASAP TIFF exists, skipping preprocess"
        elif [ ! -f "$WSI" ]; then
            log "  $SID: WSI missing, skipping"
        else
            PREPROCESS_ARGS+=("$SID")
        fi
    done

    if [ "${#PREPROCESS_ARGS[@]}" -gt 0 ]; then
        log "Preprocessing ${#PREPROCESS_ARGS[@]} samples..."
        bash "$SCRIPT_DIR/01_asap_preprocess.sh" "${PREPROCESS_ARGS[@]}"
    fi

    # Step 3: infer each sample sequentially (one GPU)
    for i in "${!BATCH_IDS[@]}"; do
        SID="${BATCH_IDS[$i]}"
        CODE="${BATCH_CODES[$i]}"
        TISSUE="${BATCH_TISSUES[$i]}"
        ASAP_TIF="$OUT_BASE/$SID/${SID}_asap.tif"
        TLS_JSON="$OUT_BASE/$SID/images/filtered/${SID}_asap_hooknettls_tls_filtered.json"
        GC_JSON="$OUT_BASE/$SID/images/filtered/${SID}_asap_hooknettls_gc_filtered.json"

        if [ ! -f "$ASAP_TIF" ]; then
            log "  $SID: no ASAP TIFF, skipping inference"
            append_result "$SID" "$CODE" "$TISSUE" 0 0 "preprocess_failed"
            continue
        fi

        T0=$(date +%s)
        log "  Inference: $SID ($CODE / $TISSUE)"

        if bash "$SCRIPT_DIR/03c_run_hooknet_native.sh" "$SID" 2>&1 | \
                tee -a "$PROJECT/outputs/batch_run.log"; then
            N_TLS=$(count_polygons "$TLS_JSON")
            N_GC=$(count_polygons "$GC_JSON")
            T1=$(date +%s)
            ELAPSED=$(( T1 - T0 ))
            log "  $SID done: TLS=$N_TLS GC=$N_GC (${ELAPSED}s)"
            append_result "$SID" "$CODE" "$TISSUE" "$N_TLS" "$N_GC" "ok"
            N_DONE=$((N_DONE + 1))
        else
            log "  $SID: inference FAILED"
            append_result "$SID" "$CODE" "$TISSUE" 0 0 "inference_failed"
        fi

        # Cleanup WSI + ASAP TIFF to recover ~1.8 GB
        if [ "$KEEP_WSI" -eq 0 ]; then
            rm -f "$DATA_DIR/wsis/${SID}.tif"
            rm -f "$ASAP_TIF"
            rm -f "$OUT_BASE/${SID}/${SID}_asap_mask.tif"
            log "  $SID: cleaned up WSI + ASAP TIFFs"
        fi

        # ETA estimate
        NOW=$(date +%s)
        ELAPSED_TOTAL=$(( NOW - START_TIME ))
        if [ "$N_DONE" -gt 0 ]; then
            PER_SAMPLE=$(( ELAPSED_TOTAL / N_DONE ))
            REMAINING=$(( TOTAL_REMAINING - IDX - i - 1 ))
            ETA=$(( REMAINING * PER_SAMPLE ))
            log "  Progress: $((N_DONE)) done / $TOTAL_REMAINING total | ~$((ETA / 60)) min remaining"
        fi
    done

    IDX="$END"
done

log "=== Inference complete ==="
log "Results: $RESULTS_TSV"
"$VENV_PY" - <<PYEOF
import pandas as pd, sys

df = pd.read_csv('$RESULTS_TSV', sep='\t')
print(f"\nTotal processed: {len(df)}")
print(f"  ok:               {(df.status=='ok').sum()}")
print(f"  inference_failed: {(df.status=='inference_failed').sum()}")
print(f"  preprocess_failed:{(df.status=='preprocess_failed').sum()}")

ok = df[df.status == 'ok']
if len(ok) == 0:
    print("No successful results yet.")
    sys.exit(0)

tls_pos = (ok.n_tls > 0).sum()
print(f"\nTLS-positive (n_tls > 0): {tls_pos} / {len(ok)} ({100*tls_pos/len(ok):.1f}%)")
print(f"  Breast IDC+ILC: {ok[ok.oncotree_code.isin(['IDC','ILC'])].shape[0]} samples")
print(f"  Melanoma SKCM+MEL: {ok[ok.oncotree_code.isin(['SKCM','MEL'])].shape[0]} samples")

print("\nTLS counts by oncotree_code (raw, cc=70 defaults -- run refilter for calibrated counts):")
print(ok.groupby('oncotree_code')[['n_tls','n_gc']].describe()[['n_tls','n_gc']].to_string())
PYEOF

fi  # end SKIP_INFER

# ============================================================
# STEP 2 -- REFILTER
# ============================================================
if [ "$RUN_REFILTER" -eq 1 ]; then
    log ""
    log "=== Step 2: Refilter (cc=${REFILTER_CC}) ==="
    log "    Examine outputs/tls_results_cc${REFILTER_CC}.tsv before proceeding to coverage."
    bash "$SCRIPT_DIR/05_refilter.sh" --tls-confidence-count "$REFILTER_CC"
fi

# ============================================================
# STEP 3 -- ST DOWNLOAD + SPOT COVERAGE
# ============================================================
if [ "$RUN_COVERAGE" -eq 1 ]; then
    if [ ! -f "$REFILTER_TSV" ]; then
        log "ERROR: $REFILTER_TSV not found -- run refilter first (--refilter or --refilter-only)"
        exit 1
    fi

    log ""
    log "=== Step 3a: Download ST data ==="
    bash "$SCRIPT_DIR/00_download_data.sh" \
        --from-tsv "$REFILTER_TSV" \
        --tls-positive-only \
        --types st \
        --dir "$DATA_DIR"

    log ""
    log "=== Step 3b: Spot coverage ==="
    bash "$SCRIPT_DIR/06_tls_spot_coverage.sh" \
        --results "$REFILTER_TSV"

    log ""
    log "=== Pipeline complete ==="
    log "    Inference:  $RESULTS_TSV"
    log "    Refiltered: $REFILTER_TSV"
    log "    Coverage:   $PROJECT/outputs/tls_spot_coverage.tsv"
fi
