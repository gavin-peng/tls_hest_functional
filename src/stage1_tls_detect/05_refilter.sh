#!/bin/bash
# 05_refilter.sh
# Re-run TLS/GC filtering on existing *_with_confidence.json files without
# re-running inference. Use this to tune confidence thresholds after a
# completed batch run.
#
# GC uses the same confidence_count as TLS by default (both structures require
# the same minimum area : GC detections in excess of TLS count are flagged as
# artifacts and excluded from training).
#
# Default confidence_count=354 corresponds to ~200 µm minimum TLS diameter (literature-grounded).
#
# Usage:
#   bash 05_refilter.sh [options]
#
# Options:
#   --tls-confidence-count N   Min TLS area sqrt-pixels (default: 70)
#   --gc-confidence-count N    Min GC area sqrt-pixels (default: same as TLS)
#   --confidence-value N       Min pixel confidence 0-255 (default: 150)
#   --results FILE             Output TSV (default: outputs/tls_results_cc<N>.tsv)
#   --dry-run                  Show what would run, change nothing

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT="$(cd "$SCRIPT_DIR/../.." && pwd)"
OUT_BASE="$PROJECT/outputs/tls_masks"
VENV_PY="$PROJECT/.venv_hooknet/bin/python3"
REPO_DIR="$PROJECT/models/pathology-hooknet-tls"
ORIG_RESULTS="$PROJECT/outputs/tls_results.tsv"

TLS_CC=354
GC_CC=""
CV=150
RESULTS_TSV=""
DRY_RUN=0

while [[ $# -gt 0 ]]; do
    case $1 in
        --tls-confidence-count) TLS_CC="$2"; shift 2 ;;
        --gc-confidence-count)  GC_CC="$2";  shift 2 ;;
        --confidence-value)     CV="$2";     shift 2 ;;
        --results)              RESULTS_TSV="$2"; shift 2 ;;
        --dry-run)              DRY_RUN=1;   shift ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# GC defaults to same threshold as TLS unless explicitly overridden
if [ -z "$GC_CC" ]; then
    GC_CC="$TLS_CC"
fi

if [ -z "$RESULTS_TSV" ]; then
    RESULTS_TSV="$PROJECT/outputs/tls_results_cc${TLS_CC}.tsv"
fi

log() { echo "[$(date '+%H:%M:%S')] $*"; }

TLS_AREA_UM2=$("$VENV_PY" -c "print(f'{int($TLS_CC)**2 * 0.25:.0f}')" 2>/dev/null || echo "?")
GC_AREA_UM2=$("$VENV_PY" -c "print(f'{int($GC_CC)**2 * 0.25:.0f}')" 2>/dev/null || echo "?")

log "Refilter settings:"
log "  TLS: confidence_count=$TLS_CC (min area ~${TLS_AREA_UM2} µm²), confidence_value=$CV"
log "  GC:  confidence_count=$GC_CC  (min area ~${GC_AREA_UM2} µm²), confidence_value=$CV"
log "  Output TSV: $RESULTS_TSV"

mapfile -t CONF_JSONS < <(find "$OUT_BASE" -name "*_with_confidence.json" | sort)
TOTAL=${#CONF_JSONS[@]}
log "Found $TOTAL samples with confidence JSON"

if [ "$DRY_RUN" -eq 1 ]; then
    log "[DRY RUN] Would refilter $TOTAL samples: cc_tls=$TLS_CC cc_gc=$GC_CC cv=$CV"
    exit 0
fi

printf 'sample_id\toncotree_code\ttissue\tn_tls\tn_gc\tstatus\tflag\n' > "$RESULTS_TSV"

# Pre-build Xenium ID set from HEST CSV (used to flag samples after filtering)
HEST_CSV="$PROJECT/data/HEST_v1_3_0.csv"

N_DONE=0
N_FAIL=0

for CONF_JSON in "${CONF_JSONS[@]}"; do
    # Path layout: .../tls_masks/{SID}/images/raw/*_with_confidence.json
    IMAGES_DIR="$(dirname "$(dirname "$CONF_JSON")")"
    SID="$(basename "$(dirname "$IMAGES_DIR")")"
    FILTERED_DIR="$IMAGES_DIR/filtered"
    BASE="${SID}_asap_hooknettls"

    TLS_FILTERED="$FILTERED_DIR/${BASE}_tls_filtered.json"
    GC_FILTERED="$FILTERED_DIR/${BASE}_gc_filtered.json"

    META=$(grep -m1 "^${SID}"$'\t' "$ORIG_RESULTS" 2>/dev/null || echo "")
    CODE=$(echo "$META" | cut -f2)
    TISSUE=$(echo "$META" | cut -f3)

    log "  $SID"

    if PYTHONPATH="$REPO_DIR" "$VENV_PY" - \
            "$CONF_JSON" "$TLS_FILTERED" "$GC_FILTERED" \
            "$TLS_CC" "$GC_CC" "$CV" <<'PYEOF' 2>&1
import sys
from pathlib import Path
from hooknettls.postprocessing import tls_filtering, gc_filtering

conf_json, tls_out, gc_out = sys.argv[1], sys.argv[2], sys.argv[3]
tls_cc, gc_cc, cv = int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])

Path(tls_out).parent.mkdir(parents=True, exist_ok=True)

tls_filtering(conf_json, tls_out, {"confidence_value": cv, "confidence_count": tls_cc})
gc_filtering(conf_json, tls_out, gc_out, {"confidence_value": cv, "confidence_count": gc_cc})
PYEOF
    then
        N_TLS=$(PYTHONPATH="$REPO_DIR" "$VENV_PY" -c "
import json; d=json.load(open('$TLS_FILTERED'))
feats=d.get('features',d) if isinstance(d,dict) else d; print(len(feats))" 2>/dev/null || echo 0)
        N_GC=$(PYTHONPATH="$REPO_DIR" "$VENV_PY" -c "
import json; d=json.load(open('$GC_FILTERED'))
feats=d.get('features',d) if isinstance(d,dict) else d; print(len(feats))" 2>/dev/null || echo 0)

        FLAG=""
        if [ "$N_GC" -gt "$N_TLS" ]; then
            FLAG="gc_exceeds_tls"
            log "    TLS=$N_TLS GC=$N_GC  *** FLAGGED: gc_exceeds_tls ***"
        else
            log "    TLS=$N_TLS GC=$N_GC"
        fi

        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$SID" "$CODE" "$TISSUE" "$N_TLS" "$N_GC" "ok" "$FLAG" >> "$RESULTS_TSV"
        N_DONE=$((N_DONE + 1))
    else
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$SID" "$CODE" "$TISSUE" 0 0 "filter_failed" "" >> "$RESULTS_TSV"
        log "    FAILED"
        N_FAIL=$((N_FAIL + 1))
    fi
done

log ""
log "=== Refilter complete: $N_DONE ok, $N_FAIL failed ==="
log "Results: $RESULTS_TSV"

# Flag Xenium samples : H&E detection and watershed run normally (step 1a/1b),
# but ST spot threshold (step 1c/1d) is not applicable; stored separately for Xenium validation.
/mnt/e/hest/.venv/bin/python3 - "$RESULTS_TSV" "$HEST_CSV" <<'PYEOF'
import pandas as pd, sys
tsv = pd.read_csv(sys.argv[1], sep='\t', dtype=str)
csv = pd.read_csv(sys.argv[2], low_memory=False)
xenium_ids = set(csv[csv['st_technology'] == 'Xenium']['id'])
n = tsv['sample_id'].isin(xenium_ids).sum()
tsv.loc[tsv['sample_id'].isin(xenium_ids), 'flag'] = 'xenium'
tsv.to_csv(sys.argv[1], sep='\t', index=False)
print(f"  Xenium flag applied to {n} samples")
PYEOF

# Detect duplicate WSI samples by polygon fingerprint.
# Identical frozenset(index, area_um2) across two samples = same physical slide.
# The later-seen sample (alphabetically) is flagged as duplicate_of the first.
/mnt/e/hest/.venv/bin/python3 - "$RESULTS_TSV" "$OUT_BASE" <<'PYEOF'
import json, sys, pandas as pd
from pathlib import Path
from shapely.geometry import Polygon

tsv_path, out_base = sys.argv[1], sys.argv[2]
tsv = pd.read_csv(tsv_path, sep='\t', dtype=str)
if 'duplicate_of' not in tsv.columns:
    tsv['duplicate_of'] = ''

HEATMAP_PX = 0.5
fingerprints = {}  # fp -> first sample_id seen

for sid in tsv['sample_id']:
    for candidate in [
        Path(out_base) / sid / 'images' / 'filtered' / f'{sid}_asap_hooknettls_tls_filtered.json',
        Path(out_base) / sid / 'filtered'        / f'{sid}_asap_hooknettls_tls_filtered.json',
    ]:
        if not candidate.exists():
            continue
        try:
            raw = json.loads(candidate.read_text())
            feats = raw.get('features', raw) if isinstance(raw, dict) else raw
            if not feats:
                break
            fp_items = []
            for f in feats:
                coords = [(x * HEATMAP_PX, y * HEATMAP_PX) for x, y in f['coordinates']]
                try:
                    area = round(Polygon(coords).area)  # µm², integer for stable hash
                    fp_items.append((f['index'], area))
                except Exception:
                    pass
            fp = frozenset(fp_items)
            if not fp:
                break
            if fp in fingerprints:
                tsv.loc[tsv['sample_id'] == sid, 'duplicate_of'] = fingerprints[fp]
            else:
                fingerprints[fp] = sid
        except Exception:
            pass
        break

n_dups = (tsv['duplicate_of'].fillna('') != '').sum()
tsv.to_csv(tsv_path, sep='\t', index=False)
print(f"  Duplicates detected: {n_dups}")
if n_dups > 0:
    dups = tsv[tsv['duplicate_of'].fillna('') != ''][['sample_id','duplicate_of']]
    for _, r in dups.iterrows():
        print(f"    {r['sample_id']} -> duplicate of {r['duplicate_of']}")
PYEOF

/mnt/e/hest/.venv/bin/python3 - "$RESULTS_TSV" <<'PYEOF'
import pandas as pd, sys

df = pd.read_csv(sys.argv[1], sep='\t')
ok = df[df.status == 'ok']
xenium   = ok[ok.flag.fillna('') == 'xenium']
gc_flag  = ok[ok.flag.fillna('') == 'gc_exceeds_tls']
clean    = ok[ok.flag.fillna('') == '']

print(f"\nTotal: {len(df)} | ok: {len(ok)} | failed: {(df.status!='ok').sum()}")
print(f"Flagged (xenium):          {len(xenium)} : Xenium validation set (separate pipeline)")
print(f"Flagged (gc_exceeds_tls): {len(gc_flag)} : excluded from training")
print(f"Clean (training-eligible): {len(clean)}")

if len(clean) > 0:
    tls_pos = (clean.n_tls > 0).sum()
    print(f"\nTLS-positive (clean only): {tls_pos} / {len(clean)} ({100*tls_pos/len(clean):.1f}%)")
    print(f"\nTLS counts by oncotree_code : clean samples only (mean / median / max):")
    print(clean.groupby('oncotree_code')[['n_tls','n_gc']].agg(['mean','median','max']).to_string())

if len(gc_flag) > 0:
    print(f"\nFlagged samples (gc_exceeds_tls):")
    print(gc_flag[['sample_id','oncotree_code','tissue','n_tls','n_gc']].to_string(index=False))
PYEOF
