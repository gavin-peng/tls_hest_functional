#!/bin/bash
# 00_make_manifest.sh
# Pre-flight metadata scan: filter all 1,276 HEST samples and assign a priority tier.
# Reads HEST_v1_3_0.csv only — no downloads, no GPU, completes in seconds.
#
# Priority tiers (written to run_manifest.tsv):
#   high      IDC, ILC, LUAD, LUSC, SKCM, MEL, COAD, READ, RCC — run first
#   medium    BLCA, STAD, HNSC, OV, UCEC, CSCC — run after high
#   low       PAAD, GBM, HCC, PRAD, EPM — run last or skip
#   xenium    Xenium / Xenium 5k — H&E detection only, then validation branch
#   visium_hd Visium HD — spot spacing 2 µm; ≥20-spot threshold needs review
#   inspect   Slide area outlier (>2 SD) or unknown cancer type — manual review first
#   skip      Quality filter failed, healthy tissue, or non-cancer disease state
#
# Usage:
#   bash 00_make_manifest.sh [--csv FILE] [--out FILE] [--min-spots N]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT="$(cd "$SCRIPT_DIR/../.." && pwd)"
HEST_CSV="$PROJECT/data/HEST_v1_3_0.csv"
OUT_TSV="$PROJECT/outputs/run_manifest.tsv"
MIN_SPOTS=1000

while [[ $# -gt 0 ]]; do
    case $1 in
        --csv)        HEST_CSV="$2"; shift 2 ;;
        --out)        OUT_TSV="$2";  shift 2 ;;
        --min-spots)  MIN_SPOTS="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

mkdir -p "$(dirname "$OUT_TSV")"

"$PROJECT/.venv/bin/python3" - "$HEST_CSV" "$OUT_TSV" "$MIN_SPOTS" <<'PYEOF'
import sys, pandas as pd, numpy as np

csv_path, out_path, min_spots = sys.argv[1], sys.argv[2], int(sys.argv[3])

df = pd.read_csv(csv_path, low_memory=False)

# ── Cancer-type priority tables ──────────────────────────────────────────────
HIGH   = {'IDC','ILC','LUAD','LUSC','SKCM','MEL','COAD','READ','KIRC','PRCC',
           'SCCRCC','RCC','COADREAD','CSCC'}
MEDIUM = {'BLCA','STAD','HNSC','HGSOC','SOC','OV','UCEC','LIHB','ACYC','CESC'}
LOW    = {'PAAD','PDAC','GBM','HCC','PRAD','EPM','LNET','ALL'}

# Organ-level inference for samples without oncotree_code.
# Only organs with a plausible cancer type are listed — anything not here and
# without an oncotree code gets skipped (non-cancer research tissue).
ORGAN_TIER = {
    'Breast':          'high',    # assume IDC
    'Lung':            'high',    # assume LUAD
    'Skin':            'high',    # assume SKCM
    'Colon':           'high',
    'Bowel':           'high',    # synonym for Colon/Rectum (COAD equivalent)
    'Rectum':          'high',
    'Kidney':          'high',
    'Bladder':         'medium',
    'Stomach':         'medium',
    'Head and Neck':   'medium',
    'Ovary':           'medium',
    'Uterus':          'medium',
    'Liver':           'low',
    'Pancreas':        'low',
    'Brain':           'low',
    'Prostate':        'low',
    'Cervix':          'medium',
}

# Disease states that indicate non-cancer / healthy (skip for TLS training)
SKIP_STATES = {'Healthy', 'Genetically modified'}

# Technology classification
XENIUM_TECHS  = {'Xenium', 'Xenium 5k'}
VIS_HD_TECHS  = {'Visium HD', "Visium HD 3'"}
GRID_TECHS    = {'Visium', 'Spatial Transcriptomics'} | VIS_HD_TECHS

# ── Slide area z-score ───────────────────────────────────────────────────────
df['area_mm2'] = (
    df['fullres_px_width']  * df['pixel_size_um_estimated'] / 1000 *
    df['fullres_px_height'] * df['pixel_size_um_estimated'] / 1000
)
area_mean = df['area_mm2'].mean()
area_std  = df['area_mm2'].std()
df['area_zscore'] = (df['area_mm2'] - area_mean) / area_std

# ── Normalised spot count ────────────────────────────────────────────────────
# Convert spots to Visium-equivalent (100 µm reference) to make the ≥1000
# threshold technology-agnostic.  inter_spot_dist missing -> use tech default.
tech_default_dist = {'Visium': 100, 'Spatial Transcriptomics': 200,
                     'Visium HD': 2, "Visium HD 3'": 2,
                     'Xenium': None, 'Xenium 5k': None}
def spot_dist(row):
    d = row.get('inter_spot_dist')
    if pd.notna(d) and d > 0:
        return float(d)
    return tech_default_dist.get(str(row.get('st_technology', '')), 100) or 100

df['inter_spot_dist_eff'] = df.apply(spot_dist, axis=1)
# Normalised = raw_spots × (dist/100)²  — cells-based techs skip this filter
df['spots_norm'] = df['spots_under_tissue'] * (df['inter_spot_dist_eff'] / 100) ** 2

# ── Assign tier ──────────────────────────────────────────────────────────────
rows = []
for _, r in df.iterrows():
    sid      = r['id']
    tech     = str(r.get('st_technology', ''))
    state    = str(r.get('disease_state', ''))
    code     = r.get('oncotree_code')        # may be NaN
    organ    = str(r.get('organ', ''))
    tissue   = str(r.get('tissue', ''))
    spots    = r['spots_under_tissue']
    norm_sp  = r['spots_norm']
    area_z   = r['area_zscore']
    area_mm2 = r['area_mm2']

    reasons = []
    species = str(r.get('species', ''))

    # ── Filter 0: species — skip non-human samples ──
    # HEST contains ~586 Mus musculus samples, some with human oncotree codes.
    # They must be excluded before any technology or cancer-type routing.
    if species != 'Homo sapiens':
        tier = 'skip'
        reasons.append(f'non_human({species})')

    # ── Filter 1: technology routing ──
    elif tech in XENIUM_TECHS:
        tier = 'xenium'

    elif (tech in VIS_HD_TECHS):
        tier = 'visium_hd'
        reasons.append('spot_threshold_review_needed')

    else:
        # ── Filter 2: disease state ──
        if state in SKIP_STATES:
            reasons.append('non_cancer')
            tier = 'skip'

        else:
            # ── Filter 2: spot quality (normalised) ──
            if pd.isna(norm_sp) or norm_sp < min_spots:
                reasons.append(f'low_spots({int(spots) if pd.notna(spots) else "?"})')
                tier = 'skip'

            else:
                # ── Filter 3: cancer type priority ──
                if pd.notna(code) and str(code) not in ('nan', 'UNKNOWN', 'TODO'):
                    code_str = str(code)
                    if code_str in HIGH:
                        tier = 'high'
                    elif code_str in MEDIUM:
                        tier = 'medium'
                    elif code_str in LOW:
                        tier = 'low'
                    else:
                        tier = 'inspect'
                        reasons.append(f'unknown_code({code_str})')
                else:
                    # Infer from organ
                    inferred = ORGAN_TIER.get(organ)
                    if inferred:
                        tier = inferred
                        reasons.append('tier_inferred_from_organ')
                    elif pd.notna(code) and str(code) not in ('nan','UNKNOWN','TODO',''):
                        # Has an unrecognised code in a cancer-plausible organ -> inspect
                        tier = 'inspect'
                        reasons.append(f'unknown_code({code})')
                    else:
                        # No oncotree code AND organ not in cancer list -> skip
                        # (spinal cord, heart, muscle, eye, etc.)
                        tier = 'skip'
                        reasons.append('no_cancer_organ')

                # ── Filter 4: slide size outlier ──
                if pd.notna(area_z) and area_z > 2.0:
                    reasons.append(f'large_slide({area_mm2:.0f}mm2,z={area_z:.1f})')
                    if tier in ('high', 'medium', 'low'):
                        tier = 'inspect'

    rows.append({
        'sample_id':       sid,
        'oncotree_code':   code if pd.notna(code) else '',
        'organ':           organ,
        'tissue':          tissue,
        'st_technology':   tech,
        'disease_state':   state,
        'priority_tier':   tier,
        'exclude_reason':  '|'.join(reasons),
        'n_spots':         int(spots) if pd.notna(spots) else 0,
        'n_spots_norm':    round(float(norm_sp), 0) if pd.notna(norm_sp) else 0,
        'nb_genes':        int(r['nb_genes']) if pd.notna(r.get('nb_genes')) else 0,
        'inter_spot_dist': r['inter_spot_dist_eff'],
        'area_mm2':        round(area_mm2, 1),
        'area_zscore':     round(area_z, 2),
        'fullres_px_w':    int(r['fullres_px_width']),
        'fullres_px_h':    int(r['fullres_px_height']),
        'pixel_size_um':   round(r['pixel_size_um_estimated'], 4),
    })

out = pd.DataFrame(rows)
out.to_csv(out_path, sep='\t', index=False)

# ── Summary ──────────────────────────────────────────────────────────────────
print(f"\nHEST-1k manifest: {len(out)} samples -> {out_path}")
print(f"\nPriority tier breakdown:")
tc = out.priority_tier.value_counts()
for tier in ['high','medium','low','xenium','visium_hd','inspect','skip']:
    n = tc.get(tier, 0)
    print(f"  {tier:12s} {n:4d}")

print(f"\nHigh-priority cancer types:")
hi = out[out.priority_tier == 'high']
code_counts = hi.oncotree_code.replace('', 'inferred').value_counts()
print(code_counts.to_string())

print(f"\nInspect samples (manual review needed):")
ins = out[out.priority_tier == 'inspect']
if len(ins) > 0:
    print(ins[['sample_id','oncotree_code','organ','exclude_reason']].to_string(index=False))

print(f"\nSkip summary — top exclude reasons:")
sk = out[out.priority_tier == 'skip']
from collections import Counter
reason_counts = Counter()
for r in sk.exclude_reason:
    for part in r.split('|'):
        reason_counts[part.split('(')[0]] += 1
for reason, n in reason_counts.most_common(10):
    print(f"  {reason}: {n}")

PYEOF

echo ""
echo "Manifest written: $OUT_TSV"
