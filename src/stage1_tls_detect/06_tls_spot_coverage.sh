#!/bin/bash
# 06_tls_spot_coverage.sh
# Spatial join: for each detected TLS polygon, count ST spots inside and compute area.
#
# Large polygons (> --large-threshold mm²) are split into instances via watershed
# on the TLS probability heatmap (heat1.tif, 0.5 µm/px):
#   1. Crop heatmap to polygon bounding box and rasterise polygon mask
#   2. Smooth with Gaussian (--sigma-um) to suppress noise
#   3. Find local maxima (TLS centres) separated by at least --min-dist-um
#   4. Watershed from those seeds : each basin = one TLS instance
#   5. Map ST spots to basins by nearest heatmap pixel
#
# Coordinate alignment:
#   TLS polygons : pixel coords in _asap.tif at 0.5 µm/px  -> µm = coord × 0.5
#   Heatmap      : same pixel space as _asap.tif (0.5 µm/px) : no conversion needed
#   ST spots     : obsm['spatial'] = (pxl_col, pxl_row) in original WSI fullres
#                  at pixel_size_um_estimated µm/px        -> µm = coord × spacing
#
# Per-instance JSON fields:
#   instance_id   sequential within sample (0-based)
#   index         original polygon index, or "N_K" for the K-th split of polygon N
#   area_um2      area in µm² (polygon area for small; watershed basin area for splits)
#   n_spots       ST spots inside the instance
#   keep          n_spots >= min_spots
#   large         True if this came from a polygon > large_threshold
#   split_from    parent polygon index (int) for split instances, else null
#
# Sample-level TSV columns:
#   n_tls         raw polygon count from filtered JSON
#   n_instances   total instances after splitting
#   n_ge20        instances with >= min_spots spots
#   n_lt20        instances with <  min_spots spots
#   n_large       original large polygons (before splitting)
#   n_split       instances produced by splitting (summed across all large polygons)
#
# Usage:
#   bash 06_tls_spot_coverage.sh [options]
#
# Options:
#   --results FILE          Refiltered results TSV (default: tls_results_cc354.tsv)
#   --min-spots N           Min spots to keep an instance for Visium/ST (default: 20)
#   --min-spots-hd N        Min bins to keep an instance for Visium HD (default: 3000)
#                           Rationale: 2µm bins vs 100µm spots; every cc=354-passing polygon
#                           has ~7,800+ bins; 3000 bins ≈ 125µm-diameter fragment floor.
#   --large-threshold MM2   Polygon area threshold in mm² for splitting (default: 1.0)
#   --sigma-um N            Gaussian smoothing sigma in µm (default: 100)
#   --min-dist-um N         Min distance between TLS centres in µm (default: 200)
#   --peak-threshold N      Min heatmap value (0-255) to count as TLS centre (default: 30)
#   --out FILE              Output TSV (default: outputs/tls_spot_coverage.tsv)
#   --dry-run               List samples that would be processed, then exit

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT="$(cd "$SCRIPT_DIR/../.." && pwd)"
OUT_BASE="$PROJECT/outputs/tls_masks"
VENV_PY="$PROJECT/.venv/bin/python3"
CSV="$PROJECT/data/HEST_v1_3_0.csv"
ST_DIR="$PROJECT/data/st"

RESULTS_TSV="$PROJECT/outputs/tls_results_cc354.tsv"
COVERAGE_TSV="$PROJECT/outputs/tls_spot_coverage.tsv"
MIN_SPOTS=20
MIN_SPOTS_HD=3000
LARGE_MM2=1.0
SIGMA_UM=100
MIN_DIST_UM=300
PEAK_THRESH=50
DRY_RUN=0

while [[ $# -gt 0 ]]; do
    case $1 in
        --results)          RESULTS_TSV="$2";  shift 2 ;;
        --min-spots)        MIN_SPOTS="$2";    shift 2 ;;
        --min-spots-hd)     MIN_SPOTS_HD="$2"; shift 2 ;;
        --large-threshold)  LARGE_MM2="$2";    shift 2 ;;
        --sigma-um)         SIGMA_UM="$2";     shift 2 ;;
        --min-dist-um)      MIN_DIST_UM="$2";  shift 2 ;;
        --peak-threshold)   PEAK_THRESH="$2";  shift 2 ;;
        --out)              COVERAGE_TSV="$2"; shift 2 ;;
        --dry-run)          DRY_RUN=1;         shift ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

log() { echo "[$(date '+%H:%M:%S')] $*"; }

mapfile -t SAMPLES < <("$VENV_PY" - "$RESULTS_TSV" <<'PYEOF'
import pandas as pd, sys
df = pd.read_csv(sys.argv[1], sep='\t')
# Include all OK, non-duplicate, TLS-positive samples regardless of technology.
# Xenium samples are included here : the technology branch happens at step 1c
# (ST spot counting), not at sample selection.
# Only gc_exceeds_tls is excluded at this stage (detection artifact, not technology issue).
ok = df[(df.status == 'ok') &
        (df.flag.fillna('') != 'gc_exceeds_tls') &
        (df.duplicate_of.fillna('') == '') &
        (df.n_tls > 0)]
for _, row in ok.iterrows():
    print(row['sample_id'])
PYEOF
)

TOTAL=${#SAMPLES[@]}
log "Clean TLS-positive unique samples: $TOTAL"
log "Split settings: large > ${LARGE_MM2} mm² | sigma=${SIGMA_UM} µm | min_dist=${MIN_DIST_UM} µm | peak_thresh=${PEAK_THRESH}"

if [ "$DRY_RUN" -eq 1 ]; then
    for SID in "${SAMPLES[@]}"; do
        H5AD="$ST_DIR/${SID}.h5ad"
        [ -f "$H5AD" ] && echo "  READY  $SID" || echo "  NO_ST  $SID"
    done
    exit 0
fi

printf 'sample_id\toncotree_code\ttissue\tst_technology\tn_tls\tn_instances\tn_ge20\tn_lt20\tn_large\tn_split\tstatus\n' \
    > "$COVERAGE_TSV"

N_PROCESSED=0
N_SKIPPED=0

for SID in "${SAMPLES[@]}"; do
    H5AD="$ST_DIR/${SID}.h5ad"
    HEATMAP="$OUT_BASE/$SID/images/${SID}_asap_hooknettls_heat1.tif"

    TLS_JSON=""
    NEW_PATH="$OUT_BASE/$SID/images/filtered/${SID}_asap_hooknettls_tls_filtered.json"
    OLD_PATH="$OUT_BASE/$SID/filtered/${SID}_asap_hooknettls_tls_filtered.json"
    if   [ -f "$NEW_PATH" ]; then TLS_JSON="$NEW_PATH"
    elif [ -f "$OLD_PATH" ]; then TLS_JSON="$OLD_PATH"
    fi

    META=$("$VENV_PY" - "$RESULTS_TSV" "$CSV" "$SID" <<'PYEOF'
import pandas as pd, sys
df  = pd.read_csv(sys.argv[1], sep='\t')
csv = pd.read_csv(sys.argv[2], low_memory=False)
row = df[df.sample_id == sys.argv[3]].iloc[0]
tech = csv[csv['id'] == sys.argv[3]]['st_technology'].iloc[0] \
       if sys.argv[3] in csv['id'].values else ''
print(row['oncotree_code'], row['tissue'], tech, row['n_tls'], sep='\t')
PYEOF
)
    CODE=$(echo "$META" | cut -f1)
    TISSUE=$(echo "$META" | cut -f2)
    TECH=$(echo "$META" | cut -f3)
    N_TLS=$(echo "$META" | cut -f4)

    if [ ! -f "$H5AD" ]; then
        log "  SKIP $SID : no h5ad"
        printf '%s\t%s\t%s\t%s\t%s\tNA\tNA\tNA\tNA\tNA\tno_st\n' \
            "$SID" "$CODE" "$TISSUE" "$TECH" "$N_TLS" >> "$COVERAGE_TSV"
        N_SKIPPED=$((N_SKIPPED + 1)); continue
    fi
    if [ -z "$TLS_JSON" ]; then
        log "  SKIP $SID : no TLS JSON"
        printf '%s\t%s\t%s\t%s\t%s\tNA\tNA\tNA\tNA\tNA\tno_tls_json\n' \
            "$SID" "$CODE" "$TISSUE" "$TECH" "$N_TLS" >> "$COVERAGE_TSV"
        N_SKIPPED=$((N_SKIPPED + 1)); continue
    fi

    log "  $SID (n_tls=$N_TLS)"

    RESULT=$("$VENV_PY" - \
        "$TLS_JSON" "$H5AD" "$CSV" "$SID" "$TECH" \
        "$OUT_BASE/$SID/tls_spot_coverage.json" \
        "$HEATMAP" \
        "$MIN_SPOTS" "$LARGE_MM2" "$SIGMA_UM" "$MIN_DIST_UM" "$PEAK_THRESH" \
        "$MIN_SPOTS_HD" <<'PYEOF'
import sys, json
import numpy as np
import anndata as ad
import pandas as pd
import tifffile
from shapely.geometry import Polygon, Point
from skimage.draw import polygon as sk_poly
from skimage.feature import peak_local_max
from skimage.segmentation import watershed
from scipy.ndimage import gaussian_filter
from pathlib import Path

tls_json_path = sys.argv[1]
h5ad_path     = sys.argv[2]
csv_path      = sys.argv[3]
sid           = sys.argv[4]
st_technology = sys.argv[5]   # e.g. "Visium", "Xenium", "Spatial Transcriptomics"
out_json      = sys.argv[6]
heatmap_path  = sys.argv[7]
min_spots     = int(sys.argv[8])
large_um2     = float(sys.argv[9]) * 1e6   # mm² -> µm²
sigma_um      = float(sys.argv[10])
min_dist_um   = float(sys.argv[11])
peak_thresh   = float(sys.argv[12])
min_spots_hd  = int(sys.argv[13])

# Step 1c branches on technology:
#   Xenium:     keep=None  : cell density incomparable to spot grids
#   Visium HD:  keep=(n >= min_spots_hd=3000) : 2µm bins; every cc=354 polygon has ~7800+ bins
#   All others: keep=(n >= min_spots=20)      : Visium/ST ~100-200µm spacing
is_xenium    = (st_technology == 'Xenium')
is_visium_hd = st_technology.startswith('Visium HD')

def keepable(n_sp):
    if is_xenium:
        return None
    if is_visium_hd:
        return n_sp >= min_spots_hd
    return n_sp >= min_spots

HEATMAP_PX = 0.5  # µm per pixel in _asap.tif / heatmap

# Pixel spacing for ST spots (same lookup regardless of technology)
csv = pd.read_csv(csv_path, low_memory=False)
spacing = float(csv[csv['id'] == sid].iloc[0]['pixel_size_um_estimated'])

# Load TLS polygons -> µm (raw coords × 0.5)
raw = json.load(open(tls_json_path))
feats = raw.get('features', raw) if isinstance(raw, dict) else raw
polygons = []
for f in feats:
    coords_um = [(x * HEATMAP_PX, y * HEATMAP_PX) for x, y in f['coordinates']]
    try:
        poly = Polygon(coords_um).buffer(0)
        if poly.is_valid and not poly.is_empty:
            polygons.append({'index': f['index'], 'polygon': poly,
                             'area_um2': poly.area})
    except Exception:
        pass

if not polygons:
    print('0\t0\t0\t0\t0')
    sys.exit(0)

# Load ST spots -> µm
adata = ad.read_h5ad(h5ad_path)
spots_um = adata.obsm['spatial'] * spacing   # (n_spots, 2): col=x, row=y

def spots_inside(poly, spots):
    minx, miny, maxx, maxy = poly.bounds
    bbox = ((spots[:, 0] >= minx) & (spots[:, 0] <= maxx) &
            (spots[:, 1] >= miny) & (spots[:, 1] <= maxy))
    cands = spots[bbox]
    inside = np.array([poly.contains(Point(x, y)) for x, y in cands], dtype=bool)
    return cands[inside]

def watershed_split(poly_um, inside_spots, heatmap):
    minx_um, miny_um, maxx_um, maxy_um = poly_um.bounds

    # Bounding box in heatmap pixel coords (col=x, row=y)
    minx_px = max(0, int(minx_um / HEATMAP_PX))
    miny_px = max(0, int(miny_um / HEATMAP_PX))
    maxx_px = min(heatmap.shape[1], int(np.ceil(maxx_um / HEATMAP_PX)) + 1)
    maxy_px = min(heatmap.shape[0], int(np.ceil(maxy_um / HEATMAP_PX)) + 1)

    crop_full = heatmap[miny_px:maxy_px, minx_px:maxx_px].astype(np.float32)
    if crop_full.size == 0:
        return None

    # Downsample 4× so watershed stays fast even on slide-spanning polygons.
    # All µm-based parameters are unchanged; only pixel arithmetic scales.
    DS = 4
    from skimage.transform import downscale_local_mean
    crop = downscale_local_mean(crop_full, (DS, DS)).astype(np.float32)
    H, W = crop.shape
    WS_PX = HEATMAP_PX * DS  # µm per pixel in downsampled space (2.0 µm/px)

    # Rasterise polygon into downsampled crop pixel space
    ext = np.array(poly_um.exterior.coords)
    col_px = ((ext[:, 0] - minx_px * HEATMAP_PX) / WS_PX).clip(0, W - 1)
    row_px = ((ext[:, 1] - miny_px * HEATMAP_PX) / WS_PX).clip(0, H - 1)
    rr, cc = sk_poly(row_px, col_px, (H, W))
    mask = np.zeros((H, W), dtype=bool)
    mask[rr, cc] = True

    # Smooth heatmap inside polygon, find local maxima
    sigma_px = sigma_um / WS_PX
    smoothed = gaussian_filter(np.where(mask, crop, 0.0), sigma=sigma_px)
    min_dist_px = max(1, int(min_dist_um / WS_PX))
    peaks = peak_local_max(smoothed, min_distance=min_dist_px,
                           threshold_abs=peak_thresh, labels=mask)

    if len(peaks) <= 1:
        return None  # unsplittable : caller keeps as single instance

    # Watershed from peak seeds
    markers = np.zeros((H, W), dtype=int)
    for i, (r, c) in enumerate(peaks, start=1):
        markers[r, c] = i
    seg = watershed(-smoothed, markers=markers, mask=mask)

    # Map inside_spots to watershed basins via nearest downsampled pixel
    instances = []
    if len(inside_spots) > 0:
        s_col = np.clip(((inside_spots[:, 0] - minx_px * HEATMAP_PX) / WS_PX).astype(int), 0, W - 1)
        s_row = np.clip(((inside_spots[:, 1] - miny_px * HEATMAP_PX) / WS_PX).astype(int), 0, H - 1)
        spot_labels = seg[s_row, s_col]
    else:
        spot_labels = np.array([], dtype=int)

    for lab in range(1, len(peaks) + 1):
        area_um2 = float((seg == lab).sum()) * (WS_PX ** 2)
        n_spots  = int((spot_labels == lab).sum())
        instances.append({'area_um2': round(area_um2, 1), 'n_spots': n_spots})

    return instances

# Load heatmap once (only if any large polygon exists)
has_large = any(p['area_um2'] > large_um2 for p in polygons)
heatmap = None
if has_large and Path(heatmap_path).exists():
    raw_tif = tifffile.imread(heatmap_path)
    # Flatten to 2D: pyramidal TIFFs may load as (levels, H, W) or (H, W)
    while raw_tif.ndim > 2:
        raw_tif = raw_tif[0]
    heatmap = raw_tif

results = []
instance_id = 0
n_large = 0
n_split = 0

for p in polygons:
    inside = spots_inside(p['polygon'], spots_um)

    if p['area_um2'] <= large_um2:
        n_sp = len(inside)
        results.append({
            'instance_id':  instance_id,
            'index':        p['index'],
            'area_um2':     round(p['area_um2'], 1),
            'n_spots':      n_sp,
            'keep':         keepable(n_sp),
            'large':        False,
            'split_from':   None,
        })
        instance_id += 1

    else:
        n_large += 1
        sub = watershed_split(p['polygon'], inside, heatmap) if heatmap is not None else None

        if sub is None or len(sub) <= 1:
            n_sp = len(inside)
            results.append({
                'instance_id':  instance_id,
                'index':        p['index'],
                'area_um2':     round(p['area_um2'], 1),
                'n_spots':      n_sp,
                'keep':         keepable(n_sp),
                'large':        True,
                'split_from':   None,
            })
            instance_id += 1
        else:
            for j, inst in enumerate(sub):
                n_sp = inst['n_spots']
                results.append({
                    'instance_id':  instance_id,
                    'index':        f"{p['index']}_{j}",
                    'area_um2':     inst['area_um2'],
                    'n_spots':      n_sp,
                    'keep':         keepable(n_sp),
                    'large':        True,
                    'split_from':   p['index'],
                })
                instance_id += 1
                n_split += 1

Path(out_json).parent.mkdir(parents=True, exist_ok=True)
json.dump(results, open(out_json, 'w'), indent=2)

n_instances = len(results)
# keep=None for Xenium : count separately so TSV columns stay comparable
n_ge  = sum(1 for r in results if r['keep'] is True)
n_lt  = sum(1 for r in results if r['keep'] is False)
n_na  = sum(1 for r in results if r['keep'] is None)
print(f'{n_instances}\t{n_ge}\t{n_lt}\t{n_na}\t{n_large}\t{n_split}')
PYEOF
)

    N_INST=$(echo "$RESULT"  | cut -f1)
    N_GE=$(echo "$RESULT"    | cut -f2)
    N_LT=$(echo "$RESULT"    | cut -f3)
    N_NA=$(echo "$RESULT"    | cut -f4)
    N_LARGE=$(echo "$RESULT" | cut -f5)
    N_SPLIT=$(echo "$RESULT" | cut -f6)

    if [ "$N_NA" -gt 0 ]; then
        log "    -> $N_INST instances ($N_LARGE large, $N_SPLIT from splits) | n_spots counted (Xenium : threshold N/A)"
    elif echo "$TECH" | grep -q "Visium HD"; then
        log "    -> $N_INST instances ($N_LARGE large, $N_SPLIT from splits) | ≥${MIN_SPOTS_HD} bins (Visium HD): $N_GE"
    else
        log "    -> $N_INST instances ($N_LARGE large, $N_SPLIT from splits) | ≥${MIN_SPOTS} spots: $N_GE"
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$SID" "$CODE" "$TISSUE" "$TECH" "$N_TLS" \
        "$N_INST" "$N_GE" "$N_LT" "$N_LARGE" "$N_SPLIT" "ok" >> "$COVERAGE_TSV"
    N_PROCESSED=$((N_PROCESSED + 1))
done

log ""
log "=== Coverage complete: $N_PROCESSED processed, $N_SKIPPED skipped ==="
log "Results: $COVERAGE_TSV"

"$PROJECT/.venv/bin/python3" - "$COVERAGE_TSV" "$MIN_SPOTS" <<'PYEOF'
import pandas as pd, sys

df  = pd.read_csv(sys.argv[1], sep='\t')
ok  = df[df.status == 'ok']
print(f"\nProcessed: {len(ok)} | Skipped: {(df.status != 'ok').sum()}")
if len(ok) == 0:
    sys.exit(0)

vis = ok[ok.st_technology != 'Xenium']
xen = ok[ok.st_technology == 'Xenium']

print(f"\n--- Visium / Spatial Transcriptomics ({len(vis)} samples) ---")
print(f"Raw TLS polygons:        {vis.n_tls.sum()}")
print(f"Instances after split:   {vis.n_instances.sum()}")
print(f"  From large polygons:   {vis.n_large.sum()} large -> {vis.n_split.sum()} split instances")
print(f"  Keep (≥{sys.argv[2]} spots): {vis.n_ge20.sum()}")
print(f"  Drop (<{sys.argv[2]} spots): {vis.n_lt20.sum()}")

if len(xen) > 0:
    print(f"\n--- Xenium ({len(xen)} samples : threshold N/A, stored for validation) ---")
    print(f"Raw TLS polygons:        {xen.n_tls.sum()}")
    print(f"Instances after split:   {xen.n_instances.sum()}")
    print(f"  From large polygons:   {xen.n_large.sum()} large -> {xen.n_split.sum()} split instances")
    print(f"  n_spots counted (Xenium cells, not comparable to Visium spots)")
PYEOF
