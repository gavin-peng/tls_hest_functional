"""
Phase 1 morphological feature extraction for Stage 3 H&E-based TLS state prediction.

Computes per-instance features from two spatial regions:
  Region A (TLS interior): shape + haematoxylin intensity from in-TLS spot patches
  Region B (300 µm ring):  haematoxylin intensity from neighbourhood spot patches

Shape features (polygon geometry, no H&E needed : all 670 instances):
  log_area_um2   from instance_features_raw_spatial.tsv
  circularity    from instance_features_raw_spatial.tsv
  convexity      from instance_features_raw_spatial.tsv
  elongation     computed here from polygon bounding rectangle

H&E intensity features (from HEST spot patches, data/patches/{sample}.h5):
  he_mean_A      mean haematoxylin intensity, TLS interior
  he_var_A       variance, TLS interior
  he_mean_B      mean haematoxylin intensity, 300 µm ring
  he_var_B       variance, ring
  he_ratio_BA    ring/interior ratio : key C4 myeloid-ring signature

Matching strategy: h5ad stores coord_x_um / coord_y_um (µm space).
  Patch h5 stores pixel coords. Convert patches to µm using pixel_size_um_estimated
  from HEST metadata, then KD-tree nearest-neighbour match (tolerance 2 µm).

Output: outputs/stage3/phase1_features.tsv
  670 rows × features + cluster + metadata
  Rows without downloaded patches have NaN for H&E features.
"""
import os, json, warnings
import numpy as np
import pandas as pd
import anndata as ad
import h5py
from pathlib import Path
from scipy.spatial import cKDTree
from scipy import sparse as sp
from shapely.geometry import Polygon, MultiPolygon
from skimage.color import rgb2hed

warnings.filterwarnings('ignore')

ROOT       = Path('/mnt/e/hest')
OUT_DIR    = ROOT / 'outputs/stage3'
PATCH_DIR  = ROOT / 'data/patches'
MASK_DIR   = ROOT / 'outputs/tls_masks'
META_CSV   = ROOT / 'data/HEST_v1_3_0.csv'
SPATIAL_H5 = ROOT / 'outputs/stage2_clustering/tls_neighborhood_spots.h5ad'
CLUST_TSV  = ROOT / 'outputs/stage2_clustering/leiden_670_r0.3.tsv'
SHAPE_TSV  = ROOT / 'outputs/stage2_clustering/instance_features_raw_spatial.tsv'

COORD_MATCH_TOL_UM = 2.0   # µm tolerance for barcode->patch matching


def extract_he_channel(patches_rgb: np.ndarray) -> np.ndarray:
    """Return haematoxylin channel values for an array of RGB patches (N,H,W,3)."""
    results = []
    for patch in patches_rgb:
        rgb = patch.astype(np.float32) / 255.0
        rgb = np.clip(rgb, 1e-6, 1.0)
        hed = rgb2hed(rgb)
        results.append(hed[:, :, 0].ravel())   # channel 0 = haematoxylin
    return np.concatenate(results) if results else np.array([])


def he_stats(he_values: np.ndarray) -> dict:
    """Mean and variance of haematoxylin pixel values."""
    if len(he_values) == 0:
        return {'mean': np.nan, 'var': np.nan}
    return {'mean': float(np.mean(he_values)), 'var': float(np.var(he_values))}


def compute_elongation(polygon_coords_px: list) -> float:
    """Major/minor axis ratio from minimum rotated bounding rectangle."""
    try:
        poly = Polygon(polygon_coords_px)
        if not poly.is_valid:
            poly = poly.buffer(0)
        rect = poly.minimum_rotated_rectangle
        coords = list(rect.exterior.coords)
        sides = []
        for i in range(len(coords) - 1):
            dx = coords[i+1][0] - coords[i][0]
            dy = coords[i+1][1] - coords[i][1]
            sides.append(np.hypot(dx, dy))
        sides = sorted(sides)
        minor = sides[0] if sides[0] > 0 else 1.0
        major = sides[2] if len(sides) > 2 else sides[-1]
        return float(major / minor)
    except Exception:
        return np.nan


def load_polygon_coords(sample_id: str, instance_id: str) -> list | None:
    """Load polygon vertex coordinates (px) for an instance from filtered JSON."""
    coverage_path = MASK_DIR / sample_id / 'tls_spot_coverage.json'
    if not coverage_path.exists():
        return None
    with open(coverage_path) as f:
        coverage = json.load(f)

    inst_short = instance_id.split('_', 1)[1]   # strip sample prefix
    entry = next((e for e in coverage if str(e.get('instance_id')) == inst_short), None)
    if entry is None:
        return None
    polygon_index = entry.get('split_from') if entry.get('split_from') is not None else entry.get('index')

    filtered_path = MASK_DIR / sample_id / 'images' / 'filtered' / \
                    f'{sample_id}_asap_hooknettls_tls_filtered.json'
    if not filtered_path.exists():
        return None
    with open(filtered_path) as f:
        polygons = json.load(f)

    for poly in polygons:
        if poly.get('index') == polygon_index:
            coords = poly.get('coordinates', [])
            if coords:
                return coords[0] if isinstance(coords[0][0], list) else coords
    return None


def build_patch_lookup(sample_id: str, pixel_size_um: float):
    """Load patch h5 and build a barcode->image dict.

    patch h5 stores top-left corner in full-res pixels.
    We join to the original ST h5ad via barcode so coord conversion is exact.
    """
    patch_path = PATCH_DIR / f'{sample_id}.h5'
    st_path    = ROOT / 'data' / 'st' / f'{sample_id}.h5ad'
    if not patch_path.exists() or not st_path.exists():
        return None

    with h5py.File(patch_path, 'r') as f:
        imgs     = f['img'][:]      # (N, 224, 224, 3) uint8
        barcodes = [b[0].decode() if isinstance(b[0], bytes) else b[0]
                    for b in f['barcode']]

    # barcode -> patch image index
    bc_to_idx = {bc: i for i, bc in enumerate(barcodes)}

    # Build KD-tree on ST coord_um (same space as h5ad coord_x_um/coord_y_um)
    st = ad.read_h5ad(str(st_path), backed='r')
    st_coords_um = st.obsm['spatial'] * pixel_size_um   # (n_spots, 2)
    st_barcodes  = st.obs_names.tolist()
    tree = cKDTree(st_coords_um)

    return {
        'tree':        tree,
        'st_barcodes': st_barcodes,
        'bc_to_idx':   bc_to_idx,
        'imgs':        imgs,
    }


def query_patches(lookup: dict, query_coords_um: np.ndarray) -> np.ndarray | None:
    """Return patch images for spots at query_coords_um via exact barcode join."""
    if lookup is None or len(query_coords_um) == 0:
        return None
    # Match h5ad coord_um -> ST barcode (distance should be ~0)
    dists, idxs = lookup['tree'].query(query_coords_um, k=1)
    matched_imgs = []
    for dist, idx in zip(dists, idxs):
        if dist > COORD_MATCH_TOL_UM:
            continue
        bc = lookup['st_barcodes'][idx]
        patch_idx = lookup['bc_to_idx'].get(bc)
        if patch_idx is not None:
            matched_imgs.append(lookup['imgs'][patch_idx])
    return np.array(matched_imgs) if matched_imgs else None


# ── Load base data ─────────────────────────────────────────────────────────────
print("Loading clustering labels and shape features...")
clust = pd.read_csv(CLUST_TSV, sep='\t')
clust = clust.rename(columns={'leiden_cluster': 'cluster'})

shape = pd.read_csv(SHAPE_TSV, sep='\t',
                    usecols=['instance_id', 'log_area_um2', 'circularity', 'convexity'])

meta_hest = pd.read_csv(META_CSV, usecols=['id', 'pixel_size_um_estimated'])
meta_hest = meta_hest.rename(columns={'id': 'sample_id'}).set_index('sample_id')

print("Loading spot neighbourhood h5ad...")
spots = ad.read_h5ad(SPATIAL_H5, backed='r')
spots_obs = spots.obs[['sample_id', 'instance_id', 'spot_type',
                        'coord_x_um', 'coord_y_um']].copy()

# ── Compute elongation from polygon JSONs ──────────────────────────────────────
print("Computing elongation from polygon JSONs...")
elongation_rows = []
for _, row in clust.iterrows():
    coords = load_polygon_coords(row['sample_id'], row['instance_id'])
    elong = compute_elongation(coords) if coords else np.nan
    elongation_rows.append({'instance_id': row['instance_id'], 'elongation': elong})
elong_df = pd.DataFrame(elongation_rows)
n_ok = elong_df['elongation'].notna().sum()
print(f"  Elongation computed for {n_ok}/{len(elong_df)} instances")

# ── Extract H&E features per instance from patches ────────────────────────────
samples_needed = clust['sample_id'].unique()
samples_with_patches = [s for s in samples_needed if (PATCH_DIR / f'{s}.h5').exists()]
samples_missing = len(samples_needed) - len(samples_with_patches)
print(f"\nSamples with patches: {len(samples_with_patches)}/{len(samples_needed)}"
      f" ({samples_missing} missing : H&E features will be NaN)")

he_rows = []
for sid in samples_with_patches:
    pixel_size = meta_hest.loc[sid, 'pixel_size_um_estimated'] \
                 if sid in meta_hest.index else 0.5
    lookup = build_patch_lookup(sid, pixel_size)
    if lookup is None:
        continue

    sample_spots = spots_obs[spots_obs['sample_id'] == sid]
    instances = clust.loc[clust['sample_id'] == sid, 'instance_id'].values

    for inst_id in instances:
        inst_spots = sample_spots[sample_spots['instance_id'] == inst_id]
        coords_A = inst_spots.loc[inst_spots['spot_type'] == 'in_tls',
                                  ['coord_x_um', 'coord_y_um']].values
        coords_B = inst_spots.loc[inst_spots['spot_type'] == 'neighborhood',
                                  ['coord_x_um', 'coord_y_um']].values

        imgs_A = query_patches(lookup, coords_A)
        imgs_B = query_patches(lookup, coords_B)

        he_A = he_stats(extract_he_channel(imgs_A) if imgs_A is not None else np.array([]))
        he_B = he_stats(extract_he_channel(imgs_B) if imgs_B is not None else np.array([]))

        ratio = (he_B['mean'] / he_A['mean']
                 if (he_A['mean'] and he_B['mean'] and he_A['mean'] > 1e-6)
                 else np.nan)

        he_rows.append({
            'instance_id': inst_id,
            'he_mean_A':   he_A['mean'],
            'he_var_A':    he_A['var'],
            'he_mean_B':   he_B['mean'],
            'he_var_B':    he_B['var'],
            'he_ratio_BA': ratio,
        })

    print(f"  {sid}: {len(instances)} instances processed", flush=True)

he_df = pd.DataFrame(he_rows) if he_rows else pd.DataFrame(
    columns=['instance_id', 'he_mean_A', 'he_var_A', 'he_mean_B', 'he_var_B', 'he_ratio_BA'])

# ── Merge everything ───────────────────────────────────────────────────────────
print("\nMerging features...")
result = (clust[['instance_id', 'sample_id', 'oncotree_code', 'st_technology', 'cluster']]
          .merge(shape, on='instance_id', how='left')
          .merge(elong_df, on='instance_id', how='left')
          .merge(he_df, on='instance_id', how='left'))

OUT_DIR.mkdir(parents=True, exist_ok=True)
out_path = OUT_DIR / 'phase1_features.tsv'
result.to_csv(out_path, sep='\t', index=False)

n_he = result['he_mean_A'].notna().sum()
print(f"\nSaved {out_path}")
print(f"  {len(result)} instances total")
print(f"  {n_he} with H&E features ({len(result)-n_he} NaN : patches not yet downloaded)")
print(f"\nFeature columns: {[c for c in result.columns if c not in ['instance_id','sample_id','oncotree_code','st_technology','cluster']]}")
