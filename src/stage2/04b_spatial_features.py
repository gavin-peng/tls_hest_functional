#!/usr/bin/env python3
"""
Stage 2 : Step 2b: Extended per-TLS instance feature matrix with spatial features.

Extends Block 1+3 gene set scores with:
  - Polygon morphology: log_area_um2, circularity, convexity
  - Radial expression gradients for each Block 1 gene set (slope of score vs normalised
    distance from TLS spot centroid; negative = centrally concentrated, e.g. organised GC)

Excluded patients:
  - BC23803 (SPA51/52/53) : 22 instances, single-patient C2 artefact (original exclusion)
  - BC23377 (SPA54/55/56) : ~53 instances, single-patient C5 artefact

Polygon coordinate system: heatmap pixel units at 0.5 µm/px.
  area_um2 from coverage JSON = authoritative (pre-computed, correct for split instances).
  Circularity and convexity derived from parent polygon coordinates (scale-invariant).
  Split instances get NaN for circularity/convexity (no watershed boundary stored as polygon).

Radial gradient: uses spot centroid (mean of in_tls spot coords) as TLS centre.
  Distance normalised by sqrt(area_um2 / π) : effective TLS radius.
  Minimum 3 in_tls spots required; NaN otherwise (imputed with 0 in z-scoring step).

Reads:
  outputs/stage2_clustering/tls_neighborhood_spots.h5ad
  outputs/tls_masks/<sample_id>/tls_spot_coverage.json
  outputs/tls_masks/<sample_id>/images/filtered/<sample_id>_asap_hooknettls_tls_filtered.json

Writes:
  outputs/stage2_clustering/instance_features_raw_spatial.tsv  (raw, unscaled)

Usage:
  python3 src/stage2_labeling/04b_spatial_features.py
  python3 src/stage2_labeling/04b_spatial_features.py --min-spots 1
"""

import argparse
import glob
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
from scipy.stats import linregress
from shapely.geometry import Polygon as ShapelyPolygon

PROJECT = Path(__file__).resolve().parents[2]

EXCLUDED_SAMPLES = {'SPA51', 'SPA52', 'SPA53', 'SPA54', 'SPA55', 'SPA56'}

GENE_SETS = {
    'GC_score':     ['BCL6', 'AICDA', 'ELL3', 'MYBL1', 'CXCR4'],
    'Tfh_score':    ['CXCL13', 'CXCR5', 'ICOS', 'IL21', 'TOX2'],
    'Plasma_score': ['IGHG1', 'IGHG3', 'IGHA1', 'JCHAIN', 'MZB1'],
    'HEV_score':    ['ACKR1', 'SELP', 'CCL21', 'GLYCAM1', 'CXCL12'],
    'Suppressive':  ['FOXP3', 'IL2RA', 'CTLA4', 'TGFB1', 'IL10'],
    'Myeloid_supp': ['ARG1', 'IL10', 'CD274', 'IDO1', 'MRC1'],
    'Cytotoxic':    ['GZMB', 'PRF1', 'IFNG', 'GNLY', 'NKG7'],
}

NEIGHBOR_GENE_SETS = {
    'Neighbor_TGFB1':            ['TGFB1'],
    'Neighbor_VEGFA':            ['VEGFA'],
    'Neighbor_CD8A':             ['CD8A'],
    'Neighbor_tumor_fraction':   ['EPCAM', 'KRT8', 'KRT18', 'KRT19', 'CDH1'],
    'Neighbor_myeloid_fraction': ['CD68', 'CSF1R', 'CD14', 'LYZ', 'S100A8'],
}

GRADIENT_COLS = [f'Grad_{k.replace("_score", "").replace("_supp", "_supp")}' for k in GENE_SETS]

META_COLS = ['instance_id', 'sample_id', 'oncotree_code', 'n_spots_in_tls', 'st_technology']
FEATURE_COLS = (
    list(GENE_SETS.keys())
    + list(NEIGHBOR_GENE_SETS.keys())
    + ['log_area_um2', 'circularity', 'convexity']
    + GRADIENT_COLS
)


def log(msg):
    print(f"[04b_spatial] {msg}", flush=True)


def normalize_log1p_cpm(X_sparse, row_indices):
    if len(row_indices) == 0:
        return np.zeros((0, X_sparse.shape[1]), dtype=np.float32)
    X = X_sparse[row_indices].toarray().astype(np.float32)
    totals = X.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1.0
    return np.log1p(X / totals * 1e6)


def gene_set_score(X_norm, gene_names, gene_set):
    cols = [gene_names[g] for g in gene_set if g in gene_names]
    if not cols or X_norm.shape[0] == 0:
        return np.nan
    return float(X_norm[:, cols].mean())


def per_spot_scores(X_norm, gene_names, gene_set):
    """Return mean score per spot (n_spots,) or None if no genes available."""
    cols = [gene_names[g] for g in gene_set if g in gene_names]
    if not cols or X_norm.shape[0] == 0:
        return None
    return X_norm[:, cols].mean(axis=1)


def polygon_morphology(coords_px):
    """
    Compute circularity and convexity from polygon vertex coordinates (pixel units).
    Both measures are scale-invariant.
    Returns (circularity, convexity) or (nan, nan) on error.
    """
    try:
        poly = ShapelyPolygon(coords_px)
        if not poly.is_valid:
            poly = poly.buffer(0)
        area = poly.area
        perim = poly.exterior.length
        convex_area = poly.convex_hull.area
        circ = (4 * np.pi * area / perim ** 2) if perim > 0 else np.nan
        conv = (area / convex_area) if convex_area > 0 else np.nan
        return float(circ), float(conv)
    except Exception:
        return np.nan, np.nan


def radial_gradient(X_norm, spot_coords_um, gene_names, gene_set, area_um2):
    """
    Slope of per-spot gene-set score regressed on normalised distance from TLS centroid.
    normalised_distance = distance_um / sqrt(area_um2 / pi)
    Returns slope (negative = centrally concentrated) or NaN if < 3 spots.
    """
    scores = per_spot_scores(X_norm, gene_names, gene_set)
    if scores is None or len(scores) < 3:
        return np.nan
    centroid = spot_coords_um.mean(axis=0)
    radius = np.sqrt(area_um2 / np.pi) if area_um2 > 0 else 1.0
    dists = np.linalg.norm(spot_coords_um - centroid, axis=1) / radius
    if np.std(dists) < 1e-9 or np.std(scores) < 1e-9:
        return 0.0
    slope, *_ = linregress(dists, scores)
    return float(slope)


def load_coverage_map(tls_masks_dir):
    """
    Returns: {sample_id: {instance_id_int: record_dict}}
    Handles two JSON formats:
      New: {instance_id, index, area_um2, n_spots, keep, large, split_from}
      Old: {index, area_um2, n_spots, keep}  : instance_id = list position, split_from = None
    """
    coverage_map = {}
    for cov_path in glob.glob(str(tls_masks_dir / '*/tls_spot_coverage.json')):
        sid = Path(cov_path).parent.name
        with open(cov_path) as f:
            records = json.load(f)
        by_iid = {}
        for pos, rec in enumerate(records):
            iid = rec.get('instance_id', pos)
            if 'split_from' not in rec:
                rec = dict(rec, split_from=None, large=False)
            by_iid[iid] = rec
        coverage_map[sid] = by_iid
    log(f"Loaded coverage JSONs for {len(coverage_map)} samples")
    return coverage_map


def load_polygon_map(tls_masks_dir):
    """
    Returns: {sample_id: {polygon_index_int: coordinates_list}}
    Skips stub files (dict format) created for segfault-skipped samples.
    """
    polygon_map = {}
    pattern = str(tls_masks_dir / '*/images/filtered/*_tls_filtered.json')
    for poly_path in glob.glob(pattern):
        sid = Path(poly_path).parts[-4]
        with open(poly_path) as f:
            polygons = json.load(f)
        if not isinstance(polygons, list):
            continue  # stub JSON (e.g. segfault_skipped samples)
        polygon_map[sid] = {int(p['index']): p['coordinates'] for p in polygons if isinstance(p, dict)}
    log(f"Loaded polygon JSONs for {len(polygon_map)} samples")
    return polygon_map


def get_polygon_for_instance(sample_id, instance_id_int, coverage_map, polygon_map):
    """
    Returns (area_um2, coords_px, is_split) or (None, None, False) if not found.
    area_um2 from coverage JSON (accurate for splits).
    coords_px from filtered polygon JSON (parent polygon for splits).
    """
    if sample_id not in coverage_map:
        return None, None, False
    cov = coverage_map[sample_id].get(instance_id_int)
    if cov is None:
        return None, None, False
    area_um2 = cov.get('area_um2')
    is_split = cov.get('split_from') is not None
    poly_idx = int(cov['split_from']) if is_split else int(cov['index'])
    coords = polygon_map.get(sample_id, {}).get(poly_idx)
    return area_um2, coords, is_split


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', default=str(PROJECT / 'outputs/stage2_clustering/tls_neighborhood_spots.h5ad'))
    parser.add_argument('--masks-dir', default=str(PROJECT / 'outputs/tls_masks'))
    parser.add_argument('--out', default=str(PROJECT / 'outputs/stage2_clustering/instance_features_raw_spatial.tsv'))
    parser.add_argument('--min-spots', type=int, default=1)
    args = parser.parse_args()

    tls_masks_dir = Path(args.masks_dir)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    log(f"Loading {args.input}")
    adata = ad.read_h5ad(args.input)
    log(f"  {adata.shape[0]} spots × {adata.shape[1]} genes")
    log(f"  {adata.obs['instance_id'].nunique()} instances across {adata.obs['sample_id'].nunique()} samples")

    gene_names = {g: i for i, g in enumerate(adata.var_names)}

    # Build instance -> row index maps
    inst_tls_rows = defaultdict(list)
    inst_neigh_rows = defaultdict(list)
    for i, (inst, stype) in enumerate(zip(adata.obs['instance_id'], adata.obs['spot_type'])):
        if stype == 'in_tls':
            inst_tls_rows[inst].append(i)
        else:
            inst_neigh_rows[inst].append(i)

    inst_meta = {}
    for inst in inst_tls_rows:
        r = adata.obs.iloc[inst_tls_rows[inst][0]]
        inst_meta[inst] = {
            'sample_id':     r['sample_id'],
            'oncotree_code': r['oncotree_code'],
            'st_technology': r['st_technology'],
        }

    # Spatial coordinates (µm) per obs row
    if 'coord_x_um' in adata.obs.columns:
        spot_x = adata.obs['coord_x_um'].values.astype(float)
        spot_y = adata.obs['coord_y_um'].values.astype(float)
    else:
        log("ERROR: coord_x_um not in obs : expected from 03_extract_tls_neighborhood.py")
        sys.exit(1)

    # Load polygon and coverage data
    coverage_map = load_coverage_map(tls_masks_dir)
    polygon_map  = load_polygon_map(tls_masks_dir)

    instances = sorted(inst_tls_rows.keys())
    log(f"Processing {len(instances)} instances (will exclude {len(EXCLUDED_SAMPLES)} sample set)")

    rows = []
    skipped = 0
    excluded = 0
    no_coverage = 0
    no_polygon = 0

    for inst_id in instances:
        sample_id = inst_meta[inst_id]['sample_id']

        if sample_id in EXCLUDED_SAMPLES:
            excluded += 1
            continue

        tls_rows  = np.array(inst_tls_rows[inst_id])
        neigh_rows = np.array(inst_neigh_rows.get(inst_id, []))
        n_tls = len(tls_rows)

        if n_tls < args.min_spots:
            skipped += 1
            continue

        # Parse instance_id_int: rightmost underscore separates sample_id from int
        inst_int = int(inst_id.rsplit('_', 1)[1])

        # Polygon and coverage lookup
        area_um2, coords_px, is_split = get_polygon_for_instance(
            sample_id, inst_int, coverage_map, polygon_map
        )

        if area_um2 is None:
            no_coverage += 1
            area_um2 = np.nan

        # --- Gene expression features ---
        X_tls   = normalize_log1p_cpm(adata.X, tls_rows)
        X_neigh = normalize_log1p_cpm(adata.X, neigh_rows) if len(neigh_rows) > 0 else None

        row = {
            'instance_id':    inst_id,
            **inst_meta[inst_id],
            'n_spots_in_tls': n_tls,
        }

        for score_name, gs in GENE_SETS.items():
            row[score_name] = gene_set_score(X_tls, gene_names, gs)

        for score_name, gs in NEIGHBOR_GENE_SETS.items():
            if X_neigh is not None:
                row[score_name] = gene_set_score(X_neigh, gene_names, gs)
            else:
                row[score_name] = np.nan

        # --- Polygon morphology ---
        row['log_area_um2'] = float(np.log(area_um2 + 1)) if not np.isnan(area_um2) else np.nan

        if coords_px is not None and not is_split:
            circ, conv = polygon_morphology(coords_px)
            row['circularity'] = circ
            row['convexity']   = conv
        else:
            row['circularity'] = np.nan
            row['convexity']   = np.nan
            if coords_px is None:
                no_polygon += 1

        # --- Radial gradients ---
        spot_coords = np.column_stack([spot_x[tls_rows], spot_y[tls_rows]])
        eff_area = area_um2 if not np.isnan(area_um2) else 1.0

        for (score_name, gs), grad_col in zip(GENE_SETS.items(), GRADIENT_COLS):
            row[grad_col] = radial_gradient(X_tls, spot_coords, gene_names, gs, eff_area)

        rows.append(row)

    log(f"Done: {len(rows)} instances included, {excluded} excluded (BC23803/BC23377), "
        f"{skipped} below min_spots, {no_coverage} missing coverage, {no_polygon} missing polygon")

    col_order = META_COLS + FEATURE_COLS
    df = pd.DataFrame(rows)[col_order]
    df.to_csv(out_path, sep='\t', index=False, float_format='%.6f')

    log(f"\nWrote {len(df)} instances × {len(FEATURE_COLS)} features -> {out_path}")
    log(f"Cohort breakdown: {df.groupby('oncotree_code')['instance_id'].count().to_dict()}")

    nan_rates = {c: f"{df[c].isna().mean()*100:.1f}%" for c in FEATURE_COLS if df[c].isna().any()}
    if nan_rates:
        log(f"NaN rates: {nan_rates}")
    else:
        log("No NaN values in feature columns")

    log("\nFeature summary (mean ± std):")
    for c in FEATURE_COLS:
        log(f"  {c:40s}  {df[c].mean():.4f} ± {df[c].std():.4f}")


if __name__ == '__main__':
    main()
