#!/usr/bin/env python3
"""
Stage 2 : Step 1 (revised): Extract ST spots inside TLS instances + 300µm neighborhood ring.

Key differences from 01_extract_tls_spots.py:
  - ALL TLS instances with n_spots ≥ 1 (not just keepable/n_ge20)
  - No spot cap
  - 300µm buffer ring spots included and labelled separately
  - obs column: spot_type = 'in_tls' | 'neighborhood'

Reads:
  outputs/tls_spot_coverage.tsv
  outputs/tls_results_cc354.tsv
  data/HEST_v1_3_0.csv
  data/metadata/{SID}.json
  data/st/{SID}.h5ad
  outputs/tls_masks/{SID}/images/filtered/*_tls_filtered.json
  outputs/tls_masks/{SID}/tls_spot_coverage.json

Writes:
  outputs/stage2_clustering/tls_neighborhood_spots.h5ad

obs columns:
  sample_id, instance_id, spot_type ('in_tls'|'neighborhood'),
  oncotree_code, tissue, st_technology, n_spots_in_tls,
  coord_x_um, coord_y_um

Usage:
  python3 src/stage2_labeling/03_extract_tls_neighborhood.py
  python3 src/stage2_labeling/03_extract_tls_neighborhood.py --samples SPA104 SPA103
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import tifffile
from shapely.geometry import Polygon, Point
from skimage.draw import polygon as sk_poly
from skimage.feature import peak_local_max
from skimage.segmentation import watershed
from scipy.ndimage import gaussian_filter

PROJECT = Path(__file__).resolve().parents[2]
HEATMAP_PX = 0.5    # µm per pixel in heatmap
LARGE_UM2  = 1e6    # 1 mm² : polygons above this are watershed-split
SIGMA_UM   = 100.0
MIN_DIST_UM = 300.0
PEAK_THRESH = 50.0
DS = 4              # heatmap downsample factor for watershed
NEIGHBORHOOD_UM = 300.0  # ring width in µm


def log(msg):
    print(f"[03_extract] {msg}", flush=True)


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def spots_inside_polygon(poly, spots_um):
    """Return boolean array len(spots_um): True if spot is inside poly."""
    minx, miny, maxx, maxy = poly.bounds
    mask = np.zeros(len(spots_um), dtype=bool)
    bbox = ((spots_um[:, 0] >= minx) & (spots_um[:, 0] <= maxx) &
            (spots_um[:, 1] >= miny) & (spots_um[:, 1] <= maxy))
    for i in np.where(bbox)[0]:
        if poly.contains(Point(spots_um[i, 0], spots_um[i, 1])):
            mask[i] = True
    return mask


def watershed_basin_labels(poly_um, inside_idx, spots_um, heatmap):
    """
    Watershed-split a large polygon. Returns (labels_for_inside_idx, n_peaks)
    where labels are 1-based basin IDs, or None if splitting yields ≤1 basin.
    """
    minx_um, miny_um, maxx_um, maxy_um = poly_um.bounds
    minx_px = max(0, int(minx_um / HEATMAP_PX))
    miny_px = max(0, int(miny_um / HEATMAP_PX))
    maxx_px = min(heatmap.shape[1], int(np.ceil(maxx_um / HEATMAP_PX)) + 1)
    maxy_px = min(heatmap.shape[0], int(np.ceil(maxy_um / HEATMAP_PX)) + 1)

    crop_full = heatmap[miny_px:maxy_px, minx_px:maxx_px].astype(np.float32)
    if crop_full.size == 0:
        return None

    from skimage.transform import downscale_local_mean
    crop = downscale_local_mean(crop_full, (DS, DS)).astype(np.float32)
    H, W = crop.shape
    WS_PX = HEATMAP_PX * DS

    ext = np.array(poly_um.exterior.coords)
    col_px = ((ext[:, 0] - minx_px * HEATMAP_PX) / WS_PX).clip(0, W - 1)
    row_px = ((ext[:, 1] - miny_px * HEATMAP_PX) / WS_PX).clip(0, H - 1)
    rr, cc = sk_poly(row_px, col_px, (H, W))
    mask = np.zeros((H, W), dtype=bool)
    mask[rr, cc] = True

    sigma_px = SIGMA_UM / WS_PX
    min_dist_px = max(1, int(MIN_DIST_UM / WS_PX))
    smoothed = gaussian_filter(np.where(mask, crop, 0.0), sigma=sigma_px)
    peaks = peak_local_max(smoothed, min_distance=min_dist_px,
                           threshold_abs=PEAK_THRESH, labels=mask)
    if len(peaks) <= 1:
        return None

    markers = np.zeros((H, W), dtype=int)
    for i, (r, c) in enumerate(peaks, start=1):
        markers[r, c] = i
    seg = watershed(-smoothed, markers=markers, mask=mask)

    if len(inside_idx) == 0:
        return np.array([], dtype=int), len(peaks)

    s = spots_um[inside_idx]
    s_col = np.clip(((s[:, 0] - minx_px * HEATMAP_PX) / WS_PX).astype(int), 0, W - 1)
    s_row = np.clip(((s[:, 1] - miny_px * HEATMAP_PX) / WS_PX).astype(int), 0, H - 1)
    labels = seg[s_row, s_col]
    return labels, len(peaks)


# ---------------------------------------------------------------------------
# Ensembl ID -> HGNC symbol conversion (same as 01_extract_tls_spots.py)
# ---------------------------------------------------------------------------

def load_ensembl_map(tsv_path):
    df = pd.read_csv(tsv_path, sep='\t', header=None,
                     names=['ensembl_id', 'hgnc_symbol'])
    df = df.dropna(subset=['hgnc_symbol'])
    df = df[df['hgnc_symbol'].str.strip() != '']
    df = df.drop_duplicates(subset='ensembl_id', keep='first')
    return dict(zip(df['ensembl_id'], df['hgnc_symbol']))


def convert_ensembl(adata, ensg_map):
    if not adata.var_names[0].startswith('ENSG'):
        return adata, False
    new_names, keep_idx, seen = [], [], set()
    for i, ensg in enumerate(adata.var_names):
        sym = ensg_map.get(ensg, '')
        if sym and sym not in seen:
            new_names.append(sym)
            keep_idx.append(i)
            seen.add(sym)
    adata_sub = adata[:, keep_idx].copy()
    adata_sub.var_names = new_names
    return adata_sub, True


# ---------------------------------------------------------------------------
# Per-sample extraction
# ---------------------------------------------------------------------------

def extract_sample(sid, cov_instances, h5ad_path, tls_json_path, heatmap_path,
                   spacing, oncotree_code, tissue, st_technology, ensg_map=None):
    """
    Returns list of dicts, one per (instance, spot_type):
      {barcodes, X, spatial_um, instance_id, spot_type, n_spots_in_tls}
    var_names is returned as second element.
    """
    adata = ad.read_h5ad(h5ad_path)

    if ensg_map is not None:
        adata, converted = convert_ensembl(adata, ensg_map)
        if converted:
            log(f"    {sid}: Ensembl->symbol converted ({adata.shape[1]} genes)")

    spots_um = adata.obsm['spatial'] * spacing  # (n_spots, 2)

    raw = json.load(open(tls_json_path))
    feats = raw.get('features', raw) if isinstance(raw, dict) else raw

    polys = {}
    for f in feats:
        coords_um = [(x * HEATMAP_PX, y * HEATMAP_PX) for x, y in f['coordinates']]
        try:
            poly = Polygon(coords_um).buffer(0)
            if poly.is_valid and not poly.is_empty:
                polys[f['index']] = poly
        except Exception:
            pass

    # Load heatmap only if any large instance needs watershed splitting
    heatmap = None
    needs_heatmap = any(
        inst.get('large') and inst.get('n_spots', 0) >= 1
        for inst in cov_instances
    )
    if needs_heatmap:
        hp = Path(heatmap_path)
        if hp.exists():
            raw_tif = tifffile.imread(str(hp))
            while raw_tif.ndim > 2:
                raw_tif = raw_tif[0]
            heatmap = raw_tif

    # Cache per parent polygon: (inside_mask, neighbor_mask, watershed_labels)
    poly_cache = {}

    def get_poly_data(poly_idx):
        if poly_idx in poly_cache:
            return poly_cache[poly_idx]
        poly = polys.get(poly_idx)
        if poly is None:
            poly_cache[poly_idx] = (None, None, None, None)
            return poly_cache[poly_idx]

        poly_buf = poly.buffer(NEIGHBORHOOD_UM)
        in_mask = spots_inside_polygon(poly, spots_um)
        buf_mask = spots_inside_polygon(poly_buf, spots_um)
        neigh_mask = buf_mask & ~in_mask
        inside_idx = np.where(in_mask)[0]

        ws_labels = None
        if heatmap is not None and poly.area > LARGE_UM2 and len(inside_idx) > 0:
            result = watershed_basin_labels(poly, inside_idx, spots_um, heatmap)
            if result is not None:
                ws_labels, _ = result

        poly_cache[poly_idx] = (in_mask, neigh_mask, inside_idx, ws_labels)
        return poly_cache[poly_idx]

    results = []
    for inst in cov_instances:
        if inst.get('n_spots', 0) < 1:
            continue

        idx_str = inst['index']
        is_split = inst.get('split_from') is not None

        if is_split:
            poly_idx = inst['split_from']
            sub_k = int(str(idx_str).split('_')[-1])  # 0-based
        else:
            poly_idx = idx_str
            sub_k = None

        in_mask, neigh_mask, inside_idx, ws_labels = get_poly_data(poly_idx)
        if in_mask is None:
            continue

        # Select in_tls spots for this instance
        if is_split and ws_labels is not None and len(inside_idx) > 0:
            basin_mask = ws_labels == (sub_k + 1)
            tls_idx = inside_idx[basin_mask]
        else:
            tls_idx = inside_idx

        if len(tls_idx) == 0:
            continue

        neigh_idx = np.where(neigh_mask)[0]
        n_tls = len(tls_idx)
        instance_name = f"{sid}_{inst['instance_id']}"

        def make_rows(idx_arr, stype):
            if len(idx_arr) == 0:
                return None
            X = adata.X[idx_arr]
            if sp.issparse(X):
                X = X.tocsr()
            barcodes = list(adata.obs_names[idx_arr])
            xy_um = spots_um[idx_arr]
            n = len(idx_arr)
            obs = pd.DataFrame({
                'sample_id':      sid,
                'instance_id':    instance_name,
                'spot_type':      stype,
                'oncotree_code':  oncotree_code,
                'tissue':         tissue,
                'st_technology':  st_technology,
                'n_spots_in_tls': n_tls,
                'coord_x_um':     xy_um[:, 0],
                'coord_y_um':     xy_um[:, 1],
            }, index=[f"{instance_name}_{stype}_{i}" for i in range(n)])
            return X.astype(np.float32), obs, xy_um

        tls_data = make_rows(tls_idx, 'in_tls')
        if tls_data is not None:
            results.append(tls_data)

        neigh_data = make_rows(neigh_idx, 'neighborhood')
        if neigh_data is not None:
            results.append(neigh_data)

    return results, adata.var_names.tolist()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--samples', nargs='+', default=None)
    parser.add_argument('--out', default=str(PROJECT / 'outputs/stage2_clustering/tls_neighborhood_spots.h5ad'))
    parser.add_argument('--ensembl-map', default=str(PROJECT / 'data/ensembl_to_symbol.tsv'))
    parser.add_argument('--species', default='Homo sapiens')
    parser.add_argument('--disease', nargs='+', default=['Cancer', 'Treated'])
    parser.add_argument('--exclude-organs', nargs='+', default=['Lymph node'])
    parser.add_argument('--exclude-samples', nargs='+', default=['MISC17'],
                        help='Exclude specific sample IDs regardless of other filters. '
                             'Default: MISC17 (bronchial BALT, organ=Lung but not tumor)')
    args = parser.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    ensg_map = None
    emap_path = Path(args.ensembl_map)
    if emap_path.exists():
        ensg_map = load_ensembl_map(emap_path)
        log(f"Loaded Ensembl map: {len(ensg_map):,} entries")
    else:
        log(f"WARNING: no Ensembl map at {emap_path}")

    cov_df = pd.read_csv(PROJECT / 'outputs/tls_spot_coverage.tsv', sep='\t')
    hest_df = pd.read_csv(PROJECT / 'data/HEST_v1_3_0.csv', low_memory=False)

    refilter_path = PROJECT / 'outputs/tls_results_cc354.tsv'
    duplicate_sids = set()
    if refilter_path.exists():
        rf = pd.read_csv(refilter_path, sep='\t')
        dup = rf[rf['duplicate_of'].notna() & (rf['duplicate_of'].astype(str).str.strip() != '')]
        duplicate_sids = set(dup['sample_id'].tolist())
        if duplicate_sids:
            log(f"Duplicates excluded: {sorted(duplicate_sids)}")

    meta_dir = PROJECT / 'data/metadata'

    # All non-Xenium samples with at least one TLS instance
    if args.samples:
        candidate_sids = args.samples
    else:
        candidate_sids = cov_df[
            (cov_df.st_technology != 'Xenium') &
            (pd.to_numeric(cov_df.n_tls, errors='coerce').fillna(0) > 0)
        ].sample_id.tolist()

    exclude_samples = set(args.exclude_samples) if args.exclude_samples else set()

    keep_sids, skip_reasons = [], {}
    for sid in candidate_sids:
        if sid in duplicate_sids:
            skip_reasons[sid] = 'duplicate'
            continue
        if sid in exclude_samples:
            skip_reasons[sid] = 'excluded_sample'
            continue
        meta_path = meta_dir / f'{sid}.json'
        if not meta_path.exists():
            skip_reasons[sid] = 'no_metadata'
            continue
        m = json.load(open(meta_path))
        if args.species and m.get('species', '') != args.species:
            skip_reasons[sid] = f"species={m.get('species', '')}"
            continue
        if args.disease and m.get('disease_state', '') not in args.disease:
            skip_reasons[sid] = f"disease={m.get('disease_state', '')}"
            continue
        if args.exclude_organs and m.get('organ', '') in args.exclude_organs:
            skip_reasons[sid] = f"organ={m.get('organ', '')}"
            continue
        keep_sids.append(sid)

    log(f"Skipped {len(skip_reasons)} samples: {dict(list(skip_reasons.items())[:10])}{'...' if len(skip_reasons)>10 else ''}")
    log(f"Processing {len(keep_sids)} samples")

    all_X, all_obs, all_spatial, all_var_per_batch = [], [], [], []

    for sid in keep_sids:
        cov_row = cov_df[cov_df.sample_id == sid]
        if cov_row.empty:
            log(f"  {sid}: not in coverage TSV, skipping")
            continue
        cr = cov_row.iloc[0]

        h5ad_path = PROJECT / 'data/st' / f'{sid}.h5ad'
        if not h5ad_path.exists():
            log(f"  {sid}: no h5ad, skipping")
            continue

        if ensg_map is None:
            import h5py
            with h5py.File(h5ad_path, 'r') as f:
                first_gene = f['var']['_index'][0].decode('utf-8') if '_index' in f['var'] \
                    else list(f['var'].keys())[0]
            if first_gene.startswith('ENSG'):
                log(f"  {sid}: Ensembl IDs, no map : skipping")
                continue

        tls_json = None
        for candidate in [
            PROJECT / 'outputs/tls_masks' / sid / 'images/filtered' / f'{sid}_asap_hooknettls_tls_filtered.json',
            PROJECT / 'outputs/tls_masks' / sid / 'filtered' / f'{sid}_asap_hooknettls_tls_filtered.json',
        ]:
            if candidate.exists():
                tls_json = candidate
                break
        if tls_json is None:
            log(f"  {sid}: no TLS JSON, skipping")
            continue

        cov_json = PROJECT / 'outputs/tls_masks' / sid / 'tls_spot_coverage.json'
        if not cov_json.exists():
            log(f"  {sid}: no coverage JSON, skipping")
            continue

        heatmap_path = PROJECT / 'outputs/tls_masks' / sid / 'images' / f'{sid}_asap_hooknettls_heat1.tif'

        hrow = hest_df[hest_df['id'] == sid]
        if hrow.empty:
            log(f"  {sid}: not in HEST CSV, skipping")
            continue
        spacing = float(hrow.iloc[0]['pixel_size_um_estimated'])

        cov_instances = json.load(open(cov_json))
        oncotree_code = str(cr['oncotree_code'])
        tissue = str(cr['tissue'])
        st_tech = str(cr['st_technology'])
        n_inst_all = len([i for i in cov_instances if i.get('n_spots', 0) >= 1])

        log(f"  {sid} ({oncotree_code}, {st_tech}): {n_inst_all} instances with ≥1 spot")
        try:
            results, var_names = extract_sample(
                sid, cov_instances, h5ad_path, tls_json, heatmap_path,
                spacing, oncotree_code, tissue, st_tech, ensg_map=ensg_map)
        except Exception as e:
            log(f"  {sid}: ERROR : {e}")
            import traceback; traceback.print_exc()
            continue

        if not results:
            log(f"  {sid}: no spots extracted")
            continue

        in_tls_count = sum(r[1]['spot_type'].iloc[0] == 'in_tls' for r in results)
        neigh_count = len(results) - in_tls_count
        log(f"    {in_tls_count} in_tls batches, {neigh_count} neighborhood batches")

        for X_batch, obs_batch, xy_batch in results:
            all_X.append(X_batch)
            all_obs.append(obs_batch)
            all_spatial.append(xy_batch)
            all_var_per_batch.append(var_names)

    if not all_X:
        log("No spots extracted : check filters and input files.")
        sys.exit(1)

    # Gene union + zero-padding
    all_genes_union = set()
    for vnames in all_var_per_batch:
        all_genes_union |= set(vnames)
    all_genes_union = sorted(all_genes_union)
    n_genes = len(all_genes_union)
    log(f"Gene union: {n_genes} genes across {len(set(map(tuple, all_var_per_batch)))} unique gene spaces")

    gene_to_new_col = {g: i for i, g in enumerate(all_genes_union)}
    aligned_X = []
    for x, vnames in zip(all_X, all_var_per_batch):
        n_spots = x.shape[0]
        old_cols = np.array([i for i, g in enumerate(vnames) if g in gene_to_new_col], dtype=np.int32)
        new_cols = np.array([gene_to_new_col[g] for g in vnames if g in gene_to_new_col], dtype=np.int32)
        if len(old_cols) > 0:
            if sp.issparse(x):
                x_sub = x[:, old_cols].tocoo()
                r, c_local, d = x_sub.row, x_sub.col, x_sub.data
                c_global = new_cols[c_local]
            else:
                x_sub = x[:, old_cols]
                nz_r, nz_c = np.nonzero(x_sub)
                d = x_sub[nz_r, nz_c]
                r, c_global = nz_r, new_cols[nz_c]
            sparse_padded = sp.csr_matrix(
                (d.astype(np.float32), (r, c_global)), shape=(n_spots, n_genes))
        else:
            sparse_padded = sp.csr_matrix((n_spots, n_genes), dtype=np.float32)
        aligned_X.append(sparse_padded)

    X_pooled = sp.vstack(aligned_X, format='csr')
    obs_pooled = pd.concat(all_obs, axis=0)
    spatial_pooled = np.vstack(all_spatial)

    adata_out = ad.AnnData(
        X=X_pooled if sp.issparse(X_pooled) else sp.csr_matrix(X_pooled),
        obs=obs_pooled,
        var=pd.DataFrame(index=all_genes_union),
    )
    adata_out.obsm['spatial_um'] = spatial_pooled

    adata_out.write_h5ad(out_path)

    n_instances = obs_pooled['instance_id'].nunique()
    n_samples = obs_pooled['sample_id'].nunique()
    n_tls_spots = (obs_pooled['spot_type'] == 'in_tls').sum()
    n_neigh_spots = (obs_pooled['spot_type'] == 'neighborhood').sum()

    log(f"\nWrote {adata_out.shape[0]} spots × {adata_out.shape[1]} genes -> {out_path}")
    log(f"Instances: {n_instances} across {n_samples} samples")
    log(f"Spot breakdown: {n_tls_spots} in_tls, {n_neigh_spots} neighborhood")
    log(f"Cohorts: {obs_pooled.groupby('oncotree_code').size().to_dict()}")


if __name__ == '__main__':
    main()
