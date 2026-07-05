#!/usr/bin/env python3
"""
Stage 2 : Step 2: Per-TLS instance feature matrix.

Reads:
  outputs/stage2_clustering/tls_neighborhood_spots.h5ad

Writes:
  outputs/stage2_clustering/instance_features.tsv

One row per TLS instance. Columns:
  Block 1 : 7 gene set scores  (in_tls spots, log1p-CPM)
  Block 3 : 5 neighborhood scores (neighborhood spots, log1p-CPM)
  Metadata: instance_id, sample_id, oncotree_code, n_spots_in_tls, st_technology

Normalization: log1p-CPM per spot, computed lazily per instance to avoid
densifying the full 64k × 40k matrix (~10 GB dense).

Moran's I (Block 2) was dropped: near-zero values (mean ≈ −0.07 to −0.10) and
57–99% NaN rates confirmed that within-TLS spatial autocorrelation at Visium
resolution measures noise, not structure (75th pct = 15 in_tls spots; too few
for meaningful spatial statistics at 100µm spot spacing).

Neighbor_tumor_fraction and Neighbor_myeloid_fraction are mean gene set scores
(proxy, not true cell fractions).

Usage:
  python3 src/stage2_labeling/04_instance_features.py
  python3 src/stage2_labeling/04_instance_features.py --min-spots 1
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

PROJECT = Path(__file__).resolve().parents[2]

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

FEATURE_COLS = list(GENE_SETS.keys()) + list(NEIGHBOR_GENE_SETS.keys())


def log(msg):
    print(f"[04_features] {msg}", flush=True)


def normalize_log1p_cpm(X_sparse, row_indices):
    """Extract rows from sparse matrix and normalize to log1p-CPM. Dense float32 output."""
    if len(row_indices) == 0:
        return np.zeros((0, X_sparse.shape[1]), dtype=np.float32)
    X = X_sparse[row_indices].toarray().astype(np.float32)
    totals = X.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1.0
    return np.log1p(X / totals * 1e6)


def gene_set_score(X_norm, gene_names, gene_set):
    """
    Mean log1p-CPM of available genes in set, averaged over spots.
    Returns NaN if no genes available or no spots.
    """
    cols = [gene_names[g] for g in gene_set if g in gene_names]
    if not cols or X_norm.shape[0] == 0:
        return np.nan
    return float(X_norm[:, cols].mean())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input',  default=str(PROJECT / 'outputs/stage2_clustering/tls_neighborhood_spots.h5ad'))
    parser.add_argument('--out',    default=str(PROJECT / 'outputs/stage2_clustering/instance_features.tsv'))
    parser.add_argument('--min-spots', type=int, default=1,
                        help='Minimum in_tls spots to include an instance (default: 1)')
    args = parser.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    log(f"Loading {args.input}")
    adata = ad.read_h5ad(args.input)
    log(f"  {adata.shape[0]} spots × {adata.shape[1]} genes")
    log(f"  {adata.obs['instance_id'].nunique()} instances across {adata.obs['sample_id'].nunique()} samples")

    gene_names = {g: i for i, g in enumerate(adata.var_names)}

    # Pre-build instance -> row indices map to avoid repeated boolean scans
    inst_tls_rows = defaultdict(list)
    inst_neigh_rows = defaultdict(list)
    for i, (inst, stype) in enumerate(zip(adata.obs['instance_id'], adata.obs['spot_type'])):
        if stype == 'in_tls':
            inst_tls_rows[inst].append(i)
        else:
            inst_neigh_rows[inst].append(i)

    # First obs row per instance for metadata
    inst_meta = {}
    for inst in inst_tls_rows:
        r = adata.obs.iloc[inst_tls_rows[inst][0]]
        inst_meta[inst] = {
            'sample_id':     r['sample_id'],
            'oncotree_code': r['oncotree_code'],
            'st_technology': r['st_technology'],
        }

    instances = sorted(inst_tls_rows.keys())
    log(f"Computing features for {len(instances)} instances...")

    rows = []
    skipped = 0

    for inst_id in instances:
        tls_rows = np.array(inst_tls_rows[inst_id])
        neigh_rows = np.array(inst_neigh_rows.get(inst_id, []))
        n_tls = len(tls_rows)

        if n_tls < args.min_spots:
            skipped += 1
            continue

        row = {
            'instance_id':    inst_id,
            **inst_meta[inst_id],
            'n_spots_in_tls': n_tls,
        }

        # Normalize lazily : only these rows, full gene width
        X_tls   = normalize_log1p_cpm(adata.X, tls_rows)
        X_neigh = normalize_log1p_cpm(adata.X, neigh_rows) if len(neigh_rows) > 0 else None

        # Block 1: gene set scores on in_tls spots
        for score_name, gs in GENE_SETS.items():
            row[score_name] = gene_set_score(X_tls, gene_names, gs)

        # Block 3: neighborhood context scores
        for score_name, gs in NEIGHBOR_GENE_SETS.items():
            if X_neigh is not None:
                row[score_name] = gene_set_score(X_neigh, gene_names, gs)
            else:
                row[score_name] = np.nan

        rows.append(row)

    if skipped:
        log(f"Skipped {skipped} instances with < {args.min_spots} in_tls spots")

    if not rows:
        log("No instances passed filters : check --min-spots and input file.")
        sys.exit(1)

    col_order = ['instance_id', 'sample_id', 'oncotree_code', 'n_spots_in_tls', 'st_technology'] + FEATURE_COLS
    df = pd.DataFrame(rows)[col_order]
    df.to_csv(out_path, sep='\t', index=False, float_format='%.6f')

    log(f"\nWrote {len(df)} instances × {len(df.columns)} columns -> {out_path}")
    log(f"Cohort breakdown: {df.groupby('oncotree_code')['instance_id'].count().to_dict()}")

    nan_rates = {c: f"{df[c].isna().mean()*100:.1f}%" for c in FEATURE_COLS if df[c].isna().any()}
    if nan_rates:
        log(f"NaN rates: {nan_rates}")
    else:
        log("No NaN values in any feature column")

    log(f"\nFeature summary (mean ± std):")
    for c in FEATURE_COLS:
        log(f"  {c:35s}  {df[c].mean():.4f} ± {df[c].std():.4f}")


if __name__ == '__main__':
    main()
