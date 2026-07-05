#!/usr/bin/env python3
"""
Stage 2 : Step 3: Batch correction on per-TLS instance feature matrix.

Reads:
  outputs/stage2_clustering/instance_features.tsv   (or --input path)

Writes:
  <stem>_scaled.tsv    (z-score only; NaN imputed with column median before scaling)
  <stem>_corrected.tsv (ComBat + z-score)

Feature columns are auto-detected: any column not in META_COLS is treated as a feature.
NaN values (e.g. circularity/convexity for split instances) are median-imputed before
z-scoring and clearly logged.

Two outputs because the dataset has a severe confound: all 566 Spatial
Transcriptomics instances are IDC; all non-IDC are Visium (+ 14 IDC Visium).
Correcting on platform = partially shifting IDC scores toward non-IDC, which
may remove real IDC biology. Both outputs are produced; the clustering step
chooses which to use.

Batch variable: st_technology (Spatial Transcriptomics vs Visium).
  - Spatial Transcriptomics (200µm spots): systematically lower scores for
    rare cell-type markers (GC, Tfh, HEV) due to spot-size dilution.
  - Visium (55µm spots): higher scores for the same markers.
  - Platform range: 0.1–1.1 z-score units across features.
  - Confound: ST ≡ IDC (566/566); Visium = all others + 14 IDC.
    The 14 IDC Visium samples are the only intra-cancer bridge.

Method: ComBat (Johnson 2007) via scanpy, batch=st_technology, no covariates.
  ComBat estimates a location+scale shift per batch feature-by-feature.
  With no covariates, biological variation is NOT protected : the correction
  assumes any mean difference between platforms is technical, which may be
  partially false given the IDC/platform confound.

Usage:
  python3 src/stage2_labeling/05_batch_correct.py
  python3 src/stage2_labeling/05_batch_correct.py --input outputs/stage2_clustering/instance_features_raw_spatial.tsv
  python3 src/stage2_labeling/05_batch_correct.py --batch-key oncotree_code
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from sklearn.preprocessing import StandardScaler

PROJECT = Path(__file__).resolve().parents[2]

META_COLS = ['instance_id', 'sample_id', 'oncotree_code', 'n_spots_in_tls', 'st_technology']


def log(msg):
    print(f"[05_batch] {msg}", flush=True)


def zscore(df, feat_cols):
    X = df[feat_cols].values.copy().astype(float)
    # Median-impute NaN columns (before winsorization so NaN don't affect percentiles)
    nan_cols = [c for c in feat_cols if df[c].isna().any()]
    if nan_cols:
        for j, c in enumerate(feat_cols):
            if c in nan_cols:
                col = X[:, j]
                median = np.nanmedian(col)
                col[np.isnan(col)] = median
                X[:, j] = col
    # Winsorize high-variance gradient features at [p2, p98] to suppress outliers from
    # noisy regression slopes on instances with only 3–5 in_tls spots
    for j, c in enumerate(feat_cols):
        if c.startswith('Grad_'):
            col = X[:, j]
            p2, p98 = np.nanpercentile(col, [2, 98])
            X[:, j] = np.clip(col, p2, p98)
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    return X_scaled, scaler


def combat_correct(X, batch_labels):
    """
    Run ComBat via scanpy on an (n_samples, n_features) array.
    Returns corrected array, same shape.
    """
    adata = ad.AnnData(
        X=X.astype(np.float32),
        obs=pd.DataFrame({'batch': batch_labels}).reset_index(drop=True),
    )
    sc.pp.combat(adata, key='batch', inplace=True)
    return adata.X


def batch_effect_report(X, labels, feat_cols, tag):
    """Log per-feature mean by batch and the max cross-batch range."""
    df = pd.DataFrame(X, columns=feat_cols)
    df['batch'] = labels
    means = df.groupby('batch')[feat_cols].mean()
    ranges = (means.max() - means.min()).sort_values(ascending=False)
    log(f"  [{tag}] Cross-batch range (z-score units) : top 5:")
    for feat, rng in ranges.head(5).items():
        vals = {b: f"{means.loc[b, feat]:+.3f}" for b in means.index}
        log(f"    {feat:<35s}  range={rng:.3f}  {vals}")
    return float(ranges.mean())


def confound_check(df, batch_key):
    """Report how strongly batch and oncotree_code are confounded."""
    xtab = pd.crosstab(df['oncotree_code'], df[batch_key])
    batches = xtab.columns.tolist()
    log(f"Confound check : oncotree_code × {batch_key}:")
    log(f"\n{xtab.to_string()}")
    # Fraction of samples where one batch dominates a cohort
    cohort_dominant = (xtab.max(axis=1) / xtab.sum(axis=1))
    fully_confounded = (cohort_dominant == 1.0).sum()
    log(f"  Cohorts fully on one platform: {fully_confounded}/{len(xtab)}")
    if fully_confounded == len(xtab):
        log("  WARNING: batch and cancer type are perfectly confounded.")
        log("  ComBat correction will be applied but interpret with caution :")
        log("  any IDC-vs-non-IDC platform shift cannot be separated from biology.")
    return fully_confounded == len(xtab)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', default=str(PROJECT / 'outputs/stage2_clustering/instance_features.tsv'))
    parser.add_argument('--out-dir', default=str(PROJECT / 'outputs/stage2_clustering'))
    parser.add_argument('--batch-key', default='st_technology',
                        help='obs column to use as batch variable (default: st_technology)')
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Derive output file names from input stem
    input_stem = Path(args.input).stem
    # Strip _raw_ infix so that instance_features_raw_spatial -> instance_features_scaled_spatial
    # For other inputs: append _scaled / _corrected suffix to the stem
    clean_stem = input_stem.replace('_raw_', '_')
    scaled_path    = out_dir / f'{clean_stem}_scaled.tsv'
    corrected_path = out_dir / f'{clean_stem}_corrected.tsv'

    log(f"Loading {args.input}")
    df = pd.read_csv(args.input, sep='\t')

    # Auto-detect feature columns: everything not in META_COLS
    FEAT_COLS = [c for c in df.columns if c not in META_COLS]
    log(f"  {len(df)} instances × {len(FEAT_COLS)} features (auto-detected)")
    log(f"  Feature columns: {FEAT_COLS}")
    log(f"  Cohorts: {df['oncotree_code'].value_counts().to_dict()}")
    log(f"  Platforms: {df['st_technology'].value_counts().to_dict()}")

    nan_rates = {c: f"{df[c].isna().mean()*100:.1f}%" for c in FEAT_COLS if df[c].isna().any()}
    if nan_rates:
        log(f"  NaN rates (will be median-imputed): {nan_rates}")

    if args.batch_key not in df.columns:
        log(f"ERROR: --batch-key '{args.batch_key}' not in columns")
        sys.exit(1)

    batch_labels = df[args.batch_key].values
    n_batches = len(np.unique(batch_labels))
    log(f"\nBatch variable: {args.batch_key} ({n_batches} levels)")

    # --- Confound check ---
    fully_confounded = confound_check(df, args.batch_key)

    # --- Z-score standardize ---
    log("\nZ-score standardizing features...")
    X_scaled, scaler = zscore(df, FEAT_COLS)
    mean_range_before = batch_effect_report(X_scaled, batch_labels, FEAT_COLS, 'before correction')

    # Save scaled (no batch correction)
    df_scaled = df[META_COLS].copy()
    df_scaled[FEAT_COLS] = X_scaled
    df_scaled.to_csv(scaled_path, sep='\t', index=False, float_format='%.6f')
    log(f"\nWrote z-scored features -> {scaled_path}")

    # --- ComBat correction ---
    if n_batches < 2:
        log("Only one batch level : skipping ComBat, writing scaled as corrected too")
        df_scaled.to_csv(corrected_path, sep='\t', index=False, float_format='%.6f')
        return

    log("\nRunning ComBat...")
    try:
        X_corrected = combat_correct(X_scaled, batch_labels)
    except Exception as e:
        log(f"ComBat failed: {e} : writing scaled as corrected")
        df_scaled.to_csv(corrected_path, sep='\t', index=False, float_format='%.6f')
        return

    mean_range_after = batch_effect_report(X_corrected, batch_labels, FEAT_COLS, 'after ComBat')
    reduction = (mean_range_before - mean_range_after) / mean_range_before * 100
    log(f"\n  Mean cross-batch range: {mean_range_before:.3f} -> {mean_range_after:.3f} "
        f"({reduction:.1f}% reduction)")

    # Save corrected
    df_corrected = df[META_COLS].copy()
    df_corrected[FEAT_COLS] = X_corrected
    df_corrected.to_csv(corrected_path, sep='\t', index=False, float_format='%.6f')
    log(f"Wrote ComBat-corrected features -> {corrected_path}")

    # --- Summary ---
    log("\nSummary:")
    log(f"  instance_features_scaled.tsv    : z-score only; preserves all inter-cohort differences")
    log(f"  instance_features_corrected.tsv : ComBat on {args.batch_key}; removes platform offset")
    if fully_confounded:
        log(f"  RECOMMENDATION: use scaled for clustering. ComBat is provided as a sensitivity check.")
        log(f"  If clusters separate purely by platform in the scaled version, revisit corrected.")
    else:
        log(f"  RECOMMENDATION: use corrected for clustering.")


if __name__ == '__main__':
    main()
