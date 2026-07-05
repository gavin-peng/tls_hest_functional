#!/usr/bin/env python3
"""
Stage 2 : Step 0: Build sample->patient map from HEST metadata.

Reads the 119 samples present in tls_neighborhood_spots.h5ad and joins
HEST_v1_3_0.csv to extract:

  patient_id       : raw 'patient' field from HEST metadata (may be NaN for
                     anonymous/singleton samples such as 10x Genomics demos)
  idc_subtype      : 'disease_comment' field (tumor subtype / diagnosis string)
  split_patient_id : globally unique patient identifier safe for train/val/test
                     splits (see note below)

Output: outputs/stage2_clustering/sample_patient_map.tsv

Note on split_patient_id disambiguation
-----------------------------------------
The raw 'patient' strings in HEST metadata are dataset-local, not globally unique.
Example: 'Patient 1' in the ZEN colorectal dataset (ZEN49, READ) and 'Patient 1'
in the MEND skin dataset (MEND37/39, SKCM) are different people.
Naively grouping by patient_id for splits would merge them into one unit, creating
incorrect cross-study leakage.

Disambiguation rules (applied in order):
  1. SPA samples whose patient matches ^BC\\d+$ (Andersson/Stahl BC##### codes)
     -> use BC##### directly; these are globally unique within that dataset.
  2. All other samples with a non-null patient string
     -> prefix with the series identifier derived from sample_id (SPA, MEND, ZEN,
        NCBI, TENX, MISC, ...) to scope the patient to its source dataset.
  3. NaN or blank patient
     -> use sample_id itself (singleton treatment for anonymous samples).

Allocatable units for GNN splits:
  Unique split_patient_id values in instance_features_725.tsv determine the units.
  Serial sections from the same patient MUST go to the same split fold.

Usage:
  python3 src/stage2_labeling/00_make_patient_map.py
  python3 src/stage2_labeling/00_make_patient_map.py --h5ad path/to/spots.h5ad
"""

import argparse
import re
from pathlib import Path

import anndata as ad
import pandas as pd

PROJECT = Path(__file__).resolve().parents[2]
H5AD_DEFAULT  = PROJECT / 'outputs/stage2_clustering/tls_neighborhood_spots.h5ad'
HEST_CSV      = PROJECT / 'data/HEST_v1_3_0.csv'
OUT_PATH      = PROJECT / 'outputs/stage2_clustering/sample_patient_map.tsv'


def log(msg):
    print(f"[00_make_patient_map] {msg}", flush=True)


def series_prefix(sample_id: str) -> str:
    """Return the alphabetic series prefix of a HEST sample ID (e.g. SPA, MEND, ZEN)."""
    m = re.match(r'^([A-Z]+)', sample_id)
    return m.group(1) if m else sample_id


def make_split_patient_id(sample_id: str, patient) -> str:
    """
    Build a globally unique patient identifier for train/val/test splits.

    SPA BC##### codes are the only globally unique patient strings in HEST.
    All other patient strings are scoped to their source dataset by prefixing
    with the series identifier. NaN patients are treated as singletons.
    """
    if pd.isna(patient) or str(patient).strip() == '':
        return sample_id  # singleton

    patient_str = str(patient).strip()
    prefix = series_prefix(sample_id)

    # SPA Andersson/Stahl BC##### codes are globally unique across both publications
    if prefix == 'SPA' and re.match(r'^BC\d+$', patient_str):
        return patient_str

    # All other patient strings: scope to source dataset
    return f"{prefix}_{patient_str}"


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    parser.add_argument('--h5ad', default=str(H5AD_DEFAULT),
                        help='Path to tls_neighborhood_spots.h5ad')
    parser.add_argument('--hest-csv', default=str(HEST_CSV))
    parser.add_argument('--out', default=str(OUT_PATH))
    args = parser.parse_args()

    h5ad_path = Path(args.h5ad)
    hest_path = Path(args.hest_csv)
    out_path  = Path(args.out)

    # ── 1. Sample IDs from h5ad (authoritative list) ─────────────────────
    log(f"Reading sample IDs from {h5ad_path.name}...")
    adata = ad.read_h5ad(h5ad_path, backed='r')
    sample_ids = sorted(adata.obs['sample_id'].unique())
    adata.file.close()
    log(f"  {len(sample_ids)} unique samples")

    # ── 2. Join HEST CSV ──────────────────────────────────────────────────
    log(f"Loading {hest_path.name}...")
    hest = pd.read_csv(hest_path, low_memory=False).set_index('id')

    rows = []
    missing = []
    for sid in sample_ids:
        if sid not in hest.index:
            missing.append(sid)
            rows.append({'sample_id': sid, 'patient_id': None, 'idc_subtype': None})
            continue
        row = hest.loc[sid]
        patient  = row.get('patient')
        subtype  = row.get('disease_comment')
        rows.append({'sample_id': sid, 'patient_id': patient, 'idc_subtype': subtype})

    if missing:
        log(f"  WARNING: {len(missing)} sample IDs not found in HEST CSV: {missing}")

    df = pd.DataFrame(rows)

    # ── 3. Build split_patient_id ─────────────────────────────────────────
    df['split_patient_id'] = df.apply(
        lambda r: make_split_patient_id(r['sample_id'], r['patient_id']), axis=1
    )

    # ── 4. Diagnostics ────────────────────────────────────────────────────
    spa_mask  = df['sample_id'].str.startswith('SPA')
    n_spa_pat = df.loc[spa_mask, 'split_patient_id'].nunique()
    n_nonspa  = (~spa_mask).sum()
    n_total   = df['split_patient_id'].nunique()

    log(f"  SPA unique patients: {n_spa_pat}")
    log(f"  Non-SPA samples (all singletons): {n_nonspa}")
    log(f"  Total allocatable split units: {n_total}")

    # Warn about raw patient strings that were disambiguated (cross-study collisions)
    raw_counts = df[df['patient_id'].notna()].groupby('patient_id')['sample_id'].apply(
        lambda s: s.apply(series_prefix).nunique()
    )
    collisions = raw_counts[raw_counts > 1]
    if len(collisions):
        log(f"  Disambiguated cross-study patient_id collisions ({len(collisions)}):")
        for raw, n_series in collisions.items():
            affected = df[df['patient_id'] == raw][['sample_id', 'split_patient_id']].values
            log(f"    '{raw}' appeared in {n_series} series -> "
                + ", ".join(f"{s}->{p}" for s, p in affected))
    else:
        log("  No cross-study patient_id collisions found.")

    # ── 5. Write ──────────────────────────────────────────────────────────
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep='\t', index=False)
    log(f"Wrote {len(df)} rows -> {out_path}")


if __name__ == '__main__':
    main()
