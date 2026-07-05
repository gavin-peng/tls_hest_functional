#!/usr/bin/env python3
"""
Stage 2 : Step 1b: Fill missing oncotree_code values in tls_neighborhood_spots.h5ad.

Some HEST samples have no oncotree_code in HEST_v1_3_0.csv but have unambiguous
organ + tissue metadata. This script fills those codes from metadata rather than
leaving them as 'nan', which would create a spurious batch in downstream correction.

Mapping (organ-based, confirmed against per-sample metadata JSONs):
  Prostate (MEND series)            -> PRAD
  Skin (MEND37, MEND39)             -> SKCM
  Brain (MEND63)                    -> GBM
  Lung Parenchyma (MISC19)          -> LUAD
  Lung Bronchi/Trachea (MISC26/29/32) -> LUSC

Reads:   outputs/stage2_clustering/tls_neighborhood_spots.h5ad
Writes:  outputs/stage2_clustering/tls_neighborhood_spots.h5ad  (in-place)

Usage:
  python3 src/stage2_labeling/03b_fix_oncotree.py
  python3 src/stage2_labeling/03b_fix_oncotree.py --dry-run
"""

import argparse
import json
from pathlib import Path

import anndata as ad
import pandas as pd

PROJECT = Path(__file__).resolve().parents[2]

# Explicit per-sample mapping derived from organ + tissue fields in metadata JSONs.
# Only covers samples with nan oncotree_code in HEST_v1_3_0.csv.
ONCOTREE_FIX = {
    # Prostate cancer (MEND series, organ=Prostate)
    'MEND37':  'SKCM',   # organ=Skin, tissue=skin
    'MEND39':  'SKCM',   # organ=Skin, tissue=skin
    'MEND61':  'PRAD',   # organ=Prostate
    'MEND63':  'GBM',    # organ=Brain
    'MEND140': 'PRAD',
    'MEND142': 'PRAD',
    'MEND149': 'PRAD',
    'MEND151': 'PRAD',
    'MEND154': 'PRAD',
    'MEND156': 'PRAD',
    'MEND157': 'PRAD',
    'MEND158': 'PRAD',
    'MEND159': 'PRAD',
    'MEND160': 'PRAD',
    'MEND161': 'PRAD',
    # Lung (MISC series, disease=Treated)
    # tissue=Parenchyma -> LUAD (peripheral); Bronchi/Trachea -> LUSC (central)
    'MISC19':  'LUAD',   # tissue=Parenchyma
    'MISC26':  'LUSC',   # tissue=Bronchi
    'MISC29':  'LUSC',   # tissue=Bronchi
    'MISC32':  'LUSC',   # tissue=Trachea
}


def log(msg):
    print(f"[03b_fix] {msg}", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', default=str(PROJECT / 'outputs/stage2_clustering/tls_neighborhood_spots.h5ad'))
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()

    path = Path(args.input)
    log(f"Loading {path}")
    adata = ad.read_h5ad(path)

    before = adata.obs['oncotree_code'].value_counts(dropna=False)
    log(f"oncotree_code before: {before[before.index == 'nan'].to_dict() or {'nan': 0}} nan spots")

    # oncotree_code is stored as Categorical : add new categories before assigning
    new_codes = set(ONCOTREE_FIX.values()) - set(adata.obs['oncotree_code'].cat.categories)
    if new_codes:
        adata.obs['oncotree_code'] = adata.obs['oncotree_code'].cat.add_categories(sorted(new_codes))

    n_fixed = 0
    for sid, new_code in ONCOTREE_FIX.items():
        mask = adata.obs['sample_id'] == sid
        n = mask.sum()
        if n == 0:
            continue
        current = adata.obs.loc[mask, 'oncotree_code'].iloc[0]
        if current != 'nan':
            log(f"  {sid}: already has oncotree_code={current}, skipping")
            continue
        log(f"  {sid}: {n} spots -> {new_code}")
        if not args.dry_run:
            adata.obs.loc[mask, 'oncotree_code'] = new_code
        n_fixed += n

    after_nan = (adata.obs['oncotree_code'] == 'nan').sum()
    log(f"\nFixed {n_fixed} spots across {len(ONCOTREE_FIX)} samples")
    log(f"Remaining nan spots: {after_nan}")

    if after_nan > 0:
        remaining = adata.obs[adata.obs['oncotree_code'] == 'nan']['sample_id'].unique()
        log(f"  Samples still nan: {sorted(remaining)}")

    if args.dry_run:
        log("Dry run : no changes written")
        return

    log(f"Writing {path}")
    adata.write_h5ad(path)

    log("\noncotree_code after fix:")
    after = adata.obs.drop_duplicates('sample_id').groupby('oncotree_code')['sample_id'].count()
    for code, n in after.sort_values(ascending=False).items():
        log(f"  {code}: {n} samples")


if __name__ == '__main__':
    main()
