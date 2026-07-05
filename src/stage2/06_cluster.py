#!/usr/bin/env python3
"""
Stage 2 : Step 4: UMAP + Leiden clustering with batch-correction diagnostics.

Reads (runs on both):
  outputs/stage2_clustering/instance_features_scaled.tsv
  outputs/stage2_clustering/instance_features_corrected.tsv

Writes per input variant ('scaled' / 'corrected'):
  outputs/stage2_clustering/leiden_{variant}_r{res}.tsv   cluster assignments
  outputs/stage2_clustering/umap_{variant}.png            UMAP panel (cluster/platform/cancer)
  outputs/stage2_clustering/cluster_summary.tsv           silhouette + IDC-Visium diagnostics

Diagnostics (run on each variant):

  1. Silhouette score for st_technology vs oncotree_code labels.
     Computed in UMAP 2-D space. If platform_silhouette > cancer_silhouette,
     batch correction is warranted. If cancer wins or they are comparable, use scaled.

  2. IDC Visium scatter.
     The 14 IDC Visium samples are ComBat's calibration anchor : if they cluster
     together (low entropy), ComBat over-leveraged a coherent subgroup. If they
     scatter across clusters (high entropy), they are biologically heterogeneous
     and the correction is epistemically weak.

Leiden sweep: resolution 0.1, 0.3, 0.5, 0.8.
Stability: ARI across 20 random seeds per resolution. Stable = mean ARI ≥ 0.8.

Usage:
  python3 src/stage2_labeling/06_cluster.py
  python3 src/stage2_labeling/06_cluster.py --variants scaled   (one variant only)
  python3 src/stage2_labeling/06_cluster.py --resolutions 0.3 0.5

  # Run on a specific features file (bypasses in-dir construction):
  python3 src/stage2_labeling/06_cluster.py \\
      --features-file outputs/stage2_clustering/instance_features_725.tsv \\
      --variant-name 725
  # Outputs: leiden_725_r{res}.tsv, umap_725.png (in --out-dir)
"""

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import scanpy as sc
import anndata as ad
from sklearn.metrics import silhouette_score, adjusted_rand_score
from scipy.stats import entropy as scipy_entropy

warnings.filterwarnings('ignore')
sc.settings.verbosity = 0

PROJECT = Path(__file__).resolve().parents[2]

FEAT_COLS = [
    'GC_score', 'Tfh_score', 'Plasma_score', 'HEV_score',
    'Suppressive', 'Myeloid_supp', 'Cytotoxic',
    'Neighbor_TGFB1', 'Neighbor_VEGFA', 'Neighbor_CD8A',
    'Neighbor_tumor_fraction', 'Neighbor_myeloid_fraction',
]

RESOLUTIONS = [0.1, 0.3, 0.5, 0.8]
N_STABILITY_SEEDS = 20
N_NEIGHBORS = 15   # KNN graph for UMAP + Leiden

# Palette for oncotree codes (≤12)
CANCER_COLORS = [
    '#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00',
    '#a65628', '#f781bf', '#999999', '#66c2a5', '#fc8d62',
    '#8da0cb', '#e78ac3',
]


def log(msg):
    print(f"[06_cluster] {msg}", flush=True)


def build_adata(df, feat_cols):
    X = df[feat_cols].values.astype(np.float32)
    obs = df.drop(columns=feat_cols).copy().reset_index(drop=True)
    adata = ad.AnnData(X=X, obs=obs)
    adata.obs_names = df['instance_id'].values
    return adata


def run_leiden(adata, resolution, seed):
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, use_rep='X', random_state=seed)
    sc.tl.leiden(adata, resolution=resolution, random_state=seed,
                 key_added=f'leiden_{resolution}_{seed}')
    return adata.obs[f'leiden_{resolution}_{seed}'].values


def leiden_sweep(adata, resolutions, n_seeds):
    """
    For each resolution, run Leiden with n_seeds seeds.
    Returns: dict res -> (labels_seed0, mean_ARI, std_ARI, n_clusters_seed0)
    """
    results = {}
    for res in resolutions:
        all_labels = [run_leiden(adata.copy(), res, s) for s in range(n_seeds)]
        aris = []
        for i in range(n_seeds):
            for j in range(i + 1, n_seeds):
                aris.append(adjusted_rand_score(all_labels[i], all_labels[j]))
        mean_ari = float(np.mean(aris)) if aris else 1.0
        std_ari  = float(np.std(aris))  if aris else 0.0
        labels0  = all_labels[0]
        n_clust  = len(np.unique(labels0))
        results[res] = (labels0, mean_ari, std_ari, n_clust)
        log(f"    r={res}: {n_clust} clusters, stability ARI={mean_ari:.3f}±{std_ari:.3f}")
    return results


def silhouette_diagnostics(umap_coords, platform_labels, cancer_labels):
    """
    Silhouette score in UMAP space for platform vs cancer labels.
    Returns (platform_sil, cancer_sil).
    """
    # Need ≥2 unique labels and ≥2 samples per label for silhouette
    def safe_sil(X, labels):
        ul = np.unique(labels)
        if len(ul) < 2:
            return np.nan
        # Filter out singleton classes (silhouette undefined for n=1)
        counts = {l: (labels == l).sum() for l in ul}
        keep = np.array([counts[l] >= 2 for l in labels])
        labels_filt, X_filt = labels[keep], X[keep]
        if len(np.unique(labels_filt)) < 2:
            return np.nan
        return float(silhouette_score(X_filt, labels_filt))

    plat_sil   = safe_sil(umap_coords, platform_labels)
    cancer_sil = safe_sil(umap_coords, cancer_labels)
    return plat_sil, cancer_sil


def idc_visium_scatter(df, cluster_labels, umap_coords):
    """
    Characterise the 14 IDC Visium samples in clustering space.
    Returns dict with cluster distribution, entropy, and concentration flag.
    """
    is_idc_vis = (df['oncotree_code'] == 'IDC') & (df['st_technology'] == 'Visium')
    n = is_idc_vis.sum()
    if n == 0:
        return {'n': 0, 'note': 'no IDC Visium samples in this variant'}

    idc_vis_clusters = cluster_labels[is_idc_vis.values]
    all_clusters     = np.unique(cluster_labels)
    counts = {str(c): int((idc_vis_clusters == c).sum()) for c in all_clusters}
    probs  = np.array([counts[c] for c in counts]) / n
    probs  = probs[probs > 0]
    ent    = float(scipy_entropy(probs, base=2))
    max_ent = np.log2(len(all_clusters))
    norm_ent = ent / max_ent if max_ent > 0 else 0.0

    # Fraction in the single most common cluster
    top_frac = float(max(probs))

    return {
        'n':               n,
        'cluster_counts':  counts,
        'entropy_bits':    round(ent, 3),
        'norm_entropy':    round(norm_ent, 3),   # 0=concentrated, 1=uniform
        'top_cluster_frac': round(top_frac, 3),
        'interpretation': 'concentrated (poor calibration anchor)'
                          if top_frac > 0.6 else 'scattered (heterogeneous anchor)',
    }


def plot_umap(adata, cluster_labels, variant, res, out_path):
    """
    3-panel UMAP: by cluster / by platform / by cancer type.
    Also marks IDC Visium samples with a distinct marker.
    """
    umap  = adata.obsm['X_umap']
    df    = adata.obs.reset_index(drop=True)
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f'UMAP : {variant}, resolution={res}', fontsize=13)

    # Panel 1: Leiden clusters
    ax = axes[0]
    unique_clusters = sorted(np.unique(cluster_labels), key=lambda x: int(x))
    cmap = plt.cm.get_cmap('tab10', len(unique_clusters))
    for i, cl in enumerate(unique_clusters):
        mask = cluster_labels == cl
        ax.scatter(umap[mask, 0], umap[mask, 1], c=[cmap(i)], s=10, alpha=0.7, label=cl)
    ax.set_title('Leiden clusters')
    ax.legend(markerscale=2, fontsize=7, loc='best')
    ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')

    # Panel 2: Platform
    ax = axes[1]
    plat_colors = {'Spatial Transcriptomics': '#2166ac', 'Visium': '#d73027'}
    for plat, color in plat_colors.items():
        mask = df['st_technology'] == plat
        ax.scatter(umap[mask.values, 0], umap[mask.values, 1],
                   c=color, s=10, alpha=0.5, label=plat)
    # Overlay IDC Visium as stars
    idc_vis = (df['oncotree_code'] == 'IDC') & (df['st_technology'] == 'Visium')
    if idc_vis.sum() > 0:
        ax.scatter(umap[idc_vis.values, 0], umap[idc_vis.values, 1],
                   c='gold', s=60, marker='*', zorder=5, label='IDC Visium (calibration anchor)')
    ax.set_title(f'Platform (★ = IDC Visium n={idc_vis.sum()})')
    ax.legend(markerscale=1.5, fontsize=7, loc='best')
    ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')

    # Panel 3: Cancer type
    ax = axes[2]
    cancers = sorted(df['oncotree_code'].unique())
    for i, code in enumerate(cancers):
        mask = df['oncotree_code'] == code
        color = CANCER_COLORS[i % len(CANCER_COLORS)]
        ax.scatter(umap[mask.values, 0], umap[mask.values, 1],
                   c=color, s=10, alpha=0.6, label=f'{code} (n={mask.sum()})')
    ax.set_title('Cancer type')
    ax.legend(markerscale=2, fontsize=6, loc='best', ncol=1)
    ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    log(f"    Saved {out_path.name}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--in-dir', default=str(PROJECT / 'outputs/stage2_clustering'))
    parser.add_argument('--out-dir', default=str(PROJECT / 'outputs/stage2_clustering'))
    parser.add_argument('--variants', nargs='+', default=['scaled', 'corrected'])
    parser.add_argument('--resolutions', nargs='+', type=float, default=RESOLUTIONS)
    parser.add_argument('--n-neighbors', type=int, default=N_NEIGHBORS)
    parser.add_argument('--features-file', default=None,
                        help='Explicit features TSV; overrides --in-dir construction. '
                             'Use with --variant-name to set output file prefix.')
    parser.add_argument('--variant-name', default=None,
                        help='Output label when --features-file is used (e.g. "725").')
    args = parser.parse_args()

    in_dir  = Path(args.in_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # If --features-file given, run only that file under the specified variant name
    if args.features_file:
        feat_path = Path(args.features_file)
        if not feat_path.exists():
            log(f"ERROR: --features-file {feat_path} not found"); return
        vname = args.variant_name or feat_path.stem.removeprefix('instance_features_') or 'custom'
        variants_to_run = [(vname, feat_path)]
    else:
        variants_to_run = [
            (v, in_dir / f'instance_features_{v}.tsv') for v in args.variants
        ]

    summary_rows = []

    for variant, tsv_path in variants_to_run:
        if not tsv_path.exists():
            log(f"Skipping {variant}: {tsv_path} not found")
            continue

        log(f"\n{'='*60}")
        log(f"Variant: {variant}")
        log(f"{'='*60}")

        df = pd.read_csv(tsv_path, sep='\t')
        # Auto-detect feature columns: any column not in META_COLS
        _meta = ['instance_id', 'sample_id', 'oncotree_code', 'n_spots_in_tls', 'st_technology']
        feat_cols = [c for c in df.columns if c not in _meta]
        adata = build_adata(df, feat_cols)

        # UMAP (fixed seed for reproducibility)
        log("  Computing UMAP...")
        sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, use_rep='X', random_state=42)
        sc.tl.umap(adata, random_state=42)
        umap_coords = adata.obsm['X_umap']

        # Silhouette diagnostics
        plat_sil, cancer_sil = silhouette_diagnostics(
            umap_coords,
            df['st_technology'].values,
            df['oncotree_code'].values,
        )
        log(f"  Silhouette (UMAP space):")
        log(f"    platform (st_technology): {plat_sil:.4f}")
        log(f"    cancer   (oncotree_code): {cancer_sil:.4f}")
        if plat_sil > cancer_sil:
            log(f"    -> Platform dominates separation. Correction warranted.")
        else:
            log(f"    -> Cancer type dominates or ties. Scaled version preferred.")

        # Leiden sweep
        log(f"  Leiden sweep (resolutions={args.resolutions}, seeds={N_STABILITY_SEEDS})...")
        leiden_results = leiden_sweep(adata, args.resolutions, N_STABILITY_SEEDS)

        # Find best resolution: highest stability with 2–5 clusters preferred
        stable_res = {r: v for r, v in leiden_results.items()
                      if v[1] >= 0.8 and 2 <= v[3] <= 6}
        if stable_res:
            best_res = max(stable_res, key=lambda r: stable_res[r][1])
        else:
            best_res = max(leiden_results, key=lambda r: leiden_results[r][1])
        log(f"  Best resolution: {best_res} ({leiden_results[best_res][3]} clusters, "
            f"ARI={leiden_results[best_res][1]:.3f})")

        # IDC Visium scatter at best resolution
        best_labels = leiden_results[best_res][0]
        idc_scatter = idc_visium_scatter(df, best_labels, umap_coords)
        log(f"  IDC Visium scatter (r={best_res}): {idc_scatter}")

        # Save cluster assignments for all resolutions
        for res, (labels, mean_ari, std_ari, n_clust) in leiden_results.items():
            out_tsv = out_dir / f'leiden_{variant}_r{res}.tsv'
            df_out = df[['instance_id', 'sample_id', 'oncotree_code',
                         'n_spots_in_tls', 'st_technology']].copy()
            df_out['leiden_cluster'] = labels
            df_out.to_csv(out_tsv, sep='\t', index=False)

        # UMAP plot at best resolution
        plot_umap(adata, best_labels, variant, best_res,
                  out_dir / f'umap_{variant}.png')

        # Summary row for each resolution
        for res, (labels, mean_ari, std_ari, n_clust) in leiden_results.items():
            idc_res = idc_visium_scatter(df, labels, umap_coords)
            summary_rows.append({
                'variant':               variant,
                'resolution':            res,
                'n_clusters':            n_clust,
                'stability_ari_mean':    round(mean_ari, 4),
                'stability_ari_std':     round(std_ari, 4),
                'platform_silhouette':   round(plat_sil, 4),
                'cancer_silhouette':     round(cancer_sil, 4),
                'platform_dominates':    plat_sil > cancer_sil,
                'idc_visium_n':          idc_res['n'],
                'idc_visium_norm_entropy': idc_res.get('norm_entropy', np.nan),
                'idc_visium_top_frac':   idc_res.get('top_cluster_frac', np.nan),
                'idc_visium_interpretation': idc_res.get('interpretation', ''),
            })

    # Write summary
    summary_path = out_dir / 'cluster_summary.tsv'
    pd.DataFrame(summary_rows).to_csv(summary_path, sep='\t', index=False)
    log(f"\nWrote cluster summary -> {summary_path}")

    # Print decision table
    log("\n" + "="*70)
    log("DECISION TABLE")
    log("="*70)
    for row in summary_rows:
        flag = "← USE THIS" if (
            not row['platform_dominates'] and
            row['stability_ari_mean'] >= 0.8 and
            2 <= row['n_clusters'] <= 5
        ) else ""
        log(f"  {row['variant']:12s}  r={row['resolution']}  "
            f"k={row['n_clusters']}  ARI={row['stability_ari_mean']:.3f}  "
            f"plat_sil={row['platform_silhouette']:.3f}  "
            f"cancer_sil={row['cancer_silhouette']:.3f}  "
            f"idc_vis_ent={row['idc_visium_norm_entropy']:.2f}  {flag}")


if __name__ == '__main__':
    main()
