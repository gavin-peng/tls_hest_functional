#!/usr/bin/env python3
"""
Assemble publication figures for manuscript.md.
Run from the paper/ directory:
    source ../.venv/bin/activate && python3 assemble_figures.py
"""

import os, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import seaborn as sns
from sklearn.preprocessing import StandardScaler
import scipy.stats as ss
import anndata as ad
import scanpy as sc
from PIL import Image

warnings.filterwarnings('ignore')
sc.settings.verbosity = 0

# ── Paths ──────────────────────────────────────────────────────────────────
ROOT    = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_DIR = os.path.join(ROOT, 'paper', 'figures')
FEAT    = os.path.join(ROOT, 'outputs/stage2_clustering/instance_features_670.tsv')
CLUST   = os.path.join(ROOT, 'outputs/stage2_clustering/leiden_670_r0.3.tsv')
H5AD    = os.path.join(ROOT, 'outputs/stage2_clustering/tls_neighborhood_spots.h5ad')
SHIAO   = os.path.join(ROOT, 'outputs/validation_geo/shiao_C1C4_pseudobulk.tsv')
TCGA_KM = os.path.join(ROOT, 'outputs/validation_tcga/TCGA_BRCA_TLS_survival.png')
MET_KM  = os.path.join(ROOT, 'outputs/validation_tcga/METABRIC_TLS_survival.png')
CESC_KM = os.path.join(ROOT, 'outputs/validation_tcga/CESC_TLS_survival.png')

os.makedirs(OUT_DIR, exist_ok=True)

# ── Global style ───────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family':    'DejaVu Sans',
    'font.size':      9,
    'axes.titlesize': 10,
    'axes.labelsize': 9,
    'xtick.labelsize':8,
    'ytick.labelsize':8,
    'axes.spines.top':   False,
    'axes.spines.right': False,
    'figure.dpi': 150,
})
DPI = 150

CLUSTER_COLORS = {
    0: '#aaaaaa',
    1: '#029E73',
    2: '#9467bd',
    3: '#ff7f0e',
    4: '#d62728',
    5: '#1f77b4',
}
CLUSTER_LABELS = {
    0: 'C0: Low-activity',
    1: 'C1: Immunogenic',
    2: 'C2: Excluded',
    3: 'C3: Angiogenic',
    4: 'C4: Myeloid-supp.',
    5: 'C5: Tumour-proximal',
}
FEATURE_NAMES = {
    'GC_score':                 'GC score',
    'Tfh_score':                'Tfh score',
    'Plasma_score':             'Plasma score',
    'HEV_score':                'HEV score',
    'Suppressive':              'Suppressive',
    'Myeloid_supp':             'Myeloid-supp.',
    'Cytotoxic':                'Cytotoxic',
    'Neighbor_TGFB1':           'Nbr. TGFB1',
    'Neighbor_VEGFA':           'Nbr. VEGFA',
    'Neighbor_CD8A':            'Nbr. CD8A',
    'Neighbor_tumor_fraction':  'Tumor frac.',
    'Neighbor_myeloid_fraction':'Myeloid frac.',
}
FEAT_COLS = list(FEATURE_NAMES.keys())

CANCER_COLORS = {
    'IDC':  '#2166ac', 'PRAD': '#e08214', 'CESC': '#4dac26',
    'COAD': '#d01c8b', 'PAAD': '#f1b6da', 'READ': '#b8e186',
    'BLCA': '#542788', 'EPM':  '#fee0b6', 'LUSC': '#7fbf7b',
    'LUAD': '#af8dc3', 'GBM':  '#8c510a', 'SKCM': '#bf812d',
}

# ── Load and merge data ─────────────────────────────────────────────────────
print("Loading features and cluster labels...")
feat  = pd.read_csv(FEAT,  sep='\t')
clust = pd.read_csv(CLUST, sep='\t')
df = feat.merge(clust[['instance_id', 'leiden_cluster']], on='instance_id')
df['cluster_label'] = df['leiden_cluster'].map(CLUSTER_LABELS)

scaler = StandardScaler()
df_z = df.copy()
df_z[FEAT_COLS] = scaler.fit_transform(df[FEAT_COLS].values)

# ── Compute UMAP ────────────────────────────────────────────────────────────
print("Computing UMAP...")
adata_umap = sc.AnnData(X=df_z[FEAT_COLS].values.astype(np.float32))
adata_umap.obs['cluster'] = df['leiden_cluster'].values.astype(str)
sc.pp.neighbors(adata_umap, n_neighbors=15, use_rep='X', random_state=42)
sc.tl.umap(adata_umap, random_state=42)
umap_coords = adata_umap.obsm['X_umap']
df['umap1'] = umap_coords[:, 0]
df['umap2'] = umap_coords[:, 1]


def add_panel_label(ax, label, x=-0.12, y=1.05, fontsize=12):
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=fontsize, fontweight='bold', va='top', ha='left')


def draw_umap(ax, df, title='', size=12, alpha=0.7):
    for cid in sorted(df['leiden_cluster'].unique()):
        sub = df[df['leiden_cluster'] == cid]
        ax.scatter(sub['umap1'], sub['umap2'],
                   c=CLUSTER_COLORS[cid], s=size, alpha=alpha,
                   label=CLUSTER_LABELS[cid], linewidths=0)
    ax.set_xlabel('UMAP 1', labelpad=2)
    ax.set_ylabel('UMAP 2', labelpad=2)
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — Study overview  (2 panels: pipeline + UMAP)
# ══════════════════════════════════════════════════════════════════════════════
print("Generating Figure 1...")

fig1 = plt.figure(figsize=(14, 4.8))
gs1  = gridspec.GridSpec(1, 2, figure=fig1, wspace=0.28,
                         width_ratios=[1.1, 1],
                         left=0.03, right=0.97, top=0.93, bottom=0.07)

# ── Panel A: pipeline (left column) + validation (right column) ──────────────
ax1a = fig1.add_subplot(gs1[0, 0])
ax1a.set_xlim(0, 10)
ax1a.set_ylim(0, 10)
ax1a.axis('off')
add_panel_label(ax1a, 'A', x=0.0, y=1.02)

# Blue pipeline: left column (x 0.3–4.6)
pipeline_steps = [
    ('H&E WSI',        '246 high-priority samples',     '#e8f4fd', '#2166ac'),
    ('HookNet-TLS',    '0.5 µm/px GPU inference',       '#cce5ff', '#2166ac'),
    ('TLS detection',  'cc≥354 → 747 instances',         '#b3d7ff', '#2166ac'),
    ('ST extraction',  'in-TLS + 300 µm ring',           '#99c9ff', '#2166ac'),
    ('Feature matrix', '12 features / instance (z)',     '#80bbff', '#2166ac'),
    ('Leiden r=0.3',   '670 instances → 6 clusters',     '#4da6ff', '#2166ac'),
]
# Green validation: right column (x 5.4–9.7), aligned with bottom 3 blue boxes
validation_steps = [
    ('Survival validation',    'TCGA-BRCA, METABRIC, TCGA-CESC', '#d4edda', '#276f3e'),
    ('Immune checkpoint\nblockade response',  'pembrolizumab, n=29 TNBC',        '#b8dfc4', '#276f3e'),
    ('Morphological transfer', 'H&E + GAT, AUROC = 0.857',        '#9cd1ae', '#276f3e'),
]

box_h = 1.2
gap   = 0.3
total_h = len(pipeline_steps) * box_h + (len(pipeline_steps) - 1) * gap  # = 8.7
y0 = (10 - total_h) / 2  # = 0.65

# Draw blue pipeline (left column)
pipe_x0, pipe_w = 0.3, 4.3
for i, (title, subtitle, color, edge) in enumerate(pipeline_steps):
    y = y0 + (len(pipeline_steps) - 1 - i) * (box_h + gap)
    rect = mpatches.FancyBboxPatch((pipe_x0, y), pipe_w, box_h,
                                   boxstyle='round,pad=0.12',
                                   facecolor=color, edgecolor=edge, linewidth=0.8)
    ax1a.add_patch(rect)
    cx = pipe_x0 + pipe_w / 2
    ax1a.text(cx, y + box_h * 0.64, title,
              ha='center', va='center', fontsize=7.5, fontweight='bold', color='#1a1a2e')
    ax1a.text(cx, y + box_h * 0.24, subtitle,
              ha='center', va='center', fontsize=6.3, color='#444444')
    if i < len(pipeline_steps) - 1:
        next_y = y0 + (len(pipeline_steps) - 2 - i) * (box_h + gap)
        ax1a.annotate('', xy=(cx, next_y + box_h + 0.01), xytext=(cx, y - 0.01),
                      arrowprops=dict(arrowstyle='->', color=edge, lw=1.1))

# Draw green validation (right column, vertically centered on bottom 3 blue boxes)
val_x0, val_w = 5.4, 4.3
n_val = len(validation_steps)
val_total_h = n_val * box_h + (n_val - 1) * gap        # 4.2
# Centre the green stack on the lower half of the blue column
blue_mid_low = y0 + (1.5 * (box_h + gap)) / 2 + box_h / 2   # ~midpoint of bottom 3
val_y0 = y0  # align bottom of green stack with bottom of blue stack

for j, (title, subtitle, color, edge) in enumerate(validation_steps):
    y = val_y0 + (n_val - 1 - j) * (box_h + gap)
    rect = mpatches.FancyBboxPatch((val_x0, y), val_w, box_h,
                                   boxstyle='round,pad=0.12',
                                   facecolor=color, edgecolor=edge, linewidth=0.8)
    ax1a.add_patch(rect)
    cx = val_x0 + val_w / 2
    ax1a.text(cx, y + box_h * 0.64, title,
              ha='center', va='center', fontsize=7.5, fontweight='bold', color='#1a1a2e')
    ax1a.text(cx, y + box_h * 0.24, subtitle,
              ha='center', va='center', fontsize=6.3, color='#444444')
    if j < n_val - 1:
        next_y = val_y0 + (n_val - 2 - j) * (box_h + gap)
        ax1a.annotate('', xy=(cx, next_y + box_h + 0.01), xytext=(cx, y - 0.01),
                      arrowprops=dict(arrowstyle='->', color=edge, lw=1.1))

# Horizontal arrow: right edge of Leiden → left edge of top green box
leiden_y   = y0  # Leiden is the bottom (i=5) blue box
leiden_cx_y = leiden_y + box_h / 2
top_green_y = val_y0 + (n_val - 1) * (box_h + gap)
top_green_cx_y = top_green_y + box_h / 2

# L-shaped connector: Leiden right → horizontal → up → top-green left
lx_start = pipe_x0 + pipe_w          # right edge of Leiden: 4.6
lx_mid   = (lx_start + val_x0) / 2  # midpoint gap: 4.95

ax1a.annotate('', xy=(val_x0 - 0.02, top_green_cx_y),
              xytext=(lx_start + 0.02, leiden_cx_y),
              arrowprops=dict(arrowstyle='->', color='#276f3e', lw=1.1,
                              connectionstyle='arc3,rad=-0.25'))

# Small label on connector
ax1a.text(lx_mid + 0.1, (leiden_cx_y + top_green_cx_y) / 2,
          'downstream\nvalidation', ha='left', va='center',
          fontsize=5.8, color='#276f3e', style='italic')

# ── Panel B: cluster UMAP ─────────────────────────────────────────────────────
ax1b = fig1.add_subplot(gs1[0, 1])
draw_umap(ax1b, df, title='TLS instances (UMAP, $n$=670)')
add_panel_label(ax1b, 'B')

legend_handles = [mpatches.Patch(facecolor=CLUSTER_COLORS[c], label=CLUSTER_LABELS[c])
                  for c in sorted(CLUSTER_LABELS)]
ax1b.legend(handles=legend_handles, fontsize=6.5, frameon=False,
            loc='lower left', ncol=1, handlelength=1.0, handleheight=0.8)

plt.savefig(os.path.join(OUT_DIR, 'fig1_overview.png'),
            dpi=DPI, bbox_inches='tight')
plt.close(fig1)
print("  Saved fig1_overview.png")


# ── cohort composition data (used in Fig 2A) ──────────────────────────────────
cancer_order = (df.groupby('oncotree_code')['instance_id']
                .count().sort_values(ascending=False).index.tolist())
cancer_counts = (df.groupby(['oncotree_code', 'leiden_cluster'])
                 .size().unstack(fill_value=0).loc[cancer_order])


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — Cluster characterisation  (3 panels: cohort + heatmap + cancer)
# ══════════════════════════════════════════════════════════════════════════════
print("Generating Figure 2...")

fig2 = plt.figure(figsize=(14, 4.8))
gs2  = gridspec.GridSpec(1, 3, figure=fig2, wspace=0.45,
                         width_ratios=[0.75, 1, 0.85],
                         left=0.05, right=0.96, top=0.90, bottom=0.10)

# ── Panel A: instances by cancer type (moved from Fig 1C) ────────────────────
ax2a = fig2.add_subplot(gs2[0, 0])
add_panel_label(ax2a, 'A')

bottoms = np.zeros(len(cancer_order))
for cid in sorted(CLUSTER_COLORS.keys()):
    if cid in cancer_counts.columns:
        vals = cancer_counts[cid].values
        ax2a.barh(cancer_order, vals, left=bottoms,
                  color=CLUSTER_COLORS[cid], label=CLUSTER_LABELS[cid])
        bottoms += vals

ax2a.set_xlabel('Number of instances')
ax2a.set_title('Instances by cancer type')
ax2a.invert_yaxis()
ax2a.legend(fontsize=5.5, frameon=False, loc='lower right',
            handlelength=0.8, handleheight=0.7)

# ── Panel B: feature heatmap ──────────────────────────────────────────────────
ax2c = fig2.add_subplot(gs2[0, 1])
add_panel_label(ax2c, 'B')

heatmap_data = (df_z.groupby(df['leiden_cluster'])[FEAT_COLS]
                .mean()
                .rename(columns=FEATURE_NAMES))

sns.heatmap(heatmap_data.T,
            ax=ax2c,
            cmap='RdBu_r', center=0,
            vmin=-1.5, vmax=1.5,
            linewidths=0.4, linecolor='white',
            xticklabels=[f'C{i}' for i in sorted(heatmap_data.index)],
            yticklabels=True,
            cbar_kws={'label': 'Mean z-score', 'shrink': 0.7})
ax2c.set_title('Feature profile per cluster')
ax2c.set_xlabel('')
ax2c.set_ylabel('')

# ── Panel C: cancer type composition per cluster ──────────────────────────────
ax2d = fig2.add_subplot(gs2[0, 2])
add_panel_label(ax2d, 'C')

clust_cancer = (df.groupby(['leiden_cluster', 'oncotree_code'])
                .size().unstack(fill_value=0))
clust_cancer_pct = clust_cancer.div(clust_cancer.sum(axis=1), axis=0) * 100
cancer_types = clust_cancer_pct.columns.tolist()

bottoms_d = np.zeros(6)
bar_x = [f'C{i}' for i in sorted(clust_cancer_pct.index)]
for ct in cancer_types:
    vals = clust_cancer_pct[ct].values
    color = CANCER_COLORS.get(ct, '#dddddd')
    ax2d.bar(bar_x, vals, bottom=bottoms_d, color=color, label=ct, width=0.65)
    bottoms_d += vals

ax2d.set_ylabel('Composition (%)')
ax2d.set_title('Cancer type per cluster')
ax2d.set_ylim(0, 105)

handles2d = [mpatches.Patch(facecolor=CANCER_COLORS.get(ct, '#dddddd'), label=ct)
             for ct in cancer_types]
ax2d.legend(handles=handles2d, fontsize=6, frameon=False,
            bbox_to_anchor=(1.01, 1), loc='upper left', ncol=1)

total_per_clust = df['leiden_cluster'].value_counts().sort_index()
for i, (bar_label, n) in enumerate(zip(bar_x, total_per_clust)):
    ax2d.text(i, 103, f'n={n}', ha='center', va='bottom', fontsize=6.5)

plt.savefig(os.path.join(OUT_DIR, 'fig2_clustering.png'),
            dpi=DPI, bbox_inches='tight')
plt.close(fig2)
print("  Saved fig2_clustering.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — Internal biological coherence of C1
# ══════════════════════════════════════════════════════════════════════════════
print("Generating Figure 3 (loading h5ad backed='r')...")

COHERENCE_GENES = ['CR2', 'LTA', 'MS4A1', 'TNFSF13B', 'LTB', 'IGHD', 'CCR7', 'CD38']

adata_full = ad.read_h5ad(H5AD, backed='r')
genes_present = [g for g in COHERENCE_GENES if g in adata_full.var_names]
print(f"  Coherence genes present: {genes_present}")

# Subset to in_tls spots and target genes only
in_tls_mask = adata_full.obs['spot_type'] == 'in_tls'
expr = adata_full[in_tls_mask, genes_present].to_memory()
adata_full.file.close()

# Normalise
sc.pp.normalize_total(expr, target_sum=1e4)
sc.pp.log1p(expr)

# Merge cluster labels
expr_df = expr.to_df()
expr_df['instance_id'] = expr.obs['instance_id'].values
cluster_map = df.set_index('instance_id')['leiden_cluster'].to_dict()
expr_df['leiden_cluster'] = expr_df['instance_id'].map(cluster_map)
expr_df = expr_df.dropna(subset=['leiden_cluster'])
expr_df['leiden_cluster'] = expr_df['leiden_cluster'].astype(int)

# Compute mean expression per cluster
mean_expr = expr_df.groupby('leiden_cluster')[genes_present].mean()
mean_all_others = (expr_df[expr_df['leiden_cluster'] != 1][genes_present].mean())

# Fold change C1 vs all others (ratio of means, pseudocount 0.01)
fc = (mean_expr.loc[1] + 0.01) / (mean_all_others + 0.01)
fc_sorted = fc.sort_values(ascending=False)

fig3, ax3 = plt.subplots(figsize=(6.5, 4.2))

colors_fc = ['#029E73' if v >= 1 else '#aaaaaa' for v in fc_sorted.values]
bars = ax3.barh(range(len(fc_sorted)), fc_sorted.values,
                color=colors_fc, edgecolor='white', linewidth=0.3)
ax3.axvline(1.0, color='#555555', lw=0.8, ls='--')

for i, (gene, val) in enumerate(fc_sorted.items()):
    ax3.text(val + 0.05, i, f'{val:.1f}×', va='center', ha='left', fontsize=8)

ax3.set_yticks(range(len(fc_sorted)))
ax3.set_yticklabels(fc_sorted.index, fontsize=9)
ax3.set_xlabel('Mean expression fold change (C1 vs all others)', labelpad=4)
ax3.set_title('C1 enrichment of TLS organiser genes\n(not included in feature set)',
              pad=8)
ax3.invert_yaxis()
ax3.set_xlim(0, fc_sorted.max() * 1.25)

# Per-cluster means for top 3 genes as small annotation
note_genes = fc_sorted.index[:3].tolist()
note_lines = []
for g in note_genes:
    c1_mean  = mean_expr.loc[1, g]
    all_mean = expr_df[genes_present].mean()[g]
    note_lines.append(f'{g}: C1={c1_mean:.3f}, all={all_mean:.3f}')
ax3.text(0.97, 0.05, '\n'.join(note_lines), transform=ax3.transAxes,
         fontsize=6.5, ha='right', va='bottom', color='#444444',
         bbox=dict(fc='#f5f5f5', ec='none', pad=3))

add_panel_label(ax3, 'A', x=-0.10, y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, 'fig3_coherence.png'),
            dpi=DPI, bbox_inches='tight')
plt.close(fig3)
print("  Saved fig3_coherence.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 4 — Survival validation
# ══════════════════════════════════════════════════════════════════════════════
print("Generating Figure 4...")

fig4 = plt.figure(figsize=(14, 9.0))
gs4  = gridspec.GridSpec(2, 2, figure=fig4, wspace=0.40, hspace=0.38,
                         left=0.04, right=0.97, top=0.94, bottom=0.06)

def embed_png(ax, path, title=''):
    img = np.array(Image.open(path).convert('RGB'))
    ax.imshow(img, aspect='auto', interpolation='lanczos')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines[:].set_visible(False)
    if title:
        ax.set_title(title, pad=4)

ax4a = fig4.add_subplot(gs4[0, 0])
embed_png(ax4a, TCGA_KM, 'TCGA-BRCA (PFI, $n$=1,211)')
add_panel_label(ax4a, 'A')

ax4b = fig4.add_subplot(gs4[0, 1])
embed_png(ax4b, MET_KM, 'METABRIC (OS, $n$=1,980)')
add_panel_label(ax4b, 'B')

ax4c = fig4.add_subplot(gs4[1, 0])
# CESC KM: only show the C1 panel (left half of the CESC figure)
cesc_img = np.array(Image.open(CESC_KM).convert('RGB'))
h, w = cesc_img.shape[:2]
ax4c.imshow(cesc_img[:, :w//2], aspect='auto', interpolation='lanczos')
ax4c.set_xticks([])
ax4c.set_yticks([])
ax4c.spines[:].set_visible(False)
ax4c.set_title('TCGA-CESC C1 (OS, $n$=304)', pad=4)
add_panel_label(ax4c, 'C')

# ── Panel D: forest plot ──────────────────────────────────────────────────────
ax4d = fig4.add_subplot(gs4[1, 1])
add_panel_label(ax4d, 'D')

forest_data = [
    # label                         HR      lo      hi      p
    ('TCGA-BRCA  C1',              0.915,  0.844,  0.992,  0.031),
    ('TCGA-BRCA  C4',              0.996,  0.886,  1.121,  0.953),
    ('METABRIC  C1',               0.855,  0.770,  0.948,  0.003),
    ('METABRIC  C4',               1.000,  0.876,  1.141,  0.999),
    ('TCGA-CESC  C1',              0.562,  0.345,  0.918,  0.021),
    ('TCGA-CESC  C4',              0.589,  0.359,  0.965,  0.036),
]

y_pos = list(range(len(forest_data)))[::-1]
colors_f = ['#029E73', '#aaaaaa', '#029E73', '#aaaaaa', '#029E73', '#aaaaaa']

ax4d.axvline(1.0, color='#999999', lw=0.8, ls='--')

for i, ((label, hr, lo, hi, p), yp, col) in enumerate(
        zip(forest_data, y_pos, colors_f)):
    ax4d.plot([lo, hi], [yp, yp], color=col, lw=2.0, solid_capstyle='round')
    ax4d.plot(hr, yp, 'o', color=col, ms=7, zorder=5)
    sig = '*' if p < 0.05 else 'ns'
    ax4d.text(max(hi, 1.0) + 0.02, yp, f' {hr:.3f} [{lo:.3f}–{hi:.3f}]  {sig}',
              va='center', fontsize=6.5, color='#222222')

ax4d.set_yticks(y_pos)
ax4d.set_yticklabels([d[0] for d in forest_data], fontsize=7.5)
ax4d.set_xlabel('Hazard ratio (per SD)', labelpad=4)
ax4d.set_title('Multivariate Cox HR\n(adj. age ± stage/PAM50)', pad=6)
ax4d.set_xlim(0.22, 1.45)
ax4d.spines['left'].set_visible(False)

plt.savefig(os.path.join(OUT_DIR, 'fig4_survival.png'),
            dpi=DPI, bbox_inches='tight')
plt.close(fig4)
print("  Saved fig4_survival.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 5 — Shiao ICB response
# ══════════════════════════════════════════════════════════════════════════════
print("Generating Figure 5...")

pb = pd.read_csv(SHIAO, sep='\t')
pre = pb[(pb['treatment'] == 'Base') & (pb['pCR'].isin(['R', 'NR']))].copy()
pre_rcb = pre.dropna(subset=['C1_score', 'RCB_num']).copy()

fig5, axes5 = plt.subplots(1, 2, figsize=(7.5, 4.2))

# ── Panel A: C1 score R vs NR boxplot ────────────────────────────────────────
ax5a = axes5[0]
groups = [pre[pre['pCR'] == 'R']['C1_score'].dropna().values,
          pre[pre['pCR'] == 'NR']['C1_score'].dropna().values]
bplot = ax5a.boxplot(groups, patch_artist=True, widths=0.45,
                     medianprops=dict(color='white', lw=2),
                     whiskerprops=dict(color='#444444'),
                     capprops=dict(color='#444444'),
                     flierprops=dict(marker='o', ms=4, alpha=0.5,
                                     markerfacecolor='#aaaaaa',
                                     markeredgecolor='none'))

box_colors = ['#029E73', '#d62728']
for patch, col in zip(bplot['boxes'], box_colors):
    patch.set_facecolor(col)
    patch.set_alpha(0.75)

# Jitter
rng = np.random.default_rng(42)
for i, (grp, col) in enumerate(zip(groups, box_colors)):
    jitter = rng.uniform(-0.12, 0.12, len(grp))
    ax5a.scatter(np.ones(len(grp)) * (i + 1) + jitter, grp,
                 c=col, s=25, alpha=0.6, zorder=5, edgecolors='white', lw=0.3)

stat, p_mwu = ss.mannwhitneyu(groups[0], groups[1], alternative='two-sided')
ax5a.set_xticks([1, 2])
ax5a.set_xticklabels([f'Responder\n(n={len(groups[0])})',
                       f'Non-responder\n(n={len(groups[1])})'], fontsize=8.5)
ax5a.set_ylabel('C1 pseudobulk score (pre-treatment)')
ax5a.set_title(f'C1 score vs pCR\n(MWU p={p_mwu:.3f})')
ax5a.axhline(0, color='#cccccc', lw=0.6, ls=':')
add_panel_label(ax5a, 'A')

# ── Panel B: C1 score vs RCB scatter ─────────────────────────────────────────
ax5b = axes5[1]
rho, p_rho = ss.spearmanr(pre_rcb['C1_score'], pre_rcb['RCB_num'])

# Jitter RCB for visibility
jitter_x = rng.uniform(-0.12, 0.12, len(pre_rcb))
sc = ax5b.scatter(pre_rcb['RCB_num'] + jitter_x, pre_rcb['C1_score'],
                  c=pre_rcb['pCR'].map({'R': '#029E73', 'NR': '#d62728'}),
                  s=40, alpha=0.75, edgecolors='white', lw=0.5, zorder=5)

# Regression line
x_fit = np.linspace(-0.1, 3.1, 50)
slope, intercept, *_ = ss.linregress(pre_rcb['RCB_num'], pre_rcb['C1_score'])
ax5b.plot(x_fit, slope * x_fit + intercept, color='#555555', lw=1.2, ls='--', zorder=3)

ax5b.set_xlabel('Residual cancer burden (RCB 0–3)')
ax5b.set_ylabel('C1 pseudobulk score (pre-treatment)')
ax5b.set_title(f'C1 score vs RCB\n(Spearman ρ={rho:.3f}, p={p_rho:.3f})')
ax5b.set_xticks([0, 1, 2, 3])
ax5b.axhline(0, color='#cccccc', lw=0.6, ls=':')

legend_patches = [mpatches.Patch(facecolor='#029E73', label='Responder (pCR)'),
                  mpatches.Patch(facecolor='#d62728', label='Non-responder')]
ax5b.legend(handles=legend_patches, fontsize=7.5, frameon=False, loc='upper right')
add_panel_label(ax5b, 'B')

plt.suptitle('Shiao et al. 2024 — TNBC ICB response (GSE246613, pre-treatment)',
             fontsize=9, y=1.01)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, 'fig5_shiao.png'),
            dpi=DPI, bbox_inches='tight')
plt.close(fig5)
print("  Saved fig5_shiao.png")

print("\nAll figures written to paper/figures/")
