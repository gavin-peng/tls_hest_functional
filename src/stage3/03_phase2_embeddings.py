"""
Phase 2 embedding extraction for Stage 3 GNN.

Extracts ResNet-50 avgpool embeddings (2048-dim) from HEST 224x224 spot patches
for all TLS instances needed for GNN training/evaluation.

For each instance, saves:
  embeddings  (N_spots, 2048)  float32
  coords_um   (N_spots, 2)     float32  : µm space, matches h5ad coord_um
  region      (N_spots,)       int8     : 0=interior (in_tls), 1=ring (neighborhood)
  platform    scalar           int8     : 0=Visium, 1=SPA

Output: outputs/stage3/phase2_embeddings.h5
  One group per instance_id. Includes all 670 instances (NaN embeddings for the
  2 instances without patches; these are Cluster 2 Visium, never in primary task).

Notes:
  - Embeddings are raw 2048-dim; PCA is applied per-fold inside the GNN script
    to avoid leakage.
  - ResNet-50 weights: IMAGENET1K_V2 (torchvision default for ResNet50).
  - Patches are ImageNet-normalised before inference (mean=[0.485,0.456,0.406],
    std=[0.229,0.224,0.225]) : standard for ImageNet-pretrained models on RGB.
  - Inference runs on GPU if available, else CPU.
  - Batch size 256 comfortably fits 6 GB VRAM.
"""
import warnings
import numpy as np
import pandas as pd
import h5py
import anndata as ad
from pathlib import Path
from scipy.spatial import cKDTree

import torch
import torch.nn as nn
from torchvision import models, transforms

warnings.filterwarnings('ignore')

ROOT      = Path('/mnt/e/hest')
OUT_DIR   = ROOT / 'outputs/stage3'
PATCH_DIR = ROOT / 'data/patches'
META_CSV  = ROOT / 'data/HEST_v1_3_0.csv'
SPATIAL_H5 = ROOT / 'outputs/stage2_clustering/tls_neighborhood_spots.h5ad'
CLUST_TSV  = ROOT / 'outputs/stage2_clustering/leiden_670_r0.3.tsv'

COORD_TOL_UM = 2.0
BATCH_SIZE   = 256

DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"Device: {DEVICE}")

# ── ResNet-50 encoder (avgpool output = 2048-dim) ─────────────────────────────
resnet = models.resnet50(weights=models.ResNet50_Weights.IMAGENET1K_V2)
resnet.eval()
# Remove the final FC layer; avgpool -> flatten gives 2048-dim
encoder = nn.Sequential(*list(resnet.children())[:-1], nn.Flatten())
encoder = encoder.to(DEVICE)

IMAGENET_NORM = transforms.Normalize(
    mean=[0.485, 0.456, 0.406],
    std=[0.229, 0.224, 0.225],
)


def embed_patches(imgs_uint8: np.ndarray) -> np.ndarray:
    """imgs_uint8: (N, 224, 224, 3) uint8 -> (N, 2048) float32 embeddings."""
    if imgs_uint8 is None or len(imgs_uint8) == 0:
        return np.empty((0, 2048), dtype=np.float32)
    results = []
    for start in range(0, len(imgs_uint8), BATCH_SIZE):
        batch = imgs_uint8[start:start + BATCH_SIZE]
        # (N, H, W, 3) uint8 -> (N, 3, H, W) float32 in [0,1]
        t = torch.from_numpy(batch).permute(0, 3, 1, 2).float() / 255.0
        t = IMAGENET_NORM(t).to(DEVICE)
        with torch.no_grad():
            emb = encoder(t).cpu().numpy()   # (N, 2048)
        results.append(emb)
    return np.concatenate(results, axis=0).astype(np.float32)


def build_patch_lookup(sample_id: str, pixel_size_um: float):
    """Same barcode-join lookup as Phase 1."""
    patch_path = PATCH_DIR / f'{sample_id}.h5'
    st_path    = ROOT / 'data' / 'st' / f'{sample_id}.h5ad'
    if not patch_path.exists() or not st_path.exists():
        return None
    with h5py.File(patch_path, 'r') as f:
        imgs     = f['img'][:]
        barcodes = [b[0].decode() if isinstance(b[0], bytes) else b[0]
                    for b in f['barcode']]
    bc_to_idx = {bc: i for i, bc in enumerate(barcodes)}
    st = ad.read_h5ad(str(st_path), backed='r')
    st_coords_um = st.obsm['spatial'] * pixel_size_um
    tree = cKDTree(st_coords_um)
    return {'tree': tree, 'st_barcodes': st.obs_names.tolist(),
            'bc_to_idx': bc_to_idx, 'imgs': imgs}


def query_patches_with_coords(lookup, query_coords_um: np.ndarray):
    """Returns (imgs, matched_coords_um) for spots within tolerance."""
    if lookup is None or len(query_coords_um) == 0:
        return None, np.empty((0, 2), dtype=np.float32)
    dists, idxs = lookup['tree'].query(query_coords_um, k=1)
    matched_imgs, matched_coords = [], []
    for dist, idx, coord in zip(dists, idxs, query_coords_um):
        if dist > COORD_TOL_UM:
            continue
        bc = lookup['st_barcodes'][idx]
        patch_idx = lookup['bc_to_idx'].get(bc)
        if patch_idx is not None:
            matched_imgs.append(lookup['imgs'][patch_idx])
            matched_coords.append(coord)
    if not matched_imgs:
        return None, np.empty((0, 2), dtype=np.float32)
    return np.array(matched_imgs), np.array(matched_coords, dtype=np.float32)


# ── Load metadata ─────────────────────────────────────────────────────────────
print("Loading metadata...")
clust = pd.read_csv(CLUST_TSV, sep='\t')
meta_hest = (pd.read_csv(META_CSV, usecols=['id', 'pixel_size_um_estimated'])
             .rename(columns={'id': 'sample_id'}).set_index('sample_id'))

spots = ad.read_h5ad(SPATIAL_H5, backed='r')
spots_obs = spots.obs[['sample_id', 'instance_id', 'spot_type',
                        'coord_x_um', 'coord_y_um']].copy()

platform_map = (clust.set_index('instance_id')['st_technology']
                .map({'Spatial Transcriptomics': 1, 'Visium': 0})
                .astype('int8'))

# ── Extract embeddings ─────────────────────────────────────────────────────────
OUT_DIR.mkdir(parents=True, exist_ok=True)
out_h5 = OUT_DIR / 'phase2_embeddings.h5'

samples_with_patches = [s for s in clust['sample_id'].unique()
                        if (PATCH_DIR / f'{s}.h5').exists()]
n_missing = clust['sample_id'].nunique() - len(samples_with_patches)
print(f"Samples with patches: {len(samples_with_patches)}/{clust['sample_id'].nunique()}"
      f" ({n_missing} missing)")

processed = 0
with h5py.File(out_h5, 'w') as hf:
    for sid in samples_with_patches:
        pixel_size = (meta_hest.loc[sid, 'pixel_size_um_estimated']
                      if sid in meta_hest.index else 0.5)
        lookup = build_patch_lookup(sid, pixel_size)
        if lookup is None:
            continue

        sample_spots = spots_obs[spots_obs['sample_id'] == sid]
        instances    = clust.loc[clust['sample_id'] == sid, 'instance_id'].values

        for inst_id in instances:
            inst_spots = sample_spots[sample_spots['instance_id'] == inst_id]
            coords_A = inst_spots.loc[inst_spots['spot_type'] == 'in_tls',
                                      ['coord_x_um', 'coord_y_um']].values
            coords_B = inst_spots.loc[inst_spots['spot_type'] == 'neighborhood',
                                      ['coord_x_um', 'coord_y_um']].values

            imgs_A, mcoords_A = query_patches_with_coords(lookup, coords_A)
            imgs_B, mcoords_B = query_patches_with_coords(lookup, coords_B)

            emb_A = embed_patches(imgs_A)
            emb_B = embed_patches(imgs_B)

            n_A, n_B = len(emb_A), len(emb_B)
            if n_A + n_B == 0:
                continue

            embeddings = np.concatenate([emb_A, emb_B], axis=0)  # (N, 2048)
            coords_um  = np.concatenate([mcoords_A, mcoords_B], axis=0)
            region     = np.array([0]*n_A + [1]*n_B, dtype=np.int8)

            grp = hf.require_group(inst_id)
            grp.create_dataset('embeddings', data=embeddings,  compression='lzf')
            grp.create_dataset('coords_um',  data=coords_um,   compression='lzf')
            grp.create_dataset('region',     data=region,      compression='lzf')
            grp.attrs['platform'] = int(platform_map.get(inst_id, 0))
            grp.attrs['n_interior'] = n_A
            grp.attrs['n_ring']     = n_B

            processed += 1

        print(f"  {sid}: {len(instances)} instances", flush=True)

print(f"\nDone. {processed} instances written to {out_h5}")
print(f"File size: {out_h5.stat().st_size / 1e6:.1f} MB")
