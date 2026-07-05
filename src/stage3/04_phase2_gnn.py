"""
Phase 2 GNN classifier : GAT on spatial patch graphs.

Reads per-instance ResNet-50 embeddings from phase2_embeddings.h5, builds
two-region spatial graphs (interior + ring nodes), and trains a 2-layer GAT
with dual region pooling.

Node features (66-dim per node):
  64-dim  PCA of ResNet-50 2048-dim embedding (PCA fitted on training fold only)
  1-dim   region indicator: 0=interior, 1=ring
  1-dim   platform indicator: 0=Visium, 1=SPA

Graph edges:
  k=6 spatial kNN within interior nodes
  k=6 spatial kNN within ring nodes
  k=3 cross-boundary: each interior node -> 3 nearest ring nodes (undirected)

Architecture:
  2-layer GAT (8 heads, hidden=64, concat -> 512 after layer 1, mean -> 64 after layer 2)
  Separate mean-pooling over interior and ring nodes -> 128-dim concatenation
  MLP head: 128 -> 64 -> 1 (binary CE for C1 vs C4)

Regularisation:
  Dropout 0.3 on attention weights
  L2 weight decay 1e-4
  Coordinate jitter augmentation (training only): ±5% of std of coords

Evaluation: same 5-metric framework as Phase 1 (see 02_phase1_classifier.py).
PCA is fitted on training instances only in each fold -- no leakage.

Output:
  outputs/stage3/phase2_gnn_results.txt
  outputs/stage3/phase2_attention_weights.h5  (per-instance node attention scores)
"""
import warnings
import numpy as np
import pandas as pd
import h5py
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score, f1_score
from sklearn.preprocessing import LabelEncoder
from scipy.spatial import cKDTree

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.data import Data, Batch
from torch_geometric.nn import GATConv, global_mean_pool

warnings.filterwarnings('ignore')

ROOT      = Path('/mnt/e/hest')
OUT_DIR   = ROOT / 'outputs/stage3'
EMB_H5    = OUT_DIR / 'phase2_embeddings.h5'
CLUST_TSV = ROOT / 'outputs/stage2_clustering/leiden_670_r0.3.tsv'

DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"Device: {DEVICE}")

PCA_DIM    = 64
KNN_INNER  = 6
KNN_RING   = 6
KNN_CROSS  = 3
NODE_DIM   = PCA_DIM + 2   # 64 + region + platform

# ── GAT model ─────────────────────────────────────────────────────────────────
class TLS_GAT(nn.Module):
    def __init__(self, in_dim=NODE_DIM, hidden=64, heads=8, dropout=0.3):
        super().__init__()
        self.dropout = dropout
        # Layer 1: concat heads -> heads*hidden
        self.gat1 = GATConv(in_dim, hidden, heads=heads, dropout=dropout,
                             concat=True)
        # Layer 2: mean over heads -> hidden
        self.gat2 = GATConv(heads * hidden, hidden, heads=1, dropout=dropout,
                             concat=False)
        # MLP after dual pooling (interior mean || ring mean)
        self.head = nn.Sequential(
            nn.Linear(2 * hidden, 64),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(64, 1),
        )

    def forward(self, data):
        x, edge_index = data.x, data.edge_index
        # Region mask , 1-indexed batch, need per-graph interior/ring masks
        region = data.region          # (total_nodes,) 0=interior 1=ring
        batch  = data.batch           # (total_nodes,)

        x = F.dropout(x, p=self.dropout, training=self.training)
        x = F.elu(self.gat1(x, edge_index))
        x = F.dropout(x, p=self.dropout, training=self.training)
        x = self.gat2(x, edge_index)   # (total_nodes, hidden)

        # Dual pooling: mean over interior, mean over ring
        interior_mask = (region == 0)
        ring_mask     = (region == 1)

        pool_int  = _masked_mean_pool(x, batch, interior_mask, data.num_graphs)
        pool_ring = _masked_mean_pool(x, batch, ring_mask,     data.num_graphs)

        graph_emb = torch.cat([pool_int, pool_ring], dim=1)  # (B, 2*hidden)
        return self.head(graph_emb).squeeze(1)               # (B,)


def _masked_mean_pool(x, batch, mask, num_graphs):
    """Mean pool over nodes where mask is True, per graph. Zero if no nodes."""
    out = torch.zeros(num_graphs, x.size(1), device=x.device)
    cnt = torch.zeros(num_graphs, 1,         device=x.device)
    sel_x     = x[mask]
    sel_batch = batch[mask]
    out.scatter_add_(0, sel_batch.unsqueeze(1).expand_as(sel_x), sel_x)
    cnt.scatter_add_(0, sel_batch.unsqueeze(1),
                     torch.ones(sel_batch.size(0), 1, device=x.device))
    cnt = cnt.clamp(min=1)
    return out / cnt


# ── Graph construction ────────────────────────────────────────────────────────
def build_graph(embeddings_pca: np.ndarray,
                coords_um: np.ndarray,
                region: np.ndarray,
                platform: int,
                augment: bool = False) -> Data:
    """Build a PyG Data object for one TLS instance.

    embeddings_pca: (N, 64)
    coords_um:      (N, 2)
    region:         (N,)  0=interior 1=ring
    platform:       int   0=Visium 1=SPA
    """
    N = len(embeddings_pca)
    if N == 0:
        return None

    if augment:
        # Random rotation around centroid + small coord jitter
        theta = np.random.uniform(0, 2 * np.pi)
        R = np.array([[np.cos(theta), -np.sin(theta)],
                      [np.sin(theta),  np.cos(theta)]], dtype=np.float32)
        centroid = coords_um.mean(axis=0, keepdims=True)
        coords_um = (coords_um - centroid) @ R.T + centroid
        jitter = np.random.randn(*coords_um.shape).astype(np.float32)
        coords_um = coords_um + jitter * 0.05 * coords_um.std(axis=0).clip(1)

    # Node features: PCA embedding + region indicator + platform indicator
    region_feat   = region.astype(np.float32).reshape(-1, 1)
    platform_feat = np.full((N, 1), platform, dtype=np.float32)
    x = np.concatenate([embeddings_pca, region_feat, platform_feat], axis=1)

    # Build edges
    int_idx  = np.where(region == 0)[0]
    ring_idx = np.where(region == 1)[0]
    edges = []

    if len(int_idx) > 1:
        edges.extend(_knn_edges(coords_um[int_idx], int_idx, KNN_INNER))
    if len(ring_idx) > 1:
        edges.extend(_knn_edges(coords_um[ring_idx], ring_idx, KNN_RING))
    if len(int_idx) > 0 and len(ring_idx) > 0:
        edges.extend(_cross_edges(coords_um[int_idx], int_idx,
                                  coords_um[ring_idx], ring_idx, KNN_CROSS))

    if edges:
        edge_index = torch.tensor(np.array(edges).T, dtype=torch.long)
    else:
        edge_index = torch.zeros((2, 0), dtype=torch.long)

    return Data(
        x=torch.tensor(x, dtype=torch.float32),
        edge_index=edge_index,
        region=torch.tensor(region, dtype=torch.long),
    )


def _knn_edges(coords, global_idx, k):
    """Undirected kNN edges within a set of nodes."""
    k_eff = min(k, len(coords) - 1)
    if k_eff < 1:
        return []
    tree = cKDTree(coords)
    dists, nbrs = tree.query(coords, k=k_eff + 1)
    edges = set()
    for i, row in enumerate(nbrs):
        for j in row[1:]:
            u, v = int(global_idx[i]), int(global_idx[j])
            edges.add((min(u, v), max(u, v)))
    return [(u, v) for u, v in edges] + [(v, u) for u, v in edges]


def _cross_edges(int_coords, int_idx, ring_coords, ring_idx, k):
    """Each interior node connects to its k nearest ring nodes (undirected)."""
    k_eff = min(k, len(ring_coords))
    if k_eff < 1:
        return []
    tree = cKDTree(ring_coords)
    _, nbrs = tree.query(int_coords, k=k_eff)
    if nbrs.ndim == 1:
        nbrs = nbrs.reshape(-1, 1)
    edges = set()
    for i, row in enumerate(nbrs):
        for j in row:
            u, v = int(int_idx[i]), int(ring_idx[j])
            edges.add((u, v))
    return [(u, v) for u, v in edges] + [(v, u) for u, v in edges]


# ── Load embeddings from h5 ───────────────────────────────────────────────────
def load_instance_data(instance_ids):
    """Returns dict: instance_id -> {embeddings, coords_um, region, platform}."""
    data = {}
    with h5py.File(EMB_H5, 'r') as hf:
        for iid in instance_ids:
            if iid not in hf:
                continue
            grp = hf[iid]
            data[iid] = {
                'embeddings': grp['embeddings'][:],
                'coords_um':  grp['coords_um'][:],
                'region':     grp['region'][:],
                'platform':   int(grp.attrs['platform']),
            }
    return data


# ── PCA fit / transform (training fold only) ─────────────────────────────────
def fit_pca(train_data: dict) -> PCA:
    """Fit PCA on all patch embeddings from training instances."""
    all_emb = np.concatenate([d['embeddings'] for d in train_data.values()], axis=0)
    pca = PCA(n_components=PCA_DIM, random_state=42)
    pca.fit(all_emb)
    return pca


def apply_pca(pca: PCA, instance_data: dict) -> dict:
    """Transform embeddings of all instances using a fitted PCA."""
    out = {}
    for iid, d in instance_data.items():
        emb_pca = pca.transform(d['embeddings']).astype(np.float32)
        out[iid] = {**d, 'embeddings_pca': emb_pca}
    return out


# ── Training / evaluation ─────────────────────────────────────────────────────
def train_epoch(model, graphs, labels, optimiser, augment=True):
    model.train()
    perm = np.random.permutation(len(graphs))
    total_loss = 0.0
    for idx in perm:
        g, y = graphs[idx], labels[idx]
        if g is None:
            continue
        if augment:
            # Rebuild graph with coordinate augmentation (jitter already applied
            # in build_graph; here we rely on the stored augment=True graphs)
            pass
        batch = Batch.from_data_list([g]).to(DEVICE)
        optimiser.zero_grad()
        logit = model(batch)
        loss  = F.binary_cross_entropy_with_logits(logit,
                    torch.tensor([y], dtype=torch.float32, device=DEVICE))
        loss.backward()
        optimiser.step()
        total_loss += loss.item()
    return total_loss / max(len(graphs), 1)


@torch.no_grad()
def predict(model, graphs):
    model.eval()
    probas = []
    for g in graphs:
        if g is None:
            probas.append(0.5)
            continue
        batch = Batch.from_data_list([g]).to(DEVICE)
        logit = model(batch)
        probas.append(torch.sigmoid(logit).item())
    return np.array(probas)


# ── LOPO-CV (mirrors Phase 1 scheme exactly) ─────────────────────────────────
def lopo_c4_auroc(df: pd.DataFrame, label_col: str,
                  c4_patients: list, all_data: dict,
                  n_epochs: int = 80) -> dict:
    rng = np.random.default_rng(42)
    neg_patients = sorted(df.loc[df[label_col] == 0, 'sample_id'].unique())
    rng.shuffle(neg_patients)
    neg_fold_assignment = {pid: i % len(c4_patients)
                           for i, pid in enumerate(neg_patients)}

    fold_scores = []
    for fold_idx, pid in enumerate(c4_patients):
        test_neg_pids = [p for p, fi in neg_fold_assignment.items() if fi == fold_idx]
        test_pids  = {pid} | set(test_neg_pids)
        test_mask  = df['sample_id'].isin(test_pids)
        train_mask = ~test_mask

        train_df = df[train_mask]
        test_df  = df[test_mask]

        if len(train_df[label_col].unique()) < 2 or len(test_df[label_col].unique()) < 2:
            continue

        train_ids = train_df['instance_id'].tolist()
        test_ids  = test_df['instance_id'].tolist()

        # Load and PCA-transform
        train_raw = {iid: all_data[iid] for iid in train_ids if iid in all_data}
        test_raw  = {iid: all_data[iid] for iid in test_ids  if iid in all_data}
        if not train_raw or not test_raw:
            continue

        pca = fit_pca(train_raw)
        train_pca = apply_pca(pca, train_raw)
        test_pca  = apply_pca(pca, test_raw)

        # Build graphs
        train_graphs = [build_graph(d['embeddings_pca'], d['coords_um'],
                                    d['region'], d['platform'], augment=True)
                        for d in train_pca.values()]
        train_labels = train_df.set_index('instance_id').loc[
            [iid for iid in train_ids if iid in train_pca], label_col].tolist()
        test_graphs  = [build_graph(d['embeddings_pca'], d['coords_um'],
                                    d['region'], d['platform'], augment=False)
                        for d in test_pca.values()]
        test_labels  = test_df.set_index('instance_id').loc[
            [iid for iid in test_ids if iid in test_pca], label_col].tolist()

        # Filter None graphs
        valid_train = [(g, y) for g, y in zip(train_graphs, train_labels) if g is not None]
        valid_test  = [(g, y) for g, y in zip(test_graphs,  test_labels)  if g is not None]
        if not valid_train or not valid_test:
            continue
        tr_graphs, tr_labels = zip(*valid_train)
        te_graphs, te_labels = zip(*valid_test)

        if len(set(te_labels)) < 2:
            continue

        # Train
        model = TLS_GAT(in_dim=NODE_DIM).to(DEVICE)
        optim = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-4)
        for epoch in range(n_epochs):
            train_epoch(model, list(tr_graphs), list(tr_labels), optim)

        probas = predict(model, list(te_graphs))
        score  = roc_auc_score(list(te_labels), probas)
        fold_scores.append(score)
        print(f"  Fold {fold_idx} (left out: {pid}, n_test={len(te_labels)}): "
              f"AUROC={score:.3f}", flush=True)

    if not fold_scores:
        return {'mean': np.nan, 'min': np.nan, 'max': np.nan, 'folds': []}
    return {'mean': float(np.mean(fold_scores)),
            'min':  float(np.min(fold_scores)),
            'max':  float(np.max(fold_scores)),
            'folds': fold_scores}


# ── Load cluster labels ───────────────────────────────────────────────────────
print("Loading cluster labels...")
clust = pd.read_csv(CLUST_TSV, sep='\t')

df_c4     = clust[clust['leiden_cluster'] == 4].copy()
df_c1     = clust[clust['leiden_cluster'] == 1].copy()
df_c1_spa = df_c1[df_c1['st_technology'] == 'Spatial Transcriptomics'].copy()
c4_patients = sorted(df_c4['sample_id'].unique())

print(f"C4 patients ({len(c4_patients)}): {c4_patients}")
print(f"C1-SPA: {len(df_c1_spa)} instances, {df_c1_spa['sample_id'].nunique()} patients")
print()

# Load all embeddings for C1-SPA + C4 instances
print("Loading embeddings from h5...")
primary_ids = set(df_c1_spa['instance_id']) | set(df_c4['instance_id'])
all_data = load_instance_data(primary_ids)
print(f"  Loaded {len(all_data)}/{len(primary_ids)} instances (others missing patches)")
print()

# ── METRIC 1: C1-SPA vs C4 AUROC (PRIMARY) ───────────────────────────────────
print("=" * 65)
print("METRIC 1: C1-SPA vs C4 AUROC (PRIMARY, biology signal)")
print("=" * 65)
df_m1 = pd.concat([df_c1_spa.assign(label=0, instance_id=df_c1_spa['instance_id']),
                   df_c4.assign(label=1, instance_id=df_c4['instance_id'])])
df_m1 = df_m1.rename(columns={'leiden_cluster': 'cluster'})
m1 = lopo_c4_auroc(df_m1, 'label', c4_patients, all_data)
print(f"  AUROC: {m1['mean']:.3f} (range [{m1['min']:.3f}, {m1['max']:.3f}])")
print()

# ── SUMMARY ───────────────────────────────────────────────────────────────────
print("=" * 65)
print("SUMMARY")
print("=" * 65)
print(f"  Phase 1 RF baseline  C1-SPA vs C4:  0.744  [0.000–1.000]")
print(f"  Phase 2 GAT primary  C1-SPA vs C4:  {m1['mean']:.3f}  [{m1['min']:.3f}–{m1['max']:.3f}]")
delta = m1['mean'] - 0.744
print(f"  Delta vs baseline: {delta:+.3f}")
print()
primary = m1['mean']
if primary > 0.744 + 0.02:
    verdict = "GNN IMPROVES over RF baseline : spatial organisation adds signal"
elif primary >= 0.744 - 0.02:
    verdict = "GNN MATCHES RF : aggregate features capture most discriminative information"
else:
    verdict = "GNN UNDERPERFORMS RF : likely overfitting at this sample size"
print(f"  Verdict: {verdict}")
