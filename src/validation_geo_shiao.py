"""
Shiao et al. Cancer Cell 2024 (GSE246613) — scRNA-seq ICB response validation (Claim 1e).

Tests whether C1 (immunogenic TLS) and C4 (suppressive) gene module scores from
pre-treatment biopsies predict pCR and RCB in TNBC patients treated with pembrolizumab
± radiotherapy.

Dataset:
  35 patients, 3 timepoints (Base / PD1 / RTPD1), human immune cells
  pCR: R (n=21) vs NR (n=8); RCB 0–3 (0=pCR, 3=extensive residual)
  File: GSE246613_immune.h5ad (342k cells × 36k genes, CSR sparse, raw counts)
  obs: cohort=patient, treatment=timepoint, pCR/RCB already embedded

Memory strategy: backed='r' → subset to target + background genes → to_memory().
Loading ~500 genes instead of 36k reduces peak RAM from ~4 GB to <200 MB.
"""
import warnings
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scipy.stats as stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')
sc.settings.verbosity = 1

DATA = "/mnt/e/hest/outputs/validation_geo"
OUT  = DATA

C1_GENES = ["CXCL13","CCL19","CCL21","LAMP3","ACKR1","MS4A1","CXCR5",
            "BCL6","IGHG1","JCHAIN","MZB1","TNFSF13B","AICDA"]
C4_GENES = ["TGFB1","IL10","ARG1","MRC1","IDO1","FOXP3","IL2RA","CTLA4","CD274"]
TARGET_GENES = list(dict.fromkeys(C1_GENES + C4_GENES))  # deduped, order-preserving

# ── 1. Memory-mapped load — no X in RAM yet ───────────────────────────────
print("Memory-mapping immune h5ad (backed='r')...")
adata_full = ad.read_h5ad(f"{DATA}/GSE246613_immune.h5ad", backed='r')
print(f"  Shape: {adata_full.shape}")
print(f"  obs columns: {list(adata_full.obs.columns)}")

# ── 2. Identify column names (obs uses cohort/treatment, not Patient_Number) ─
pat_col = next((c for c in ["cohort","Patient_Number","patient","patient_id","batch"]
                if c in adata_full.obs.columns), None)
tp_col  = next((c for c in ["treatment","timepoint","time_point","visit"]
                if c in adata_full.obs.columns), None)
print(f"  patient col={pat_col!r}  timepoint col={tp_col!r}")

# ── 3. Score genes — use cache to skip expensive gene loading on restart ──
import os
SCORED_CACHE = f"{OUT}/shiao_scored_cells.csv"

rng = np.random.default_rng(42)

if os.path.exists(SCORED_CACHE):
    print(f"\n  Cached scores found — skipping gene matrix load.")
    scored = pd.read_csv(SCORED_CACHE, index_col=0)
    # Only need obs metadata; close backed file, build lightweight obs df
    obs_df = adata_full.obs[[pat_col, tp_col, "pCR", "RCB"]].copy() if (pat_col and tp_col) \
             else adata_full.obs[["pCR","RCB"]].copy()
    adata_full.file.close()
    obs_df["C1_score"] = scored["C1_score"].values
    obs_df["C4_score"] = scored["C4_score"].values
    c1_avail = C1_GENES  # all present (confirmed on first run)
    c4_avail = C4_GENES
else:
    target_set = set(TARGET_GENES)
    c1_avail = [g for g in C1_GENES if g in adata_full.var_names]
    c4_avail = [g for g in C4_GENES if g in adata_full.var_names]
    print(f"\n  C1 genes found: {len(c1_avail)}/{len(C1_GENES)}: {c1_avail}")
    print(f"  C4 genes found: {len(c4_avail)}/{len(C4_GENES)}: {c4_avail}")
    if len(c1_avail) == 0:
        print(f"  WARNING: sample var_names: {list(adata_full.var_names[:20])}")

    non_target = [g for g in adata_full.var_names if g not in target_set]
    bg_genes = rng.choice(non_target, size=min(500, len(non_target)),
                          replace=False).tolist()
    genes_to_load = [g for g in TARGET_GENES if g in adata_full.var_names] + bg_genes
    print(f"\n  Loading {len(genes_to_load)} genes out of {adata_full.n_vars}...")
    expr = adata_full[:, genes_to_load].to_memory()
    adata_full.file.close()
    print(f"  Loaded shape: {expr.shape}")

    # ── 4. Normalise ──────────────────────────────────────────────────────
    print("  Normalising: total→1e4, log1p...")
    sc.pp.normalize_total(expr, target_sum=1e4)
    sc.pp.log1p(expr)

    if len(c1_avail) > 0:
        sc.tl.score_genes(expr, c1_avail, score_name="C1_score", random_state=0)
    else:
        expr.obs["C1_score"] = np.nan
    if len(c4_avail) > 0:
        sc.tl.score_genes(expr, c4_avail, score_name="C4_score", random_state=0)
    else:
        expr.obs["C4_score"] = np.nan

    expr.obs[["C1_score","C4_score"]].to_csv(SCORED_CACHE)
    print(f"\n  Saved cell scores → {SCORED_CACHE}")
    obs_df = expr.obs[[pat_col, tp_col, "pCR", "RCB", "C1_score", "C4_score"]].copy() \
             if (pat_col and tp_col) else \
             expr.obs[["pCR","RCB","C1_score","C4_score"]].copy()

# ── 5. Build analysis obs (pCR/RCB already in obs_df) ────────────────────
adata_obs = obs_df.copy()

# Optionally enrich with clinical TSV (adds numeric RCB if not already present)
try:
    clin = pd.read_csv(f"{DATA}/shiao_patient_clinical.tsv", sep="\t")
    clin["RCB_num"] = pd.to_numeric(clin["RCB"], errors="coerce")
    print(f"\n  Clinical TSV: {len(clin)} patients, "
          f"pCR={clin.pCR.value_counts().to_dict()}")
    if pat_col and "Patient_Number" in clin.columns:
        adata_obs = adata_obs.merge(
            clin[["Patient_Number","RCB_num"]],
            left_on=pat_col, right_on="Patient_Number", how="left"
        )
except FileNotFoundError:
    print("  Clinical TSV not found — using RCB from obs directly")
    adata_obs["RCB_num"] = pd.to_numeric(adata_obs["RCB"], errors="coerce")

# ── 7. Pseudo-bulk: mean C1/C4 per patient × timepoint ───────────────────
group_cols = [c for c in [pat_col, tp_col] if c is not None]
print(f"\n  Grouping by: {group_cols}")

pb = adata_obs.groupby(group_cols)[["C1_score","C4_score"]].mean().reset_index()
pb["n_cells"] = adata_obs.groupby(group_cols).size().values

# Carry over pCR / RCB_num (one value per patient, constant within group)
for col in ["pCR","RCB_num"]:
    if col in adata_obs.columns:
        pb = pb.merge(
            adata_obs.groupby(group_cols)[col].first().reset_index(),
            on=group_cols, how="left"
        )

print(f"  Pseudo-bulk rows: {len(pb)}")
cols_show = group_cols + ["C1_score","C4_score","pCR"] + (["RCB_num"] if "RCB_num" in pb else [])
print(pb[cols_show].head(10).to_string())

pb.to_csv(f"{OUT}/shiao_C1C4_pseudobulk.tsv", sep="\t", index=False)
print(f"\n  Saved shiao_C1C4_pseudobulk.tsv")

# ── 8. Pre-treatment subset ───────────────────────────────────────────────
if tp_col and tp_col in pb.columns:
    pre_vals = [v for v in pb[tp_col].unique()
                if any(k in str(v).lower() for k in ["base","pre","a_p","tp1","_0","0h"])]
    if not pre_vals:
        pre_vals = [pb[tp_col].value_counts().index[0]]
    print(f"\n  Pre-treatment values: {pre_vals}")
    pb_pre = pb[pb[tp_col].isin(pre_vals)].copy()
else:
    pb_pre = pb.copy()

pb_pre_known = pb_pre[pb_pre["pCR"].isin(["R","NR"])].copy()
print(f"  Pre-treatment samples with known pCR: {len(pb_pre_known)} "
      f"(R={(pb_pre_known.pCR=='R').sum()}, NR={(pb_pre_known.pCR=='NR').sum()})")

# ── 9. Statistical tests ─────────────────────────────────────────────────
print("\n--- Pre-treatment C1/C4 vs pCR ---")
for score in ["C1_score","C4_score"]:
    r  = pb_pre_known[pb_pre_known["pCR"]=="R"][score].dropna()
    nr = pb_pre_known[pb_pre_known["pCR"]=="NR"][score].dropna()
    if len(r) > 1 and len(nr) > 1:
        stat, p = stats.mannwhitneyu(r, nr, alternative="two-sided")
        print(f"  {score}: R median={r.median():.3f}  NR median={nr.median():.3f}  "
              f"MWU p={p:.4f}")
    else:
        print(f"  {score}: insufficient data (R={len(r)}, NR={len(nr)})")

if "RCB_num" in pb_pre.columns:
    pb_rcb = pb_pre.dropna(subset=["RCB_num"]).copy()
    if len(pb_rcb) > 5:
        print("\n--- Pre-treatment score vs RCB (Spearman) ---")
        for score in ["C1_score","C4_score"]:
            d = pb_rcb[[score,"RCB_num"]].dropna()
            if len(d) > 5:
                r, p = stats.spearmanr(d[score], d["RCB_num"])
                print(f"  {score} ~ RCB: rho={r:.3f} p={p:.4f}")
else:
    pb_rcb = pd.DataFrame()

# ── 10. Figures ───────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 3, figsize=(15, 9))

for ax, score in zip(axes[0, :2], ["C1_score","C4_score"]):
    d = pb_pre_known.dropna(subset=[score])
    r_vals  = d[d["pCR"]=="R"][score]
    nr_vals = d[d["pCR"]=="NR"][score]
    ax.boxplot([r_vals, nr_vals],
               labels=[f"R (n={len(r_vals)})", f"NR (n={len(nr_vals)})"],
               patch_artist=True,
               boxprops=dict(facecolor="#2ca02c", alpha=0.7))
    title = f"{score} — pre-treatment"
    if len(r_vals) > 1 and len(nr_vals) > 1:
        _, p = stats.mannwhitneyu(r_vals, nr_vals, alternative="two-sided")
        title += f"\nMWU p={p:.3f}"
    ax.set_title(title, fontsize=10)
    ax.set_ylabel("Score")

ax = axes[0, 2]
if len(pb_rcb) > 3:
    d = pb_rcb.dropna(subset=["C1_score","RCB_num"])
    ax.scatter(d["RCB_num"] + rng.standard_normal(len(d))*0.05,
               d["C1_score"], alpha=0.7, color="#2ca02c")
    r, p = stats.spearmanr(d["C1_score"], d["RCB_num"])
    ax.set_title(f"C1 vs RCB (pre-tx)\nSpearman rho={r:.2f} p={p:.3f}", fontsize=10)
    ax.set_xlabel("RCB score")
    ax.set_ylabel("C1_score")
    ax.set_xticks([0,1,2,3])
else:
    ax.set_visible(False)

if tp_col and tp_col in pb.columns:
    timepoints = sorted(pb[tp_col].unique())
    for i, tp in enumerate(timepoints[:3]):
        ax = axes[1, i]
        sub = pb[(pb[tp_col]==tp) & pb["pCR"].isin(["R","NR"])].copy()
        r_vals  = sub[sub["pCR"]=="R"]["C1_score"].dropna()
        nr_vals = sub[sub["pCR"]=="NR"]["C1_score"].dropna()
        if len(r_vals) > 0 or len(nr_vals) > 0:
            ax.boxplot([r_vals  if len(r_vals)  else [np.nan],
                        nr_vals if len(nr_vals) else [np.nan]],
                       labels=[f"R (n={len(r_vals)})", f"NR (n={len(nr_vals)})"],
                       patch_artist=True)
            title = f"C1_score — {tp}"
            if len(r_vals) > 1 and len(nr_vals) > 1:
                _, p = stats.mannwhitneyu(r_vals, nr_vals, alternative="two-sided")
                title += f"\np={p:.3f}"
            ax.set_title(title, fontsize=10)
            ax.set_ylabel("C1_score")
        else:
            ax.set_visible(False)
else:
    for ax in axes[1]:
        ax.set_visible(False)

plt.suptitle("Shiao et al. 2024 — C1/C4 TLS scores vs ICB response (TNBC, n=35)",
             fontsize=12, y=1.01)
plt.tight_layout()
plt.savefig(f"{OUT}/shiao_TLS_ICB.png", dpi=150, bbox_inches="tight")
plt.close()
print("\nSaved shiao_TLS_ICB.png")

# ── 11. Summary ──────────────────────────────────────────────────────────
print("\n=== SHIAO GSE246613 VALIDATION SUMMARY ===")
print(f"n={len(pb_pre_known)} pre-treatment samples with pCR known")
for score in ["C1_score","C4_score"]:
    r  = pb_pre_known[pb_pre_known["pCR"]=="R"][score].dropna()
    nr = pb_pre_known[pb_pre_known["pCR"]=="NR"][score].dropna()
    if len(r) > 1 and len(nr) > 1:
        _, p = stats.mannwhitneyu(r, nr, alternative="two-sided")
        print(f"  {score}: R={r.median():.3f} NR={nr.median():.3f} MWU p={p:.4f}")
