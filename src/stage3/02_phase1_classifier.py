"""
Phase 1 morphological classifier : Random Forest evaluation.

Implements the five-metric evaluation framework from implementation_notes.txt:
  1. C1-SPA vs C4 AUROC  (n≈182, SPA only, LOPO-CV over C4 patients) : PRIMARY
  2. C1 vs C4 AUROC full (n≈238, mixed platform, LOPO-CV over C4 patients) : secondary
  3. C1-Visium vs C1-SPA AUROC (n=168, C1 only, LOPO-CV* over C1-SPA patients) : platform ctrl
  4. C1 vs C0 AUROC      (n≈359, mixed, 3-fold patient CV) : ST necessity
  5. 6-class macro F1    (n=670, 3-fold patient CV) : context

All CV splits are patient-level (split by sample_id) to avoid data leakage.
Missing features (circularity/convexity for split instances, rare H&E NaN) are
median-imputed within each training fold and applied to test fold.

Decision gates:
  primary AUROC > 0.70 -> proceed to Phase 2 GNN
  primary AUROC 0.60–0.70 -> cautious proceed
  primary AUROC < 0.60 -> stop; confirms ST necessity

Output: outputs/stage3/phase1_classifier_results.txt
        outputs/stage3/phase1_feature_importance.tsv
"""
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, f1_score
from sklearn.preprocessing import LabelEncoder
from sklearn.impute import SimpleImputer

warnings.filterwarnings('ignore')

ROOT    = Path('/mnt/e/hest')
OUT_DIR = ROOT / 'outputs/stage3'

FEAT_COLS = [
    'log_area_um2', 'elongation', 'circularity', 'convexity',
    'he_mean_A', 'he_var_A', 'he_mean_B', 'he_var_B', 'he_ratio_BA',
]

RF_PARAMS = dict(n_estimators=500, max_features='sqrt',
                 class_weight='balanced', random_state=42, n_jobs=-1)


def lopo_c4_auroc(df: pd.DataFrame, label_col: str, c4_patients: list) -> dict:
    """LOPO-CV over C4 patients for binary C1 vs C4 classification.

    Each fold i leaves out ALL instances from one C4 patient (test positives).
    To ensure the test set has both classes, we also hold out a proportional
    sample of negative-class (C1) patients: ceil(n_c1_patients / n_c4_patients)
    patients per fold, assigned by round-robin with a fixed seed.

    Returns dict with mean, min, max, and per-fold AUROC scores.
    """
    rng = np.random.default_rng(42)

    # Patients with label=0 (C1 / negative class)
    neg_patients = sorted(df.loc[df[label_col] == 0, 'sample_id'].unique())
    rng.shuffle(neg_patients)

    n_c4 = len(c4_patients)
    n_neg = len(neg_patients)
    # Assign neg patients round-robin to folds so each fold holds out some negatives
    neg_fold_assignment = {pid: i % n_c4 for i, pid in enumerate(neg_patients)}

    fold_scores = []
    for fold_idx, pid in enumerate(c4_patients):
        # test: this C4 patient + C1 patients assigned to this fold
        test_neg_pids = [p for p, fi in neg_fold_assignment.items() if fi == fold_idx]
        test_pids = {pid} | set(test_neg_pids)

        test_mask  = df['sample_id'].isin(test_pids)
        train_mask = ~test_mask

        X_train = df.loc[train_mask, FEAT_COLS].values
        y_train = df.loc[train_mask, label_col].values
        X_test  = df.loc[test_mask, FEAT_COLS].values
        y_test  = df.loc[test_mask, label_col].values

        if len(np.unique(y_test)) < 2 or len(np.unique(y_train)) < 2:
            continue

        imp = SimpleImputer(strategy='median')
        X_train = imp.fit_transform(X_train)
        X_test  = imp.transform(X_test)

        rf = RandomForestClassifier(**RF_PARAMS)
        rf.fit(X_train, y_train)
        pos_idx = list(rf.classes_).index(1)
        proba = rf.predict_proba(X_test)[:, pos_idx]
        fold_scores.append(roc_auc_score(y_test, proba))

    if not fold_scores:
        return {'mean': np.nan, 'min': np.nan, 'max': np.nan, 'folds': []}
    return {
        'mean': float(np.mean(fold_scores)),
        'min':  float(np.min(fold_scores)),
        'max':  float(np.max(fold_scores)),
        'folds': fold_scores,
    }


def kfold_patient_auroc(df: pd.DataFrame, label_col: str, k: int = 3) -> dict:
    """3-fold patient-level CV AUROC (for C1 vs C0 and other binary tasks)."""
    patients = sorted(df['sample_id'].unique())
    rng = np.random.default_rng(42)
    rng.shuffle(patients)
    folds = np.array_split(patients, k)

    fold_scores = []
    for i, test_patients in enumerate(folds):
        test_mask  = df['sample_id'].isin(test_patients)
        train_mask = ~test_mask

        X_train = df.loc[train_mask, FEAT_COLS].values
        y_train = df.loc[train_mask, label_col].values
        X_test  = df.loc[test_mask, FEAT_COLS].values
        y_test  = df.loc[test_mask, label_col].values

        if len(np.unique(y_test)) < 2 or len(np.unique(y_train)) < 2:
            continue

        imp = SimpleImputer(strategy='median')
        X_train = imp.fit_transform(X_train)
        X_test  = imp.transform(X_test)

        rf = RandomForestClassifier(**RF_PARAMS)
        rf.fit(X_train, y_train)
        pos_idx = list(rf.classes_).index(1)
        proba = rf.predict_proba(X_test)[:, pos_idx]
        fold_scores.append(roc_auc_score(y_test, proba))

    return {
        'mean': float(np.mean(fold_scores)) if fold_scores else np.nan,
        'folds': fold_scores,
    }


def kfold_patient_f1(df: pd.DataFrame, k: int = 3) -> dict:
    """3-fold patient-level CV 6-class macro F1."""
    le = LabelEncoder()
    df = df.copy()
    df['label_enc'] = le.fit_transform(df['cluster'])

    patients = sorted(df['sample_id'].unique())
    rng = np.random.default_rng(42)
    rng.shuffle(patients)
    folds = np.array_split(patients, k)

    fold_scores = []
    for test_patients in folds:
        test_mask  = df['sample_id'].isin(test_patients)
        train_mask = ~test_mask

        X_train = df.loc[train_mask, FEAT_COLS].values
        y_train = df.loc[train_mask, 'label_enc'].values
        X_test  = df.loc[test_mask, FEAT_COLS].values
        y_test  = df.loc[test_mask, 'label_enc'].values

        imp = SimpleImputer(strategy='median')
        X_train = imp.fit_transform(X_train)
        X_test  = imp.transform(X_test)

        rf = RandomForestClassifier(**RF_PARAMS)
        rf.fit(X_train, y_train)
        pred = rf.predict(X_test)
        fold_scores.append(f1_score(y_test, pred, average='macro'))

    return {
        'mean': float(np.mean(fold_scores)) if fold_scores else np.nan,
        'folds': fold_scores,
    }


def platform_kfold_auroc(df_c1: pd.DataFrame, k: int = 5) -> dict:
    """Metric 3: C1-Visium vs C1-SPA AUROC via patient-level k-fold CV.

    Label: SPA=1, Visium=0. Patients (sample_ids) are assigned to k folds;
    each fold trains on k-1 folds of patients (both SPA + Visium) and tests
    on the held-out patients. Reports mean AUROC ± [min, max].

    AUROC near 0.5 -> platform not visible from H&E features (good).
    AUROC near 1.0 -> platform is a confounder (caveat required).
    """
    df = df_c1.copy()
    df['label'] = (df['st_technology'] == 'Spatial Transcriptomics').astype(int)

    patients = sorted(df['sample_id'].unique())
    rng = np.random.default_rng(42)
    rng.shuffle(patients)
    folds = np.array_split(patients, k)

    fold_scores = []
    for test_patients in folds:
        test_mask  = df['sample_id'].isin(test_patients)
        train_mask = ~test_mask

        X_train = df.loc[train_mask, FEAT_COLS].values
        y_train = df.loc[train_mask, 'label'].values
        X_test  = df.loc[test_mask, FEAT_COLS].values
        y_test  = df.loc[test_mask, 'label'].values

        if len(np.unique(y_test)) < 2 or len(np.unique(y_train)) < 2:
            continue

        imp = SimpleImputer(strategy='median')
        X_train = imp.fit_transform(X_train)
        X_test  = imp.transform(X_test)

        rf = RandomForestClassifier(**RF_PARAMS)
        rf.fit(X_train, y_train)
        pos_idx = list(rf.classes_).index(1)
        proba = rf.predict_proba(X_test)[:, pos_idx]
        fold_scores.append(roc_auc_score(y_test, proba))

    if not fold_scores:
        return {'mean': np.nan, 'min': np.nan, 'max': np.nan, 'folds': []}
    return {
        'mean': float(np.mean(fold_scores)),
        'min':  float(np.min(fold_scores)),
        'max':  float(np.max(fold_scores)),
        'folds': fold_scores,
    }


# ── Load data ─────────────────────────────────────────────────────────────────
print("Loading features...")
df = pd.read_csv(OUT_DIR / 'phase1_features.tsv', sep='\t')
print(f"  {len(df)} instances, {df[FEAT_COLS].notna().all(axis=1).sum()} with all features complete")
print(f"  Feature NaN counts: {df[FEAT_COLS].isna().sum().to_dict()}")
print()

# Print per-cluster sample counts
print("Cluster composition:")
tech_counts = df.groupby(['cluster', 'st_technology'])['instance_id'].count().unstack(fill_value=0)
print(tech_counts)
print()

# ── Define subsets ─────────────────────────────────────────────────────────────
df_c4      = df[df['cluster'] == 4].copy()
df_c1      = df[df['cluster'] == 1].copy()
df_c1_spa  = df_c1[df_c1['st_technology'] == 'Spatial Transcriptomics'].copy()
df_c1_vis  = df_c1[df_c1['st_technology'] == 'Visium'].copy()
df_c0      = df[df['cluster'] == 0].copy()

c4_patients = sorted(df_c4['sample_id'].unique())
print(f"C4 patients ({len(c4_patients)}): {c4_patients}")
print(f"C1 patients: {df_c1['sample_id'].nunique()} total, "
      f"{df_c1_spa['sample_id'].nunique()} SPA, "
      f"{df_c1_vis['sample_id'].nunique()} Visium")
print()

# ── METRIC 1: C1-SPA vs C4 AUROC (PRIMARY) ───────────────────────────────────
print("=" * 65)
print("METRIC 1: C1-SPA vs C4 AUROC (PRIMARY, biology signal)")
print("=" * 65)
df_m1 = pd.concat([df_c1_spa.assign(label=0), df_c4.assign(label=1)])
print(f"  Instances: {len(df_m1)} (C1-SPA={len(df_c1_spa)}, C4={len(df_c4)})")
m1 = lopo_c4_auroc(df_m1, 'label', c4_patients)
print(f"  AUROC: {m1['mean']:.3f} (range [{m1['min']:.3f}, {m1['max']:.3f}])")
print(f"  Per-fold: {[f'{s:.3f}' for s in m1['folds']]}")
print()

# ── METRIC 2: C1 vs C4 full (secondary) ──────────────────────────────────────
print("=" * 65)
print("METRIC 2: C1 vs C4 AUROC (full, mixed platform)")
print("=" * 65)
df_m2 = pd.concat([df_c1.assign(label=0), df_c4.assign(label=1)])
print(f"  Instances: {len(df_m2)} (C1={len(df_c1)}, C4={len(df_c4)})")
m2 = lopo_c4_auroc(df_m2, 'label', c4_patients)
print(f"  AUROC: {m2['mean']:.3f} (range [{m2['min']:.3f}, {m2['max']:.3f}])")
print(f"  Per-fold: {[f'{s:.3f}' for s in m2['folds']]}")
print()

# ── METRIC 3: C1-Visium vs C1-SPA (platform control) ─────────────────────────
print("=" * 65)
print("METRIC 3: C1-Visium vs C1-SPA AUROC (platform confound control)")
print("=" * 65)
print(f"  Instances: {len(df_c1)} (C1-SPA={len(df_c1_spa)}, C1-Visium={len(df_c1_vis)})")
m3 = platform_kfold_auroc(df_c1)
print(f"  AUROC: {m3['mean']:.3f} (range [{m3['min']:.3f}, {m3['max']:.3f}])")
print(f"  Interpretation: near 0.5 = platform invisible; near 1.0 = platform confounds")
print()

# ── METRIC 4: C1 vs C0 AUROC (ST necessity check) ────────────────────────────
print("=" * 65)
print("METRIC 4: C1 vs C0 AUROC (ST necessity : expect ~0.5)")
print("=" * 65)
df_m4 = pd.concat([df_c0.assign(label=0), df_c1.assign(label=1)])
print(f"  Instances: {len(df_m4)} (C0={len(df_c0)}, C1={len(df_c1)})")
m4 = kfold_patient_auroc(df_m4, 'label', k=3)
print(f"  AUROC: {m4['mean']:.3f}")
print(f"  Per-fold: {[f'{s:.3f}' for s in m4['folds']]}")
print()

# ── METRIC 5: 6-class macro F1 (context) ─────────────────────────────────────
print("=" * 65)
print("METRIC 5: 6-class macro F1 (context)")
print("=" * 65)
print(f"  Instances: {len(df)}")
m5 = kfold_patient_f1(df, k=3)
print(f"  Macro F1: {m5['mean']:.3f}")
print(f"  Per-fold: {[f'{s:.3f}' for s in m5['folds']]}")
print()

# ── Feature importance (train on full C1 vs C4 SPA) ───────────────────────────
print("Computing feature importance (full data, no CV)...")
df_fi = pd.concat([df_c1_spa.assign(label=0), df_c4.assign(label=1)])
imp_all = SimpleImputer(strategy='median')
X_fi = imp_all.fit_transform(df_fi[FEAT_COLS].values)
y_fi = df_fi['label'].values
rf_fi = RandomForestClassifier(**RF_PARAMS)
rf_fi.fit(X_fi, y_fi)
importance_df = pd.DataFrame({'feature': FEAT_COLS,
                               'importance': rf_fi.feature_importances_})
importance_df = importance_df.sort_values('importance', ascending=False)
print(importance_df.to_string(index=False))
importance_df.to_csv(OUT_DIR / 'phase1_feature_importance.tsv', sep='\t', index=False)
print()

# ── Summary and decision gate ─────────────────────────────────────────────────
print("=" * 65)
print("SUMMARY")
print("=" * 65)
print(f"  Primary   C1-SPA vs C4 AUROC:       {m1['mean']:.3f}  [{m1['min']:.3f}–{m1['max']:.3f}]")
print(f"  Secondary C1 vs C4 full AUROC:       {m2['mean']:.3f}  [{m2['min']:.3f}–{m2['max']:.3f}]")
print(f"  Platform  C1-Visium vs C1-SPA AUROC: {m3['mean']:.3f}  [{m3['min']:.3f}–{m3['max']:.3f}]")
print(f"  ST-nec.   C1 vs C0 AUROC:            {m4['mean']:.3f}")
print(f"  Context   6-class macro F1:          {m5['mean']:.3f}")
print()
primary = m1['mean']
if primary > 0.70:
    verdict = "PROCEED to Phase 2 GNN (primary AUROC > 0.70)"
elif primary >= 0.60:
    verdict = "CAUTIOUS PROCEED to Phase 2 (0.60 ≤ AUROC ≤ 0.70)"
else:
    verdict = "STOP : negative result; H&E morphology cannot predict C4 without ST labels"
print(f"  Decision: {verdict}")
