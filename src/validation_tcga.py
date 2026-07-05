"""
TCGA-BRCA TLS cluster survival validation (Claim 1b).

Input files (all in outputs/validation_tcga/):
  tcga_RSEM_gene_tpm.gz      — UCSC Xena TOIL RSEM log2(TPM+0.001), pan-cancer
  BRCA_expression_modules.tsv — pre-extracted BRCA samples × 33 genes (z-scored)
  BRCA_clinical.tsv          — TCGA-BRCA clinical (from Xena)
  BRCA_survival.txt          — TCGA-BRCA OS/DFS/PFI (from Xena)

If BRCA_expression_modules.tsv is absent, the script re-extracts it from the raw
TOIL matrix using data/ensembl_to_symbol.tsv.
"""
import gzip, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

warnings.filterwarnings('ignore')

DATA = "/mnt/e/hest/outputs/validation_tcga"
ENSEMBL_MAP = "/mnt/e/hest/data/ensembl_to_symbol.tsv"
OUT = DATA

# ── Gene module definitions ────────────────────────────────────────────────
C1_GENES = ["CXCL13","CCL19","CCL21","LAMP3","ACKR1","MS4A1","CXCR5",
            "BCL6","IGHG1","JCHAIN","MZB1","TNFSF13B","AICDA"]
C4_GENES = ["TGFB1","IL10","ARG1","MRC1","IDO1","FOXP3","IL2RA","CTLA4","CD274"]
ALL_GENES = sorted(set(C1_GENES + C4_GENES +
    ["ICOS","IL21","SELP","SELL","CXCR4","CXCL12","CR2","IGHA1","TOX2","MYBL1","ELL3"]))

# Manual Ensembl IDs for genes missing from symbol mapping (IG loci / pseudogenes)
MANUAL_ENSEMBL = {
    "IGHG1":  "ENSG00000211896",
    "IGHA1":  "ENSG00000211895",
    "JCHAIN": "ENSG00000132465",
}

# ── 1. Load or extract expression ─────────────────────────────────────────
expr_path = f"{DATA}/BRCA_expression_modules.tsv"
try:
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    print(f"Loaded pre-extracted expression: {expr.shape}")
except FileNotFoundError:
    print("BRCA_expression_modules.tsv not found — extracting from raw TOIL matrix...")

    # Load Ensembl→symbol map
    emap = pd.read_csv(ENSEMBL_MAP, sep="\t")
    sym2ens = dict(zip(emap["gene_name"], emap["gene_id"]))
    sym2ens.update(MANUAL_ENSEMBL)

    target_ens = {sym2ens[g]: g for g in ALL_GENES if g in sym2ens}
    print(f"  Ensembl IDs to extract: {len(target_ens)}/{len(ALL_GENES)}")

    # Stream the 707MB matrix
    with gzip.open(f"{DATA}/tcga_RSEM_gene_tpm.gz", "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        # header[0] = "sample" (gene ID column); rest = TCGA sample IDs
        brca_idx = [i for i in range(1, len(header)) if header[i].startswith("TCGA-")]
        brca_ids = [header[i] for i in brca_idx]
        found, rows = {}, []
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            ens_id = parts[0].split(".")[0]  # strip version suffix
            if ens_id in target_ens:
                vals = [float(parts[i]) for i in brca_idx]
                found[target_ens[ens_id]] = vals

    print(f"  Genes found in matrix: {sorted(found.keys())}")
    raw = pd.DataFrame(found, index=brca_ids)

    # z-score across samples for each gene
    expr = (raw - raw.mean()) / raw.std()
    expr.index.name = "sample"
    expr.to_csv(expr_path, sep="\t")
    print(f"  Saved {expr_path}")

# ── 2. Load clinical + survival ────────────────────────────────────────────
clin = pd.read_csv(f"{DATA}/BRCA_clinical.tsv", sep="\t", index_col=0, low_memory=False)
surv = pd.read_csv(f"{DATA}/BRCA_survival.txt", sep="\t", index_col=0)

# Use PAM50 from RNA-seq (most reliable)
pam50 = clin[["PAM50Call_RNAseq","pathologic_stage",
              "age_at_initial_pathologic_diagnosis"]].copy()
pam50.columns = ["PAM50","stage","age"]
pam50["age"] = pd.to_numeric(pam50["age"], errors="coerce")

# Recode stage to numeric (I=1…IV=4)
stage_map = {"Stage I":1,"Stage IA":1,"Stage IB":1,
             "Stage II":2,"Stage IIA":2,"Stage IIB":2,
             "Stage III":3,"Stage IIIA":3,"Stage IIIB":3,"Stage IIIC":3,
             "Stage IV":4,"Stage X":np.nan}
pam50["stage_n"] = pam50["stage"].map(stage_map)

# PFI (Progression-Free Interval) — preferred endpoint for BRCA immunological validation.
# More events than OS (173 vs 198) and captures immune-relevant recurrence.
# OS has only ~16% event rate in TCGA-BRCA (good-prognosis cohort; confounded by
# competing causes of death); DSS has only 107 events.
os_df = surv[["PFI","PFI.time"]].copy()
os_df.columns = ["EVENT","MONTHS"]
os_df["MONTHS"] = os_df["MONTHS"] / 30.44  # days → months

# Module scores
c1_avail = [g for g in C1_GENES if g in expr.columns]
c4_avail = [g for g in C4_GENES if g in expr.columns]
print(f"C1 genes: {len(c1_avail)}/{len(C1_GENES)}  C4 genes: {len(c4_avail)}/{len(C4_GENES)}")

expr["C1_score"] = expr[c1_avail].mean(axis=1)
expr["C4_score"] = expr[c4_avail].mean(axis=1)

# Join
df = expr[["C1_score","C4_score"]].join(os_df).join(pam50)
df = df.dropna(subset=["EVENT","MONTHS"])
print(f"Joined: {df.shape}, events={int(df['EVENT'].sum())}")
print(f"PAM50: {df['PAM50'].value_counts().to_dict()}")

# ── 3. KM helper ─────────────────────────────────────────────────────────────
def km_q(ax, df, score_col, title, q_hi=0.75, q_lo=0.25):
    d = df[[score_col,"MONTHS","EVENT"]].dropna()
    hi = d[d[score_col] >= d[score_col].quantile(q_hi)]
    lo = d[d[score_col] <= d[score_col].quantile(q_lo)]
    lr = logrank_test(hi["MONTHS"], lo["MONTHS"], hi["EVENT"], lo["EVENT"])
    p = lr.p_value
    kmf = KaplanMeierFitter()
    for grp, color, lbl in [
        (hi, "#d62728", f"High (n={len(hi)})"),
        (lo, "#1f77b4", f"Low (n={len(lo)})")
    ]:
        kmf.fit(grp["MONTHS"], grp["EVENT"], label=lbl)
        kmf.plot_survival_function(ax=ax, ci_show=False, color=color)
    p_str = f"p={p:.2e}" if p >= 1e-4 else "p<1e-4"
    ax.set_title(f"{title}\n{p_str}", fontsize=10)
    ax.set_xlabel("Months")
    ax.set_ylabel("PFI")
    ax.legend(fontsize=7)
    return p

# ── 4. Main KM figure (6 panels) ────────────────────────────────────────────
fig, axes = plt.subplots(2, 3, figsize=(15, 9))
axes = axes.flatten()

panels = [
    (df,                            "C1_score", "C1 — All TCGA-BRCA"),
    (df,                            "C4_score", "C4 — All TCGA-BRCA"),
    (df[df["PAM50"]=="LumB"],       "C1_score", "C1 — LumB"),
    (df[df["PAM50"]=="LumA"],       "C1_score", "C1 — LumA"),
    (df[df["PAM50"]=="Basal"],      "C1_score", "C1 — Basal"),
    (df[df["PAM50"]=="LumB"],       "C4_score", "C4 — LumB"),
]

p_values = {}
for ax, (sub, score, title) in zip(axes, panels):
    p = km_q(ax, sub, score, title)
    p_values[title] = p

plt.suptitle("TCGA-BRCA TLS Cluster Survival (PFI endpoint)", fontsize=13, y=1.02)
plt.tight_layout()
plt.savefig(f"{OUT}/TCGA_BRCA_TLS_survival.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved TCGA_BRCA_TLS_survival.png")

# ── 5. C4 PAM50-stratified KM (all subtypes) ──────────────────────────────
subtypes = [s for s in ["LumA","LumB","Her2","Basal","Normal"]
            if s in df["PAM50"].values]
fig, axes = plt.subplots(1, len(subtypes), figsize=(4*len(subtypes), 5))
if len(subtypes) == 1:
    axes = [axes]
for ax, sub in zip(axes, subtypes):
    km_q(ax, df[df["PAM50"]==sub], "C4_score", f"C4 — {sub}")
plt.suptitle("C4 score by PAM50 subtype (TCGA-BRCA)", fontsize=12, y=1.03)
plt.tight_layout()
plt.savefig(f"{OUT}/TCGA_BRCA_C4_stratified.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved TCGA_BRCA_C4_stratified.png")

# ── 6. Cox multivariate ────────────────────────────────────────────────────
print("\nRunning Cox multivariate...")
cox_df = df[["MONTHS","EVENT","C1_score","C4_score","age","stage_n","PAM50"]].copy()
cox_df = cox_df.dropna()
print(f"Cox N={len(cox_df)}, events={int(cox_df['EVENT'].sum())} (PFI endpoint)")

pam_dummies = pd.get_dummies(cox_df["PAM50"], prefix="PAM50", drop_first=False)
ref = "PAM50_LumA"
pam_cols = [c for c in pam_dummies.columns if c != ref]
cox_df = pd.concat([cox_df, pam_dummies[pam_cols]], axis=1)

def run_cox(df, score_col):
    cph = CoxPHFitter(penalizer=0.1)
    cols = ["MONTHS","EVENT", score_col, "age","stage_n"] + pam_cols
    sub = df[cols].dropna()
    cph.fit(sub, duration_col="MONTHS", event_col="EVENT")
    row = cph.summary.loc[score_col]
    return np.exp(row["coef"]), row["p"], row["exp(coef) lower 95%"], row["exp(coef) upper 95%"]

hr1, p1, lo1, hi1 = run_cox(cox_df, "C1_score")
hr4, p4, lo4, hi4 = run_cox(cox_df, "C4_score")

# ── 7. Summary ──────────────────────────────────────────────────────────────
print("\n=== TCGA-BRCA SURVIVAL SUMMARY ===")
print(f"N={len(df)}, events={int(df['EVENT'].sum())}")
print("\nKM log-rank p-values:")
for k, v in p_values.items():
    sig = "***" if v < 0.001 else ("**" if v < 0.01 else ("*" if v < 0.05 else "ns"))
    print(f"  {k:<30s}  p={v:.5f}  {sig}")
print(f"\nCox multivariate (adj age+stage+PAM50, ref=LumA):")
print(f"  C1  HR={hr1:.3f} [{lo1:.3f}–{hi1:.3f}] p={p1:.4f}")
print(f"  C4  HR={hr4:.3f} [{lo4:.3f}–{hi4:.3f}] p={p4:.4f}")
