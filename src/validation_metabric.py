"""METABRIC TLS cluster survival validation (Claim 1c)."""
import requests, time, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

warnings.filterwarnings('ignore')

API = "https://www.cbioportal.org/api"
STUDY = "brca_metabric"
OUT = "/mnt/e/hest/outputs/validation_tcga"

# ── Gene module definitions (same split as TCGA analysis) ──────────────────
C1_GENES = ["CXCL13","CCL19","CCL21","LAMP3","ACKR1","MS4A1","CXCR5",
            "BCL6","IGHG1","JCHAIN","MZB1","TNFSF13B","AICDA"]
C4_GENES = ["TGFB1","IL10","ARG1","MRC1","IDO1","FOXP3","IL2RA","CTLA4","CD274"]

# ── 1. Load expression ──────────────────────────────────────────────────────
expr = pd.read_csv(f"{OUT}/METABRIC_expr_zscores.tsv", sep="\t", index_col=0)
print(f"Expression: {expr.shape}")

c1_avail = [g for g in C1_GENES if g in expr.columns]
c4_avail = [g for g in C4_GENES if g in expr.columns]
print(f"C1 genes available: {len(c1_avail)}/{len(C1_GENES)}: {c1_avail}")
print(f"C4 genes available: {len(c4_avail)}/{len(C4_GENES)}: {c4_avail}")

expr["C1_score"] = expr[c1_avail].mean(axis=1)
expr["C4_score"] = expr[c4_avail].mean(axis=1)

# ── 2. Fetch clinical data fresh from cBioPortal ───────────────────────────
def fetch_clinical(attr, dtype="PATIENT"):
    url = f"{API}/studies/{STUDY}/clinical-data"
    params = {"attributeId": attr, "clinicalDataType": dtype,
              "pageSize": 100000, "pageNumber": 0}
    r = requests.get(url, params=params, timeout=60)
    r.raise_for_status()
    data = r.json()
    if not data:
        print(f"  WARNING: {attr} returned 0 rows")
        return pd.Series(dtype="object", name=attr)
    df = pd.DataFrame(data)
    s = df.set_index("patientId")["value"]
    s.name = attr
    print(f"  {attr}: {len(s)} patients")
    return s

print("\nFetching clinical attributes...")
os_months = fetch_clinical("OS_MONTHS")
os_status = fetch_clinical("OS_STATUS")
age       = fetch_clinical("AGE_AT_DIAGNOSIS")
pam50     = fetch_clinical("CLAUDIN_SUBTYPE")
npi       = fetch_clinical("NPI")

clin = pd.concat([os_months, os_status, age, pam50, npi], axis=1)
print(f"Clinical joined: {clin.shape}")

# Parse OS
clin["OS_MONTHS"] = pd.to_numeric(clin["OS_MONTHS"], errors="coerce")
clin["EVENT"] = (clin["OS_STATUS"] == "1:DECEASED").astype(int)
clin["AGE"] = pd.to_numeric(clin["AGE_AT_DIAGNOSIS"], errors="coerce")
clin["NPI"] = pd.to_numeric(clin["NPI"], errors="coerce")

# Remap METABRIC sample IDs: clinical uses patientId (MB-XXXX); expr index may differ
# cBioPortal METABRIC patient IDs match expr index directly
df = expr.join(clin, how="inner")
print(f"After join: {df.shape}, events={df['EVENT'].sum()}")
print(f"PAM50: {df['CLAUDIN_SUBTYPE'].value_counts().to_dict()}")

# ── 3. KM helper ─────────────────────────────────────────────────────────────
def km_q(ax, df, score_col, title, q_hi=0.75, q_lo=0.25):
    d = df[[score_col, "OS_MONTHS", "EVENT"]].dropna()
    hi = d[d[score_col] >= d[score_col].quantile(q_hi)]
    lo = d[d[score_col] <= d[score_col].quantile(q_lo)]
    lr = logrank_test(hi["OS_MONTHS"], lo["OS_MONTHS"],
                      hi["EVENT"], lo["EVENT"])
    p = lr.p_value
    kmf = KaplanMeierFitter()
    for grp, color, lbl in [
        (hi, "#d62728", f"High (n={len(hi)})"),
        (lo, "#1f77b4", f"Low (n={len(lo)})")
    ]:
        kmf.fit(grp["OS_MONTHS"], grp["EVENT"], label=lbl)
        kmf.plot_survival_function(ax=ax, ci_show=False, color=color)
    p_str = f"p={p:.2e}" if p >= 1e-4 else f"p<1e-4"
    ax.set_title(f"{title}\n{p_str}", fontsize=10)
    ax.set_xlabel("Months")
    ax.set_ylabel("OS")
    ax.legend(fontsize=7)
    return p

# ── 4. 6-panel KM figure ────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 3, figsize=(15, 9))
axes = axes.flatten()

panels = [
    (df,                                    "C1_score", "C1 — All METABRIC"),
    (df,                                    "C4_score", "C4 — All METABRIC"),
    (df[df["CLAUDIN_SUBTYPE"]=="LumB"],     "C1_score", "C1 — LumB"),
    (df[df["CLAUDIN_SUBTYPE"]=="LumB"],     "C4_score", "C4 — LumB"),
    (df[df["CLAUDIN_SUBTYPE"]=="LumA"],     "C4_score", "C4 — LumA"),
    (df[df["CLAUDIN_SUBTYPE"]=="Basal"],    "C4_score", "C4 — Basal"),
]

p_values = {}
for ax, (sub, score, title) in zip(axes, panels):
    p = km_q(ax, sub, score, title)
    p_values[title] = p
    print(f"  {title}: p={p:.5f}")

plt.suptitle("METABRIC TLS Cluster Survival (n=1,979)", fontsize=13, y=1.02)
plt.tight_layout()
plt.savefig(f"{OUT}/METABRIC_TLS_survival.png", dpi=150, bbox_inches="tight")
plt.close()
print("\nSaved METABRIC_TLS_survival.png")

# ── 5. Cox multivariate (age + NPI + PAM50) ─────────────────────────────────
print("\nRunning Cox multivariate...")

pam50_order = ["LumA","LumB","Her2","Basal","Normal","claudin-low"]
cox_df = df[["OS_MONTHS","EVENT","C1_score","C4_score","AGE","NPI","CLAUDIN_SUBTYPE"]].copy()
cox_df = cox_df.dropna()
print(f"Cox N={len(cox_df)}, events={cox_df['EVENT'].sum()}")

# PAM50 dummies, drop LumA as reference
pam_dummies = pd.get_dummies(cox_df["CLAUDIN_SUBTYPE"], prefix="PAM50", drop_first=False)
ref = "PAM50_LumA"
pam_cols = [c for c in pam_dummies.columns if c != ref]
cox_df = pd.concat([cox_df, pam_dummies[pam_cols]], axis=1)

def run_cox(df, score_col):
    cph = CoxPHFitter(penalizer=0.1)
    cols = ["OS_MONTHS","EVENT", score_col, "AGE","NPI"] + pam_cols
    sub = df[cols].dropna()
    cph.fit(sub, duration_col="OS_MONTHS", event_col="EVENT")
    row = cph.summary.loc[score_col]
    return np.exp(row["coef"]), row["p"], row["exp(coef) lower 95%"], row["exp(coef) upper 95%"]

hr1, p1, lo1, hi1 = run_cox(cox_df, "C1_score")
hr4, p4, lo4, hi4 = run_cox(cox_df, "C4_score")

print(f"\nCox results (penalizer=0.1, ref=LumA):")
print(f"  C1_score: HR={hr1:.3f} [{lo1:.3f}–{hi1:.3f}] p={p1:.4f}")
print(f"  C4_score: HR={hr4:.3f} [{lo4:.3f}–{hi4:.3f}] p={p4:.4f}")

# ── 6. Summary ──────────────────────────────────────────────────────────────
print("\n=== METABRIC SURVIVAL SUMMARY ===")
print(f"N={len(df)}, events={int(df['EVENT'].sum())}")
print("\nKM log-rank p-values:")
for k, v in p_values.items():
    sig = "***" if v < 0.001 else ("**" if v < 0.01 else ("*" if v < 0.05 else "ns"))
    print(f"  {k:<30s}  p={v:.5f}  {sig}")
print(f"\nCox multivariate (adj age+NPI+PAM50):")
print(f"  C1  HR={hr1:.3f} [{lo1:.3f}–{hi1:.3f}] p={p1:.4f}")
print(f"  C4  HR={hr4:.3f} [{lo4:.3f}–{hi4:.3f}] p={p4:.4f}")
