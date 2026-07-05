"""
TCGA-CESC TLS cluster survival validation.

Tests whether C1 (immunogenic TLS) gene module scores predict OS in cervical
squamous cell carcinoma (TCGA, Firehose Legacy; n=304 RNA-seq samples).
CESC is represented in the HEST-1k cohort (n=12 ST instances).

Fetches expression z-scores and clinical data from cBioPortal.
Cox PH adjusted for age (FIGO staging not available in this cohort).

Output:
    outputs/validation_tcga/cesc_results.txt
    outputs/validation_tcga/CESC_TLS_survival.png
"""
import requests, time, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

warnings.filterwarnings('ignore')

API   = "https://www.cbioportal.org/api"
STUDY = "cesc_tcga"
PROFILE = "cesc_tcga_rna_seq_v2_mrna_median_Zscores"
OUT   = "/mnt/e/hest/outputs/validation_tcga"

C1_GENES = ["CXCL13","CCL19","CCL21","LAMP3","ACKR1","MS4A1","CXCR5",
            "BCL6","IGHG1","JCHAIN","MZB1","TNFSF13B","AICDA"]
C4_GENES = ["TGFB1","IL10","ARG1","MRC1","IDO1","FOXP3","IL2RA","CTLA4","CD274"]


def fetch_genes(genes, profile, study):
    r = requests.post(f"{API}/genes/fetch",
                      params={"geneIdType": "HUGO_GENE_SYMBOL"},
                      json=genes, timeout=30)
    r.raise_for_status()
    sym_to_entrez = {g["hugoGeneSymbol"]: g["entrezGeneId"] for g in r.json()}
    entrez_to_sym = {v: k for k, v in sym_to_entrez.items()}
    print(f"  Resolved {len(sym_to_entrez)}/{len(genes)} genes")

    r2 = requests.post(
        f"{API}/molecular-profiles/{profile}/molecular-data/fetch",
        params={"projection": "SUMMARY"},
        json={"entrezGeneIds": list(sym_to_entrez.values()),
              "sampleListId": f"{study}_all"},
        timeout=120)
    r2.raise_for_status()
    rows = [{"patientId": d["patientId"],
             "gene": entrez_to_sym[d["entrezGeneId"]],
             "value": d["value"]} for d in r2.json()]
    df = pd.DataFrame(rows).groupby(["patientId","gene"])["value"].first().unstack()
    df.columns.name = None
    return df


def fetch_clinical(attr, study=STUDY):
    r = requests.get(f"{API}/studies/{study}/clinical-data",
                     params={"attributeId": attr, "clinicalDataType": "PATIENT",
                             "pageSize": 100000}, timeout=60)
    r.raise_for_status()
    s = pd.Series({d["patientId"]: d["value"] for d in r.json()}, name=attr)
    print(f"  {attr}: {len(s)} records")
    return s


# ── 1. Expression ──────────────────────────────────────────────────────────
print("Fetching expression...")
expr = fetch_genes(sorted(set(C1_GENES + C4_GENES)), PROFILE, STUDY)

c1_avail = [g for g in C1_GENES if g in expr.columns]
c4_avail = [g for g in C4_GENES if g in expr.columns]
print(f"C1 genes: {len(c1_avail)}/{len(C1_GENES)}  missing: {set(C1_GENES)-set(c1_avail)}")
print(f"C4 genes: {len(c4_avail)}/{len(C4_GENES)}")

expr["C1_score"] = expr[c1_avail].mean(axis=1)
expr["C4_score"] = expr[c4_avail].mean(axis=1)

# ── 2. Clinical ────────────────────────────────────────────────────────────
print("\nFetching clinical data...")
df = expr[["C1_score","C4_score"]].join(
    [fetch_clinical("OS_MONTHS"), fetch_clinical("OS_STATUS"), fetch_clinical("AGE")])

df["EVENT"] = df["OS_STATUS"].str.startswith("1").astype(float)
df["MONTHS"] = pd.to_numeric(df["OS_MONTHS"], errors="coerce")
df["AGE_N"]  = pd.to_numeric(df["AGE"], errors="coerce")
df = df.dropna(subset=["EVENT","MONTHS","C1_score","AGE_N"])
print(f"\nAnalysis set: n={len(df)}, events={int(df['EVENT'].sum())}")

# ── 3. KM plot ─────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

def km_panel(ax, score_col, title, color_hi="#029E73", color_lo="#aaaaaa"):
    d = df[[score_col,"MONTHS","EVENT"]].dropna()
    hi = d[d[score_col] >= d[score_col].quantile(0.75)]
    lo = d[d[score_col] <= d[score_col].quantile(0.25)]
    lr = logrank_test(hi["MONTHS"], lo["MONTHS"], hi["EVENT"], lo["EVENT"])
    p = lr.p_value
    kmf = KaplanMeierFitter()
    for grp, col, lbl in [(hi, color_hi, f"High (n={len(hi)})"),
                          (lo, color_lo, f"Low (n={len(lo)})")]:
        kmf.fit(grp["MONTHS"], grp["EVENT"], label=lbl)
        kmf.plot_survival_function(ax=ax, ci_show=True, color=col)
    p_str = f"p={p:.3f}" if p >= 0.001 else f"p={p:.2e}"
    ax.set_title(f"{title}\n(log-rank {p_str})", fontsize=10)
    ax.set_xlabel("Months")
    ax.set_ylabel("Overall survival")
    ax.legend(fontsize=8)
    ax.spines[["top","right"]].set_visible(False)

km_panel(axes[0], "C1_score", "C1 (Immunogenic TLS) — TCGA-CESC")
km_panel(axes[1], "C4_score", "C4 (Myeloid-suppressive TLS) — TCGA-CESC",
         color_hi="#d62728")

plt.suptitle(f"TCGA Cervical SCC (n={len(df)}, OS endpoint)", fontsize=10, y=1.01)
plt.tight_layout()
plt.savefig(f"{OUT}/CESC_TLS_survival.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved CESC_TLS_survival.png")

# ── 4. Multivariate Cox ────────────────────────────────────────────────────
results = {}
for score_col, label in [("C1_score","C1"), ("C4_score","C4")]:
    cox_df = df[["MONTHS","EVENT",score_col,"AGE_N"]].dropna()
    cox_df = cox_df.rename(columns={score_col:"score"})
    cph = CoxPHFitter()
    cph.fit(cox_df, duration_col="MONTHS", event_col="EVENT",
            formula="score + AGE_N")
    row = cph.summary.loc["score"]
    hr = np.exp(row["coef"])
    lo = np.exp(row["coef lower 95%"])
    hi = np.exp(row["coef upper 95%"])
    p  = row["p"]
    results[label] = (hr, lo, hi, p, len(cox_df))
    print(f"{label}: HR={hr:.3f} [{lo:.3f}–{hi:.3f}]  p={p:.4f}  n={len(cox_df)}")

# ── 5. Save results ────────────────────────────────────────────────────────
lines = [
    "TCGA-CESC TLS survival validation",
    f"Study: {STUDY}  Profile: {PROFILE}",
    f"n={len(df)}  events={int(df['EVENT'].sum())}  endpoint=OS  adj=age",
    "",
]
for label, (hr, lo, hi, p, n) in results.items():
    sig = "*" if p < 0.05 else "ns"
    lines.append(f"{label}:  HR={hr:.3f} [{lo:.3f}–{hi:.3f}]  p={p:.4f}  n={n}  {sig}")

with open(f"{OUT}/cesc_results.txt", "w") as f:
    f.write("\n".join(lines) + "\n")

print("\n".join(lines))
print(f"\nDone — results in {OUT}/cesc_results.txt")
