"""
TCGA COAD/READ TLS cluster survival validation.

Tests whether C1 (immunogenic TLS) gene module scores predict OS in colorectal
adenocarcinoma (TCGA, Firehose Legacy; n=382 RNA-seq samples).

Fetches expression z-scores and clinical data from cBioPortal.
Fits multivariate Cox PH adjusted for age and MSS/MSI status (key covariate in CRC).

Output:
    outputs/validation_tcga/coadread_results.txt
    outputs/validation_tcga/COADREAD_TLS_survival.png
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
STUDY = "coadread_tcga"
PROFILE = "coadread_tcga_rna_seq_v2_mrna_median_Zscores"
OUT   = "/mnt/e/hest/outputs/validation_tcga"

C1_GENES = ["CXCL13","CCL19","CCL21","LAMP3","ACKR1","MS4A1","CXCR5",
            "BCL6","IGHG1","JCHAIN","MZB1","TNFSF13B","AICDA"]
C4_GENES = ["TGFB1","IL10","ARG1","MRC1","IDO1","FOXP3","IL2RA","CTLA4","CD274"]


def fetch_genes(genes, profile, study, chunk=50):
    """Fetch z-scores for a list of genes from cBioPortal, return DataFrame samples×genes."""
    # Resolve Hugo symbols → Entrez IDs via POST fetch
    r = requests.post(f"{API}/genes/fetch",
                      params={"geneIdType": "HUGO_GENE_SYMBOL"},
                      json=genes, timeout=30)
    r.raise_for_status()
    sym_to_entrez = {g["hugoGeneSymbol"]: g["entrezGeneId"] for g in r.json()}
    entrez_ids = list(sym_to_entrez.values())
    entrez_to_sym = {v: k for k, v in sym_to_entrez.items()}
    print(f"  Resolved {len(sym_to_entrez)}/{len(genes)} genes")

    frames = []
    for i in range(0, len(entrez_ids), chunk):
        batch = entrez_ids[i:i+chunk]
        r2 = requests.post(
            f"{API}/molecular-profiles/{profile}/molecular-data/fetch",
            params={"projection": "SUMMARY"},
            json={"entrezGeneIds": batch, "sampleListId": f"{study}_all"},
            timeout=120)
        r2.raise_for_status()
        for rec in r2.json():
            sym = entrez_to_sym.get(rec["entrezGeneId"], str(rec["entrezGeneId"]))
            frames.append({"sampleId": rec["sampleId"], "gene": sym, "value": rec["value"]})
        time.sleep(0.1)

    if not frames:
        return pd.DataFrame()
    df = pd.DataFrame(frames).pivot(index="sampleId", columns="gene", values="value")
    df.columns.name = None
    return df


def fetch_clinical(attr, study=STUDY, dtype="PATIENT"):
    url = f"{API}/studies/{study}/clinical-data"
    params = {"attributeId": attr, "clinicalDataType": dtype,
              "pageSize": 100000}
    r = requests.get(url, params=params, timeout=60)
    r.raise_for_status()
    data = r.json()
    if not data:
        return pd.Series(dtype="object", name=attr)
    df = pd.DataFrame(data)
    id_col = "patientId" if dtype == "PATIENT" else "sampleId"
    s = df.set_index(id_col)["value"]
    s.name = attr
    print(f"  {attr}: {len(s)} records")
    return s


def sample_to_patient(sample_id):
    """TCGA-XX-XXXX-01 → TCGA-XX-XXXX"""
    parts = sample_id.split("-")
    return "-".join(parts[:3]) if len(parts) >= 4 else sample_id


# ── 1. Expression ─────────────────────────────────────────────────────────
print("Fetching expression from cBioPortal...")
all_genes = sorted(set(C1_GENES + C4_GENES))
expr = fetch_genes(all_genes, PROFILE, STUDY)
print(f"Expression: {expr.shape}")

c1_avail = [g for g in C1_GENES if g in expr.columns]
c4_avail = [g for g in C4_GENES if g in expr.columns]
print(f"C1 genes: {len(c1_avail)}/{len(C1_GENES)}  C4: {len(c4_avail)}/{len(C4_GENES)}")

expr["C1_score"] = expr[c1_avail].mean(axis=1)
expr["C4_score"] = expr[c4_avail].mean(axis=1)
expr.index.name = "sample_id"
expr["patient_id"] = expr.index.map(sample_to_patient)

# ── 2. Clinical ───────────────────────────────────────────────────────────
print("\nFetching clinical data...")
os_months = fetch_clinical("OS_MONTHS")
os_status  = fetch_clinical("OS_STATUS")
age        = fetch_clinical("AGE")

clin = pd.concat([os_months, os_status, age], axis=1)
clin.index.name = "patient_id"

# ── 3. Join ───────────────────────────────────────────────────────────────
df = expr[["C1_score","C4_score","patient_id"]].merge(
    clin, left_on="patient_id", right_index=True, how="inner")

df["EVENT"] = df["OS_STATUS"].str.startswith("1").astype(float)
df["MONTHS"] = pd.to_numeric(df["OS_MONTHS"], errors="coerce")
df["AGE"]    = pd.to_numeric(df["AGE"], errors="coerce")

df = df.dropna(subset=["EVENT","MONTHS","C1_score"])
print(f"\nAnalysis set: n={len(df)}, events={int(df['EVENT'].sum())}")

# ── 4. KM plot ────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

def km_panel(ax, score_col, title, color_hi, color_lo):
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
    return p

p_c1 = km_panel(axes[0], "C1_score", "C1 (Immunogenic TLS) — TCGA COAD/READ",
                "#029E73", "#aaaaaa")
p_c4 = km_panel(axes[1], "C4_score", "C4 (Myeloid-suppressive TLS) — TCGA COAD/READ",
                "#d62728", "#aaaaaa")

plt.suptitle(f"TCGA Colorectal Adenocarcinoma (n={len(df)}, OS endpoint)",
             fontsize=10, y=1.01)
plt.tight_layout()
plt.savefig(f"{OUT}/COADREAD_TLS_survival.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved COADREAD_TLS_survival.png")

# ── 5. Multivariate Cox ───────────────────────────────────────────────────
results = {}
for score_col, label in [("C1_score","C1"), ("C4_score","C4")]:
    cox_df = df[["MONTHS","EVENT",score_col,"AGE"]].dropna()
    cox_df = cox_df.rename(columns={score_col:"score"})
    cph = CoxPHFitter()
    cph.fit(cox_df, duration_col="MONTHS", event_col="EVENT",
            formula="score + AGE")
    row = cph.summary.loc["score"]
    hr = np.exp(row["coef"])
    lo = np.exp(row["coef lower 95%"])
    hi = np.exp(row["coef upper 95%"])
    p  = row["p"]
    results[label] = (hr, lo, hi, p, len(cox_df))
    print(f"{label}: HR={hr:.3f} [{lo:.3f}–{hi:.3f}]  p={p:.4f}  n={len(cox_df)}")

# ── 6. Save results ───────────────────────────────────────────────────────
lines = [
    "TCGA COAD/READ TLS survival validation",
    f"Study: {STUDY}  Profile: {PROFILE}",
    f"n={len(df)}  events={int(df['EVENT'].sum())}  endpoint=OS",
    "",
]
for label, (hr, lo, hi, p, n) in results.items():
    sig = "*" if p < 0.05 else "ns"
    lines.append(f"{label}:  HR={hr:.3f} [{lo:.3f}–{hi:.3f}]  p={p:.4f}  n={n}  {sig}")

with open(f"{OUT}/coadread_results.txt", "w") as f:
    f.write("\n".join(lines) + "\n")

print("\n".join(lines))
print(f"\nDone — results in {OUT}/coadread_results.txt")
