from paths import DATA_DIR, RESULTS_DIR, FIG_DIR, PLOTS_DIR

import argparse
import gc
import os

import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr
from scipy import sparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "--allow-heuristic-doublet-fallback",
    action="store_true",
    help="Allow a labeled gene-count heuristic if Scrublet fails or cannot set a threshold.",
)
args = parser.parse_args()

IN_FP = os.path.join(RESULTS_DIR, "GSE178318_CRC_LM_raw.h5ad")
OUT_FP = os.path.join(RESULTS_DIR, "GSE178318_CRC_LM_qc.h5ad")

print(f"[*] Loading raw CRC+LM object: {IN_FP}")
adata = sc.read_h5ad(IN_FP)
if "patient" not in adata.obs.columns:
    raise KeyError("The AnnData object must contain a patient column for per-patient doublet detection.")
if adata.obs["patient"].isna().any():
    raise ValueError("Missing patient IDs; refusing to pool patients for doublet detection.")

if "gene_symbol" in adata.var.columns:
    mito_mask = adata.var["gene_symbol"].astype(str).str.upper().str.startswith("MT-")
else:
    mito_mask = np.array([False] * adata.n_vars)
adata.var["mt"] = mito_mask

print("[*] Calculating basic QC metrics...")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

valid_cells = (
    (adata.obs["n_genes_by_counts"] >= 200)
    & (adata.obs["n_genes_by_counts"] <= 6000)
    & (adata.obs["pct_counts_mt"] <= 20.0)
)
before_qc = adata.n_obs
adata = adata[valid_cells].copy()
print(f"[+] Retained after QC: {adata.n_obs:,} cells (filtered: {before_qc - adata.n_obs:,})")
gc.collect()

doublet_scores = pd.Series(0.0, index=adata.obs_names, dtype=float)
predicted_doublets = pd.Series(False, index=adata.obs_names, dtype=bool)
doublet_method = pd.Series("scrublet", index=adata.obs_names, dtype=object)

print("[*] Running doublet detection separately for each patient...")
for patient in sorted(adata.obs["patient"].astype(str).unique()):
    patient_mask = adata.obs["patient"].astype(str) == patient
    print(f"  -> Patient {patient}: {int(patient_mask.sum()):,} cells")
    sub_x = adata[patient_mask].X
    sub_x = sub_x.tocsr() if sparse.issparse(sub_x) else sparse.csr_matrix(sub_x)
    try:
        scrub = scr.Scrublet(sub_x, expected_doublet_rate=0.05)
        scores, calls = scrub.scrub_doublets(
            min_counts=2,
            min_cells=3,
            min_gene_variability_pctl=85,
            verbose=False,
        )
        doublet_scores.loc[patient_mask] = scores
        if scrub.threshold_ is None or calls is None:
            if not args.allow_heuristic_doublet_fallback:
                raise RuntimeError("Scrublet did not return a threshold or calls.")
            cutoff = np.nanpercentile(scores, 95)
            predicted_doublets.loc[patient_mask] = scores > cutoff
            doublet_method.loc[patient_mask] = "scrublet_95th_percentile_threshold"
        else:
            predicted_doublets.loc[patient_mask] = calls
    except Exception as error:
        if not args.allow_heuristic_doublet_fallback:
            raise RuntimeError(
                f"Scrublet failed for patient {patient!r}. To permit the labeled gene-count fallback, rerun with --allow-heuristic-doublet-fallback."
            ) from error
        print(f"     [!] Using explicitly enabled heuristic fallback for {patient}: {error}")
        high_counts = (
            adata.obs.loc[patient_mask, "n_genes_by_counts"]
            > adata.obs.loc[patient_mask, "n_genes_by_counts"].quantile(0.95)
        )
        predicted_doublets.loc[patient_mask] = high_counts.to_numpy()
        doublet_method.loc[patient_mask] = "gene_count_top_5_percent"
    del sub_x
    gc.collect()

adata.obs["doublet_score"] = doublet_scores
adata.obs["predicted_doublet"] = predicted_doublets
adata.obs["doublet_method"] = doublet_method
print("[+] Doublet calls by method:")
print(adata.obs["doublet_method"].value_counts(dropna=False))

n_doublets = int(predicted_doublets.sum())
print(f"[+] Predicted doublets: {n_doublets:,} ({n_doublets / max(adata.n_obs, 1) * 100:.2f}%)")
adata = adata[~adata.obs["predicted_doublet"]].copy()
print(f"[+] Final clean cells: {adata.n_obs:,}")

print(f"[*] Saving to: {OUT_FP}")
adata.write_h5ad(OUT_FP, compression="gzip")
print("[+] Step 05 complete.")
