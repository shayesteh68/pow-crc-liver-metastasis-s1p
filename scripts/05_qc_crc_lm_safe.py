from paths import DATA_DIR, RESULTS_DIR, FIG_DIR, PLOTS_DIR

import os
import gc
import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr
from scipy import sparse

IN_FP = os.path.join(RESULTS_DIR, "GSE178318_CRC_LM_raw.h5ad")
OUT_FP = os.path.join(RESULTS_DIR, "GSE178318_CRC_LM_qc.h5ad")

print(f"[*] Loading raw CRC+LM object: {IN_FP}")
adata = sc.read_h5ad(IN_FP)

# 1. Mito calculation
if "gene_symbol" in adata.var.columns:
    mito_mask = adata.var["gene_symbol"].astype(str).str.upper().str.startswith("MT-")
else:
    mito_mask = np.array([False] * adata.n_vars)

adata.var["mt"] = mito_mask

print("[*] Calculating basic QC metrics...")
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True
)

# 2. Standard Threshold Filtering
min_genes = 200
max_genes = 6000
max_mt = 20.0

valid_cells = (
    (adata.obs["n_genes_by_counts"] >= min_genes) &
    (adata.obs["n_genes_by_counts"] <= max_genes) &
    (adata.obs["pct_counts_mt"] <= max_mt)
)

before_qc = adata.n_obs
adata = adata[valid_cells].copy()
after_qc = adata.n_obs
print(f"[+] Retained after QC: {after_qc:,} cells (Filtered: {before_qc - after_qc:,})")

gc.collect()

# 3. Ultra Memory-Efficient Doublet Detection (Sparse input to Scrublet)
print("\n[*] Running Memory-Optimized Doublet Detection per patient...")
doublet_scores = pd.Series(0.0, index=adata.obs_names, dtype=float)
predicted_doublets = pd.Series(False, index=adata.obs_names, dtype=bool)

for pt in sorted(adata.obs["patient"].unique()):
    pt_idx = adata.obs["patient"] == pt
    n_pt = pt_idx.sum()
    print(f"  -> Patient {pt}: {n_pt:,} cells...")
    
    # Slice sparse matrix directly without converting whole matrix to dense
    sub_X = adata[pt_idx].X
    if not sparse.issparse(sub_X):
        sub_X = sparse.csr_matrix(sub_X)
    else:
        sub_X = sub_X.tocsr()
        
    # Scrublet accepts scipy.sparse.csr_matrix natively! No toarray() needed.
    try:
        scrub = scr.Scrublet(sub_X, expected_doublet_rate=0.05)
        d_scores, d_preds = scrub.scrub_doublets(
            min_counts=2, 
            min_cells=3, 
            min_gene_variability_pctl=85, 
            verbose=False
        )
        doublet_scores.loc[pt_idx] = d_scores
        if scrub.threshold_ is not None and d_preds is not None:
            predicted_doublets.loc[pt_idx] = d_preds
        else:
            thresh = np.nanpercentile(d_scores, 95)
            predicted_doublets.loc[pt_idx] = d_scores > thresh
    except Exception as e:
        print(f"     [!] Fallback for {pt}: {e}")
        # High gene count as doublet proxy if simulation fails
        high_counts = adata.obs.loc[pt_idx, "n_genes_by_counts"] > adata.obs.loc[pt_idx, "n_genes_by_counts"].quantile(0.95)
        predicted_doublets.loc[pt_idx] = high_counts

    del sub_X
    gc.collect()

adata.obs["doublet_score"] = doublet_scores
adata.obs["predicted_doublet"] = predicted_doublets

n_doublets = predicted_doublets.sum()
print(f"\n[+] Total Predicted Doublets: {n_doublets:,} ({n_doublets/adata.n_obs*100:.2f}%)")

# Remove doublets
adata = adata[~adata.obs["predicted_doublet"]].copy()
print(f"[+] Final Clean Cells: {adata.n_obs:,}")

print(f"[*] Saving to: {OUT_FP}")
adata.write_h5ad(OUT_FP, compression="gzip")
print("[✓] Step 05 complete!")
