import os
import numpy as np
import scanpy as sc
import scrublet as scr
from scipy import sparse

RESULTS_DIR = "/home/digital_home/targeted_researcher_discovery/pow_yeatman_colorectal/results"
IN_FP = os.path.join(RESULTS_DIR, "GSE178318_CRC_LM_raw.h5ad")
OUT_FP = os.path.join(RESULTS_DIR, "GSE178318_CRC_LM_qc.h5ad")

print(f"[*] Reading: {IN_FP}")
adata = sc.read_h5ad(IN_FP)

# --- QC metrics ---
# Detect mitochondrial genes (human: MT- prefix in gene symbol)
if "gene_symbol" in adata.var.columns:
    mito_mask = adata.var["gene_symbol"].astype(str).str.upper().str.startswith("MT-")
else:
    # fallback: var_names (gene_id) won't have MT-, so this is conservative
    mito_mask = np.array([False] * adata.n_vars)

adata.var["mt"] = mito_mask

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True
)

print("[*] QC summary (before filtering)")
print("    cells:", adata.n_obs)
print("    genes:", adata.n_vars)
print("    n_genes_by_counts: min/median/max =",
      int(adata.obs["n_genes_by_counts"].min()),
      int(adata.obs["n_genes_by_counts"].median()),
      int(adata.obs["n_genes_by_counts"].max()))
print("    total_counts: min/median/max =",
      int(adata.obs["total_counts"].min()),
      int(adata.obs["total_counts"].median()),
      int(adata.obs["total_counts"].max()))
print("    pct_counts_mt: min/median/max =",
      float(adata.obs["pct_counts_mt"].min()),
      float(adata.obs["pct_counts_mt"].median()),
      float(adata.obs["pct_counts_mt"].max()))

# --- Conservative filtering (can be tuned later) ---
# Keep cells with reasonable complexity; avoid harsh cutoffs now.
min_genes = 200
max_genes = 8000
max_mt = 25.0

before = adata.n_obs
adata = adata[
    (adata.obs["n_genes_by_counts"] >= min_genes) &
    (adata.obs["n_genes_by_counts"] <= max_genes) &
    (adata.obs["pct_counts_mt"] <= max_mt)
].copy()
after = adata.n_obs
print(f"[*] Filtered cells: {before} -> {after} (removed {before-after})")

# --- Doublet detection with Scrublet ---
# Scrublet expects dense counts; we do per-batch chunks to avoid RAM blowups.
# We'll run on a random subset for scoring threshold estimation, then score all in chunks.
print("[*] Doublet detection (Scrublet)")
np.random.seed(0)

X = adata.X
if sparse.issparse(X):
    X_csr = X.tocsr()
else:
    X_csr = sparse.csr_matrix(X)

n = adata.n_obs
subset_n = min(50000, n)  # cap for threshold estimation
subset_idx = np.random.choice(n, size=subset_n, replace=False)

scrub = scr.Scrublet(X_csr[subset_idx].toarray(), expected_doublet_rate=0.06)
subset_scores, subset_preds = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85)

# Use Scrublet's automatically set threshold
thr = scrub.threshold_
print(f"    estimated threshold = {thr:.4f}")
print(f"    subset doublet fraction = {subset_preds.mean():.4f} on n={subset_n}")

# Score all cells in chunks to limit memory
chunk = 20000
all_scores = np.zeros(n, dtype=np.float32)

for start in range(0, n, chunk):
    end = min(n, start + chunk)
    X_chunk = X_csr[start:end].toarray()
    # Reuse scrublet pipeline by instantiating a new object per chunk with same threshold logic
    s = scr.Scrublet(X_chunk, expected_doublet_rate=0.06)
    scores, _ = s.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85)
    all_scores[start:end] = scores
    print(f"    scored cells {start}:{end}")

adata.obs["doublet_score"] = all_scores
adata.obs["predicted_doublet"] = adata.obs["doublet_score"] >= thr

dbl_before = adata.n_obs
adata = adata[~adata.obs["predicted_doublet"]].copy()
dbl_after = adata.n_obs
print(f"[*] Removed predicted doublets: {dbl_before} -> {dbl_after} (removed {dbl_before-dbl_after})")

print(f"[*] Writing QC'ed object: {OUT_FP}")
adata.write_h5ad(OUT_FP, compression="gzip")
print("[+] Done.")
