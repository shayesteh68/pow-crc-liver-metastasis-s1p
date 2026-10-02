import os
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse

RESULTS_DIR = "/home/digital_home/targeted_researcher_discovery/pow_yeatman_colorectal/results"
IN_FP = os.path.join(RESULTS_DIR, "GSE178318_CRC_LM_clustered.h5ad")

print(f"[*] Loading dataset: {IN_FP}")
adata = sc.read_h5ad(IN_FP)

# Check where the expression matrix resides
if adata.raw is not None and adata.raw.X is not None:
    print("[*] Found expression matrix in adata.raw.X")
    X = adata.raw.X
    genes = adata.raw.var_names.to_numpy()
elif adata.X is not None:
    print("[*] Found expression matrix in adata.X")
    X = adata.X
    genes = adata.var_names.to_numpy()
elif len(adata.layers) > 0:
    first_layer = list(adata.layers.keys())[0]
    print(f"[*] Found expression matrix in adata.layers['{first_layer}']")
    X = adata.layers[first_layer]
    genes = adata.var_names.to_numpy()
else:
    raise ValueError("No expression matrix found in adata.X, adata.raw, or adata.layers!")

if not sparse.isspmatrix_csr(X):
    X = sparse.csr_matrix(X)

cluster_col = "leiden_res05"
clusters = sorted(adata.obs[cluster_col].unique(), key=lambda x: int(x))

print("[*] Computing cluster markers via sparse differential means...")
overall_mean = np.array(X.mean(axis=0)).flatten()

print("\n" + "="*80)
print("TOP 6 DISTINCTIVE MARKERS PER CLUSTER")
print("="*80)

summary_rows = []
for clus in clusters:
    cell_mask = (adata.obs[cluster_col] == clus).to_numpy()
    cell_count = int(cell_mask.sum())
    pct = (cell_count / adata.n_obs) * 100
    
    clus_mean = np.array(X[cell_mask, :].mean(axis=0)).flatten()
    diff = clus_mean - overall_mean
    top_indices = np.argsort(diff)[::-1][:6]
    top_markers = genes[top_indices].tolist()
    top_genes_str = ", ".join(top_markers)
    
    print(f"Cluster {clus:>2} ({cell_count:>5} cells | {pct:4.1f}%): {top_genes_str}")
    summary_rows.append({
        "cluster": clus,
        "cell_count": cell_count,
        "pct": pct,
        "top_genes": top_genes_str
    })

out_csv = os.path.join(RESULTS_DIR, "cluster_marker_summary.csv")
pd.DataFrame(summary_rows).to_csv(out_csv, index=False)
print("="*80)
print(f"[✓] Summary saved to: {out_csv}")
