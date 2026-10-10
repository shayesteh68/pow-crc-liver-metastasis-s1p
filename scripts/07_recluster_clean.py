from paths import DATA_DIR, RESULTS_DIR, FIG_DIR, PLOTS_DIR

import os
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse

IN_FP = os.path.join(RESULTS_DIR, "GSE178318_CRC_LM_processed.h5ad")
OUT_FP = os.path.join(RESULTS_DIR, "GSE178318_CRC_LM_clustered.h5ad")

print(f"[*] Reading clean processed dataset: {IN_FP}")
adata = sc.read_h5ad(IN_FP)
print(f"    Loaded shape: {adata.shape}")

# Ensure gene symbols are var_names and clean up duplicate column names
if 'gene_symbol' in adata.var.columns:
    adata.var_names = adata.var['gene_symbol'].astype(str)
    adata.var_names_make_unique()
    adata.var.drop(columns=['gene_symbol'], inplace=True)
adata.var.index.name = None

print("[*] Running PCA (30 components)...")
sc.tl.pca(adata, n_comps=30, svd_solver="arpack")

print("[*] Computing Neighborhood Graph (n_neighbors=15, n_pcs=30)...")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)

print("[*] Running Leiden clustering (resolution=0.5, flavor='igraph')...")
sc.tl.leiden(adata, resolution=0.5, key_added="leiden_res05", flavor="igraph", n_iterations=2, directed=False)

print("[*] Computing UMAP coordinates...")
sc.tl.umap(adata, min_dist=0.3)

print(f"[*] Saving verified clustered object to: {OUT_FP}")
adata.write_h5ad(OUT_FP, compression="gzip")
print("[✓] Clustered object saved successfully.")

# Extract markers cleanly
print("\n[*] Computing distinctive cluster markers...")
X = adata.X if sparse.isspmatrix_csr(adata.X) else adata.X.tocsr()
overall_mean = np.array(X.mean(axis=0)).flatten()
genes = adata.var_names.to_numpy()
cluster_col = "leiden_res05"
clusters = sorted(adata.obs[cluster_col].unique(), key=lambda x: int(x))

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

pd.DataFrame(summary_rows).to_csv(os.path.join(RESULTS_DIR, "cluster_marker_summary.csv"), index=False)
print("="*80)
print(f"[✓] Marker summary saved to: {os.path.join(RESULTS_DIR, 'cluster_marker_summary.csv')}")
