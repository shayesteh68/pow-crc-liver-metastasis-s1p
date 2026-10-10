from paths import DATA_DIR, RESULTS_DIR, FIG_DIR, PLOTS_DIR

import os
import gc
import scanpy as sc

IN_FP = os.path.join(RESULTS_DIR, "GSE178318_CRC_LM_processed.h5ad")
OUT_FP = os.path.join(RESULTS_DIR, "GSE178318_CRC_LM_clustered.h5ad")

print(f"[*] Loading processed dataset: {IN_FP}")
adata = sc.read_h5ad(IN_FP)

# Run PCA on the 3000 Highly Variable Genes (memory efficient)
print("[*] Running PCA (50 components)...")
sc.tl.pca(adata, svd_solver="arpack", use_highly_variable=True, n_comps=50)

# Compute Neighborhood Graph
print("[*] Computing Neighborhood Graph (n_neighbors=15, n_pcs=30)...")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)

# Fast igraph Leiden Clustering
print("[*] Running Leiden clustering (resolution=0.5, flavor='igraph')...")
sc.tl.leiden(adata, resolution=0.5, key_added="leiden_res05", flavor="igraph", directed=False, n_iterations=2)

# Compute UMAP for 2D visualization
print("[*] Computing UMAP coordinates...")
sc.tl.umap(adata, min_dist=0.3)

gc.collect()

print(f"[*] Saving clustered object: {OUT_FP}")
adata.write_h5ad(OUT_FP, compression="gzip")
print("[✓] Step 07 completed successfully!")
