import scanpy as sc
import os

RESULTS_DIR = "/home/digital_home/targeted_researcher_discovery/pow_yeatman_colorectal/results"
IN_FP = os.path.join(RESULTS_DIR, "GSE178318_CRC_LM_qc.h5ad")
OUT_FP = os.path.join(RESULTS_DIR, "GSE178318_CRC_LM_processed.h5ad")

print(f"[*] Loading clean data: {IN_FP}")
adata = sc.read_h5ad(IN_FP)

# Normalization (Total count normalization)
print("[*] Normalizing data...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Highly Variable Genes
print("[*] Identifying Highly Variable Genes (HVG)...")
sc.pp.highly_variable_genes(adata, n_top_genes=3000, subset=False, flavor="seurat")

# Save
adata.write_h5ad(OUT_FP, compression="gzip")
print(f"[✓] Normalization complete. File saved: {OUT_FP}")
