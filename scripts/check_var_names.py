import scanpy as sc
import os

RESULTS_DIR = "/home/digital_home/targeted_researcher_discovery/pow_yeatman_colorectal/results"
IN_FP = os.path.join(RESULTS_DIR, "GSE178318_CRC_LM_clustered.h5ad")

adata = sc.read_h5ad(IN_FP)

print("--- Data Summary ---")
print(f"Total cells: {adata.n_obs}")
print(f"Total genes: {adata.n_vars}")
print("\nFirst 10 var_names (Index):")
print(adata.var_names[:10].tolist())

print("\nColumns in adata.var:")
print(adata.var.columns.tolist())

print("\nSample of adata.var head:")
print(adata.var.head(5))

# Check some canonical markers in var_names vs var columns
test_genes = ["EPCAM", "CD3D", "CD68", "COL1A1", "PECAM1", "SPHK1", "S1PR1"]
print("\nDirect matching check in var_names:")
for g in test_genes:
    print(f"  {g} in var_names: {g in adata.var_names}")

if not any(g in adata.var_names for g in test_genes):
    for col in adata.var.columns:
        matches = [g for g in test_genes if g in adata.var[col].values]
        print(f"  Matches in column '{col}': {matches}")
