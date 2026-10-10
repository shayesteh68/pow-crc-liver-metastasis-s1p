from paths import DATA_DIR, RESULTS_DIR, FIG_DIR, PLOTS_DIR

import scanpy as sc
import os

IN_FP = os.path.join(RESULTS_DIR, "GSE178318_CRC_LM_clustered.h5ad")

adata = sc.read_h5ad(IN_FP)

# Fix: Drop the redundant column that causes the saving error
if 'gene_symbol' in adata.var.columns:
    adata.var = adata.var.drop(columns=['gene_symbol'])

# Save successfully
adata.write_h5ad(IN_FP, compression="gzip")
print("[✓] File saved successfully after cleaning index conflict.")
