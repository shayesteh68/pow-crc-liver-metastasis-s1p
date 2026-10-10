from paths import DATA_DIR, RESULTS_DIR, FIG_DIR, PLOTS_DIR

import os
import gzip
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import io
from scipy import sparse


mtx_fp = os.path.join(DATA_DIR, "GSE178318_matrix.mtx.gz")
genes_fp = os.path.join(DATA_DIR, "GSE178318_genes.tsv.gz")
barcodes_fp = os.path.join(DATA_DIR, "GSE178318_barcodes.tsv.gz")

print("[*] Reading matrix.mtx.gz (this can take a few minutes)...")
with gzip.open(mtx_fp, "rb") as f:
    X = io.mmread(f).tocsr()  # genes x cells (GEO MTX often is genes x cells)
print(f"[+] Raw matrix shape: {X.shape} (rows x cols)")

print("[*] Reading genes.tsv.gz ...")
genes = pd.read_csv(genes_fp, sep="\t", header=None)
# Expected: gene_id, gene_symbol
if genes.shape[1] == 1:
    genes.columns = ["gene_id"]
    genes["gene_symbol"] = genes["gene_id"]
else:
    genes = genes.iloc[:, :2]
    genes.columns = ["gene_id", "gene_symbol"]

print("[*] Reading barcodes.tsv.gz ...")
barcodes = pd.read_csv(barcodes_fp, sep="\t", header=None, names=["barcode_raw"])

# --- Detect orientation and build AnnData as cells x genes ---
n_genes_file = genes.shape[0]
n_cells_file = barcodes.shape[0]

# We expect final AnnData shape: n_cells x n_genes
# If X is genes x cells -> transpose
if X.shape[0] == n_genes_file and X.shape[1] == n_cells_file:
    print("[+] Detected orientation: genes x cells. Transposing to cells x genes...")
    X = X.transpose().tocsr()
elif X.shape[0] == n_cells_file and X.shape[1] == n_genes_file:
    print("[+] Detected orientation: cells x genes. Keeping as-is...")
else:
    raise RuntimeError(
        f"Matrix shape {X.shape} does not match genes ({n_genes_file}) and barcodes ({n_cells_file})."
    )

print(f"[+] Final X shape (cells x genes): {X.shape}")

print("[*] Parsing barcode fields -> cell_barcode, patient, tissue")
# Example: AAAC..._COL07_CRC
parts = barcodes["barcode_raw"].str.split("_", expand=True)
if parts.shape[1] < 3:
    raise RuntimeError("Barcode format is not BARCODE_PATIENT_TISSUE; cannot parse tissue/patient.")

barcodes["cell_barcode"] = parts[0].astype(str)
barcodes["patient"] = parts[1].astype(str)
barcodes["tissue"] = parts.iloc[:, 2:].agg("_".join, axis=1).astype(str)

# Construct AnnData
adata = sc.AnnData(X=X)
adata.obs = barcodes.set_index("barcode_raw")
adata.var = genes.set_index("gene_id")
adata.var["gene_symbol"] = genes["gene_symbol"].values

# Keep gene symbols as a convenient alias (without destroying gene_id index)
adata.var_names_make_unique()

# Basic sanity checks
print("[*] Sanity checks")
print("    cells:", adata.n_obs)
print("    genes:", adata.n_vars)
print("    tissues:", adata.obs["tissue"].value_counts().to_dict())
print("    patients:", adata.obs["patient"].nunique())

# Write full object (optional but useful)
full_out = os.path.join(RESULTS_DIR, "GSE178318_full_raw.h5ad")
print(f"[*] Writing: {full_out}")
adata.write_h5ad(full_out, compression="gzip")

# Split tumor tissues vs PBMC
tumor_mask = adata.obs["tissue"].isin(["CRC", "LM"])
pbmc_mask = adata.obs["tissue"].isin(["PBMC"])

adata_tumor = adata[tumor_mask].copy()
adata_pbmc = adata[pbmc_mask].copy()

tumor_out = os.path.join(RESULTS_DIR, "GSE178318_CRC_LM_raw.h5ad")
pbmc_out = os.path.join(RESULTS_DIR, "GSE178318_PBMC_raw.h5ad")

print(f"[*] Writing: {tumor_out}")
adata_tumor.write_h5ad(tumor_out, compression="gzip")

print(f"[*] Writing: {pbmc_out}")
adata_pbmc.write_h5ad(pbmc_out, compression="gzip")

print("[+] Done.")
