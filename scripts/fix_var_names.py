import sys
import os
import pandas as pd
import scanpy as sc

sys.path.append(os.path.abspath(os.path.dirname(__file__)))
import paths

genes_file = os.path.join(paths.DATA_DIR, "GSE178318_genes.tsv.gz")
clustered_file = os.path.join(paths.RESULTS_DIR, "GSE178318_CRC_LM_clustered.h5ad")

print(f"[*] Reading genes reference: {genes_file}")
genes_df = pd.read_csv(genes_file, sep="\t", header=None)

if genes_df.shape[1] == 1:
    genes_df.columns = ["gene_symbol"]
    genes_df["gene_id"] = genes_df["gene_symbol"]
elif genes_df.shape[1] >= 2:
    genes_df = genes_df.iloc[:, :2]
    genes_df.columns = ["gene_id", "gene_symbol"]

print(f"[*] Reference genes count: {len(genes_df)}")

print(f"[*] Loading AnnData: {clustered_file}")
adata = sc.read_h5ad(clustered_file)
print(f"[*] AnnData gene count: {adata.n_vars}")

if len(genes_df) != adata.n_vars:
    raise ValueError(f"Gene count mismatch: reference has {len(genes_df)}, adata has {adata.n_vars}")

adata.var = pd.DataFrame(index=genes_df["gene_symbol"].astype(str).values)
adata.var["gene_id"] = genes_df["gene_id"].astype(str).values
adata.var["gene_symbol"] = genes_df["gene_symbol"].astype(str).values
adata.var_names_make_unique()
adata.var.index.name = None

print("[*] Verifying critical marker genes:")
test_genes = ["EPCAM", "CD3D", "CD68", "COL1A1", "PECAM1", "SPHK1", "S1PR1"]
for g in test_genes:
    print(f"    {g} present: {g in adata.var_names}")

print(f"[*] Overwriting repaired dataset: {clustered_file}")
adata.write_h5ad(clustered_file)
print("[✓] Successfully fixed and updated clustered AnnData!")
