import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from paths import FIG_DIR, RESULTS_DIR

input_file = os.path.join(RESULTS_DIR, "GSE178318_CRC_LM_annotated.h5ad")
s1p_genes = [
    "SPHK1", "SPHK2", "SGPL1", "SGPP1", "SGPP2",
    "S1PR1", "S1PR2", "S1PR3", "S1PR4", "S1PR5",
    "CERK", "ASAH1", "PLPP1", "PLPP2", "PLPP3",
]

print(f"[*] Loading annotated dataset: {input_file}")
adata = sc.read_h5ad(input_file)
if adata.X is None:
    raise ValueError("adata.X is None in the annotated file.")
for column in ("cell_type", "tissue"):
    if column not in adata.obs.columns:
        raise ValueError(f"Annotated dataset is missing adata.obs[{column!r}]; run step 09 first.")

matrix = adata.X
if sparse.issparse(matrix):
    matrix = matrix.tocsr()
else:
    matrix = sparse.csr_matrix(np.asarray(matrix))

genes_present = [gene for gene in s1p_genes if gene in adata.var_names]
missing = [gene for gene in s1p_genes if gene not in adata.var_names]
print(f"[*] S1P genes present: {len(genes_present)}/{len(s1p_genes)}")
if missing:
    print("[!] Missing genes (not in var_names): " + ", ".join(missing))
if not genes_present:
    raise ValueError("None of the S1P genes were found in adata.var_names; confirm var_names are gene symbols.")

cell_types = adata.obs["cell_type"].astype("string")
tissues = adata.obs["tissue"].astype("string")
if cell_types.isna().any() or tissues.isna().any():
    raise ValueError("cell_type and tissue must not contain missing values.")
group_labels = cell_types.astype(str) + " | " + tissues.astype(str)
groups = sorted(group_labels.unique().tolist())

gene_indices = [adata.var_names.get_loc(gene) for gene in genes_present]
selected = matrix[:, gene_indices]
means = []
counts = []
for group in groups:
    mask = group_labels.to_numpy() == group
    counts.append(int(mask.sum()))
    means.append(np.asarray(selected[mask, :].mean(axis=0)).ravel())

mean_df = pd.DataFrame(means, index=groups, columns=genes_present)
count_series = pd.Series(counts, index=groups, name="n_cells")

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)
mean_path = os.path.join(RESULTS_DIR, "s1p_mean_by_celltype_tissue.csv")
mean_df.to_csv(mean_path)
count_path = os.path.join(RESULTS_DIR, "s1p_group_cellcounts.csv")
count_series.to_csv(count_path)
print(f"[*] Saved descriptive mean matrix: {mean_path}")
print(f"[*] Saved cell counts: {count_path}")

# This figure summarizes per-cell means; it does not perform patient-level inference.
values = mean_df.to_numpy()
vmin, vmax = float(np.nanmin(values)), float(np.nanmax(values))
normalized = (values - vmin) / (vmax - vmin + 1e-9)
plt.figure(figsize=(max(10, 0.6 * len(genes_present)), max(6, 0.35 * len(groups))))
for row, group in enumerate(groups):
    for column, gene in enumerate(genes_present):
        plt.scatter(
            column,
            row,
            s=20 + 180 * normalized[row, column],
            c=[values[row, column]],
            cmap="viridis",
            vmin=vmin,
            vmax=vmax,
        )
plt.yticks(range(len(groups)), groups)
plt.xticks(range(len(genes_present)), genes_present, rotation=45, ha="right")
plt.gca().invert_yaxis()
plt.colorbar(label="Mean expression (log-normalized)")
plt.title("S1P pathway genes: mean expression by cell type and tissue")
plt.tight_layout()

figure_path = os.path.join(FIG_DIR, "s1p_dotplot_celltype_tissue.png")
plt.savefig(figure_path, dpi=300, bbox_inches="tight")
plt.close()
print(f"[*] Saved figure: {figure_path}")
