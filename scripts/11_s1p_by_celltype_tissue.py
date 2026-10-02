import os
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt

RESULTS_DIR = "/home/digital_home/targeted_researcher_discovery/pow_yeatman_colorectal/results"
FIG_DIR = "/home/digital_home/targeted_researcher_discovery/pow_yeatman_colorectal/figures"
os.makedirs(FIG_DIR, exist_ok=True)

IN_FP = os.path.join(RESULTS_DIR, "GSE178318_CRC_LM_annotated.h5ad")

S1P_GENES = [
    "SPHK1","SPHK2","SGPL1","SGPP1","SGPP2",
    "S1PR1","S1PR2","S1PR3","S1PR4","S1PR5",
    "CERK","ASAH1","PLPP1","PLPP2","PLPP3"
]

print(f"[*] Loading annotated dataset: {IN_FP}")
adata = sc.read_h5ad(IN_FP)

# Ensure expression matrix exists
X = adata.X
if X is None:
    raise ValueError("adata.X is None in annotated file. This should not happen—check the file.")
if not sparse.isspmatrix_csr(X):
    X = X.tocsr()

# Intersect genes
genes_present = [g for g in S1P_GENES if g in adata.var_names]
missing = [g for g in S1P_GENES if g not in adata.var_names]

print(f"[*] S1P genes present: {len(genes_present)}/{len(S1P_GENES)}")
if missing:
    print("[!] Missing genes (not in var_names): " + ", ".join(missing))

if len(genes_present) == 0:
    raise ValueError("None of the S1P genes were found in adata.var_names. Confirm var_names are gene symbols.")

celltype_col = "cell_type_major"
tissue_col = "tissue"

# Prepare group labels
adata.obs["group"] = (adata.obs[celltype_col].astype(str) + " | " + adata.obs[tissue_col].astype(str)).astype("category")
groups = list(adata.obs["group"].cat.categories)

# Extract submatrix for genes
gene_idx = [adata.var_names.get_loc(g) for g in genes_present]
Xg = X[:, gene_idx]

# Compute mean per group (sparse-safe)
means = []
counts = []
for grp in groups:
    mask = (adata.obs["group"] == grp).to_numpy()
    counts.append(int(mask.sum()))
    if counts[-1] == 0:
        means.append(np.zeros(len(genes_present), dtype=float))
    else:
        means.append(np.array(Xg[mask, :].mean(axis=0)).ravel())

mean_df = pd.DataFrame(means, index=groups, columns=genes_present)
count_s = pd.Series(counts, index=groups, name="n_cells")

out_csv = os.path.join(RESULTS_DIR, "s1p_mean_by_celltype_tissue.csv")
mean_df.to_csv(out_csv)
count_csv = os.path.join(RESULTS_DIR, "s1p_group_cellcounts.csv")
count_s.to_csv(count_csv)

print(f"[✓] Saved mean matrix: {out_csv}")
print(f"[✓] Saved cell counts: {count_csv}")

# Simple dotplot-like figure (size encodes mean; color encodes mean)
# Normalize for visualization only
vals = mean_df.to_numpy()
vmin, vmax = np.nanmin(vals), np.nanmax(vals)
norm = (vals - vmin) / (vmax - vmin + 1e-9)

plt.figure(figsize=(max(10, 0.6*len(genes_present)), max(6, 0.35*len(groups))))
for i, grp in enumerate(groups):
    for j, g in enumerate(genes_present):
        size = 20 + 180 * norm[i, j]
        plt.scatter(j, i, s=size, c=vals[i, j], cmap="viridis", vmin=vmin, vmax=vmax)
plt.yticks(range(len(groups)), groups)
plt.xticks(range(len(genes_present)), genes_present, rotation=45, ha="right")
plt.gca().invert_yaxis()
plt.colorbar(label="Mean expression (log-normalized)")
plt.title("S1P pathway genes: mean expression by cell type and tissue")
plt.tight_layout()

fig_fp = os.path.join(FIG_DIR, "s1p_dotplot_celltype_tissue.png")
plt.savefig(fig_fp, dpi=300, bbox_inches="tight")
plt.close()
print(f"[✓] Saved figure: {fig_fp}")
