import sys
import os
import scanpy as sc
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath(os.path.dirname(__file__)))
import paths

sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=150, facecolor="white")

clustered_file = os.path.join(paths.RESULTS_DIR, "GSE178318_CRC_LM_clustered.h5ad")
print(f"[*] Loading dataset: {clustered_file}")
adata = sc.read_h5ad(clustered_file)

cluster_annotations = {
    "0": "Epithelial / Tumor Cells",
    "1": "Epithelial / Tumor Cells",
    "2": "Myeloid / Macrophages",
    "3": "CAFs / Stromal Cells",
    "4": "B Cells",
    "5": "CD8+ T Cells (Cytotoxic)",
    "6": "CD4+ T / Exhausted T Cells",
    "7": "Naive / Memory T Cells",
    "8": "Activated T Cells",
    "9": "NK Cells",
    "10": "Stress-response Cells",
    "11": "Plasma Cells (IgA+)",
    "12": "Plasma Cells (IgA+)",
    "13": "Plasma Cells (General)",
    "14": "Plasma Cells (IgG+)",
    "15": "Mast Cells",
    "16": "Plasmacytoid Dendritic Cells",
    "17": "Dendritic / Antigen-Presenting Cells",
    "18": "Stress-response Cells",
    "19": "Epithelial / Tumor Cells"
}

print("[*] Assigning cell type annotations...")
adata.obs["cell_type"] = adata.obs["leiden_res05"].astype(str).map(cluster_annotations).astype("category")

print("\n--- Cell Type Distribution ---")
print(adata.obs["cell_type"].value_counts())

if "X_umap" not in adata.obsm:
    print("[*] X_umap not found in obsm. Computing neighbors and UMAP...")
    if "X_pca" not in adata.obsm:
        print("[*] Computing PCA...")
        sc.tl.pca(adata, svd_solver="arpack", n_comps=50)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata)

print("[*] Generating UMAP plots with cell annotations...")
os.makedirs(paths.PLOTS_DIR, exist_ok=True)

fig, ax = plt.subplots(figsize=(10, 8))
sc.pl.umap(
    adata,
    color="cell_type",
    title="Cell Types Annotation (GSE178318 CRC & LM)",
    legend_loc="on data",
    legend_fontsize=8,
    frameon=False,
    show=False,
    ax=ax
)
umap_annot_path = os.path.join(paths.PLOTS_DIR, "umap_cell_types.png")
fig.savefig(umap_annot_path, bbox_inches="tight", dpi=300)
plt.close(fig)
print(f"[✓] Saved UMAP plot: {umap_annot_path}")

annotated_file = os.path.join(paths.RESULTS_DIR, "GSE178318_CRC_LM_annotated.h5ad")
adata.var.index.name = None
print(f"[*] Saving annotated AnnData: {annotated_file}")
adata.write_h5ad(annotated_file)
print("[✓] Step 09 completed successfully!")
