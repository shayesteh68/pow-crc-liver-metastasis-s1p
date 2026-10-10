from paths import DATA_DIR, RESULTS_DIR, FIG_DIR, PLOTS_DIR

import os
import scanpy as sc
import matplotlib.pyplot as plt


IN_FP = os.path.join(RESULTS_DIR, "GSE178318_CRC_LM_clustered.h5ad")
OUT_FP = os.path.join(RESULTS_DIR, "GSE178318_CRC_LM_annotated.h5ad")

print(f"[*] Loading clustered data: {IN_FP}")
adata = sc.read_h5ad(IN_FP)

annotation_dict = {
    "0": "Epithelial_Tumor",
    "1": "Myeloid_TAM",
    "2": "CAF_Fibroblast",
    "3": "Epithelial_Secretory",
    "4": "Plasma_Cell",
    "5": "Inflammatory_Monocyte",
    "6": "B_Cell",
    "7": "T_Cell_Activated",
    "8": "Plasma_Cell_IgA",
    "9": "T_Cell_Helper",
    "10": "NK_Cytotoxic_T",
    "11": "T_Cell_CXCL13",
    "12": "Stress_Response_T_NK",
    "13": "Stress_Response_General",
    "14": "Plasma_Cell_IgA",
    "15": "Mast_Cell",
    "16": "Plasmacytoid_DC",
    "17": "Plasma_Cell_IgG"
}

# Major lineage grouping for high-level visualization
major_dict = {
    "0": "Epithelial",
    "1": "Myeloid",
    "2": "Stromal_CAF",
    "3": "Epithelial",
    "4": "B_Plasma",
    "5": "Myeloid",
    "6": "B_Plasma",
    "7": "T_NK",
    "8": "B_Plasma",
    "9": "T_NK",
    "10": "T_NK",
    "11": "T_NK",
    "12": "Stress_Response",
    "13": "Stress_Response",
    "14": "B_Plasma",
    "15": "Mast",
    "16": "Myeloid",
    "17": "B_Plasma"
}

adata.obs["cell_type_detailed"] = adata.obs["leiden_res05"].map(annotation_dict).astype("category")
adata.obs["cell_type_major"] = adata.obs["leiden_res05"].map(major_dict).astype("category")

print("[*] Generating UMAP plots...")
sc.settings.set_figure_params(dpi=300, facecolor="white")

fig, axes = plt.subplots(1, 2, figsize=(18, 7))
sc.pl.umap(adata, color="cell_type_major", title="Major Lineages", ax=axes[0], show=False, frameon=False)
sc.pl.umap(adata, color="tissue", title="Tissue Distribution (CRC Primary vs Liver Met)", ax=axes[1], show=False, frameon=False)
plt.tight_layout()
umap_fp = os.path.join(FIG_DIR, "umap_celltypes_tissue.png")
fig.savefig(umap_fp, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"[✓] UMAP plot saved to: {umap_fp}")

print(f"[*] Saving annotated AnnData to: {OUT_FP}")
adata.write_h5ad(OUT_FP, compression="gzip")
print("[✓] Step 10 Annotation completed successfully!")
