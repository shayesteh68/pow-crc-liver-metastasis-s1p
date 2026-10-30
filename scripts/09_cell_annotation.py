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
    # Reviewed against the marker table from the seeded GSE178318 run on 2026-10-28.
    # Leiden IDs are implementation labels; this exact-ID guard must remain in place.
    "0": "Epithelial / Tumor Cells",
    "1": "Epithelial / Tumor Cells",
    "2": "Myeloid / Macrophages",
    "3": "Stress-associated Cytotoxic Lymphocytes (provisional)",
    "4": "Cytotoxic NK Cells",
    "5": "Myeloid Antigen-presenting Cells",
    "6": "Plasma Cells (IgG-enriched)",
    "7": "B Cells",
    "8": "Naive / Memory T Cells",
    "9": "GZMK+ Cytotoxic T/NK Cells",
    "10": "Plasma Cells (IgA-enriched)",
    "11": "Activated / Regulatory T Cells",
    "12": "Stress-associated Lymphocytes (provisional)",
    "13": "CXCL13+ Cytotoxic Lymphocytes (provisional)",
    "14": "Plasma Cells (IgA-enriched)",
    "15": "Mast Cells",
    "16": "Plasmacytoid Dendritic Cells",
    "17": "Fibroblasts / Stromal Cells",
    "18": "Epithelial / Tumor Cells",
    "19": "Mixed IgA-high Cytotoxic Lymphocytes (unresolved)",
}

cluster_annotation_evidence = {
    "0": "KRT8, KRT18, KRT19, EPCAM",
    "1": "CLDN3, TFF3, ELF3, KRT8, EPCAM",
    "2": "TYROBP, FCER1G, LYZ, CD68, AIF1",
    "3": "HSPA1A, HSPA1B, HSP90AA1 with CCL5, NKG7, GNLY",
    "4": "NKG7, GNLY, KLRD1, GZMB, PRF1",
    "5": "HLA-DRA, LST1, LYZ, CD74, HLA-DPB1",
    "6": "IGHG1, MZB1, IGHG4, JCHAIN, XBP1",
    "7": "CD79A, CD74, HLA-DRA",
    "8": "IL7R, TSC22D3, CXCR4, ZFP36L2",
    "9": "GZMK, CCL4, CCL5, GZMA, NKG7",
    "10": "IGHA1, IGHA2, JCHAIN, immunoglobulin light chains",
    "11": "IL2RA, CTLA4, TIGIT, TNFRSF4, TNFRSF18",
    "12": "HSPA1A, HSPA1B, DNAJB1, HSP90AA1",
    "13": "CXCL13, GZMB, CD3D, RGS1",
    "14": "JCHAIN, IGHA1, MZB1, IGHA2, TNFRSF17",
    "15": "TPSB2, TPSAB1, CPA3, MS4A2",
    "16": "GZMB, LILRA4, PLD4, TCF4, IRF7",
    "17": "IGFBP7, CALD1, SPARC, COL6A2, MGP",
    "18": "ADIRF, KRT18, KRT8, KRT19, MSLN",
    "19": "IGHA1, IGHA2, JCHAIN with NKG7, CCL5, GZMK; mixed profile",
}

cluster_annotation_status = {
    cluster: "provisional_manual_marker_reviewed" for cluster in cluster_annotations
}
for cluster in ("3", "12", "13", "19"):
    cluster_annotation_status[cluster] = "provisional_mixed_or_stress_profile"

if "leiden_res05" not in adata.obs.columns:
    raise KeyError("leiden_res05 is missing; run clustering first.")
if adata.obs["leiden_res05"].isna().any():
    raise ValueError("leiden_res05 contains missing cluster IDs.")
observed_clusters = set(adata.obs["leiden_res05"].astype(str).unique())
expected_clusters = set(cluster_annotations)
missing_clusters = sorted(expected_clusters - observed_clusters, key=int)
unexpected_clusters = sorted(observed_clusters - expected_clusters, key=lambda value: (not value.isdigit(), value))
if missing_clusters or unexpected_clusters:
    raise ValueError(
        "Cluster IDs differ from the reviewed marker map. Missing clusters: "
        + ", ".join(missing_clusters)
        + "; unexpected clusters: "
        + ", ".join(unexpected_clusters)
    )

print("[*] Assigning cell type annotations...")
adata.obs["cell_type"] = adata.obs["leiden_res05"].astype(str).map(cluster_annotations).astype("category")
adata.obs["annotation_status"] = (
    adata.obs["leiden_res05"].astype(str).map(cluster_annotation_status).astype("category")
)
adata.uns["cluster_annotation_evidence"] = cluster_annotation_evidence
adata.uns["cluster_annotation_status"] = cluster_annotation_status

print("\n--- Cell Type Distribution ---")
print(adata.obs["cell_type"].value_counts())

if "X_umap" not in adata.obsm:
    print("[*] X_umap not found in obsm. Computing neighbors and UMAP...")
    if "X_pca" not in adata.obsm:
        print("[*] Computing PCA...")
        sc.tl.pca(adata, svd_solver="arpack", n_comps=50, random_state=0)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=0)
    sc.tl.umap(adata, random_state=0)

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
