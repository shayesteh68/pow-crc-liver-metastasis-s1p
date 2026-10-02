import os
import scanpy as sc
import matplotlib.pyplot as plt

RESULTS_DIR = "/home/digital_home/targeted_researcher_discovery/pow_yeatman_colorectal/results"
PLOTS_DIR = "/home/digital_home/targeted_researcher_discovery/pow_yeatman_colorectal/results/plots"
os.makedirs(PLOTS_DIR, exist_ok=True)

IN_FP = os.path.join(RESULTS_DIR, "GSE178318_CRC_LM_clustered.h5ad")

print(f"[*] Loading clustered dataset: {IN_FP}")
adata = sc.read_h5ad(IN_FP)

# Map gene symbols to var_names for standard downstream access
print("[*] Setting var_names to gene symbols...")
adata.var["gene_id"] = adata.var_names
adata.var_names = adata.var["gene_symbol"].astype(str)
adata.var_names_make_unique()

# Define Canonical Lineage Markers
marker_dict = {
    "Tumor_Epithelial": ["EPCAM", "KRT8", "KRT18"],
    "T_NK": ["CD3D", "CD3E", "CD8A", "NKG7"],
    "B_Plasma": ["MS4A1", "CD79A", "MZB1"],
    "TAMs_Myeloid": ["CD68", "CD163", "CSF1R", "MARCO"],
    "CAFs_Stromal": ["COL1A1", "COL1A2", "ACTA2", "FAP", "DCN"],
    "Endothelial": ["PECAM1", "VWF", "CDH5"],
    "Mast": ["KIT", "TPSAB1"]
}

# Filter only markers that exist in our dataset
valid_markers = {k: [g for g in v if g in adata.var_names] for k, v in marker_dict.items()}

print("[*] Generating Marker DotPlot...")
sc.settings.set_figure_params(dpi=150, facecolor="white")
dp = sc.pl.dotplot(
    adata,
    valid_markers,
    groupby="leiden_res05",
    dendrogram=False,
    standard_scale="var",
    return_fig=True,
    show=False
)
dp.savefig(os.path.join(PLOTS_DIR, "lineage_markers_dotplot.png"), bbox_inches="tight")

print("[*] Finding Top Differential Markers per Leiden Cluster...")
sc.tl.rank_genes_groups(adata, groupby="leiden_res05", method="wilcoxon", key_added="rank_leiden")
markers_df = sc.get.rank_genes_groups_df(adata, group=None, key="rank_leiden")

# Print summary table of top 5 markers for each cluster
print("\n--- Top Markers per Leiden Cluster ---")
top_markers = markers_df.groupby("group").head(5)
for clus, df_c in top_markers.groupby("group"):
    genes = ", ".join(df_c["names"].tolist())
    print(f"Cluster {clus}: {genes}")

# Save the updated object with gene symbols as var_names
adata.write_h5ad(IN_FP, compression="gzip")

print(f"\n[✓] Dotplot saved to: {os.path.join(PLOTS_DIR, 'lineage_markers_dotplot.png')}")
print("[✓] Step 08 complete!")
