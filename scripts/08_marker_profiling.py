import sys
import os
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath(os.path.dirname(__file__)))
import paths

sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=150, facecolor="white")

clustered_file = os.path.join(paths.RESULTS_DIR, "GSE178318_CRC_LM_clustered.h5ad")
print(f"[*] Loading clustered dataset: {clustered_file}")
adata = sc.read_h5ad(clustered_file)

lineage_markers = {
    "Epithelial/Tumor": ["EPCAM", "KRT8", "KRT18"],
    "T/NK Cells": ["CD3D", "CD3E", "CD8A", "NKG7"],
    "B Cells": ["CD79A", "MS4A1", "CD19"],
    "Plasma Cells": ["SDC1", "MZB1", "IGHG1"],
    "Myeloid/Macrophage": ["CD68", "CD14", "LYZ", "FCGR3A"],
    "Fibroblasts": ["COL1A1", "DCN", "ACTA2"],
    "Endothelial": ["PECAM1", "VWF", "CDH5"],
    "S1P Axis": ["SPHK1", "SPHK2", "S1PR1", "S1PR2", "S1PR3", "S1PR4", "S1PR5", "SGPP1", "SGPL1"]
}

available_markers = {}
for lineage, genes in lineage_markers.items():
    found = [g for g in genes if g in adata.var_names]
    if found:
        available_markers[lineage] = found

print("[*] Generating Lineage Marker DotPlot...")
os.makedirs(paths.PLOTS_DIR, exist_ok=True)
dp = sc.pl.dotplot(
    adata,
    var_names=available_markers,
    groupby="leiden_res05",
    standard_scale="var",
    return_fig=True,
    show=False
)
dotplot_path = os.path.join(paths.PLOTS_DIR, "lineage_markers_dotplot.png")
dp.savefig(dotplot_path, bbox_inches="tight")
print(f"[✓] Saved DotPlot: {dotplot_path}")

print("[*] Computing Differential Expression across Leiden Clusters...")
sc.tl.rank_genes_groups(adata, groupby="leiden_res05", method="wilcoxon", key_added="rank_leiden")

markers_df = sc.get.rank_genes_groups_df(adata, group=None, key="rank_leiden")
markers_csv_path = os.path.join(paths.RESULTS_DIR, "leiden_cluster_markers.csv")
markers_df.to_csv(markers_csv_path, index=False)
print(f"[✓] Saved Cluster Markers Table: {markers_csv_path}")

print("\n--- Top 5 Markers per Leiden Cluster ---")
top_markers = markers_df.groupby("group").head(5)
for clus, df_c in top_markers.groupby("group"):
    genes_str = ", ".join(df_c["names"].tolist())
    print(f"Cluster {clus}: {genes_str}")

adata.var.index.name = None
adata.write_h5ad(clustered_file)
print(f"[✓] Successfully updated: {clustered_file}")
