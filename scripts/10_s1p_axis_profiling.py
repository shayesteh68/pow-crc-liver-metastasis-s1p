import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


sys.path.append(os.path.abspath(os.path.dirname(__file__)))
import paths

sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=150, facecolor="white")

annotated_file = os.path.join(paths.RESULTS_DIR, "GSE178318_CRC_LM_annotated.h5ad")
print(f"[*] Loading annotated dataset: {annotated_file}")
adata = sc.read_h5ad(annotated_file)

if "cell_type" not in adata.obs.columns:
    raise ValueError("Annotated dataset is missing adata.obs['cell_type']; run step 09 first.")
if adata.X is None:
    raise ValueError("Annotated dataset has no expression matrix in adata.X.")

s1p_genes = ["SPHK1", "SPHK2", "S1PR1", "S1PR2", "S1PR3", "S1PR4", "S1PR5", "SGPP1", "SGPL1"]
available_s1p = [gene for gene in s1p_genes if gene in adata.var_names]
print(f"[*] S1P axis genes available: {available_s1p}")
if not available_s1p:
    raise ValueError("None of the S1P genes were found in adata.var_names; confirm feature names are gene symbols.")

os.makedirs(paths.RESULTS_DIR, exist_ok=True)
os.makedirs(paths.PLOTS_DIR, exist_ok=True)

print("[*] Generating S1P axis DotPlot by cell type...")
dotplot = sc.pl.dotplot(
    adata,
    var_names=available_s1p,
    groupby="cell_type",
    standard_scale="var",
    return_fig=True,
    show=False,
)
dotplot_path = os.path.join(paths.PLOTS_DIR, "s1p_axis_celltypes_dotplot.png")
dotplot.savefig(dotplot_path, bbox_inches="tight")
print(f"[*] Saved S1P cell type DotPlot: {dotplot_path}")

if "tissue" in adata.obs.columns:
    print("[*] Generating S1P axis DotPlot by tissue (CRC vs LM)...")
    tissue_plot = sc.pl.dotplot(
        adata,
        var_names=available_s1p,
        groupby="tissue",
        standard_scale="var",
        return_fig=True,
        show=False,
    )
    tissue_path = os.path.join(paths.PLOTS_DIR, "s1p_axis_tissue_dotplot.png")
    tissue_plot.savefig(tissue_path, bbox_inches="tight")
    print(f"[*] Saved S1P tissue DotPlot: {tissue_path}")

print("[*] Generating UMAP feature plots for S1P key genes...")
core_genes = [gene for gene in ["SPHK1", "S1PR1", "S1PR2", "S1PR3"] if gene in adata.var_names]
if core_genes and "X_umap" in adata.obsm:
    fig, axes = plt.subplots(1, len(core_genes), figsize=(5 * len(core_genes), 4))
    if len(core_genes) == 1:
        axes = [axes]
    for ax, gene in zip(axes, core_genes):
        sc.pl.umap(
            adata,
            color=gene,
            cmap="Reds",
            frameon=False,
            show=False,
            ax=ax,
            title=f"Expression: {gene}",
        )
    umap_path = os.path.join(paths.PLOTS_DIR, "umap_s1p_expression.png")
    fig.savefig(umap_path, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"[*] Saved S1P UMAP expression plot: {umap_path}")
else:
    print("[!] Skipping S1P UMAP plots: no selected genes or no X_umap coordinates.")

print("[*] Calculating descriptive mean S1P expression per cell type...")
labels = adata.obs["cell_type"].astype("string")
if labels.isna().any():
    raise ValueError("cell_type contains missing values; resolve annotations before profiling.")
groups = sorted(labels.unique().tolist())
gene_indices = [adata.var_names.get_loc(gene) for gene in available_s1p]
expression = adata.X[:, gene_indices]
means = []
for group in groups:
    mask = labels.to_numpy() == group
    means.append(np.asarray(expression[mask, :].mean(axis=0)).ravel())
mean_df = pd.DataFrame(means, index=groups, columns=available_s1p)
mean_path = os.path.join(paths.RESULTS_DIR, "s1p_axis_mean_expression.csv")
mean_df.to_csv(mean_path)
print(f"[*] Saved descriptive S1P mean expression table: {mean_path}")
print("[*] Step 10 completed successfully.")
