import sys
import os
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath(os.path.dirname(__file__)))
import paths

sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=150, facecolor="white")

annotated_file = os.path.join(paths.RESULTS_DIR, "GSE178318_CRC_LM_annotated.h5ad")
print(f"[*] Loading annotated dataset: {annotated_file}")
adata = sc.read_h5ad(annotated_file)

s1p_genes = ["SPHK1", "SPHK2", "S1PR1", "S1PR2", "S1PR3", "S1PR4", "S1PR5", "SGPP1", "SGPL1"]
available_s1p = [g for g in s1p_genes if g in adata.var_names]
print(f"[*] S1P axis genes available: {available_s1p}")

os.makedirs(paths.PLOTS_DIR, exist_ok=True)

print("[*] Generating S1P axis DotPlot by cell type...")
dp = sc.pl.dotplot(
    adata,
    var_names=available_s1p,
    groupby="cell_type",
    standard_scale="var",
    return_fig=True,
    show=False
)
s1p_dotplot_path = os.path.join(paths.PLOTS_DIR, "s1p_axis_celltypes_dotplot.png")
dp.savefig(s1p_dotplot_path, bbox_inches="tight")
print(f"[✓] Saved S1P cell types DotPlot: {s1p_dotplot_path}")

if "tissue" in adata.obs.columns:
    print("[*] Generating S1P axis DotPlot by tissue (CRC vs LM)...")
    dp_tissue = sc.pl.dotplot(
        adata,
        var_names=available_s1p,
        groupby="tissue",
        standard_scale="var",
        return_fig=True,
        show=False
    )
    s1p_tissue_path = os.path.join(paths.PLOTS_DIR, "s1p_axis_tissue_dotplot.png")
    dp_tissue.savefig(s1p_tissue_path, bbox_inches="tight")
    print(f"[✓] Saved S1P tissue DotPlot: {s1p_tissue_path}")

print("[*] Generating UMAP feature plots for S1P key genes...")
s1p_core = [g for g in ["SPHK1", "S1PR1", "S1PR2", "S1PR3"] if g in adata.var_names]
if s1p_core:
    fig, axes = plt.subplots(1, len(s1p_core), figsize=(5 * len(s1p_core), 4))
    if len(s1p_core) == 1:
        axes = [axes]
    for ax, g in zip(axes, s1p_core):
        sc.pl.umap(adata, color=g, cmap="Reds", frameon=False, show=False, ax=ax, title=f"Expression: {g}")
    umap_s1p_path = os.path.join(paths.PLOTS_DIR, "umap_s1p_expression.png")
    fig.savefig(umap_s1p_path, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"[✓] Saved S1P UMAP expression plot: {umap_s1p_path}")

print("[*] Calculating mean S1P gene expression per cell type...")
mean_df = pd.DataFrame(index=adata.obs["cell_type"].cat.categories)
for gene in available_s1p:
    idx = adata.var_names.get_loc(gene)
    vals = adata.X[:, idx]
    if hasattr(vals, "toarray"):
        vals = vals.toarray().flatten()
    else:
        vals = vals.flatten()
    mean_df[gene] = adata.obs.groupby("cell_type", observed=False).apply(lambda df: vals[df.index.to_numpy() if isinstance(df.index[0], int) else [adata.obs_names.get_loc(x) for x in df.index]].mean())

mean_csv_path = os.path.join(paths.RESULTS_DIR, "s1p_axis_mean_expression.csv")
mean_df.to_csv(mean_csv_path)
print(f"[✓] Saved S1P mean expression table: {mean_csv_path}")

print("[✓] Step 10 completed successfully!")
