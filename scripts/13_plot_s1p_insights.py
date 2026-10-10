from paths import DATA_DIR, RESULTS_DIR, FIG_DIR, PLOTS_DIR

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


DIFF_CSV = os.path.join(RESULTS_DIR, "s1p_detailed_differential_summary.csv")
df_diff = pd.read_csv(DIFF_CSV)

pivot_diff = df_diff.pivot(index="cell_type", columns="gene", values="diff_LM_minus_CRC")

plt.figure(figsize=(12, 6), dpi=300)
sns.set_theme(style="white")

cmap = sns.diverging_palette(240, 10, as_cmap=True)

ax = sns.heatmap(
    pivot_diff,
    cmap=cmap,
    center=0,
    annot=True,
    fmt=".3f",
    linewidths=0.75,
    cbar_kws={"label": "Mean Difference (Liver Met - Primary CRC)"}
)

plt.title("S1P Axis Remodeling: Liver Metastasis vs Primary CRC Across Lineages", fontsize=13, weight="bold", pad=15)
plt.xlabel("S1P Pathway Genes", fontsize=11, weight="bold")
plt.ylabel("Cell Lineage", fontsize=11, weight="bold")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()

out_fig = os.path.join(FIG_DIR, "s1p_differential_heatmap.png")
plt.savefig(out_fig, dpi=300)
plt.close()

print(f"[✓] Differential heatmap successfully saved to: {out_fig}")
