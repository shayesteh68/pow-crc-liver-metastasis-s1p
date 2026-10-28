from paths import DATA_DIR, RESULTS_DIR, FIG_DIR, PLOTS_DIR

import os
import pandas as pd

CSV_FP = os.path.join(RESULTS_DIR, "s1p_mean_by_celltype_tissue.csv")
COUNT_FP = os.path.join(RESULTS_DIR, "s1p_group_cellcounts.csv")

df = pd.read_csv(CSV_FP, index_col=0)
counts = pd.read_csv(COUNT_FP, index_col=0)
print("Descriptive cell-level mean-expression differences only; no patient-level test or inferential p-values are computed.")

cell_types = sorted({str(idx).split(" | ")[0] for idx in df.index if " | " in str(idx)})
records = []

for cell_type in cell_types:
    crc_idx = f"{cell_type} | CRC"
    lm_idx = f"{cell_type} | LM"
    if crc_idx not in df.index or lm_idx not in df.index:
        continue

    crc_row = df.loc[crc_idx]
    lm_row = df.loc[lm_idx]
    n_crc = int(counts.loc[crc_idx].iloc[0]) if crc_idx in counts.index else 0
    n_lm = int(counts.loc[lm_idx].iloc[0]) if lm_idx in counts.index else 0
    difference = lm_row - crc_row

    for gene in df.columns:
        records.append({
            "cell_type": cell_type,
            "gene": gene,
            "mean_CRC": round(float(crc_row[gene]), 4),
            "mean_LM": round(float(lm_row[gene]), 4),
            "diff_LM_minus_CRC": round(float(difference[gene]), 4),
            "n_cells_CRC": n_crc,
            "n_cells_LM": n_lm,
            "interpretation": "descriptive_cell_mean_no_inference",
        })

    up = difference.sort_values(ascending=False).head(4)
    down = difference.sort_values(ascending=True).head(3)
    up_text = ", ".join(f"{gene} (+{value:.4f})" for gene, value in up.items() if value > 0.001)
    down_text = ", ".join(f"{gene} ({value:.4f})" for gene, value in down.items() if value < -0.001)
    print(f"\n[Cell type: {cell_type}] CRC={n_crc} cells | LM={n_lm} cells")
    print(f"  Higher in LM (descriptive): {up_text or 'None / minimal'}")
    print(f"  Lower in LM (descriptive): {down_text or 'None / minimal'}")

out_fp = os.path.join(RESULTS_DIR, "s1p_detailed_differential_summary.csv")
pd.DataFrame(records).to_csv(out_fp, index=False)
print(f"\n[+] Descriptive table saved to: {out_fp}")
