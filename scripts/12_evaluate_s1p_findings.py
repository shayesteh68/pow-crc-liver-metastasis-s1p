from paths import DATA_DIR, RESULTS_DIR, FIG_DIR, PLOTS_DIR

import os
import pandas as pd
import numpy as np

CSV_FP = os.path.join(RESULTS_DIR, "s1p_mean_by_celltype_tissue.csv")
COUNT_FP = os.path.join(RESULTS_DIR, "s1p_group_cellcounts.csv")

df = pd.read_csv(CSV_FP, index_col=0)
counts = pd.read_csv(COUNT_FP, index_col=0)

print("\n" + "="*85)
print("S1P PATHWAY QUANTITATIVE PROFILING (CRC PRIMARY vs LIVER METASTASIS [LM])")
print("="*85)

cell_types = sorted(list(set([str(idx).split(" | ")[0] for idx in df.index if " | " in str(idx)])))

detailed_records = []

for ct in cell_types:
    crc_idx = f"{ct} | CRC"
    lm_idx = f"{ct} | LM"
    
    if crc_idx in df.index and lm_idx in df.index:
        p_row = df.loc[crc_idx]
        m_row = df.loc[lm_idx]
        
        n_p = counts.loc[crc_idx].values[0] if crc_idx in counts.index else 0
        n_m = counts.loc[lm_idx].values[0] if lm_idx in counts.index else 0
        
        diff = m_row - p_row
        
        # Collect for table
        for gene in df.columns:
            detailed_records.append({
                "cell_type": ct,
                "gene": gene,
                "mean_CRC": round(p_row[gene], 4),
                "mean_LM": round(m_row[gene], 4),
                "diff_LM_minus_CRC": round(diff[gene], 4)
            })
            
        top_up = diff.sort_values(ascending=False).head(4)
        top_down = diff.sort_values(ascending=True).head(3)
        
        up_str = ", ".join([f"{g} (+{v:.4f})" for g, v in top_up.items() if v > 0.001])
        down_str = ", ".join([f"{g} ({v:.4f})" for g, v in top_down.items() if v < -0.001])
        
        print(f"\n[Cell Type: {ct:<15}] (CRC: {n_p:>5} cells | LM: {n_m:>5} cells)")
        print(f"  ▲ Upregulated in Liver Met: {up_str if up_str else 'None / Minimal'}")
        print(f"  ▼ Downregulated in Liver Met: {down_str if down_str else 'None / Minimal'}")

out_diff_csv = os.path.join(RESULTS_DIR, "s1p_detailed_differential_summary.csv")
pd.DataFrame(detailed_records).to_csv(out_diff_csv, index=False)
print("\n" + "="*85)
print(f"[✓] Complete differential table saved to: {out_diff_csv}")
print("="*85 + "\n")
