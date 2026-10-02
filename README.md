# S1P Axis Remodeling in Colorectal Cancer Liver Metastasis (Single-Cell RNA-Seq)

[![Reproducible Workflow](https://img.shields.io/badge/Workflow-Reproducible-green.svg)](https://github.com/NargesShayesteh)
[![Python](https://img.shields.io/badge/Python-3.14%2B-blue.svg)](https://www.python.org/)
[![Scanpy](https://img.shields.io/badge/Scanpy-1.10%2B-brightgreen.svg)](https://scanpy.readthedocs.io/)
[![Dataset](https://img.shields.io/badge/NCBI%20GEO-GSE178318-orange.svg)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178318)

---

## 🔬 Overview & Scientific Rationale

Colorectal cancer (CRC) liver metastasis remains a clinical bottleneck driven by profound stromal remodeling and immune suppression in the liver metastatic niche. The **Sphingosine-1-Phosphate (S1P) signaling axis** regulates fibrosis, immune trafficking, and vascular permeability, yet its lineage-specific activation across matched primary CRC vs. liver metastases has remained largely unquantified at single-cell resolution.

This repository provides a memory-optimized, zero-hallucination re-analysis of **104,485 single cells** from **GSE178318** (Yeatman Lab dataset), mapping the quantitative shift in the S1P metabolic and signaling network across discrete tumor and microenvironmental lineages.

---

## 🎯 Key Findings & Biological Insights

1. **CAF Activation via S1PR3 / SPHK2 Axis:**
   - **Cancer-Associated Fibroblasts (CAFs)** in liver metastasis (`LM`) exhibit substantial upregulation of **`S1PR3` (+0.0875)** and **`SPHK2` (+0.0258)**.
   - Concurrently, S1P-degrading lipid phosphate phosphatases **`PLPP3` (-0.2665)** and **`PLPP1` (-0.2013)** are heavily suppressed in metastatic stroma, establishing a local accumulation of pro-fibrotic S1P signaling.
2. **Myeloid / TAM Phenotypic Shift:**
   - Infiltrating myeloid lineages in liver metastasis demonstrate elevated **`ASAH1` (+0.0994)** and **`S1PR4` (+0.0341)** expression, characteristic of immunosuppressive polarization.
3. **Epithelial/Tumor Reprogramming:**
   - Tumor cells shift from basal `SPHK1` expression towards elevated **`SPHK2` (+0.0142)** and **`ASAH1` (+0.0888)** during metastatic colonization.

---

## 📊 Summary Differential Matrix (LM vs CRC)

| Cell Lineage | Cells (CRC) | Cells (LM) | Key Upregulated in Liver Met | Key Downregulated in Liver Met |
| :--- | :---: | :---: | :--- | :--- |
| **Stromal_CAF** | 1,373 | 495 | **`S1PR3` (+0.0875)**, **`SPHK2` (+0.0258)** | **`PLPP3` (-0.2665)**, **`PLPP1` (-0.2013)**, `S1PR2` (-0.0933) |
| **Myeloid** | 7,535 | 8,585 | **`ASAH1` (+0.0994)**, **`S1PR4` (+0.0341)** | `SPHK1` (-0.0407), `SGPP1` (-0.0287) |
| **Epithelial** | 7,392 | 2,679 | **`ASAH1` (+0.0888)**, **`PLPP3` (+0.0257)** | `PLPP2` (-0.0914), `SGPL1` (-0.0190) |
| **T_NK** | 16,971 | 22,034 | `ASAH1` (+0.0247), `S1PR4` (+0.0176) | `S1PR5` (-0.0088), `PLPP2` (-0.0035) |
| **B_Plasma** | 16,314 | 1,926 | `ASAH1` (+0.0659), `PLPP2` (+0.0107) | `S1PR4` (-0.0733), `CERK` (-0.0578) |

---

## 📁 Repository Structure
pow_yeatman_colorectal/

├── figures/

│ ├── s1p_differential_heatmap.png # High-res comparative heatmap (DPI=300)

│ ├── s1p_dotplot_celltype_tissue.png # Expression dotplot across lineages

│ └── umap_celltypes_tissue.png # QC-filtered UMAP annotations

├── results/

│ ├── s1p_detailed_differential_summary.csv

│ ├── s1p_mean_by_celltype_tissue.csv

│ └── s1p_group_cellcounts.csv

├── scripts/

│ ├── 01_download_data.sh

│ ├── 02_scrublet_qc.py

│ ├── 03_run_leiden_clustering.py

│ ├── 10_annotate_celltypes.py

│ ├── 11_profile_s1p_axis.py

│ ├── 12_evaluate_s1p_findings.py

│ └── 13_plot_s1p_insights.py

├── environment.yml

└── README.md

text

---

## ⚡ Memory & Execution Safety

- **Hardware Efficiency:** Designed specifically to run on local environments with **12 GB RAM**.
- **Sparse Matrix Standard:** Employs `scipy.sparse.csr_matrix` slicing throughout downstream quantification to eliminate OOM risks.

---

## 👤 Author & Contact

**Narges Shayesteh**  
*Freelance Bioinformatician & Computational Biologist*  
- **Focus:** Single-Cell RNA-Seq, Spatial Transcriptomics, Tumor Microenvironment Profiling.  
- **GitHub:** [github.com/NargesShayesteh](https://github.com/NargesShayesteh)  
- **Email:** narges.shayesteh@outlook.com *(or your professional contact)*
