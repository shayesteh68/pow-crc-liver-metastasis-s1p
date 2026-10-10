# Single-Cell Profiling of the S1P Metabolic and Signaling Axis in Colorectal Cancer Liver Metastasis

[![Python](https://img.shields.io/badge/Python-3.10-blue.svg)](https://www.python.org/)
[![Scanpy](https://img.shields.io/badge/Scanpy-1.9+-brightgreen.svg)](https://scanpy.readthedocs.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

This repository provides an end-to-end, reproducible single-cell RNA-sequencing (scRNA-seq) workflow investigating the **Sphingosine-1-Phosphate (S1P) metabolic and signaling axis** across Primary Colorectal Cancer (CRC) and Colorectal Liver Metastasis (CRLM).

Utilizing the public single-cell atlas **GSE178318** (140,281 initial cells), this pipeline performs rigorous quality control, normalization, high-resolution Leiden clustering, lineage-informed annotation, and comprehensive expression mapping of S1P kinases, receptors, and degradation enzymes.

## Key Biological Insights

- **S1P Production Machinery (`SPHK1` vs `SPHK2`):**
  - `SPHK1` is stromal/myeloid-enriched: Dendritic / APC (0.3427), CAFs (0.2705), Macrophages (0.2139).
  - `SPHK2` shows baseline expression in Epithelial / Tumor and Stress-response compartments.
- **Lineage-Restricted S1P Receptors:**
  - `S1PR4` is enriched across T cells, B cells, and Plasma cells.
  - `S1PR5` is an exclusive NK-cell marker (0.3569).
  - `S1PR3` / `S1PR2` dominate the CAF/stromal compartment (ECM remodeling).
- **S1P Rheostat Degradation:**
  - `SGPL1` (S1P lyase) is strongly expressed in Dendritic / APC (0.2176) and Macrophages (0.1371), regulating local S1P gradients.

## Repository Structure

pow-crc-liver-metastasis-s1p/
├── data/                               # Raw datasets (not tracked)
├── results/
│   ├── plots/                          # Publication-ready figures
│   ├── leiden_cluster_markers.csv      # Top markers across clusters
│   └── s1p_axis_mean_expression.csv    # S1P expression per lineage
├── scripts/
│   ├── paths.py                        # Centralized path resolution
│   ├── 01_fetch_data.py                # Dataset download automation
│   ├── 05_qc_crc_lm_safe.py            # Memory-efficient QC and filtering
│   ├── 06_preprocess_normalize.py      # Normalization and HVG selection
│   ├── 07_dimred_clustering.py         # PCA, graph, Leiden clustering
│   ├── 08_marker_profiling.py          # Differential gene expression
│   ├── 09_cell_annotation.py           # Marker-based lineage annotation
│   ├── 10_s1p_axis_profiling.py        # S1P profiling and visualization
│   └── fix_var_names.py                # Gene symbol alignment routine
├── environment.yml                     # Conda environment definition
├── LICENSE                             # MIT License
└── README.md                           # Project documentation

## Installation & Environment Setup

git clone https://github.com/shayesteh68/pow-crc-liver-metastasis-s1p.git
cd pow-crc-liver-metastasis-s1p
conda env create -f environment.yml
conda activate pow_s1p_crc

## Pipeline Execution

1. Preprocessing and quality control:

python scripts/05_qc_crc_lm_safe.py
python scripts/06_preprocess_normalize.py
python scripts/07_dimred_clustering.py

2. Annotation and marker profiling:

python scripts/fix_var_names.py
python scripts/08_marker_profiling.py
python scripts/09_cell_annotation.py

3. S1P axis interrogation:

python scripts/10_s1p_axis_profiling.py

## Citation & Data Source

- **Dataset:** GSE178318 (NCBI Gene Expression Omnibus)
- **Reference:** Single-cell transcriptomic atlas of colorectal cancer and liver metastases.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
