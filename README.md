# S1P axis profiling in colorectal liver metastasis

This repository contains a Scanpy workflow for the GSE178318 single-cell atlas. The scripts download GEO inputs, create AnnData objects, run cell-level QC and clustering, annotate cell populations, and summarize S1P-axis expression.

## Environment

Create and activate the environment with:

- conda env create -f environment.yml
- conda activate pow_s1p_crc

## Reproducible workflow

Run from the repository root, in this order:

1. python scripts/01_fetch_data.py
2. python scripts/01_fetch_metadata.py
3. python scripts/02_inspect_data.py
4. python scripts/03_summarize_samples.py
5. python scripts/04_build_anndata_and_split.py
6. python scripts/05_qc_crc_lm_safe.py
7. python scripts/06_preprocess_normalize.py
8. python scripts/07_dimred_clustering.py
9. python scripts/08_marker_profiling.py
10. python scripts/09_cell_annotation.py
11. python scripts/10_s1p_axis_profiling.py
12. python scripts/11_s1p_by_celltype_tissue.py
13. python scripts/12_evaluate_s1p_findings.py
14. python scripts/13_plot_s1p_insights.py

The QC script uses Scrublet per patient. If Scrublet fails or cannot provide a threshold, the script stops. The optional heuristic based on the highest 5% of gene counts is available only when explicitly requested with --allow-heuristic-doublet-fallback; fallback calls are recorded in the output object.

The S1P group comparison reports per-cell mean-expression differences by cell type and tissue. It does not compute patient-level pseudobulk statistics or inferential p-values. Treat these values as descriptive; cells from one patient are not independent biological replicates.

The cluster-to-cell-type table in script 09 is a manual annotation. Check the marker profiles and the current clustering before treating those labels as biological conclusions.

The current `results/leiden_cluster_markers.csv` snapshot contains clusters 0–19 and supports the broad lineages used in script 09. The committed `results/cluster_marker_summary.csv` contains only clusters 0–17 and conflicts with that marker table, so it is a stale snapshot from a different clustering run; do not interpret the two files together. Script 09 now requires an exact cluster-ID match before applying the manual map.

Scripts 05_qc_crc_lm.py, 07_recluster_clean.py, 08_fix_and_save.py, 09_inspect_all_clusters.py, and 10_annotate_celltypes.py are exploratory/legacy alternatives and are not part of the documented pipeline.

The MIT LICENSE applies to repository code and does not change the terms for GEO data. The GSE178318 GEO record directs users to contact the principal investigator for data-use terms. The committed output tables are snapshots; regenerate them from inputs before interpreting them as current results.
