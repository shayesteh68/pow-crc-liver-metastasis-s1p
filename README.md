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

The QC script runs Scrublet independently per patient/tissue library and records scores, thresholds, call candidates, score quantiles, and observed-versus-simulated histograms. Automatic Scrublet thresholds require visual review; this workflow does not remove cells based on them by default. To apply manually reviewed thresholds, copy `doublet_thresholds.template.tsv`, replace every `NA` with a reviewed numeric threshold for each library after reviewing its histogram, and pass the completed file with `--thresholds-tsv`. The optional gene-count heuristic is available only with `--allow-heuristic-doublet-fallback`; heuristic candidates are also not removed unless reviewed thresholds are supplied.

The S1P group comparison reports per-cell mean-expression differences by cell type and tissue. It does not compute patient-level pseudobulk statistics or inferential p-values. Treat these values as descriptive; cells from one patient are not independent biological replicates.

The cluster-to-cell-type table in script 09 is a provisional manual annotation reviewed against marker profiles. Script 07 fixes the random states and `environment.yml` pins the analysis packages. Script 09 checks that the observed Leiden IDs match the reviewed mapping and stops if they differ, rather than silently assigning labels to changed clusters. Re-review the marker table whenever inputs, preprocessing, or package versions change.

Committed result tables and figures are historical snapshots; some come from different clustering runs and must not be interpreted together as current results. Regenerate the complete workflow from the inputs before using them for biological conclusions.

Scripts 05_qc_crc_lm.py, 07_recluster_clean.py, 08_fix_and_save.py, 09_inspect_all_clusters.py, and 10_annotate_celltypes.py are exploratory/legacy alternatives and are not part of the documented pipeline.

The MIT LICENSE applies to repository code and does not change the terms for GEO data. The GSE178318 GEO record directs users to contact the principal investigator for data-use terms. The committed output tables are snapshots; regenerate them from inputs before interpreting them as current results.
