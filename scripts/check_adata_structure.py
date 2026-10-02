import os
import scanpy as sc
import anndata as ad

RESULTS_DIR = "/home/digital_home/targeted_researcher_discovery/pow_yeatman_colorectal/results"
clustered_fp = os.path.join(RESULTS_DIR, "GSE178318_CRC_LM_clustered.h5ad")
processed_fp = os.path.join(RESULTS_DIR, "GSE178318_CRC_LM_processed.h5ad")
qc_fp = os.path.join(RESULTS_DIR, "GSE178318_CRC_LM_qc.h5ad")

print("=== Checking Available Files ===")
for name, fp in [("clustered", clustered_fp), ("processed", processed_fp), ("qc", qc_fp)]:
    if os.path.exists(fp):
        size_gb = os.path.getsize(fp) / (1024**3)
        print(f"[File] {name}: {fp} ({size_gb:.2f} GB)")
        try:
            a = sc.read_h5ad(fp, backed='r')
            print(f"       shape: {a.shape}")
            print(f"       X is None: {a.X is None}")
            print(f"       raw is None: {a.raw is None}")
            print(f"       layers: {list(a.layers.keys())}")
            print(f"       obs columns: {list(a.obs.columns)}")
            print(f"       obsm keys: {list(a.obsm.keys())}")
        except Exception as e:
            print(f"       Error reading {name}: {e}")
    else:
        print(f"[File] {name}: NOT FOUND")
