from paths import DATA_DIR, RESULTS_DIR, FIG_DIR, PLOTS_DIR

import os
import subprocess


base_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178318/suppl/"
files = [
    "GSE178318_matrix.mtx.gz",
    "GSE178318_barcodes.tsv.gz",
    "GSE178318_genes.tsv.gz"
]

for f in files:
    target_path = os.path.join(DATA_DIR, f)
    file_url = base_url + f
    print(f"[*] Downloading {f}...")
    cmd = f"curl -L -s -o '{target_path}' '{file_url}'"
    subprocess.run(cmd, shell=True, check=True)
    size_mb = os.path.getsize(target_path) / (1024 * 1024)
    print(f"[+] {f} -> {size_mb:.2f} MB")

