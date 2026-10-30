from paths import DATA_DIR, RESULTS_DIR, FIG_DIR, PLOTS_DIR

import gzip
import os
from collections import Counter


genes_file = os.path.join(DATA_DIR, "GSE178318_genes.tsv.gz")
barcodes_file = os.path.join(DATA_DIR, "GSE178318_barcodes.tsv.gz")

print("=== 1. Checking S1P Axis Target Genes ===")
target_genes = ["SPHK1", "SPHK2", "S1PR1", "S1PR2", "S1PR3", "S1PR4", "S1PR5", "SPNS2", "CD68", "CD163", "ACTA2", "COL1A1", "EPCAM", "CD3D"]
found_genes = {}

with gzip.open(genes_file, 'rt') as f:
    for idx, line in enumerate(f):
        parts = line.strip().split('\t')
        gene_id = parts[0]
        gene_symbol = parts[1] if len(parts) > 1 else parts[0]
        for tg in target_genes:
            if tg.upper() in [gene_id.upper(), gene_symbol.upper()]:
                found_genes[tg] = (idx + 1, gene_id, gene_symbol)

for tg in target_genes:
    if tg in found_genes:
        print(f"[FOUND] {tg:<8} -> Index: {found_genes[tg][0]}, ID: {found_genes[tg][1]}, Symbol: {found_genes[tg][2]}")
    else:
        print(f"[MISSING] {tg}")

print("\n=== 2. Checking Barcode Structure & Sample Tags (First 15 Barcodes) ===")
with gzip.open(barcodes_file, 'rt') as f:
    sample_tags = Counter()
    malformed = 0
    for i, line in enumerate(f):
        bc = line.strip()
        if i < 15:
            print(f"  Barcode {i+1}: {bc}")
        # GEO uses <cell-barcode>_<patient>_<tissue>; the first field is not a sample ID.
        parts = bc.rsplit("_", 2)
        if len(parts) == 3 and parts[1] and parts[2]:
            sample_tags[(parts[1], parts[2])] += 1
        else:
            malformed += 1

if sample_tags:
    patients = sorted({patient for patient, _ in sample_tags})
    tissues = sorted({tissue for _, tissue in sample_tags})
    print(f"\n[+] Parsed {len(patients)} patient IDs and {len(tissues)} tissue tags:")
    for (patient, tissue), count in sorted(sample_tags.items()):
        print(f"  {patient} / {tissue}: {count:,} barcodes")
    print(f"[+] Malformed/unlabeled barcodes: {malformed:,}")
else:
    print("\n[-] No patient/tissue tags found in barcode names.")

