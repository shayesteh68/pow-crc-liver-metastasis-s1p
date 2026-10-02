import gzip
import os
from collections import Counter

DATA_DIR = "/home/digital_home/targeted_researcher_discovery/pow_yeatman_colorectal/data"
barcodes_file = os.path.join(DATA_DIR, "GSE178318_barcodes.tsv.gz")

suffix_counter = Counter()
patient_counter = Counter()
total_cells = 0

with gzip.open(barcodes_file, 'rt') as f:
    for line in f:
        total_cells += 1
        parts = line.strip().split('_')
        # Structure is usually: CELLBARCODE_PATIENT_TISSUE or CELLBARCODE_SAMPLE
        if len(parts) >= 3:
            patient = parts[1]
            tissue = "_".join(parts[2:])
            suffix_counter[tissue] += 1
            patient_counter[patient] += 1
        elif len(parts) == 2:
            patient_counter[parts[1]] += 1
            suffix_counter[parts[1]] += 1

print(f"Total Cells in Dataset: {total_cells:,}")
print("\n--- Breakdown by Tissue / Condition ---")
for tissue, count in suffix_counter.most_common():
    print(f"  {tissue:<20}: {count:>7,} cells ({count/total_cells*100:.2f}%)")

print(f"\n--- Total Patients / Donors: {len(patient_counter)} ---")
for pt, count in sorted(patient_counter.items())[:15]:
    print(f"  Patient {pt:<10}: {count:>6,} cells")
if len(patient_counter) > 15:
    print(f"  ... and {len(patient_counter)-15} more patients")

