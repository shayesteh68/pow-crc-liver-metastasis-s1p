from paths import DATA_DIR

import gzip
import os
import urllib.request

GEO_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178318/matrix/GSE178318_series_matrix.txt.gz"
DESTINATION = os.path.join(DATA_DIR, "GSE178318_series_matrix.txt.gz")


def validate_gzip(path):
    with gzip.open(path, "rb") as compressed:
        while compressed.read(1024 * 1024):
            pass


if os.path.exists(DESTINATION) and os.path.getsize(DESTINATION) > 0:
    validate_gzip(DESTINATION)
    print(f"[=] Valid metadata file already exists: {DESTINATION}")
else:
    temporary = DESTINATION + ".part"
    request = urllib.request.Request(GEO_URL, headers={"User-Agent": "Mozilla/5.0"})
    try:
        print(f"[*] Downloading series matrix from: {GEO_URL}")
        with urllib.request.urlopen(request, timeout=120) as response, open(temporary, "wb") as output:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                output.write(chunk)
        validate_gzip(temporary)
        os.replace(temporary, DESTINATION)
        print(f"[+] Saved valid metadata: {DESTINATION} ({os.path.getsize(DESTINATION) / (1024 * 1024):.2f} MB)")
    finally:
        if os.path.exists(temporary):
            os.remove(temporary)
