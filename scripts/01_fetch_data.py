from paths import DATA_DIR

import gzip
import os
import urllib.request

BASE_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178318/suppl/"
FILES = [
    "GSE178318_matrix.mtx.gz",
    "GSE178318_barcodes.tsv.gz",
    "GSE178318_genes.tsv.gz",
]


def validate_gzip(path):
    with gzip.open(path, "rb") as compressed:
        while compressed.read(1024 * 1024):
            pass


def download_one(filename):
    destination = os.path.join(DATA_DIR, filename)
    if os.path.exists(destination) and os.path.getsize(destination) > 0:
        try:
            validate_gzip(destination)
        except (OSError, EOFError) as error:
            raise RuntimeError(f"Existing file is not a valid gzip; inspect it before replacing: {destination}") from error
        print(f"[=] Valid local input exists: {filename}")
        return

    temporary = destination + ".part"
    url = BASE_URL + filename
    print(f"[*] Downloading {filename}...")
    try:
        request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(request, timeout=120) as response, open(temporary, "wb") as output:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                output.write(chunk)
        validate_gzip(temporary)
        os.replace(temporary, destination)
        size_mb = os.path.getsize(destination) / (1024 * 1024)
        print(f"[+] {filename} -> {size_mb:.2f} MB")
    finally:
        if os.path.exists(temporary):
            os.remove(temporary)


for filename in FILES:
    download_one(filename)
