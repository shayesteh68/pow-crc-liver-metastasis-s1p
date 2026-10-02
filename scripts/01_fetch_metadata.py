import os
import urllib.request
import re

DATA_DIR = "/home/digital_home/targeted_researcher_discovery/pow_yeatman_colorectal/data"
os.makedirs(DATA_DIR, exist_ok=True)

# GSE178318 metadata url
geo_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178318/matrix/GSE178318_series_matrix.txt.gz"
target_file = os.path.join(DATA_DIR, "GSE178318_series_matrix.txt.gz")

print(f"[*] Downloading series matrix from: {geo_url}")
try:
    headers = {'User-Agent': 'Mozilla/5.0'}
    req = urllib.request.Request(geo_url, headers=headers)
    with urllib.request.urlopen(req, timeout=60) as response, open(target_file, 'wb') as out_file:
        out_file.write(response.read())
    print(f"[+] Download complete: {target_file}")
    file_size_mb = os.path.getsize(target_file) / (1024 * 1024)
    print(f"[+] File size: {file_size_mb:.2f} MB")
except Exception as e:
    print(f"[-] Download error: {e}")

