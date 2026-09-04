"""
download_data.py
----------------
Programmatically downloads the required datasets from the NHGRI-EBI GWAS Catalog:
1. Full GWAS Catalog association TSV (zipped): gwas-catalog-associations-full.zip
2. Study Ancestry TSV: gwas-catalog-ancestry.tsv
"""

import os
import sys
import logging
import requests

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

ASSOC_URL = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations-full.zip"
ANCESTRY_URL = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-ancestry.tsv"

def download_file(url: str, dest_path: str, chunk_size: int = 1024 * 1024) -> None:
    """Download file from url to dest_path with progress logging."""
    if os.path.exists(dest_path) and os.path.getsize(dest_path) > 0:
        logger.info(f"File already exists: {dest_path} ({os.path.getsize(dest_path) / (1024*1024):.2f} MB)")
        return

    logger.info(f"Downloading {url} -> {dest_path}...")
    response = requests.get(url, stream=True, timeout=60)
    response.raise_for_status()

    total_size = int(response.headers.get('content-length', 0))
    downloaded = 0

    with open(dest_path, 'wb') as f:
        for chunk in response.iter_content(chunk_size=chunk_size):
            if chunk:
                f.write(chunk)
                downloaded += len(chunk)
                if total_size > 0:
                    pct = (downloaded / total_size) * 100
                    if downloaded % (5 * 1024 * 1024) < chunk_size:
                        logger.info(f"Progress: {downloaded/(1024*1024):.2f}/{total_size/(1024*1024):.2f} MB ({pct:.1f}%)")

    logger.info(f"Download complete: {dest_path} ({os.path.getsize(dest_path) / (1024*1024):.2f} MB)")

def main(data_dir: str = None) -> tuple[str, str]:
    """Ensure data directory exists and datasets are downloaded."""
    if data_dir is None:
        data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
    os.makedirs(data_dir, exist_ok=True)

    assoc_path = os.path.join(data_dir, "gwas_catalog_associations.zip")
    ancestry_path = os.path.join(data_dir, "gwas_catalog_ancestry.tsv")

    download_file(ASSOC_URL, assoc_path)
    download_file(ANCESTRY_URL, ancestry_path)

    return assoc_path, ancestry_path

if __name__ == "__main__":
    main()
