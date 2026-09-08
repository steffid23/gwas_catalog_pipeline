#!/usr/bin/env python3
"""
bin/download_data.py
--------------------
Programmatically downloads datasets from NHGRI-EBI GWAS Catalog.
"""

import os
import sys
import argparse
import logging
import requests

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

ASSOC_URL = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations-full.zip"
ANCESTRY_URL = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-ancestry.tsv"

def download_file(url: str, dest_path: str, chunk_size: int = 1024 * 1024) -> None:
    if os.path.exists(dest_path) and os.path.getsize(dest_path) > 0:
        logger.info(f"File already exists: {dest_path}")
        return

    logger.info(f"Downloading {url} -> {dest_path}...")
    response = requests.get(url, stream=True, timeout=60)
    response.raise_for_status()

    with open(dest_path, 'wb') as f:
        for chunk in response.iter_content(chunk_size=chunk_size):
            if chunk:
                f.write(chunk)
    logger.info(f"Downloaded: {dest_path}")

def main():
    parser = argparse.ArgumentParser(description="Download GWAS Catalog datasets.")
    parser.add_argument("--data_dir", type=str, default="data", help="Directory to store datasets")
    args = parser.parse_args()

    os.makedirs(args.data_dir, exist_ok=True)
    assoc_path = os.path.join(args.data_dir, "gwas_catalog_associations.zip")
    ancestry_path = os.path.join(args.data_dir, "gwas_catalog_ancestry.tsv")

    download_file(ASSOC_URL, assoc_path)
    download_file(ANCESTRY_URL, ancestry_path)

if __name__ == "__main__":
    main()
