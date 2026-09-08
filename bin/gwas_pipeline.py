#!/usr/bin/env python3
"""
bin/gwas_pipeline.py
--------------------
CLI module for executing the GWAS Catalog ETL, population classification,
low-frequency filtering, and gene-disease mapping pipeline in Nextflow.
"""

import os
import sys
import json
import zipfile
import argparse
import logging
import pandas as pd
import numpy as np

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

def process_gwas(assoc_zip_path: str, ancestry_path: str, results_dir: str):
    os.makedirs(results_dir, exist_ok=True)

    logger.info(f"Loading association zip: {assoc_zip_path}")
    with zipfile.ZipFile(assoc_zip_path) as z:
        tsv_name = z.namelist()[0]
        df_assoc = pd.read_csv(z.open(tsv_name), sep='\t', low_memory=False)

    logger.info(f"Loading ancestry TSV: {ancestry_path}")
    df_ancestry = pd.read_csv(ancestry_path, sep='\t', low_memory=False)

    # 1. Parse Risk Allele Frequency (RAF)
    logger.info("Parsing Risk Allele Frequency (RAF)...")
    df_assoc['RAF_num'] = pd.to_numeric(df_assoc['RISK ALLELE FREQUENCY'].astype(str).str.strip(), errors='coerce')
    df_assoc['is_target_raf'] = (df_assoc['RAF_num'] >= 0.005) & (df_assoc['RAF_num'] <= 0.020)

    # 2. Population Mapping (European vs South Asian)
    logger.info("Mapping European (EUR) and South Asian (SAS) populations...")
    pmid_anc_map = df_ancestry.groupby('PUBMEDID')['BROAD ANCESTRAL CATEGORY'].apply(
        lambda s: set([c.strip() for cat in s.dropna() for c in str(cat).split(',')])
    ).to_dict()

    def check_eur(row):
        cats = pmid_anc_map.get(row['PUBMEDID'], set())
        if any('European' in c for c in cats): return True
        return 'European' in str(row['INITIAL SAMPLE SIZE']) or 'European' in str(row['REPLICATION SAMPLE SIZE'])

    def check_sas(row):
        cats = pmid_anc_map.get(row['PUBMEDID'], set())
        if any('South Asian' in c for c in cats): return True
        return 'South Asian' in str(row['INITIAL SAMPLE SIZE']) or 'South Asian' in str(row['REPLICATION SAMPLE SIZE'])

    df_assoc['is_EUR'] = df_assoc.apply(check_eur, axis=1)
    df_assoc['is_SAS'] = df_assoc.apply(check_sas, axis=1)

    # 3. Analysis 1: Low-Frequency Variant Cohort
    logger.info("Executing Analysis 1: Filtering low-frequency variants (0.5% - 2.0% RAF)...")
    df_low_freq = df_assoc[df_assoc['is_target_raf'] & (df_assoc['is_EUR'] | df_assoc['is_SAS'])].copy()
    df_low_freq['gene_clean'] = df_low_freq['MAPPED_GENE'].fillna(df_low_freq['REPORTED GENE(S)'])

    def parse_genes(gene_series):
        gene_set = set()
        for g_str in gene_series.dropna():
            for g in str(g_str).replace(' - ', ',').replace(';', ',').split(','):
                g_clean = g.strip()
                if g_clean and g_clean != 'NR':
                    gene_set.add(g_clean)
        return gene_set

    low_freq_eur_genes = parse_genes(df_low_freq[df_low_freq['is_EUR']]['gene_clean'])
    low_freq_sas_genes = parse_genes(df_low_freq[df_low_freq['is_SAS']]['gene_clean'])
    low_freq_shared_genes = low_freq_eur_genes.intersection(low_freq_sas_genes)
    low_freq_eur_only_genes = low_freq_eur_genes - low_freq_sas_genes
    low_freq_sas_only_genes = low_freq_sas_genes - low_freq_eur_genes

    low_freq_out = os.path.join(results_dir, "low_freq_variants_eur_sas.csv")
    cols_to_export = [
        'DATE', 'PUBMEDID', 'FIRST AUTHOR', 'STUDY', 'DISEASE/TRAIT',
        'CHR_ID', 'CHR_POS', 'SNPS', 'STRONGEST SNP-RISK ALLELE',
        'REPORTED GENE(S)', 'MAPPED_GENE', 'RISK ALLELE FREQUENCY',
        'RAF_num', 'is_EUR', 'is_SAS'
    ]
    df_low_freq[cols_to_export].to_csv(low_freq_out, index=False)

    # 4. Analysis 2: Full Catalog Gene-Disease Mapping
    logger.info("Executing Analysis 2: Mapping genes, traits, and populations...")
    df_pop = df_assoc[df_assoc['is_EUR'] | df_assoc['is_SAS']].copy()
    df_pop['gene_clean'] = df_pop['MAPPED_GENE'].fillna(df_pop['REPORTED GENE(S)'])

    gene_map, trait_gene_map = {}, {}
    for idx, row in df_pop.iterrows():
        g_str, trait = row['gene_clean'], row['DISEASE/TRAIT']
        if pd.isna(g_str) or str(g_str).strip() == 'NR': continue
        genes = [g.strip() for g in str(g_str).replace(' - ', ',').replace(';', ',').split(',') if g.strip() and g.strip() != 'NR']
        is_e, is_s = row['is_EUR'], row['is_SAS']

        for g in genes:
            if g not in gene_map: gene_map[g] = {'EUR': False, 'SAS': False, 'traits': set(), 'snps': set()}
            if is_e: gene_map[g]['EUR'] = True
            if is_s: gene_map[g]['SAS'] = True
            if pd.notna(trait): gene_map[g]['traits'].add(trait)
            if pd.notna(row['SNPS']): gene_map[g]['snps'].add(row['SNPS'])

            if pd.notna(trait):
                if trait not in trait_gene_map: trait_gene_map[trait] = {'EUR_genes': set(), 'SAS_genes': set()}
                if is_e: trait_gene_map[trait]['EUR_genes'].add(g)
                if is_s: trait_gene_map[trait]['SAS_genes'].add(g)

    gene_rows = []
    for g, info in gene_map.items():
        anc_cat = "Both (Shared)" if (info['EUR'] and info['SAS']) else ("European-Only" if info['EUR'] else "South Asian-Only")
        gene_rows.append({
            "Gene": g, "Ancestry_Category": anc_cat,
            "EUR_Associated": info['EUR'], "SAS_Associated": info['SAS'],
            "Associated_Traits_Count": len(info['traits']),
            "Associated_Traits_List": "; ".join(sorted(list(info['traits'])[:10])),
            "Associated_Variants_Count": len(info['snps'])
        })

    df_gene_mapping = pd.DataFrame(gene_rows)
    gene_map_out = os.path.join(results_dir, "gene_population_trait_mapping.csv")
    df_gene_mapping.to_csv(gene_map_out, index=False)

    trait_rows = []
    for trait, info in trait_gene_map.items():
        eur_g, sas_g = info['EUR_genes'], info['SAS_genes']
        trait_rows.append({
            "Disease/Trait": trait, "Total Mapped Genes": len(eur_g.union(sas_g)),
            "EUR-Only Genes": len(eur_g - sas_g), "Shared Genes": len(eur_g.intersection(sas_g)),
            "SAS-Only Genes": len(sas_g - eur_g)
        })

    df_trait_summary = pd.DataFrame(trait_rows).sort_values(by="Total Mapped Genes", ascending=False)
    trait_summary_out = os.path.join(results_dir, "population_trait_summary.csv")
    df_trait_summary.to_csv(trait_summary_out, index=False)

    metrics = {
        "dataset_statistics": {
            "total_gwas_catalog_records": int(len(df_assoc)),
            "records_with_numeric_raf": int(df_assoc['RAF_num'].notna().sum()),
            "records_with_raf_0.5_to_2.0_pct": int(df_assoc['is_target_raf'].sum()),
            "total_european_records": int(df_assoc['is_EUR'].sum()),
            "total_south_asian_records": int(df_assoc['is_SAS'].sum()),
            "total_multi_ancestry_both_records": int((df_assoc['is_EUR'] & df_assoc['is_SAS']).sum())
        },
        "analysis_1_low_frequency_cohort": {
            "filtered_records_count": int(len(df_low_freq)),
            "filtered_eur_records": int(df_low_freq['is_EUR'].sum()),
            "filtered_sas_records": int(df_low_freq['is_SAS'].sum()),
            "filtered_both_records": int((df_low_freq['is_EUR'] & df_low_freq['is_SAS']).sum()),
            "unique_variants_count": int(df_low_freq['SNPS'].nunique()),
            "unique_diseases_count": int(df_low_freq['DISEASE/TRAIT'].nunique()),
            "total_unique_eur_genes": int(len(low_freq_eur_genes)),
            "total_unique_sas_genes": int(len(low_freq_sas_genes)),
            "low_freq_shared_genes": int(len(low_freq_shared_genes)),
            "low_freq_eur_only_genes": int(len(low_freq_eur_only_genes)),
            "low_freq_sas_only_genes": int(len(low_freq_sas_only_genes))
        },
        "analysis_2_full_catalog_gene_disease_mapping": {
            "total_unique_mapped_genes": int(len(df_gene_mapping)),
            "total_unique_diseases_traits": int(len(df_trait_summary)),
            "full_catalog_eur_genes": int((df_gene_mapping['EUR_Associated']).sum()),
            "full_catalog_sas_genes": int((df_gene_mapping['SAS_Associated']).sum()),
            "full_catalog_shared_genes": int((df_gene_mapping['Ancestry_Category'] == 'Both (Shared)').sum()),
            "full_catalog_eur_only_genes": int((df_gene_mapping['Ancestry_Category'] == 'European-Only').sum()),
            "full_catalog_sas_only_genes": int((df_gene_mapping['Ancestry_Category'] == 'South Asian-Only').sum())
        }
    }

    metrics_out = os.path.join(results_dir, "pipeline_metrics_summary.json")
    with open(metrics_out, 'w') as f:
        json.dump(metrics, f, indent=2)
    logger.info(f"Pipeline complete. Outputs saved in {results_dir}")

def main():
    parser = argparse.ArgumentParser(description="Run GWAS Pipeline ETL")
    parser.add_argument("--assoc_zip", type=str, required=True, help="Path to associations ZIP")
    parser.add_argument("--ancestry_tsv", type=str, required=True, help="Path to ancestry TSV")
    parser.add_argument("--results_dir", type=str, default="results", help="Directory for output tables")
    args = parser.parse_args()

    process_gwas(args.assoc_zip, args.ancestry_tsv, args.results_dir)

if __name__ == "__main__":
    main()
