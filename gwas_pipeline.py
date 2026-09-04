"""
gwas_pipeline.py
----------------
Core ETL and Analysis engine for the GWAS Catalog Analysis Pipeline.
Performs data loading, column documentation, allele-frequency parsing,
population classification (EUR vs SAS), low-frequency filtering (0.5%-2.0%),
and gene-disease-population mapping.
"""

import os
import json
import zipfile
import logging
import pandas as pd
import numpy as np

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

COLUMN_DOCUMENTATION = {
    "variant_mutation": ["SNPS", "STRONGEST SNP-RISK ALLELE", "CONTEXT", "SNP_ID_CURRENT"],
    "gene": ["MAPPED_GENE", "REPORTED GENE(S)", "SNP_GENE_IDS", "UPSTREAM_GENE_ID", "DOWNSTREAM_GENE_ID"],
    "allele_frequency": ["RISK ALLELE FREQUENCY"],
    "population": ["INITIAL SAMPLE SIZE", "REPLICATION SAMPLE SIZE", "BROAD ANCESTRAL CATEGORY"],
    "disease_trait": ["DISEASE/TRAIT"],
    "chromosome_position": ["CHR_ID", "CHR_POS", "REGION"],
    "study_information": ["PUBMEDID", "FIRST AUTHOR", "DATE", "JOURNAL", "STUDY", "LINK"]
}

def load_and_inspect_data(assoc_zip_path: str, ancestry_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load association zip and ancestry TSV files, logging column statistics."""
    logger.info(f"Loading association zip: {assoc_zip_path}")
    with zipfile.ZipFile(assoc_zip_path) as z:
        tsv_name = z.namelist()[0]
        logger.info(f"Extracting file: {tsv_name}")
        df_assoc = pd.read_csv(z.open(tsv_name), sep='\t', low_memory=False)

    logger.info(f"Loading ancestry TSV: {ancestry_path}")
    df_ancestry = pd.read_csv(ancestry_path, sep='\t', low_memory=False)

    logger.info(f"Loaded Associations shape: {df_assoc.shape}")
    logger.info(f"Loaded Ancestry shape: {df_ancestry.shape}")

    # Inspect columns
    print("\n" + "="*80)
    print("GWAS CATALOG ASSOCIATION DATASET - COLUMN SCHEMATICS & CATEGORIZATION")
    print("="*80)
    for category, cols in COLUMN_DOCUMENTATION.items():
        present = [c for c in cols if c in df_assoc.columns or c in df_ancestry.columns]
        print(f"\n  [{category.upper().replace('_', ' ')}]")
        for c in present:
            dtype = df_assoc[c].dtype if c in df_assoc.columns else df_ancestry[c].dtype
            null_cnt = df_assoc[c].isna().sum() if c in df_assoc.columns else df_ancestry[c].isna().sum()
            total = len(df_assoc) if c in df_assoc.columns else len(df_ancestry)
            pct_valid = ((total - null_cnt) / total) * 100
            print(f"   - {c:<30} | Type: {str(dtype):<10} | Non-Null: {total-null_cnt:,}/{total:,} ({pct_valid:.1f}%)")
    print("="*80 + "\n")

    return df_assoc, df_ancestry

def process_pipeline(df_assoc: pd.DataFrame, df_ancestry: pd.DataFrame, results_dir: str) -> dict:
    """Execute complete analysis pipeline and export structured outputs."""
    os.makedirs(results_dir, exist_ok=True)

    # 1. Parse Risk Allele Frequency (RAF)
    logger.info("Parsing Risk Allele Frequency (RAF)...")
    df_assoc['RAF_num'] = pd.to_numeric(df_assoc['RISK ALLELE FREQUENCY'].astype(str).str.strip(), errors='coerce')
    
    total_records = len(df_assoc)
    valid_raf_count = df_assoc['RAF_num'].notna().sum()
    
    # Low frequency range: 0.5% to 2.0% inclusive (0.005 <= RAF <= 0.020)
    df_assoc['is_target_raf'] = (df_assoc['RAF_num'] >= 0.005) & (df_assoc['RAF_num'] <= 0.020)
    total_target_raf = df_assoc['is_target_raf'].sum()

    # 2. Population Mapping (European vs South Asian)
    logger.info("Mapping European (EUR) and South Asian (SAS) populations...")
    
    # Map PMID -> Broad Ancestral Category set from ancestry dataset
    pmid_anc_map = df_ancestry.groupby('PUBMEDID')['BROAD ANCESTRAL CATEGORY'].apply(
        lambda s: set([c.strip() for cat in s.dropna() for c in str(cat).split(',')])
    ).to_dict()

    def check_eur(row):
        cats = pmid_anc_map.get(row['PUBMEDID'], set())
        if any('European' in c for c in cats):
            return True
        init_str = str(row['INITIAL SAMPLE SIZE'])
        rep_str = str(row['REPLICATION SAMPLE SIZE'])
        return 'European' in init_str or 'European' in rep_str

    def check_sas(row):
        cats = pmid_anc_map.get(row['PUBMEDID'], set())
        if any('South Asian' in c for c in cats):
            return True
        init_str = str(row['INITIAL SAMPLE SIZE'])
        rep_str = str(row['REPLICATION SAMPLE SIZE'])
        return 'South Asian' in init_str or 'South Asian' in rep_str

    df_assoc['is_EUR'] = df_assoc.apply(check_eur, axis=1)
    df_assoc['is_SAS'] = df_assoc.apply(check_sas, axis=1)

    total_eur = df_assoc['is_EUR'].sum()
    total_sas = df_assoc['is_SAS'].sum()
    total_both = (df_assoc['is_EUR'] & df_assoc['is_SAS']).sum()

    # 3. Analysis 1: Low-Frequency Variant Cohort
    logger.info("Executing Analysis 1: Filtering low-frequency variants (0.5% - 2.0% RAF)...")
    df_low_freq = df_assoc[df_assoc['is_target_raf'] & (df_assoc['is_EUR'] | df_assoc['is_SAS'])].copy()

    # Define clean gene column
    df_low_freq['gene_clean'] = df_low_freq['MAPPED_GENE'].fillna(df_low_freq['REPORTED GENE(S)'])
    
    # Helper to parse multi-gene strings
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

    # Export Analysis 1 result table
    low_freq_out = os.path.join(results_dir, "low_freq_variants_eur_sas.csv")
    cols_to_export = [
        'DATE', 'PUBMEDID', 'FIRST AUTHOR', 'STUDY', 'DISEASE/TRAIT',
        'CHR_ID', 'CHR_POS', 'SNPS', 'STRONGEST SNP-RISK ALLELE',
        'REPORTED GENE(S)', 'MAPPED_GENE', 'RISK ALLELE FREQUENCY',
        'RAF_num', 'is_EUR', 'is_SAS'
    ]
    df_low_freq[cols_to_export].to_csv(low_freq_out, index=False)
    logger.info(f"Saved Analysis 1 results to {low_freq_out} ({len(df_low_freq):,} records)")

    # 4. Analysis 2: Unique Gene-Disease-Population Mapping Across Catalog
    logger.info("Executing Analysis 2: Mapping genes, traits, and populations across full catalog...")
    df_pop = df_assoc[df_assoc['is_EUR'] | df_assoc['is_SAS']].copy()
    df_pop['gene_clean'] = df_pop['MAPPED_GENE'].fillna(df_pop['REPORTED GENE(S)'])

    # Build gene -> population & traits dictionary
    gene_map = {}
    trait_gene_map = {} # trait -> {EUR_genes, SAS_genes}

    for idx, row in df_pop.iterrows():
        g_str = row['gene_clean']
        trait = row['DISEASE/TRAIT']
        if pd.isna(g_str) or str(g_str).strip() == 'NR':
            continue
        
        genes = [g.strip() for g in str(g_str).replace(' - ', ',').replace(';', ',').split(',') if g.strip() and g.strip() != 'NR']
        is_e = row['is_EUR']
        is_s = row['is_SAS']

        for g in genes:
            if g not in gene_map:
                gene_map[g] = {'EUR': False, 'SAS': False, 'traits': set(), 'snps': set()}
            if is_e:
                gene_map[g]['EUR'] = True
            if is_s:
                gene_map[g]['SAS'] = True
            if pd.notna(trait):
                gene_map[g]['traits'].add(trait)
            if pd.notna(row['SNPS']):
                gene_map[g]['snps'].add(row['SNPS'])

            # Trait mapping
            if pd.notna(trait):
                if trait not in trait_gene_map:
                    trait_gene_map[trait] = {'EUR_genes': set(), 'SAS_genes': set()}
                if is_e:
                    trait_gene_map[trait]['EUR_genes'].add(g)
                if is_s:
                    trait_gene_map[trait]['SAS_genes'].add(g)

    # Classify genes
    gene_rows = []
    for g, info in gene_map.items():
        anc_cat = "Both (Shared)" if (info['EUR'] and info['SAS']) else ("European-Only" if info['EUR'] else "South Asian-Only")
        gene_rows.append({
            "Gene": g,
            "Ancestry_Category": anc_cat,
            "EUR_Associated": info['EUR'],
            "SAS_Associated": info['SAS'],
            "Associated_Traits_Count": len(info['traits']),
            "Associated_Traits_List": "; ".join(sorted(list(info['traits'])[:10])), # cap for CSV
            "Associated_Variants_Count": len(info['snps'])
        })

    df_gene_mapping = pd.DataFrame(gene_rows)
    gene_map_out = os.path.join(results_dir, "gene_population_trait_mapping.csv")
    df_gene_mapping.to_csv(gene_map_out, index=False)
    logger.info(f"Saved Gene Mapping results to {gene_map_out} ({len(df_gene_mapping):,} unique genes)")

    # Trait summary table
    trait_rows = []
    for trait, info in trait_gene_map.items():
        eur_g = info['EUR_genes']
        sas_g = info['SAS_genes']
        all_g = eur_g.union(sas_g)
        shared_g = eur_g.intersection(sas_g)
        eur_only_g = eur_g - sas_g
        sas_only_g = sas_g - eur_g

        trait_rows.append({
            "Disease/Trait": trait,
            "Total Mapped Genes": len(all_g),
            "EUR-Only Genes": len(eur_only_g),
            "Shared Genes": len(shared_g),
            "SAS-Only Genes": len(sas_only_g),
            "Total EUR Genes": len(eur_g),
            "Total SAS Genes": len(sas_g)
        })

    df_trait_summary = pd.DataFrame(trait_rows).sort_values(by="Total Mapped Genes", ascending=False)
    trait_summary_out = os.path.join(results_dir, "population_trait_summary.csv")
    df_trait_summary.to_csv(trait_summary_out, index=False)
    logger.info(f"Saved Trait Summary results to {trait_summary_out} ({len(df_trait_summary):,} unique traits)")

    # 5. Compile Metrics JSON Summary
    metrics = {
        "dataset_statistics": {
            "total_gwas_catalog_records": int(total_records),
            "records_with_numeric_raf": int(valid_raf_count),
            "records_with_raf_0.5_to_2.0_pct": int(total_target_raf),
            "total_european_records": int(total_eur),
            "total_south_asian_records": int(total_sas),
            "total_multi_ancestry_both_records": int(total_both)
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
    logger.info(f"Saved pipeline metrics summary to {metrics_out}")

    return {
        "df_low_freq": df_low_freq,
        "df_gene_mapping": df_gene_mapping,
        "df_trait_summary": df_trait_summary,
        "metrics": metrics
    }
