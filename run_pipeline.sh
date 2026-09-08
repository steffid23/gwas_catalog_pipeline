#!/usr/bin/env bash

# ==============================================================================
# NHGRI-EBI GWAS Catalog Pipeline - Pure Bash Shell Implementation
# ==============================================================================

set -euo pipefail

DATA_DIR="data_shell"
RESULTS_DIR="results_shell"
PLOTS_DIR="plots_shell"

mkdir -p "$DATA_DIR" "$RESULTS_DIR" "$PLOTS_DIR"

ASSOC_ZIP="$DATA_DIR/gwas_catalog_associations.zip"
ANCESTRY_TSV="$DATA_DIR/gwas_catalog_ancestry.tsv"

ASSOC_URL="https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations-full.zip"
ANCESTRY_URL="https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-ancestry.tsv"

echo "========================================================================"
echo "Step 1: Programmatic Data Acquisition via Curl"
echo "========================================================================"

if [ ! -f "$ASSOC_ZIP" ]; then
    echo "Downloading GWAS Catalog Associations..."
    curl -sSL -o "$ASSOC_ZIP" "$ASSOC_URL"
else
    echo "Associations ZIP already exists: $ASSOC_ZIP"
fi

if [ ! -f "$ANCESTRY_TSV" ]; then
    echo "Downloading GWAS Catalog Ancestry..."
    curl -sSL -o "$ANCESTRY_TSV" "$ANCESTRY_URL"
else
    echo "Ancestry TSV already exists: $ANCESTRY_TSV"
fi

echo ""
echo "========================================================================"
echo "Step 2: Inspecting Schema & Headers via Unzip & Awk"
echo "========================================================================"

unzip -p "$ASSOC_ZIP" | head -n 1 | tr '\t' '\n' | awk '{printf "%2d: %s\n", NR, $0}' > "$RESULTS_DIR/schema_headers.txt"
cat "$RESULTS_DIR/schema_headers.txt"

echo ""
echo "========================================================================"
echo "Step 3: Filtering Low-Frequency Variants (0.5% - 2.0% RAF) & Population Mapping"
echo "========================================================================"

python3 -c "
import zipfile, pandas as pd, numpy as np

with zipfile.ZipFile('$ASSOC_ZIP') as z:
    df = pd.read_csv(z.open(z.namelist()[0]), sep='\t', low_memory=False)

df_anc = pd.read_csv('$ANCESTRY_TSV', sep='\t', low_memory=False)

# Convert RAF and handle NR
df['RAF_num'] = pd.to_numeric(df['RISK ALLELE FREQUENCY'].astype(str).str.strip(), errors='coerce')
df['is_target_raf'] = (df['RAF_num'] >= 0.005) & (df['RAF_num'] <= 0.020)

# PubMed lookup
pmid_anc = df_anc.groupby('PUBMEDID')['BROAD ANCESTRAL CATEGORY'].apply(
    lambda s: set([c.strip() for cat in s.dropna() for c in str(cat).split(',')])
).to_dict()

df['is_EUR'] = df.apply(lambda r: any('European' in c for c in pmid_anc.get(r['PUBMEDID'], set())) or 'European' in str(r['INITIAL SAMPLE SIZE']) or 'European' in str(r['REPLICATION SAMPLE SIZE']), axis=1)
df['is_SAS'] = df.apply(lambda r: any('South Asian' in c for c in pmid_anc.get(r['PUBMEDID'], set())) or 'South Asian' in str(r['INITIAL SAMPLE SIZE']) or 'South Asian' in str(r['REPLICATION SAMPLE SIZE']), axis=1)

# Analysis 1 Filtered Output
df_filt = df[df['is_target_raf'] & (df['is_EUR'] | df['is_SAS'])].copy()
cols = ['DATE', 'PUBMEDID', 'FIRST AUTHOR', 'STUDY', 'DISEASE/TRAIT', 'CHR_ID', 'CHR_POS', 'SNPS', 'STRONGEST SNP-RISK ALLELE', 'REPORTED GENE(S)', 'MAPPED_GENE', 'RISK ALLELE FREQUENCY', 'RAF_num', 'is_EUR', 'is_SAS']
df_filt[cols].to_csv('$RESULTS_DIR/low_freq_variants_eur_sas.csv', index=False)

# Analysis 2 Gene & Trait Mapping
df_pop = df[df['is_EUR'] | df['is_SAS']].copy()
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
    anc_cat = 'Both (Shared)' if (info['EUR'] and info['SAS']) else ('European-Only' if info['EUR'] else 'South Asian-Only')
    gene_rows.append({'Gene': g, 'Ancestry_Category': anc_cat, 'EUR_Associated': info['EUR'], 'SAS_Associated': info['SAS'], 'Associated_Traits_Count': len(info['traits']), 'Associated_Variants_Count': len(info['snps'])})

pd.DataFrame(gene_rows).to_csv('$RESULTS_DIR/gene_population_trait_mapping.csv', index=False)

trait_rows = []
for trait, info in trait_gene_map.items():
    eur_g, sas_g = info['EUR_genes'], info['SAS_genes']
    trait_rows.append({'Disease/Trait': trait, 'Total Mapped Genes': len(eur_g.union(sas_g)), 'EUR-Only Genes': len(eur_only_g := eur_g - sas_g), 'Shared Genes': len(eur_g.intersection(sas_g)), 'SAS-Only Genes': len(sas_g - eur_g)})

pd.DataFrame(trait_rows).sort_values(by='Total Mapped Genes', ascending=False).to_csv('$RESULTS_DIR/population_trait_summary.csv', index=False)
"

echo "Results exported to $RESULTS_DIR"
