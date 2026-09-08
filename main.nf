nextflow.enable.dsl=2

/*
========================================================================================
    NHGRI-EBI GWAS Catalog Pipeline - Nextflow DSL2 Workflow
========================================================================================
*/

params.assoc_url    = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations-full.zip"
params.ancestry_url = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-ancestry.tsv"
params.outdir       = "results_nextflow"
params.plots_dir    = "plots_nextflow"

/*
 * Process 1: Download Datasets
 */
process DOWNLOAD_DATA {
    tag "downloading_gwas_data"
    publishDir "${params.outdir}/raw_data", mode: 'copy'

    output:
    path "gwas_catalog_associations.zip", emit: assoc_zip
    path "gwas_catalog_ancestry.tsv",     emit: ancestry_tsv

    script:
    """
    echo "Downloading GWAS Catalog Associations..."
    curl -sSL -o gwas_catalog_associations.zip "${params.assoc_url}"

    echo "Downloading GWAS Catalog Ancestry..."
    curl -sSL -o gwas_catalog_ancestry.tsv "${params.ancestry_url}"
    """
}

/*
 * Process 2: Inspect Column Schema
 */
process INSPECT_SCHEMA {
    tag "inspect_schema"
    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    path assoc_zip

    output:
    path "column_schema_report.txt", emit: schema_report

    script:
    """
    unzip -p ${assoc_zip} | head -n 1 | tr '\t' '\n' | awk '{print NR": " \$0}' > column_schema_report.txt
    echo "=== GWAS CATALOG HEADER COLUMNS ==="
    cat column_schema_report.txt
    """
}

/*
 * Process 3: Filter Low-Frequency Variants (0.5% - 2.0% RAF) & Population Scoping
 */
process FILTER_LOW_FREQ_VARIANTS {
    tag "filter_variants"
    publishDir "${params.outdir}/analysis_1", mode: 'copy'

    input:
    path assoc_zip
    path ancestry_tsv

    output:
    path "low_freq_variants_eur_sas.csv", emit: low_freq_csv
    path "low_freq_summary.txt",          emit: summary_txt

    script:
    """
    python3 -c "
import zipfile, pandas as pd, numpy as np

with zipfile.ZipFile('${assoc_zip}') as z:
    df_assoc = pd.read_csv(z.open(z.namelist()[0]), sep='\\t', low_memory=False)

df_ancestry = pd.read_csv('${ancestry_tsv}', sep='\\t', low_memory=False)

# Parse numeric RAF
df_assoc['RAF_num'] = pd.to_numeric(df_assoc['RISK ALLELE FREQUENCY'].astype(str).str.strip(), errors='coerce')
df_assoc['is_target_raf'] = (df_assoc['RAF_num'] >= 0.005) & (df_assoc['RAF_num'] <= 0.020)

# Build PMID ancestry lookup
pmid_anc_map = df_ancestry.groupby('PUBMEDID')['BROAD ANCESTRAL CATEGORY'].apply(
    lambda s: set([c.strip() for cat in s.dropna() for c in str(cat).split(',')])
).to_dict()

def check_eur(row):
    cats = pmid_anc_map.get(row['PUBMEDID'], set())
    return any('European' in c for c in cats) or 'European' in str(row['INITIAL SAMPLE SIZE']) or 'European' in str(row['REPLICATION SAMPLE SIZE'])

def check_sas(row):
    cats = pmid_anc_map.get(row['PUBMEDID'], set())
    return any('South Asian' in c for c in cats) or 'South Asian' in str(row['INITIAL SAMPLE SIZE']) or 'South Asian' in str(row['REPLICATION SAMPLE SIZE'])

df_assoc['is_EUR'] = df_assoc.apply(check_eur, axis=1)
df_assoc['is_SAS'] = df_assoc.apply(check_sas, axis=1)

# Filter low frequency cohort
df_low_freq = df_assoc[df_assoc['is_target_raf'] & (df_assoc['is_EUR'] | df_assoc['is_SAS'])].copy()

cols = ['DATE', 'PUBMEDID', 'FIRST AUTHOR', 'STUDY', 'DISEASE/TRAIT', 'CHR_ID', 'CHR_POS', 'SNPS', 'STRONGEST SNP-RISK ALLELE', 'REPORTED GENE(S)', 'MAPPED_GENE', 'RISK ALLELE FREQUENCY', 'RAF_num', 'is_EUR', 'is_SAS']
df_low_freq[cols].to_csv('low_freq_variants_eur_sas.csv', index=False)

with open('low_freq_summary.txt', 'w') as f:
    f.write(f'Filtered Records: {len(df_low_freq)}\\n')
    f.write(f'EUR Records: {df_low_freq[\"is_EUR\"].sum()}\\n')
    f.write(f'SAS Records: {df_low_freq[\"is_SAS\"].sum()}\\n')
"
    """
}

/*
 * Process 4: Gene & Disease Population Mapping
 */
process MAP_GENES_AND_TRAITS {
    tag "gene_trait_mapping"
    publishDir "${params.outdir}/analysis_2", mode: 'copy'

    input:
    path assoc_zip
    path ancestry_tsv

    output:
    path "gene_population_trait_mapping.csv", emit: gene_map_csv
    path "population_trait_summary.csv",      emit: trait_summary_csv
    path "pipeline_metrics_summary.json",     emit: metrics_json

    script:
    """
    python3 -c "
import zipfile, json, pandas as pd, numpy as np

with zipfile.ZipFile('${assoc_zip}') as z:
    df = pd.read_csv(z.open(z.namelist()[0]), sep='\\t', low_memory=False)

df_anc = pd.read_csv('${ancestry_tsv}', sep='\\t', low_memory=False)
df['RAF_num'] = pd.to_numeric(df['RISK ALLELE FREQUENCY'].astype(str).str.strip(), errors='coerce')
df['is_target_raf'] = (df['RAF_num'] >= 0.005) & (df['RAF_num'] <= 0.020)

pmid_anc = df_anc.groupby('PUBMEDID')['BROAD ANCESTRAL CATEGORY'].apply(
    lambda s: set([c.strip() for cat in s.dropna() for c in str(cat).split(',')])
).to_dict()

df['is_EUR'] = df.apply(lambda r: any('European' in c for c in pmid_anc.get(r['PUBMEDID'], set())) or 'European' in str(r['INITIAL SAMPLE SIZE']) or 'European' in str(r['REPLICATION SAMPLE SIZE']), axis=1)
df['is_SAS'] = df.apply(lambda r: any('South Asian' in c for c in pmid_anc.get(r['PUBMEDID'], set())) or 'South Asian' in str(r['INITIAL SAMPLE SIZE']) or 'South Asian' in str(r['REPLICATION SAMPLE SIZE']), axis=1)

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

pd.DataFrame(gene_rows).to_csv('gene_population_trait_mapping.csv', index=False)

trait_rows = []
for trait, info in trait_gene_map.items():
    eur_g, sas_g = info['EUR_genes'], info['SAS_genes']
    trait_rows.append({'Disease/Trait': trait, 'Total Mapped Genes': len(eur_g.union(sas_g)), 'EUR-Only Genes': len(eur_g - sas_g), 'Shared Genes': len(eur_g.intersection(sas_g)), 'SAS-Only Genes': len(sas_g - eur_g)})

pd.DataFrame(trait_rows).sort_values(by='Total Mapped Genes', ascending=False).to_csv('population_trait_summary.csv', index=False)

metrics = {'total_records': len(df), 'mapped_genes': len(gene_rows), 'traits_count': len(trait_rows)}
with open('pipeline_metrics_summary.json', 'w') as f: json.dump(metrics, f, indent=2)
"
    """
}

/*
 * Process 5: Generate Publication Figures
 */
process GENERATE_PLOTS {
    tag "generate_plots"
    publishDir "${params.plots_dir}", mode: 'copy'

    input:
    path low_freq_csv
    path trait_summary_csv
    path metrics_json

    output:
    path "allele_freq_distribution.png",         emit: plot1
    path "gene_disease_population_landscape.png", emit: plot2

    script:
    """
    python3 -c "
import pandas as pd, numpy as np, matplotlib.pyplot as plt, seaborn as sns

df_low_freq = pd.read_csv('${low_freq_csv}')
df_trait = pd.read_csv('${trait_summary_csv}')

# Plot 1: RAF distribution
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
eur_raf = df_low_freq[df_low_freq['is_EUR']]['RAF_num'].dropna() * 100
sas_raf = df_low_freq[df_low_freq['is_SAS']]['RAF_num'].dropna() * 100

sns.histplot(eur_raf, kde=True, color='#1f77b4', label=f'European (n={len(eur_raf):,})', ax=axes[0])
sns.histplot(sas_raf, kde=True, color='#ff7f0e', label=f'South Asian (n={len(sas_raf):,})', ax=axes[0])
axes[0].set_title('A. Risk Allele Frequency Distribution (0.5% - 2.0%)')
axes[0].legend()

plot_data = pd.DataFrame({'Population': ['European']*len(eur_raf) + ['South Asian']*len(sas_raf), 'RAF (%)': list(eur_raf) + list(sas_raf)})
sns.boxplot(x='Population', y='RAF (%)', data=plot_data, palette=['#1f77b4', '#ff7f0e'], ax=axes[1])
axes[1].set_title('B. Population RAF Summary Statistics')

plt.tight_layout()
fig.savefig('allele_freq_distribution.png', dpi=300)

# Plot 2: Heatmap
fig2, ax = plt.subplots(figsize=(10, 6))
top_traits = df_trait.head(15).set_index('Disease/Trait')[['EUR-Only Genes', 'Shared Genes', 'SAS-Only Genes']]
sns.heatmap(top_traits, annot=True, fmt='d', cmap='YlGnBu', ax=ax)
ax.set_title('Gene-Disease Population Landscape (Top 15 Traits)')
plt.tight_layout()
fig2.savefig('gene_disease_population_landscape.png', dpi=300)
"
    """
}

/*
 * Workflow Orchestration
 */
workflow {
    DOWNLOAD_DATA()
    INSPECT_SCHEMA(DOWNLOAD_DATA.out.assoc_zip)
    FILTER_LOW_FREQ_VARIANTS(DOWNLOAD_DATA.out.assoc_zip, DOWNLOAD_DATA.out.ancestry_tsv)
    MAP_GENES_AND_TRAITS(DOWNLOAD_DATA.out.assoc_zip, DOWNLOAD_DATA.out.ancestry_tsv)
    GENERATE_PLOTS(
        FILTER_LOW_FREQ_VARIANTS.out.low_freq_csv,
        MAP_GENES_AND_TRAITS.out.trait_summary_csv,
        MAP_GENES_AND_TRAITS.out.metrics_json
    )
}
