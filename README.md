# NHGRI-EBI GWAS Catalog Pipeline (Python & Nextflow DSL2)

![Python](https://img.shields.io/badge/Python-3.10%2B-blue?logo=python)
![Nextflow](https://img.shields.io/badge/Nextflow-DSL2-green?logo=nextflow)
![License](https://img.shields.io/badge/License-MIT-brightgreen)

A reproducible, production-ready **Python and Nextflow DSL2 pipeline** designed to query, process, filter, and visualize cross-ancestry association data from the official [NHGRI-EBI GWAS Catalog](https://www.ebi.ac.uk/gwas/).

This project analyzes over **1.19 million GWAS catalog records** to evaluate **low-frequency variants (0.5% – 2.0% Risk Allele Frequency)** across **European (EUR)** and **South Asian (SAS)** populations, and maps the cross-ancestry gene-disease landscape.

---

## Quick Start

### Option 1: Run with Nextflow DSL2
```bash
nextflow run main.nf
```

### Option 2: Run with Python Engine directly
```bash
python run_pipeline.py
```

### Option 3: Run with Docker / Singularity
```bash
nextflow run main.nf -profile docker
```

---

## Project Objectives

1. **Programmatic Data Acquisition**: Automatically fetch the latest complete GWAS Catalog association and ancestry datasets using Python HTTP stream tools.
2. **Nextflow Workflow Orchestration**: Execute pipeline processes modularly using Nextflow DSL2 (`main.nf`) calling Python CLI scripts in `bin/`.
3. **Schema & Column Documentation**: Explicitly document and categorize catalog fields spanning variants, genes, allele frequencies, ancestral populations, diseases/traits, genomic positions, and study metadata.
4. **Low-Frequency Variant Filtering (Analysis 1)**: Isolate associations with Risk Allele Frequency (RAF) between **0.5% and 2.0% inclusive** ($0.005 \le \text{RAF} \le 0.020$) specifically for European and South Asian populations.
5. **Gene-Disease-Population Mapping (Analysis 2)**: Identify unique genes associated with traits/diseases and map them to population cohorts, delineating **European-specific**, **South Asian-specific**, and **Shared (Both)** genes.
6. **Publication-Quality Visualization**: Generate 300 DPI figures illustrating allele frequency distributions and population-gene-disease landscapes.

---

## Data Source & Automated Download

The pipeline ingests data directly from the official EBI FTP repository:
- **Association Dataset**: `gwas-catalog-associations-full.zip` (~66.2 MB compressed, containing `gwas-catalog-download-associations-v1.0-full.tsv`).
- **Ancestry Dataset**: `gwas-catalog-ancestry.tsv` (~49.2 MB, containing study-level curated sample and ancestry metadata).

Data download is handled programmatically via `bin/download_data.py` with HTTP streaming, size verification, and caching to avoid redundant network calls.

---

## Documented Column Schematic

| Category | Primary Columns | Description |
| :--- | :--- | :--- |
| **Variant / Mutation** | `SNPS`, `STRONGEST SNP-RISK ALLELE`, `CONTEXT`, `SNP_ID_CURRENT` | Variant rsIDs, risk allele details (e.g., `rs17723514-A`), functional context (e.g., `3_prime_UTR_variant`). |
| **Gene Information** | `MAPPED_GENE`, `REPORTED GENE(S)`, `SNP_GENE_IDS`, `UPSTREAM_GENE_ID`, `DOWNSTREAM_GENE_ID` | Curated gene symbols, reported genes, Ensembl IDs, and flanking gene distances. |
| **Allele Frequency** | `RISK ALLELE FREQUENCY` | Reported frequency of the risk allele in the discovery/replication sample. |
| **Population / Ancestry** | `INITIAL SAMPLE SIZE`, `REPLICATION SAMPLE SIZE`, `BROAD ANCESTRAL CATEGORY` | Cohort sample descriptions and curated broad ancestry classifications (linked via `PUBMEDID`). |
| **Disease / Trait** | `DISEASE/TRAIT` | Phenotype or clinical outcome investigated in the GWAS. |
| **Genomic Position** | `CHR_ID`, `CHR_POS`, `REGION` | Chromosome number, base-pair coordinate, and cytogenetic band (e.g., `18q22.1`). |
| **Study Metadata** | `PUBMEDID`, `FIRST AUTHOR`, `DATE`, `JOURNAL`, `STUDY`, `LINK` | Publication identification, study title, journal, and PubMed URL. |

---

## Methodology & Filtering Decisions

### 1. Risk Allele Frequency (RAF) Standardizing
- Non-numeric placeholders (e.g., `'NR'`, blanks) are explicitly identified and converted to `NaN`.
- Valid numerical RAF values are converted to float.
- Low-frequency variants are defined as $0.005 \le \text{RAF} \le 0.020$.

### 2. Population Scoping & Dual-Layer Matching
To prevent population misclassification, a dual-layer strategy is employed:
- **Layer 1 (Curated Ancestry)**: Merges PubMed IDs (`PUBMEDID`) with `gwas-catalog-ancestry.tsv` to match study-level `BROAD ANCESTRAL CATEGORY` tags (`European`, `South Asian`).
- **Layer 2 (Text Parsing)**: Searches `INITIAL SAMPLE SIZE` and `REPLICATION SAMPLE SIZE` text fields for `'European'` and `'South Asian'` substrings.

### 3. Multi-Gene Normalization
When `MAPPED_GENE` or `REPORTED GENE(S)` contains multiple comma-, semicolon-, or hyphen-separated gene symbols, the pipeline parses and standardizes each distinct gene symbol.

---

## Empirical Findings & Result Summary

All metrics are calculated directly from the downloaded GWAS Catalog dataset:

### Global Dataset Overview
- **Total GWAS Catalog Records**: `1,191,572`
- **Records with Valid Numeric RAF**: `575,408`
- **Records with RAF in [0.5%, 2.0%] Range**: `30,959`
- **Total European Associated Records**: `1,101,993`
- **Total South Asian Associated Records**: `267,018`
- **Total Multi-Ancestry (Both EUR & SAS) Records**: `262,153`

### Analysis 1: Low-Frequency Variant Cohort ($0.5\% \le \text{RAF} \le 2.0\%$)
- **Filtered Records**: `28,679`
  - *European Records*: `28,629`
  - *South Asian Records*: `2,990`
  - *Both (Multi-Ancestry) Records*: `2,940`
- **Unique Variants (rsIDs)**: `19,698`
- **Unique Diseases / Traits**: `5,337`
- **Low-Frequency Gene Distribution**:
  - *Total EUR Genes*: `11,969`
  - *Total SAS Genes*: `1,451`
  - *Shared Genes (Both EUR & SAS)*: `1,432`
  - *EUR-Only Genes*: `10,537`
  - *SAS-Only Genes*: `19`

### Analysis 2: Full Catalog Gene-Trait-Population Landscape
- **Total Unique Mapped Genes**: `30,238`
- **Total Unique Diseases / Traits**: `40,201`
- **Ancestral Gene Breakdown**:
  - *European-Associated Genes*: `30,187`
  - *South Asian-Associated Genes*: `20,977`
  - *Shared (Both EUR & SAS) Genes*: `20,926`
  - *EUR-Specific Genes*: `9,261`
  - *SAS-Specific Genes*: `51`

---

## Publication Figures

### Figure 1: Allele-Frequency Distribution & Comparison (`plots/allele_freq_distribution.png`)
Dual-panel plot showing:
- **Panel A**: Overlapping KDE density curves and histograms comparing European (blue) and South Asian (orange) risk allele frequency distributions in the low-frequency window (0.5% – 2.0%), annotated with median markers.
- **Panel B**: Boxplot showing distribution quartiles, interquartile ranges, and mean markers for each population cohort.

### Figure 2: Gene-Disease-Population Landscape (`plots/gene_disease_population_landscape.png`)
Dual-panel plot showing:
- **Panel A**: Stacked bar chart of unique genes categorized by ancestral specificity (EUR-Specific, Shared, SAS-Specific) on a logarithmic scale.
- **Panel B**: Heatmap showing gene counts mapped to the top 15 disease/trait categories across ancestral classifications.

---

## Repository Structure

```
gwas_catalog_pipeline/
├── bin/                                 # Executable Python CLI scripts invoked by Nextflow
│   ├── download_data.py
│   ├── gwas_pipeline.py
│   └── plot_generator.py
├── main.nf                              # Nextflow DSL2 workflow
├── nextflow.config                      # Nextflow configuration parameters
├── run_pipeline.py                      # Standalone Python entrypoint
├── run_pipeline.sh                      # Standalone Bash entrypoint
├── data/                                # Raw dataset storage
├── results/                             # Processed analytical outputs
│   ├── low_freq_variants_eur_sas.csv     # Filtered variant table (0.5% - 2.0% RAF)
│   ├── gene_population_trait_mapping.csv # Mapped gene table with ancestry breakdown
│   ├── population_trait_summary.csv      # Disease/trait summary by population
│   └── pipeline_metrics_summary.json     # Empirical counts & JSON metrics
├── plots/                                # Publication-grade figures
│   ├── allele_freq_distribution.png
│   └── gene_disease_population_landscape.png
├── requirements.txt                     # Python dependencies
└── README.md                            # Documentation
```

---

## Reproduction Instructions

### 1. Environment Setup
Clone the repository and ensure Python 3.10+ and Nextflow are installed:
```bash
git clone https://github.com/steffid23/gwas_catalog_pipeline.git
cd gwas_catalog_pipeline
pip install -r requirements.txt
```

### 2. Execute Nextflow Workflow
```bash
nextflow run main.nf
```

### 3. Execute Standalone Python Pipeline
```bash
python run_pipeline.py
```
