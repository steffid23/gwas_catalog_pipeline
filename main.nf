nextflow.enable.dsl=2

/*
========================================================================================
    NHGRI-EBI GWAS Catalog Pipeline - Nextflow DSL2 Workflow (Python Engine)
========================================================================================
*/

params.assoc_url    = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations-full.zip"
params.ancestry_url = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-ancestry.tsv"
params.outdir       = "results_nextflow"
params.plots_dir    = "plots_nextflow"

/*
 * Process 1: Download Datasets using bin/download_data.py
 */
process DOWNLOAD_DATA {
    tag "download_data"
    publishDir "${params.outdir}/raw_data", mode: 'copy'

    output:
    path "data/gwas_catalog_associations.zip", emit: assoc_zip
    path "data/gwas_catalog_ancestry.tsv",     emit: ancestry_tsv

    script:
    """
    python3 ${projectDir}/bin/download_data.py --data_dir data
    """
}

/*
 * Process 2: Execute GWAS Pipeline ETL using bin/gwas_pipeline.py
 */
process RUN_GWAS_PIPELINE {
    tag "run_gwas_pipeline"
    publishDir "${params.outdir}/processed_results", mode: 'copy'

    input:
    path assoc_zip
    path ancestry_tsv

    output:
    path "results/low_freq_variants_eur_sas.csv",     emit: low_freq_csv
    path "results/gene_population_trait_mapping.csv", emit: gene_map_csv
    path "results/population_trait_summary.csv",      emit: trait_summary_csv
    path "results/pipeline_metrics_summary.json",     emit: metrics_json

    script:
    """
    python3 ${projectDir}/bin/gwas_pipeline.py \\
        --assoc_zip ${assoc_zip} \\
        --ancestry_tsv ${ancestry_tsv} \\
        --results_dir results
    """
}

/*
 * Process 3: Generate Publication Plots using bin/plot_generator.py
 */
process GENERATE_PLOTS {
    tag "generate_plots"
    publishDir "${params.plots_dir}", mode: 'copy'

    input:
    path low_freq_csv
    path trait_summary_csv
    path metrics_json

    output:
    path "plots/allele_freq_distribution.png",         emit: plot1
    path "plots/gene_disease_population_landscape.png", emit: plot2

    script:
    """
    python3 ${projectDir}/bin/plot_generator.py \\
        --low_freq_csv ${low_freq_csv} \\
        --trait_summary_csv ${trait_summary_csv} \\
        --metrics_json ${metrics_json} \\
        --plots_dir plots
    """
}

/*
 * Workflow Orchestration
 */
workflow {
    DOWNLOAD_DATA()
    RUN_GWAS_PIPELINE(
        DOWNLOAD_DATA.out.assoc_zip,
        DOWNLOAD_DATA.out.ancestry_tsv
    )
    GENERATE_PLOTS(
        RUN_GWAS_PIPELINE.out.low_freq_csv,
        RUN_GWAS_PIPELINE.out.trait_summary_csv,
        RUN_GWAS_PIPELINE.out.metrics_json
    )
}
