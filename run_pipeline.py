"""
run_pipeline.py
---------------
Master entry point for the NHGRI-EBI GWAS Catalog Analysis Pipeline.
Orchestrates downloading, data processing, statistical analysis, output exporting,
and publication figure generation.
"""

import os
import sys
import logging
import json

from download_data import main as download_datasets
from gwas_pipeline import load_and_inspect_data, process_pipeline
from plot_generator import generate_allele_freq_plot, generate_gene_disease_landscape_plot

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(base_dir, "data")
    results_dir = os.path.join(base_dir, "results")
    plots_dir = os.path.join(base_dir, "plots")

    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(plots_dir, exist_ok=True)

    logger.info("Step 1: Programmatic Data Acquisition...")
    assoc_zip_path, ancestry_path = download_datasets(data_dir)

    logger.info("Step 2: Loading & Inspecting GWAS Catalog Datasets...")
    df_assoc, df_ancestry = load_and_inspect_data(assoc_zip_path, ancestry_path)

    logger.info("Step 3: Running Pipeline ETL & Statistical Analyses...")
    pipeline_results = process_pipeline(df_assoc, df_ancestry, results_dir)

    df_low_freq = pipeline_results["df_low_freq"]
    df_gene_mapping = pipeline_results["df_gene_mapping"]
    df_trait_summary = pipeline_results["df_trait_summary"]
    metrics = pipeline_results["metrics"]

    logger.info("Step 4: Generating Publication-Quality Figures...")
    fig1_path = generate_allele_freq_plot(df_low_freq, plots_dir)
    fig2_path = generate_gene_disease_landscape_plot(
        metrics["analysis_1_low_frequency_cohort"],
        df_trait_summary,
        plots_dir
    )

    print("\n" + "="*80)
    print("GWAS PIPELINE COMPLETE - SUMMARY METRICS & RESULTS OVERVIEW")
    print("="*80)
    print(json.dumps(metrics, indent=2))
    print("="*80)
    print(f"\nResult Tables exported to: {results_dir}")
    print(f"Publication Plots exported to: {plots_dir}")
    print(f"  - Figure 1: {fig1_path}")
    print(f"  - Figure 2: {fig2_path}\n")

if __name__ == "__main__":
    main()
