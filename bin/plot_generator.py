#!/usr/bin/env python3
"""
bin/plot_generator.py
---------------------
CLI module for generating publication plots in Nextflow.
"""

import os
import json
import argparse
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

COLOR_EUR = "#1f77b4"
COLOR_SAS = "#ff7f0e"
COLOR_BOTH = "#2ca02c"

def generate_plots(low_freq_csv: str, trait_summary_csv: str, metrics_json: str, plots_dir: str):
    os.makedirs(plots_dir, exist_ok=True)

    df_filt = pd.read_csv(low_freq_csv)
    df_trait = pd.read_csv(trait_summary_csv)
    with open(metrics_json) as f:
        metrics = json.load(f)

    gene_summary = metrics["analysis_1_low_frequency_cohort"]

    # 1. Plot 1: Allele Frequency Distribution
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={'width_ratios': [2, 1]})

    eur_raf_pct = df_filt[df_filt['is_EUR']]['RAF_num'].dropna() * 100
    sas_raf_pct = df_filt[df_filt['is_SAS']]['RAF_num'].dropna() * 100

    ax1 = axes[0]
    sns.histplot(eur_raf_pct, kde=True, color=COLOR_EUR, label=f"European (n={len(eur_raf_pct):,})", bins=30, stat="density", alpha=0.3, ax=ax1)
    sns.histplot(sas_raf_pct, kde=True, color=COLOR_SAS, label=f"South Asian (n={len(sas_raf_pct):,})", bins=30, stat="density", alpha=0.3, ax=ax1)

    ax1.set_title("A. Risk Allele Frequency Distribution (0.5% - 2.0%)", fontweight='bold', pad=12)
    ax1.set_xlabel("Risk Allele Frequency (%)", fontweight='bold')
    ax1.set_ylabel("Density", fontweight='bold')
    ax1.set_xlim(0.45, 2.05)

    eur_median = eur_raf_pct.median()
    sas_median = sas_raf_pct.median()
    ax1.axvline(eur_median, color=COLOR_EUR, linestyle='--', linewidth=1.5, label=f"EUR Median: {eur_median:.3f}%")
    ax1.axvline(sas_median, color=COLOR_SAS, linestyle=':', linewidth=2.0, label=f"SAS Median: {sas_median:.3f}%")
    ax1.legend(loc='upper right')

    ax2 = axes[1]
    plot_data = pd.DataFrame({
        "Population": ["European"] * len(eur_raf_pct) + ["South Asian"] * len(sas_raf_pct),
        "RAF (%)": list(eur_raf_pct) + list(sas_raf_pct)
    })

    sns.boxplot(x="Population", y="RAF (%)", data=plot_data, hue="Population", palette=[COLOR_EUR, COLOR_SAS], legend=False, width=0.4, ax=ax2, boxprops=dict(alpha=0.8))
    ax2.set_title("B. Population RAF Summary Statistics", fontweight='bold', pad=12)
    ax2.set_xlabel("Population Cohort", fontweight='bold')
    ax2.set_ylabel("Risk Allele Frequency (%)", fontweight='bold')
    ax2.set_ylim(0.45, 2.05)

    fig.suptitle("Allele-Frequency Distribution & Comparison: European vs. South Asian GWAS Variants", fontsize=15, fontweight='bold', y=0.98)
    fig.subplots_adjust(top=0.88, wspace=0.25)

    out_path1 = os.path.join(plots_dir, "allele_freq_distribution.png")
    fig.savefig(out_path1, dpi=300, bbox_inches='tight')
    plt.close(fig)

    # 2. Plot 2: Landscape Plot
    fig2 = plt.figure(figsize=(15, 7.5))
    gs = fig2.add_gridspec(1, 2, width_ratios=[1, 1.8], wspace=0.35)

    ax1_f2 = fig2.add_subplot(gs[0])
    ax2_f2 = fig2.add_subplot(gs[1])

    cats = ["EUR-Specific", "Shared (Both)", "SAS-Specific"]
    low_freq_counts = [
        gene_summary.get("low_freq_eur_only_genes", 0),
        gene_summary.get("low_freq_shared_genes", 0),
        gene_summary.get("low_freq_sas_only_genes", 0)
    ]

    bars = ax1_f2.bar(cats, low_freq_counts, color=[COLOR_EUR, COLOR_BOTH, COLOR_SAS], alpha=0.85, edgecolor='black', linewidth=1)
    ax1_f2.set_title("A. Low-Frequency Gene Ancestry Breakdown", fontweight='bold', pad=12)
    ax1_f2.set_ylabel("Number of Unique Genes", fontweight='bold')
    ax1_f2.set_yscale('log')
    ax1_f2.set_ylim(1, max(low_freq_counts) * 5)

    for bar in bars:
        height = bar.get_height()
        ax1_f2.text(bar.get_x() + bar.get_width()/2., height * 1.25, f"{height:,}", ha='center', va='bottom', fontsize=10, fontweight='bold')

    if not df_trait.empty:
        top_traits = df_trait.head(15).set_index("Disease/Trait")[["EUR-Only Genes", "Shared Genes", "SAS-Only Genes"]]
        sns.heatmap(top_traits, annot=True, fmt="d", cmap="YlGnBu", ax=ax2_f2, cbar_kws={'label': 'Gene Count'}, linewidths=0.5, linecolor='gray')
        ax2_f2.set_title("B. Gene-Trait Mapping Across Ancestral Categories (Top 15 Traits)", fontweight='bold', pad=12)
        ax2_f2.set_xlabel("Ancestral Classification", fontweight='bold')
        ax2_f2.set_ylabel("Disease / Trait", fontweight='bold')

    fig2.suptitle("Gene-Disease-Population Landscape in GWAS Catalog Associations", fontsize=16, fontweight='bold', y=0.98)
    fig2.subplots_adjust(top=0.90, left=0.08, right=0.95, bottom=0.10)

    out_path2 = os.path.join(plots_dir, "gene_disease_population_landscape.png")
    fig2.savefig(out_path2, dpi=300, bbox_inches='tight')
    plt.close(fig2)
    logger.info(f"Saved publication figures to {plots_dir}")

def main():
    parser = argparse.ArgumentParser(description="Generate GWAS Publication Figures")
    parser.add_argument("--low_freq_csv", type=str, required=True)
    parser.add_argument("--trait_summary_csv", type=str, required=True)
    parser.add_argument("--metrics_json", type=str, required=True)
    parser.add_argument("--plots_dir", type=str, default="plots")
    args = parser.parse_args()

    generate_plots(args.low_freq_csv, args.trait_summary_csv, args.metrics_json, args.plots_dir)

if __name__ == "__main__":
    main()
