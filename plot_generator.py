"""
plot_generator.py
-----------------
Generates publication-quality figures for the GWAS Catalog analysis pipeline:
Figure 1: Allele-frequency distribution comparison (EUR vs SAS)
Figure 2: Gene-disease relationship across population categories (Heatmap & Bar Plot)
"""

import os
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# Set publication plotting style
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 14,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 16,
    'pdf.fonttype': 42,
    'ps.fonttype': 42
})

COLOR_EUR = "#1f77b4"      # Steel Blue
COLOR_SAS = "#ff7f0e"      # Muted Orange
COLOR_BOTH = "#2ca02c"     # Forest Green

def generate_allele_freq_plot(df_filt: pd.DataFrame, output_dir: str) -> str:
    """
    Figure 1: Allele-frequency distribution comparison between EUR and SAS.
    Combines KDE density distribution with boxplot statistics.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={'width_ratios': [2, 1]})

    # Extract RAF values for EUR and SAS in target range (0.005 to 0.020)
    eur_raf = df_filt[df_filt['is_EUR']]['RAF_num'].dropna()
    sas_raf = df_filt[df_filt['is_SAS']]['RAF_num'].dropna()

    # Convert to percentage for readable axis (0.5% - 2.0%)
    eur_raf_pct = eur_raf * 100
    sas_raf_pct = sas_raf * 100

    # Panel A: KDE + Histogram Overlay
    ax1 = axes[0]
    sns.histplot(eur_raf_pct, kde=True, color=COLOR_EUR, label=f"European (n={len(eur_raf):,})",
                 bins=30, stat="density", alpha=0.3, ax=ax1)
    sns.histplot(sas_raf_pct, kde=True, color=COLOR_SAS, label=f"South Asian (n={len(sas_raf):,})",
                 bins=30, stat="density", alpha=0.3, ax=ax1)

    ax1.set_title("A. Risk Allele Frequency Distribution (0.5% - 2.0%)", fontweight='bold', pad=12)
    ax1.set_xlabel("Risk Allele Frequency (%)", fontweight='bold')
    ax1.set_ylabel("Density", fontweight='bold')
    ax1.set_xlim(0.45, 2.05)

    # Annotate summary statistics
    eur_median = eur_raf_pct.median()
    sas_median = sas_raf_pct.median()
    ax1.axvline(eur_median, color=COLOR_EUR, linestyle='--', linewidth=1.5, label=f"EUR Median: {eur_median:.3f}%")
    ax1.axvline(sas_median, color=COLOR_SAS, linestyle=':', linewidth=2.0, label=f"SAS Median: {sas_median:.3f}%")
    ax1.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.9)

    # Panel B: Boxplot comparison
    ax2 = axes[1]
    plot_data = pd.DataFrame({
        "Population": ["European"] * len(eur_raf_pct) + ["South Asian"] * len(sas_raf_pct),
        "RAF (%)": list(eur_raf_pct) + list(sas_raf_pct)
    })

    sns.boxplot(x="Population", y="RAF (%)", data=plot_data, hue="Population",
                palette=[COLOR_EUR, COLOR_SAS], legend=False,
                width=0.4, ax=ax2, boxprops=dict(alpha=0.8), fliersize=2)
    
    ax2.set_title("B. Population RAF Summary Statistics", fontweight='bold', pad=12)
    ax2.set_xlabel("Population Cohort", fontweight='bold')
    ax2.set_ylabel("Risk Allele Frequency (%)", fontweight='bold')
    ax2.set_ylim(0.45, 2.05)

    # Add mean markers
    means = plot_data.groupby("Population")["RAF (%)"].mean()
    for idx, pop in enumerate(["European", "South Asian"]):
        mean_val = means[pop]
        ax2.scatter(idx, mean_val, color='darkred', marker='D', s=40, zorder=5)
        ax2.text(idx + 0.08, mean_val, f"Mean: {mean_val:.3f}%", fontsize=9, verticalalignment='center')

    fig.suptitle("Allele-Frequency Distribution & Comparison: European vs. South Asian GWAS Variants",
                 fontsize=15, fontweight='bold', y=0.98)
    fig.subplots_adjust(top=0.88, wspace=0.25)

    out_path = os.path.join(output_dir, "allele_freq_distribution.png")
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Saved Figure 1 to {out_path}")
    return out_path

def generate_gene_disease_landscape_plot(gene_summary: dict, trait_summary_df: pd.DataFrame, output_dir: str) -> str:
    """
    Figure 2: Relationship between genes/diseases and populations.
    Panel A: Stacked bar chart of unique genes by ancestral classification (EUR-Only, SAS-Only, Both).
    Panel B: Heatmap showing gene counts for top disease/trait categories across population classifications.
    """
    fig = plt.figure(figsize=(15, 7.5))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.8], wspace=0.35)

    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])

    # Panel A: Bar plot of gene populations (Low-frequency cohort)
    cats = ["EUR-Specific", "Shared (Both)", "SAS-Specific"]
    low_freq_counts = [
        gene_summary.get("low_freq_eur_only_genes", 0),
        gene_summary.get("low_freq_shared_genes", 0),
        gene_summary.get("low_freq_sas_only_genes", 0)
    ]
    
    bars = ax1.bar(cats, low_freq_counts, color=[COLOR_EUR, COLOR_BOTH, COLOR_SAS], alpha=0.85, edgecolor='black', linewidth=1)
    ax1.set_title("A. Low-Frequency Gene Ancestry Breakdown", fontweight='bold', pad=12)
    ax1.set_ylabel("Number of Unique Genes", fontweight='bold')
    ax1.set_yscale('log')
    ax1.set_ylim(1, max(low_freq_counts) * 5)
    
    # Add count labels above bars
    for bar in bars:
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height * 1.25,
                 f"{height:,}", ha='center', va='bottom', fontsize=10, fontweight='bold')

    # Panel B: Heatmap of top traits vs population categories
    if not trait_summary_df.empty:
        # Take top 15 traits by total gene count
        top_traits = trait_summary_df.head(15).set_index("Disease/Trait")
        heatmap_data = top_traits[["EUR-Only Genes", "Shared Genes", "SAS-Only Genes"]]

        sns.heatmap(heatmap_data, annot=True, fmt="d", cmap="YlGnBu", ax=ax2, cbar_kws={'label': 'Gene Count'},
                    linewidths=0.5, linecolor='gray')
        ax2.set_title("B. Gene-Trait Mapping Across Ancestral Categories (Top 15 Traits)", fontweight='bold', pad=12)
        ax2.set_xlabel("Ancestral Classification", fontweight='bold')
        ax2.set_ylabel("Disease / Trait", fontweight='bold')
        ax2.set_xticklabels(["EUR-Specific", "Shared (Both)", "SAS-Specific"], rotation=0)

    fig.suptitle("Gene-Disease-Population Landscape in GWAS Catalog Associations",
                 fontsize=16, fontweight='bold', y=0.98)
    fig.subplots_adjust(top=0.90, left=0.08, right=0.95, bottom=0.10)

    out_path = os.path.join(output_dir, "gene_disease_population_landscape.png")
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Saved Figure 2 to {out_path}")
    return out_path
