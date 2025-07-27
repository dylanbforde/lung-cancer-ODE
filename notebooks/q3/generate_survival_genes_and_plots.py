#!/usr/bin/env python3
"""
Script to generate gene lists for survival analysis and save visualizations
from beeline GENIE3 and PPCOR outputs
"""

from beeline_analysis import BeelineAnalyzer
import pandas as pd


def main():
    # Initialize analyzer
    analyzer = BeelineAnalyzer()

    # Load and preprocess data
    print("Loading data...")
    analyzer.load_data()
    analyzer.preprocess_data()

    # Calculate statistics
    print("Calculating statistics...")
    stats = analyzer.calculate_statistics()

    # Print summary
    print(analyzer.generate_summary_report())

    # Generate and save gene lists for survival analysis
    print("Generating gene lists for survival analysis...")
    gene_lists = analyzer.save_gene_lists_for_survival()

    # Print gene list summaries
    print("\n=== GENE LISTS FOR SURVIVAL ANALYSIS ===")
    for list_name, genes in gene_lists.items():
        print(f"{list_name}: {len(genes)} genes")
        if len(genes) <= 20:
            print(f"  Genes: {', '.join(genes)}")
        else:
            print(f"  Top 10: {', '.join(genes[:10])}")
        print()

    # Save visualizations
    print("Generating and saving visualizations...")
    analyzer.save_visualizations()

    # Find and display common interactions
    common_interactions = analyzer.find_common_interactions(1000, 1000)
    print(f"\nCommon interactions in top 1000: {len(common_interactions)}")

    if len(common_interactions) > 0:
        print("\nTop 10 common interactions:")
        top_common = common_interactions.sort_values(
            "GENIE3_importance", ascending=False
        ).head(10)
        for _, row in top_common.iterrows():
            print(
                f"  {row['Gene1']} <-> {row['Gene2']}: "
                f"GENIE3={row['GENIE3_importance']:.6f}, "
                f"PPCOR={row['PPCOR_correlation']:.6f}"
            )

    print("\nAnalysis complete!")
    return analyzer, gene_lists


if __name__ == "__main__":
    analyzer, gene_lists = main()

