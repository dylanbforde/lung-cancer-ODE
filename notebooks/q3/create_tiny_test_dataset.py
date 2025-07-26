#!/usr/bin/env python3
"""
Create a tiny test dataset for faster Beeline algorithm testing.
Reduces both genes and cells significantly.
"""
import pandas as pd
import numpy as np
import os

def create_tiny_dataset():
    # Input paths
    input_dir = "/home/dylan/Code/lung-cancer-ODE/Beeline/inputs/q3_data/LungCancer_small"
    output_dir = "/home/dylan/Code/lung-cancer-ODE/Beeline/inputs/q3_data/LungCancer_tiny"
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Load expression data
    print("Loading expression data...")
    expr_df = pd.read_csv(f"{input_dir}/ExpressionData.csv", index_col=0)
    print(f"Original size: {expr_df.shape[0]} genes x {expr_df.shape[1]} cells")
    
    # Select highly variable genes (top 500)
    gene_vars = expr_df.var(axis=1).sort_values(ascending=False)
    top_genes = gene_vars.head(500).index
    
    # Select subset of cells (first 50)
    selected_cells = expr_df.columns[:50]
    
    # Create tiny expression dataset
    tiny_expr = expr_df.loc[top_genes, selected_cells]
    print(f"Tiny size: {tiny_expr.shape[0]} genes x {tiny_expr.shape[1]} cells")
    
    # Save expression data
    tiny_expr.to_csv(f"{output_dir}/ExpressionData.csv")
    print(f"Saved: {output_dir}/ExpressionData.csv")
    
    # Copy and filter pseudotime data
    print("Processing pseudotime data...")
    pseudo_df = pd.read_csv(f"{input_dir}/PseudoTime.csv")  # Has header, comma-separated
    
    # Filter pseudotime to match selected cells
    pseudo_filtered = pseudo_df[pseudo_df['CellID'].isin(selected_cells)]
    pseudo_filtered.to_csv(f"{output_dir}/PseudoTime.csv", header=True, index=False)
    print(f"Saved: {output_dir}/PseudoTime.csv")
    
    # Create dummy true edges file (since we don't have real regulatory network)
    print("Creating dummy true edges...")
    np.random.seed(42)
    n_edges = min(100, len(top_genes) * 2)  # Create some dummy edges
    selected_genes = list(top_genes)
    
    edges = []
    for _ in range(n_edges):
        source = np.random.choice(selected_genes)
        target = np.random.choice(selected_genes)
        if source != target:  # No self-loops
            edges.append([source, target])
    
    edges_df = pd.DataFrame(edges, columns=['Gene1', 'Gene2'])
    edges_df = edges_df.drop_duplicates()
    edges_df.to_csv(f"{output_dir}/trueEdges.csv", header=False, index=False)
    print(f"Saved: {output_dir}/trueEdges.csv with {len(edges_df)} edges")
    
    print(f"\nTiny test dataset created in: {output_dir}")
    print("You can now test Beeline algorithms much faster!")

if __name__ == "__main__":
    create_tiny_dataset()