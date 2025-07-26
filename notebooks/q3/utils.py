import pandas as pd
import scanpy as sc
import os

def check_cell_ids(expression_data_path, pseudotime_path):
    """
    Checks and prints the first 5 cell IDs from expression data and pseudotime files.
    """
    try:
        e = pd.read_csv(expression_data_path, sep='	', nrows=0).columns.tolist()
        # Assuming the first column is 'GeneID' and not a cell ID
        if e and e[0] == 'GeneID':
            e = e[1:]
        p = pd.read_csv(pseudotime_path, sep='	', header=None, nrows=5)[0].tolist()
        print(f"First 5 expression data cell IDs: {e[:5]}")
        print(f"First 5 pseudotime cell IDs: {p[:5]}")
    except FileNotFoundError as fnf_error:
        print(f"Error: One of the input files not found: {fnf_error}")
    except Exception as e:
        print(f"An error occurred while checking cell IDs: {e}")

def check_h5ad(h5ad_path):
    """
    Loads and prints information about an h5ad file.
    """
    try:
        adata = sc.read_h5ad(h5ad_path, backed='r')
        print(f"Successfully loaded h5ad. Shape: {adata.shape}")
        print(f"Number of genes (obs): {len(adata.obs_names)}")
        print(f"Number of cells (var): {len(adata.var_names)}")

        if not adata.obs_names.empty:
            print(f"First 5 obs_names: {adata.obs_names[:5].tolist()}")
        else:
            print("adata.obs_names is empty.")

        if not adata.var_names.empty:
            print(f"First 5 var_names: {adata.var_names[:5].tolist()}")
        else:
            print("adata.var_names is empty.")

    except FileNotFoundError:
        print(f"Error: h5ad file not found at {h5ad_path}")
    except Exception as e:
        print(f"Error loading h5ad: {e}")

def remove_first_line(infile_path, outfile_path):
    """
    Removes the first line from a file and writes the rest to a new file.
    """
    try:
        with open(infile_path, 'r') as fin:
            fin.readline() # Skip the first line
            with open(outfile_path, 'w') as fout:
                for line in fin:
                    fout.write(line)
        print(f"Removed first line from {infile_path} and saved to {outfile_path}")
    except FileNotFoundError:
        print(f"Error: Input file not found at {infile_path}")
    except Exception as e:
        print(f"An error occurred: {e}")
