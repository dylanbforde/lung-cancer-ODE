import pandas as pd
import os

input_expr_csv = "/home/dylan/Code/lung-cancer-ODE/Beeline/inputs/q3_data/LungCancer/ExpressionData.csv"
output_pseudotime_csv = "/home/dylan/Code/lung-cancer-ODE/Beeline/inputs/q3_data/LungCancer/PseudoTime.csv"

if not os.path.exists(input_expr_csv):
    print(f"Error: {input_expr_csv} not found. Cannot create PseudoTime.csv.")
else:
    print(f"Reading cell IDs from {input_expr_csv}...")
    # Read only the header to get cell IDs
    df_expr = pd.read_csv(input_expr_csv, sep='\t', nrows=0)
    cell_ids = df_expr.columns.tolist()

    # Remove the first column name (GeneID) which is not a cell ID
    if cell_ids and cell_ids[0] == 'GeneID':
        cell_ids = cell_ids[1:]

    if not cell_ids:
        print("No cell IDs found in ExpressionData.csv. Cannot create PseudoTime.csv.")
    else:
        # Create a dummy pseudotime DataFrame with CellID as a regular column
        df_pseudotime = pd.DataFrame({'CellID': cell_ids, 'PseudoTime': range(len(cell_ids))})
        df_pseudotime.to_csv(output_pseudotime_csv, index=False) # Do not write DataFrame index
        print(f"Dummy PseudoTime.csv created at {output_pseudotime_csv} with {len(cell_ids)} cells.")