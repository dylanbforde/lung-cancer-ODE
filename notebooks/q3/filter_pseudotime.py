import csv

expression_data_header_file = "/home/dylan/Code/lung-cancer-ODE/Beeline/inputs/q3_data/LungCancer_small/ExpressionData.csv"
pseudotime_input_file = "/home/dylan/Code/lung-cancer-ODE/Beeline/inputs/q3_data/LungCancer/PseudoTime.csv"
pseudotime_output_file = "/home/dylan/Code/lung-cancer-ODE/Beeline/inputs/q3_data/LungCancer_small/PseudoTime.csv"

# Get valid cell IDs from the header of the smaller ExpressionData.csv
with open(expression_data_header_file, 'r') as f:
    reader = csv.reader(f)
    header = next(reader)
    # The first column is 'GeneID', so skip it to get actual cell IDs
    valid_cell_ids = set(header[1:])

# Filter PseudoTime.csv
with open(pseudotime_input_file, 'r') as infile, open(pseudotime_output_file, 'w', newline='') as outfile:
    reader = csv.reader(infile)
    writer = csv.writer(outfile)

    header_pseudotime = next(reader) # Read PseudoTime.csv header
    writer.writerow(header_pseudotime) # Write PseudoTime.csv header

    for row in reader:
        if row and row[0] in valid_cell_ids:
            writer.writerow(row)
