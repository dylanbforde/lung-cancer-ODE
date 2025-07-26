import csv

input_file = "/home/dylan/Code/lung-cancer-ODE/Beeline/inputs/q3_data/LungCancer/ExpressionData.csv"
output_file = "/home/dylan/Code/lung-cancer-ODE/Beeline/inputs/q3_data/LungCancer_small/ExpressionData.csv"
num_columns = 100

with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
    reader = csv.reader(infile)
    writer = csv.writer(outfile)

    for i, row in enumerate(reader):
        if i == 0:  # Write header row
            writer.writerow(row[:num_columns])
        else:
            writer.writerow(row[:num_columns])
