# file: scripts/03_write_pseudotime_for_beeline.py
import pandas as pd, argparse

ap = argparse.ArgumentParser()
ap.add_argument("--pseudotime_csv", required=True, help="output of 02_infer_pseudotime.py")
ap.add_argument("--outfile",        required=True, help="BEELINE‑style PseudoTime.txt")
args = ap.parse_args()

df = pd.read_csv(args.pseudotime_csv, index_col=0)

# Make sure pseudotime column is numeric & sort
df.iloc[:,0] = pd.to_numeric(df.iloc[:,0], errors="coerce")
df = df.sort_values(df.columns[0])

df_out = pd.DataFrame({'CellID': df.index, 'PseudoTime': df.iloc[:,0]})

df_out.to_csv(args.outfile,
          sep='	',
          header=False,
          index=False,
          columns=['CellID', 'PseudoTime'])