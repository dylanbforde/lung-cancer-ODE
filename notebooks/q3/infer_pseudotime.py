# file: scripts/02_infer_pseudotime.py
import scanpy as sc
import pandas as pd
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--h5ad', required=True,
                    help='Path to the processed .h5ad from txt_to_h5ad_stream.py')
parser.add_argument('--out_pseudotime', required=True,
                    help='CSV to write <CellID, PseudoTime>')
args = parser.parse_args()

print("[INFO] Loading data (backed=False so we can modify)")
adata = sc.read_h5ad(args.h5ad)

print("[INFO] Normalising and selecting HVGs …")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(
    adata, min_mean=0.0125, max_mean=3, min_disp=0.5, subset=True
)

print("[INFO] PCA → neighbours → UMAP → diffusion → DPT")
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.diffmap(adata)
sc.tl.dpt(adata)                     # diffusion pseudotime

print("[INFO] Writing pseudotime CSV")
pd.DataFrame({
    'Cell': adata.obs_names,
    'PseudoTime': adata.obs['dpt_pseudotime'].values
}).to_csv(args.out_pseudotime, index=False)

print("[DONE] Pseudotime saved to", args.out_pseudotime)