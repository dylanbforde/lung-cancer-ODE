# file: txt_to_h5ad_stream.py
import argparse, gc, sys, h5py, numpy as np, scanpy as sc
from pathlib import Path
import pandas as pd

def parse():
    p = argparse.ArgumentParser()
    p.add_argument("--txt", required=True)
    p.add_argument("--h5ad", required=True)
    p.add_argument("--csv_out", required=True, help="Path to output ExpressionData.csv")
    p.add_argument("--pseudotime_txt", required=True, help="Path to PseudoTime.txt for cell ordering")
    p.add_argument("--buffer", type=int, default=200)   # genes per flush
    p.add_argument("--dtype", default="uint32")
    return p.parse_args()

def main():
    args = parse()
    txt   = Path(args.txt).expanduser()
    out   = Path(args.h5ad).expanduser()
    assert txt.exists()

    print(f"[INFO] streaming {txt.name} → {out.name} (buf={args.buffer})")

    # ── open input & parse header ────────────────────────────────────────────
    fh = txt.open()
    header = fh.readline().rstrip("\n").split("\t")
    if header[0].lower() in {"gene", "genes", "gene_id", "symbol"}:
        header = header[1:]                       # drop gene column label
    cells = np.array(header, dtype="U")
    n_cells = len(cells)

    # Remove 'Index' if present in cells (from original header)
    if cells[0].lower() == "index":
        cells = cells[1:]
        n_cells = len(cells)

    # Read PseudoTime.txt for ordered cell IDs
    ptime_df = pd.read_csv(args.pseudotime_txt, sep='\t', header=None, usecols=[0])
    ordered_cells = ptime_df.iloc[:,0].tolist()

    # Reorder the cells array based on pseudotime order
    # Create a mapping from original cell ID to its index
    cell_id_to_idx = {cell: i for i, cell in enumerate(cells)}
    # Create a new array of counts in the ordered sequence
    ordered_indices = np.array([cell_id_to_idx[cell] for cell in ordered_cells])

    cells = np.array(ordered_cells, dtype="U") # Update cells to be in ordered sequence
    n_cells = len(cells)

    # Initialize CSV output
    csv_out = Path(args.csv_out).expanduser()
    csv_fh = csv_out.open("w")
    csv_fh.write("GeneID," + ",".join(cells) + "\n") # Write CSV header

    # ── prepare HDF5 container ──────────────────────────────────────────────
    with h5py.File(out, "w") as h5:
        # 1️⃣  root‑level attrs required by AnnData
        h5.attrs["encoding-type"]      = "anndata"
        h5.attrs["encoding-version"]   = "0.1.0"
        h5.attrs["n_obs"]              = 0
        h5.attrs["n_vars"]             = n_cells

        # 2️⃣  /X group – holds CSR pieces
        gX = h5.create_group("X")
        gX.attrs["encoding-type"]    = "csr_matrix"
        gX.attrs["encoding-version"] = "0.1.0"

        gX.create_dataset("data",   (0,), maxshape=(None,), dtype=args.dtype,   chunks=True)
        gX.create_dataset("indices",(0,), maxshape=(None,), dtype="uint32",     chunks=True)
        gX.create_dataset("indptr", (1,), maxshape=(None,), dtype="uint32",     chunks=True)
        gX["indptr"][0] = 0

        # 3️⃣  /obs and /var
        g_obs = h5.create_group("obs")
        g_var = h5.create_group("var")

        g_obs.create_dataset("_index", (0,), maxshape=(None,), dtype="S32")  # gene names
        g_var.create_dataset("_index", data=cells.astype("S32"))             # cell barcodes

        # keep handles for later
        ds_data   = gX["data"]
        ds_idx    = gX["indices"]
        ds_indptr = gX["indptr"]
        ds_obs    = g_obs["_index"]

        nnz, buf_rows, data_collected, idx_collected, csv_rows = 0, [], [], [], []
        indptr = [0]

        for r, line in enumerate(fh, 1):
            gene, rest = line.split("\t", 1)
            counts = np.fromstring(rest, sep="\t", dtype=args.dtype)
            # Reorder counts based on ordered_cells
            counts = counts[ordered_indices]
            nz = counts.nonzero()[0]
            if nz.size:
                data_collected.append(counts[nz])
                idx_collected.append(nz.astype("uint32"))
                nnz += nz.size
            indptr.append(nnz)
            buf_rows.append(gene.encode())  # bytes for fixed‑length dtype
            counts_str = ','.join(map(str, counts))
            csv_rows.append(f"{gene},{counts_str}") # For CSV output

            # flush every buffer rows
            if len(buf_rows) >= args.buffer:
                flush(h5, ds_data, ds_idx, ds_indptr, ds_obs, buf_rows, data_collected, idx_collected, indptr)
                csv_flush(csv_fh, csv_rows)
                buf_rows, data_collected, idx_collected, indptr, csv_rows = [], [], [], [nnz], []

        # final flush
        flush(h5, ds_data, ds_idx, ds_indptr, ds_obs, buf_rows, data_collected, idx_collected, indptr, final=True)
        csv_flush(csv_fh, csv_rows, final=True)

    # Debugging: Inspect the HDF5 file directly before Scanpy tries to read it
    try:
        with h5py.File(out, "r") as f_debug:
            print(f"[DEBUG] HDF5 file root keys: {list(f_debug.keys())}")
            if "obs" in f_debug:
                print(f"[DEBUG] 'obs' group keys: {list(f_debug['obs'].keys())}")
                if "_index" in f_debug["obs"]:
                    print(f"[DEBUG] 'obs/_index' shape: {f_debug['obs']['_index'].shape}")
            else:
                print("[DEBUG] 'obs' group NOT found at root level.")
    except Exception as e:
        print(f"[DEBUG] Error inspecting HDF5 file: {e}")

    # ── wrap in AnnData stub so Scanpy understands it ───────────────────────
    ad = sc.read_h5ad(out, backed="r")
    print("[DONE]", ad)
    csv_fh.close()

def csv_flush(csv_fh, rows, final=False):
    for row in rows:
        csv_fh.write(row + "\n")
    if final:
        csv_fh.close()

def flush(h5, ds_data, ds_idx, ds_indptr, ds_obs,
          rows, datas, idxs, indptr, final=False):
    if not rows:
        return

    # grow datasets
    new_nnz   = sum(len(d) for d in datas)
    new_rows  = len(rows)

    ds_data.resize(  (ds_data.shape[0]   + new_nnz,) )
    ds_idx.resize(   (ds_idx.shape[0]    + new_nnz,) )
    ds_indptr.resize((ds_indptr.shape[0] + len(indptr) - 1,) )
    ds_obs.resize(   (ds_obs.shape[0]    + new_rows,) )

    # writable views
    d_view  = ds_data   [-new_nnz:]
    i_view  = ds_idx    [-new_nnz:]
    ip_view = ds_indptr[-(len(indptr) - 1):]
    obs_view = ds_obs   [-new_rows:]

    # copy CSR pieces
    offset = 0
    for d, ix in zip(datas, idxs):
        nn = len(d)
        d_view[offset:offset+nn] = d
        i_view[offset:offset+nn] = ix
        offset += nn
    ip_view[:] = indptr[1:]                   # skip duplicate first entry
    obs_view[:] = rows                        # gene names

    # update global attrs
    h5.attrs["n_obs"] = ds_obs.shape[0]

    print(f"  · flushed {new_rows} genes "
          f"(rows total {ds_obs.shape[0]}, nnz {ds_data.shape[0]})")
    gc.collect()

if __name__ == "__main__":
    main()