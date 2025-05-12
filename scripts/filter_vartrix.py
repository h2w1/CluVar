import os
import argparse
import numpy as np
import pandas as pd
from scipy.io import mmread, mmwrite
from scipy.sparse import csr_matrix

# 1. Load barcodes and variants
def load_metadata(barcodes_path, vcf_path):
    if not os.path.exists(barcodes_path) or not os.path.exists(vcf_path):
        raise FileNotFoundError("Barcodes or VCF file not found.")

    barcodes = pd.read_csv(barcodes_path, header=None, sep="\t").iloc[:, 0].astype(str)
    vcf = pd.read_csv(vcf_path, comment="#", header=None, sep="\t",
                      usecols=[0, 1, 3, 4],
                      names=["chrom", "pos", "ref", "alt"])
    variants = vcf.chrom.astype(str) + "_" + vcf.pos.astype(str) + "_" + vcf.ref + "_" + vcf.alt  # ref_alt로 수정
    return barcodes, variants

# 2. Load sparse matrix safely
def load_sparse_matrix_without_zeros(mtx_path: str) -> csr_matrix:
    if not os.path.exists(mtx_path):
        raise FileNotFoundError(f"Matrix file not found: {mtx_path}")
    print(f"Loading matrix from {mtx_path}...")
    mat = mmread(mtx_path).tocsr()
    mat.sum_duplicates()
    if (mat.data == 0).any():
        print("Zero entries found — cleaning...")
        mat.eliminate_zeros()
    print(f"Loaded matrix: {mat.shape[0]} rows × {mat.shape[1]} columns with {mat.nnz} nonzero entries.")
    return mat

# 3. Replace genotype values
def replace_value_sparse_csr(matrix: csr_matrix, old_val: int = 3, new_val: int = 2) -> csr_matrix:
    modified = matrix.copy()
    np.place(modified.data, modified.data == old_val, new_val)
    return modified

# 4. Calculate genotype counts and proportions
def genotype_counts_and_proportions(mat: csr_matrix):
    rows, cols = mat.shape
    total = rows * cols
    if total == 0:
        raise ValueError("Matrix has zero size (no rows or columns).")

    cnt_2 = np.sum(mat.data == 2)
    cnt_1 = np.sum(mat.data == 1)
    cnt_0 = total - mat.nnz

    counts = {2: cnt_2, 1: cnt_1, 0: cnt_0}
    props = {k: v / total for k, v in counts.items()}

    return {'counts': counts, 'props': props}

# 5. Summarize counts
def summarize_counts(label, stats, shape=None):
    c, p = stats['counts'], stats['props']
    if shape:
        print(f"[{label}] Shape: {shape[0]} rows × {shape[1]} columns")
    print(f"    2 (only alt or common): {c[2]:,} ({p[2]:.4f}),"
          f" 1 (only ref): {c[1]:,} ({p[1]:.4f}),"
          f" 0 (missing): {c[0]:,} ({p[0]:.4f})")

# 6. Select top-n variants
def filter_variants(mat: csr_matrix, top_n: int, variant_list: pd.Series):
    counts = mat.getnnz(axis=1)
    idx = np.argsort(counts)[::-1][:top_n]
    selected_variants = variant_list.iloc[idx].values  # idx에 맞춰 variants를 선택
    return mat[idx, :], selected_variants

# 7. Select cells based on missing ratio
def filter_cells(mat: csr_matrix, cell_thr: float):
    n_vars = mat.shape[0]
    counts = mat.getnnz(axis=0)
    missing_ratio = 1 - (counts / n_vars)
    keep = np.where(missing_ratio <= cell_thr)[0]
    return mat[:, keep], keep

# 8. Main pipeline
def main():
    parser = argparse.ArgumentParser(
        description="Build & filter custom genotype matrix with counts"
    )
    parser.add_argument("--barcodes", required=True)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--vartrix-mtx", required=True)
    parser.add_argument("--top-n", type=int, default=300)
    parser.add_argument("--cell-thr", type=float, default=0.95)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    # Prepare output directory
    outdir = args.out_dir
    os.makedirs(outdir, exist_ok=True)

    # 0) Load metadata
    barcodes, variants = load_metadata(args.barcodes, args.vcf)

    # 1) Load sparse matrix
    var_csr = load_sparse_matrix_without_zeros(args.vartrix_mtx)

    # 2) Build custom CSR
    custom = replace_value_sparse_csr(var_csr)

    # 3) Raw stats
    stats_raw = genotype_counts_and_proportions(custom)
    summarize_counts("custom-raw", stats_raw, shape=custom.shape)

    # 4) Filtering
    custom_top, selected_variants = filter_variants(custom, args.top_n, variants)
    custom_filt, cell_idx = filter_cells(custom_top, args.cell_thr)

    # 5) Filtered stats
    stats_filt = genotype_counts_and_proportions(custom_filt)
    summarize_counts("custom-filtered", stats_filt, shape=custom_filt.shape)

    # 6) Save outputs
    # 6) Save outputs
    print("Saving outputs...")
    mmwrite(os.path.join(outdir, "custom_filtered.mtx"), custom_filt)
    np.savetxt(os.path.join(outdir, "filtered_variants.idx"), selected_variants, fmt="%s")
    np.savetxt(os.path.join(outdir, "filtered_cells.idx"), barcodes.iloc[cell_idx], fmt="%s")
    print("All done.")

if __name__ == "__main__":
    main()
