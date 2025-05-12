import pandas as pd
from scipy.io import mmread
import argparse

def convert_mtx_to_tsv(mtx_path, output_path):
    print(f"Reading matrix from: {mtx_path}")
    mat = mmread(mtx_path).tocoo()
    df = pd.DataFrame.sparse.from_spmatrix(mat).sparse.to_dense().T

    print(f"Saving to temporary TSV...")
    temp_path = output_path + ".temp"
    df.to_csv(temp_path, sep='\t', index=False, header=False, float_format="%.1f")

    print(f"Post-processing to ensure final column is well-terminated...")
    with open(temp_path, "r") as fin, open(output_path, "w") as fout:
        for line in fin:
            line = line.rstrip('\n')  # 기존 줄 끝 개행 제거
            fout.write(line + '\n')   # 강제로 한 번 개행 추가

    print("✅ Done.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mtx", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    convert_mtx_to_tsv(args.mtx, args.out)
