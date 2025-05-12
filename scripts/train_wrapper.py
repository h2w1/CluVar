# scripts/train_wrapper.py

import argparse
import json
import subprocess

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--hp', required=True, help="Path to best_hyperparams.json file")
    parser.add_argument('--data_type', required=True, help="Dataset type: sim or real")
    parser.add_argument('--num_clusters', type=int, required=True, help="Number of clusters")
    parser.add_argument('--data_path', type=str, required=True, help="Path to input filtered_matrix.tsv")
    parser.add_argument('--out_path', type=str, required=True, help="Path to output directory")
    args = parser.parse_args()

    # best_hyperparams.json 읽기
    with open(args.hp) as f:
        hp = json.load(f)

    lr = hp["lr"]
    init_type = hp["init"]
    act_type = hp["activation"]

    # main.py 실행
    cmd = [
        "python", "scripts/main.py",
        "--data_type", args.data_type,
        "--learning_rate", str(lr),
        "--num_clusters", str(args.num_clusters),
        "--init_type", init_type,
        "--activation_type", act_type,
        "--data_path", args.data_path,
        "--out_path", args.out_path
    ]
    subprocess.run(cmd, check=True)

if __name__ == "__main__":
    main()
