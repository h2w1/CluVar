
import torch
import os
import pandas as pd
import argparse
from model import AutoEncoder
from dataset import load_dataset
from train import train_autoencoder, run_clustering, save_decoder_outputs
from utils import set_seed, evaluate_clustering, save_results, initialize_weights
import numpy as np
def main(args):
    set_seed(42)
    device = 'cuda' if torch.cuda.is_available() else 'cpu'

    # 데이터 경로 설정
    data_path = f"data/{args.data_type}_ALT_c{args.num_clusters}.tsv" if args.data_type == 'sim' else f"data/{args.data_type}_ALT.tsv"
    target_path = f"data/sim_clades_c{args.num_clusters}.tsv" if args.data_type == 'sim' else None

    # 데이터 로딩
    dataloader, raw_data = load_dataset(data_path, batch_size=args.batch_size)
    full_tensor = torch.tensor(raw_data.values, dtype=torch.float32)

    # 모델 학습
    model = AutoEncoder(input_dim=args.num_feature, latent_dim=args.latent_dim)
    initialize_weights(model, init_type=args.init_type)
    model = train_autoencoder(model, dataloader, device, lr=args.learning_rate, epochs=args.epochs)

    # 디코더 출력 저장
    save_decoder_outputs(model, full_tensor, output_dir='outputs/decoder_outputs')

    # 클러스터링
    pred_labels = run_clustering(model, full_tensor, n_clusters=args.num_clusters)

    # 결과 저장
    if args.data_type == 'sim':
        true_labels = pd.read_csv(target_path, sep='\t', header=None).values.flatten()
        ari = evaluate_clustering(true_labels, pred_labels)
        print(f"ARI score: {ari:.4f}")
        save_results("outputs/sim_pred_labels.csv", pred_labels)
    else:
        save_results("outputs/real_pred_labels.csv", pred_labels)
        print(f"Saved cluster predictions of {args.data_type} dataset to {target_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_type", type=str, default="sim", help="Type of data: 'sim' or 'real'")
    parser.add_argument("--num_feature", type=int, default=300, help="Number of input features")
    parser.add_argument("--learning_rate", type=float, default=0.001, help="Learning rate for optimizer")
    parser.add_argument("--epochs", type=int, default=250, help="Number of training epochs")
    parser.add_argument("--num_clusters", type=int, default=7, help="Number of clusters for BGMM")
    parser.add_argument("--latent_dim", type=int, default=10, help="Latent dimension size")
    parser.add_argument("--batch_size", type=int, default=64, help="Batch size")
    parser.add_argument("--init_type", type=str, default="Xavier",
                        help="Weight initialization method: Normal | Xavier | Xavieru | He | Kaiming_Uniform")

    args = parser.parse_args()
    main(args)
