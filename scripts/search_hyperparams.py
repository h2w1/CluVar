import torch
import itertools
import pandas as pd
import argparse
import time
import json
from model import AutoEncoder
from dataset import load_dataset
from train import train_autoencoder, run_clustering
from utils import set_seed, initialize_weights, evaluate_clustering

# Hyperparameter grid
init_list = ['He', 'Xavieru', 'Xavier', 'Normal']
activation_list = ['ELU', 'leaky_relu']
lr_list = [0.001, 0.0005]
param_combinations = list(itertools.product(init_list, activation_list, lr_list))

# Constants
BATCH_SIZE = 64
LATENT_DIM = 10
EPOCHS = 250
DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'

def main(args):
    set_seed(42)

    if args.data_type == 'sim':
        data_path = f"data/sim_ALT_c{args.num_clusters}.tsv"
        target_path = f"data/sim_clades_c{args.num_clusters}.tsv"
        true_labels = pd.read_csv(target_path, sep='\t', header=None).values.flatten()
    else:
        data_path = args.input_tsv
        true_labels = None

    dataloader, raw_data = load_dataset(data_path, batch_size=BATCH_SIZE)
    full_tensor = torch.tensor(raw_data.values, dtype=torch.float32)

    # 🔥 핵심 수정: input_dim 자동 추출
    input_dim = full_tensor.shape[1]
    print(f"Input dim: {input_dim}")
    print(f"Using device: {DEVICE} (CUDA available: {torch.cuda.is_available()})")
    if DEVICE == "cuda":
        print(f"CUDA Device: {torch.cuda.get_device_name(0)}")

    results = []
    start_time = time.time()

    for i, (init_type, act_type, lr) in enumerate(param_combinations):
        print(f"[{i+1}/{len(param_combinations)}] Trying: Init={init_type}, Act={act_type}, LR={lr}")

        model = AutoEncoder(input_dim=input_dim, latent_dim=LATENT_DIM, activation_type=act_type)
        initialize_weights(model, init_type=init_type)
        model = train_autoencoder(model, dataloader, DEVICE, lr=lr, epochs=EPOCHS, verbose=False)

        # ✅ 명시적으로 GPU에 모델 다시 올리기
        model.to(DEVICE)
        model.eval()
        with torch.no_grad():
            recon, _ = model(full_tensor.to(DEVICE))
            mask = (full_tensor != -1).float().to(DEVICE)
            loss = torch.nn.functional.mse_loss(recon * mask, full_tensor.to(DEVICE) * mask, reduction='sum') / mask.sum()
            loss_val = loss.item()

        if args.data_type == 'sim':
            pred_labels = run_clustering(model, full_tensor, n_clusters=args.num_clusters)
            ari_val = evaluate_clustering(true_labels, pred_labels)
            results.append(((init_type, act_type, lr), loss_val, ari_val))
        else:
            results.append(((init_type, act_type, lr), loss_val))

    elapsed = time.time() - start_time

    print("\n🎯 Recommended Top 3 Hyperparameter Combinations:")
    if args.data_type == 'sim':
        results.sort(key=lambda x: (-x[2], x[1]))
        for rank, ((init_type, act_type, lr), loss_val, ari_val) in enumerate(results[:3], 1):
            print(f"{rank}. Init: {init_type}, Activation: {act_type}, LR: {lr}, Loss: {loss_val:.4f}, ARI: {ari_val:.4f}")
        best = results[0]
        best_result = {
            "init": best[0][0],
            "activation": best[0][1],
            "lr": best[0][2],
            "loss": best[1],
            "ari": best[2]
        }
    else:
        results.sort(key=lambda x: x[1])
        for rank, ((init_type, act_type, lr), loss_val) in enumerate(results[:3], 1):
            print(f"{rank}. Init: {init_type}, Activation: {act_type}, LR: {lr}, Loss: {loss_val:.4f}")
        best = results[0]
        best_result = {
            "init": best[0][0],
            "activation": best[0][1],
            "lr": best[0][2],
            "loss": best[1]
        }

    # Save best to JSON
    with open(args.output_json, "w") as f:
        json.dump(best_result, f, indent=2)
    print(f"\n✅ Saved best hyperparams to best_hyperparams.json")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_type", type=str, default="sim", help="'sim' or 'real'")
    parser.add_argument("--input_tsv", type=str, help="Path to real ALT tsv input (real only)")
    parser.add_argument("--num_clusters", type=int, default=7, help="Used only in 'sim'")
    parser.add_argument("--output_json", type=str, default="best_hyperparams.json", help="Path to save best hyperparameters JSON")

    args = parser.parse_args()
    main(args)
