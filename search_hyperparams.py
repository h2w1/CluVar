import torch
from model import AutoEncoder
from dataset import load_dataset
from train import train_autoencoder
from utils import set_seed, initialize_weights, evaluate_clustering
import itertools
import pandas as pd
import argparse
import time

# Parameter search space
init_list = ['He', 'Xavieru','Xavier', 'Normal']
activation_list = ['ELU', 'leaky_relu']
lr_list = [0.001, 0.0005]
param_combinations = list(itertools.product(init_list, activation_list, lr_list))

BATCH_SIZE = 64
INPUT_DIM = 300
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
        data_path = f"data/{args.data_type}_ALT.tsv"
        true_labels = None

    dataloader, raw_data = load_dataset(data_path, batch_size=BATCH_SIZE)
    full_tensor = torch.tensor(raw_data.values, dtype=torch.float32)

    results = []
    start_time = time.time()

    for i, (init_type, act_type, lr) in enumerate(param_combinations):
        print(f"[{i+1}/{len(param_combinations)}] Searching...")

        model = AutoEncoder(input_dim=INPUT_DIM, latent_dim=LATENT_DIM, activation_type=act_type)
        initialize_weights(model, init_type=init_type)
        model = train_autoencoder(model, dataloader, DEVICE, lr=lr, epochs=EPOCHS, verbose=False)

        model.eval()
        with torch.no_grad():
            recon, _ = model(full_tensor.to(DEVICE))
            mask = (full_tensor != -1).float().to(DEVICE)
            loss = torch.nn.functional.mse_loss(recon * mask, full_tensor.to(DEVICE) * mask, reduction='sum') / mask.sum()
            loss_val = loss.item()

        if args.data_type == 'sim':
            from train import run_clustering
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
            print(f"{rank}. Init: {init_type}, Activation: {act_type}, LR: {lr}")
    else:
        results.sort(key=lambda x: x[1])
        for rank, ((init_type, act_type, lr), loss_val) in enumerate(results[:3], 1):
            print(f"{rank}. Init: {init_type}, Activation: {act_type}, LR: {lr}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_type", type=str, default="sim", help="'sim' or 'real'")
    parser.add_argument("--num_clusters", type=int, default=7, help="Number of clusters (sim only)")
    args = parser.parse_args()
    main(args)
