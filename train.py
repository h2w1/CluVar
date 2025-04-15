
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.mixture import BayesianGaussianMixture
from utils import evaluate_clustering, save_results
from tqdm.auto import tqdm
import numpy as np

def train_autoencoder(model, dataloader, device, lr=0.0005, epochs=250):
    model.to(device)
    optimizer = optim.Adam(model.parameters(), lr=lr)
    criterion = nn.MSELoss(reduction='none')
    model.train()
    
    for epoch in tqdm(range(epochs), desc="Training AE"):
        avg_loss  = 0
        epoch_loss = 0
        for batch,x in enumerate(dataloader):
            x = x.float().to(device)
            mask=(x!= -1).float()
            x_hat, _ = model(x)
            optimizer.zero_grad()
            loss = criterion(x_hat, x)
            loss = loss * mask
            loss = loss.sum() / mask.sum()

            loss.backward()
            optimizer.step()

            epoch_loss += loss.item()
        avg_loss = epoch_loss / len(dataloader)
        # Optional: print or log epoch loss here
    return model


def run_clustering(model, full_data_tensor, n_clusters=5):
    model.eval()
    with torch.no_grad():
        device = next(model.parameters()).device

        # 인코더 출력
        latent_feature = model.encoder(full_data_tensor.to(device)).cpu().numpy()

        # 디코더 출력
        _, decoder_outputs = model(full_data_tensor.to(device))
        decoder_output_0 = decoder_outputs[0].cpu().detach().numpy()  # 첫 번째 디코더 레이어 출력

        # 결합: latent + decoder layer 1
        combined_features = np.concatenate((latent_feature, decoder_output_0), axis=1)

    # 클러스터링
    bgmm = BayesianGaussianMixture(n_components=n_clusters, max_iter=500, n_init=10, random_state=42)
    pred_labels = bgmm.fit_predict(combined_features)
    return pred_labels


def save_decoder_outputs(model, full_data_tensor, output_dir='outputs/decoder_outputs'):
    model.eval()
    import os
    import numpy as np
    os.makedirs(output_dir, exist_ok=True)
    with torch.no_grad():
        _, decoder_activations = model(full_data_tensor.to(next(model.parameters()).device))
    for i, layer_out in enumerate(decoder_activations):
        np.save(f"{output_dir}/decoder_layer_{i+1}.npy", layer_out.cpu().numpy())
