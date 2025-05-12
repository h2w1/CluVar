import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
from sklearn.mixture import BayesianGaussianMixture
from tqdm import tqdm

def train_autoencoder(model, dataloader, device, lr=0.0005, epochs=250, verbose=True):
    model = model.to(device)        # GPU로 모델 이동
    model.train()
    optimizer = optim.Adam(model.parameters(), lr=lr)
    criterion = nn.MSELoss(reduction='none')
    epoch_range = tqdm(range(epochs), desc="Training AE") if verbose else range(epochs)

    for epoch in epoch_range:
        epoch_loss = 0
        for batch, x in enumerate(dataloader):
            x = x.float().to(device)         # 배치도 GPU로 이동
            mask = (x != -1).float().to(device)

            optimizer.zero_grad()
            x_hat, _ = model(x)              # 모델은 이미 GPU에 있음
            loss = criterion(x_hat, x)
            loss = (loss * mask).sum() / mask.sum()
            loss.backward()
            optimizer.step()

            epoch_loss += loss.item()
        avg_loss = epoch_loss / len(dataloader)
        # Optionally log avg_loss
    return model


def run_clustering(model, full_data_tensor, n_clusters=5):
    model.eval()
    with torch.no_grad():
        device = next(model.parameters()).device  # 모델이 올라가 있는 디바이스 확인
        model = model.to(device)                  # 혹시라도 CPU에 있다면 GPU로 보내기
        input_tensor = full_data_tensor.to(device)

        # Encoder output
        latent_feature = model.encoder(input_tensor).detach().cpu().numpy()
        print("[DEBUG] latent_feature shape:", latent_feature.shape)

        # Decoder layer output
        _, decoder_outputs = model(input_tensor)
        decoder_output_0 = decoder_outputs[0].detach().cpu().numpy()
        print("[DEBUG] decoder_output_0 shape:", decoder_output_0.shape)

        # Concatenate latent + decoder output
        combined_features = np.concatenate((latent_feature, decoder_output_0), axis=1)
        print("[DEBUG] combined_features shape:", combined_features.shape)

    # Clustering
    bgmm = BayesianGaussianMixture(n_components=n_clusters, max_iter=500, n_init=10, random_state=42)
    pred_labels = bgmm.fit_predict(combined_features)
    print("[DEBUG] pred_labels shape:", len(pred_labels))
    return pred_labels




import os
import numpy as np
import torch

def save_decoder_outputs(model, full_tensor, output_dir, device):
    """
    Save the reconstructed output from the decoder.
    """
    os.makedirs(output_dir, exist_ok=True)  # 디렉터리 생성

    model = model.to(device)
    model.eval()
    with torch.no_grad():
        full_tensor = full_tensor.to(device)
        x_hat, _ = model(full_tensor)

    decoded = x_hat.detach().cpu().numpy()

    # 🔥 저장 파일 경로는 output_dir 기준으로 정확히 만듬
    output_file = os.path.join(output_dir, "decoded_matrix.tsv")

    # 저장
    np.savetxt(output_file, decoded, delimiter="\t")

    print(f"[INFO] Decoder output saved to: {output_file}")  # 정확한 전체 경로 출력
