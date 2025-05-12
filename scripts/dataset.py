import pandas as pd
import torch
from torch.utils.data import DataLoader
import numpy as np

def load_dataset(file_path, batch_size=64):
    """
    Load dataset from TSV file and return a PyTorch DataLoader and DataFrame.

    Args:
        file_path (str): Path to the input .tsv file
        batch_size (int): Batch size for DataLoader

    Returns:
        DataLoader: PyTorch DataLoader for the dataset
        DataFrame: Original pandas DataFrame (after processing)
    """

    is_sim = 'sim' in file_path.lower()

    if is_sim:
        # 🔵 Simulated data
        data = pd.read_csv(file_path, sep='\t', header=None, engine='python', dtype=str)
        data = data.apply(pd.to_numeric, errors='coerce')  # convert to float, invalid → NaN
        data = data.fillna(-1)
        data = data.T  # Transpose so rows = cells, cols = variants
    else:
        # 🟢 Real data
        data = pd.read_csv(file_path, sep='\t', header=None, engine='python', dtype=str)
        data = data.apply(pd.to_numeric, errors='coerce')
        data = data.fillna(-1)
        # (index_col=0 삭제: 지금 구조는 첫 열이 인덱스가 아니니까)
        print("[DEBUG] Loaded raw_data shape:", data.shape)
    # 🔄 VALUE MAPPING: 0 → -1, 1 → 0, 2 → 1
    data = data.astype(np.int32)
    data = data.replace({0: -1, 1: 0, 2: 1})

    # 🔥 Convert to float32 tensor for PyTorch
    dataset = data.to_numpy().astype(np.float32)
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=False, pin_memory=True)

    return loader, data

