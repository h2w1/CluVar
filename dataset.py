import pandas as pd
import torch
from torch.utils.data import DataLoader
import numpy as np

def load_dataset(file_path, batch_size=64):
    # Determine data type from file name
    is_sim = 'sim' in file_path.lower()

    if is_sim:
        # For simulated data
        data = pd.read_csv(file_path, sep='\t', header=None, dtype=str)
        data = data.apply(pd.to_numeric, errors='coerce')  # convert to float, invalid → NaN
        data = data.fillna(-1)
        data = data.T
        if 'ALT' in file_path:
            data[(data != 0) & (data != -1)] = 1
    else:
        # For real data
        data = pd.read_csv(file_path, sep='\t', index_col=0, dtype=str)
        data = data.apply(pd.to_numeric, errors='coerce')
        data = data.fillna(-1)
        if 'ALT' in file_path:
            data[(data != 0) & (data != -1)] = 1

    dataset = data.to_numpy().astype(np.float32)
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=False, pin_memory=True)
    return loader, data
