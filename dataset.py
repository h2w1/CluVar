
import pandas as pd
import torch
from torch.utils.data import TensorDataset, DataLoader
import numpy as np
def load_dataset(file_path, batch_size=64):

    data = pd.read_csv(file_path, sep='\t', header=None)
    data = data.fillna(-1)  # Replace NaNs with -1
    data = data.T  # Transpose to make shape (samples, features)
    if 'ALT' in file_path:
        data[(data != 0) & (data != -1)] = 1
    dataset = data.to_numpy()
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=False, pin_memory=True)
    return loader, data
