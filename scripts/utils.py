
import numpy as np
import torch
import random
from sklearn.metrics import adjusted_rand_score
import os
import torch.nn as nn

def set_seed(seed=42):
    np.random.seed(seed)
    torch.manual_seed(seed)
    random.seed(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)

def evaluate_clustering(true_labels, pred_labels):
    return adjusted_rand_score(true_labels, pred_labels)

def save_results(output_path, labels):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        for label in labels:
            f.write(f"{label}\n")


def initialize_weights(net, init_type='Xavier'):
    for m in net.modules():
        if isinstance(m, (nn.Conv2d, nn.ConvTranspose2d, nn.Linear)):
            if init_type == 'He':
                nn.init.kaiming_normal_(m.weight, nonlinearity='leaky_relu')
            elif init_type == 'Xavieru':
                nn.init.xavier_uniform_(m.weight)
            elif init_type == 'Xavier':
                nn.init.xavier_normal_(m.weight)
            elif init_type == 'Normal':
                m.weight.data.normal_(0, 0.02)
            elif init_type == 'Kaiming_Uniform':
                nn.init.kaiming_uniform_(m.weight, nonlinearity='leaky_relu')
            else:
                raise ValueError(f"Unsupported initialization type: {init_type}")
            if m.bias is not None:
                nn.init.constant_(m.bias, 0)
