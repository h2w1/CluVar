
import torch.nn as nn
import torch

class AutoEncoder(nn.Module):
    def __init__(self, input_dim=300, latent_dim=10, activation_type='leaky_relu'):
        super(AutoEncoder, self).__init__()

        if activation_type == 'ELU':
            self.activation = nn.ELU()
        elif activation_type == 'leaky_relu':
            self.activation = nn.LeakyReLU()
        else:
            raise ValueError(f"Unsupported activation_type: {activation_type}")

        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 200),
            nn.BatchNorm1d(200),
            self.activation,
            nn.Dropout(),
            nn.Linear(200, 100),
            nn.BatchNorm1d(100),
            self.activation,
            nn.Dropout(),
            nn.Linear(100, latent_dim)
        )

        self.decoder_layer1 = nn.Linear(latent_dim, 100)
        self.bn1 = nn.BatchNorm1d(100)
        self.dropout1 = nn.Dropout()

        self.decoder_layer2 = nn.Linear(100, 200)
        self.bn2 = nn.BatchNorm1d(200)
        self.dropout2 = nn.Dropout()
        self.decoder_output = nn.Linear(200, input_dim)

    def forward(self, x):
        z = self.encoder(x)
        d1 = self.activation(self.bn1(self.decoder_layer1(z)))
        d1 = self.dropout1(d1)
        d2 = self.activation(self.bn2(self.decoder_layer2(d1)))
        d2 = self.dropout2(d2)
        x_hat = self.decoder_output(d2)
        return x_hat, [d1, d2]

