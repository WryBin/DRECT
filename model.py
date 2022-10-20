import copy
import math

import torch
import torch.nn as nn
import torch.nn.functional as F

import pytorch_lightning as pl

import config

DEVICE = torch.device(f"cuda:{config.gpu}")

def clones(module, N):
    """Cloning model blocks"""
    return nn.ModuleList([copy.deepcopy(module) for _ in range(N)])

class LayerNorm(nn.Module):
    def __init__(self, features, eps=1e-6):
        super(LayerNorm, self).__init__()
        # Initialize α to 1, and β to 0
        self.a_2 = nn.Parameter(torch.ones(features))
        self.b_2 = nn.Parameter(torch.zeros(features))
        # Smoothing items
        self.eps = eps

    def forward(self, x):
        # Calculate the mean and variance by the last dimension
        mean = x.mean(-1, keepdim=True)
        std = x.std(-1, keepdim=True)

        # Return the result of Layer Norm
        return self.a_2 * (x - mean) / torch.sqrt(std ** 2 + self.eps) + self.b_2


class SublayerConnection(nn.Module):
    """
    The role of SublayerConnection is to connect Multi-Head Attention and Feed Forward layers together.
    Only after each layer is output, we have to do Layer Norm first and then connect the residuals.
    sublayer is lambda function
    """
    def __init__(self, size, dropout):
        super(SublayerConnection, self).__init__()
        self.norm = LayerNorm(size)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x, sublayer):
        return x + self.dropout(sublayer(self.norm(x)))

class Encoder(nn.Module):
    # layer = EncoderLayer
    def __init__(self, layer, N):
        super(Encoder, self).__init__()
        # Clone the encoder layer
        self.layers = clones(layer, N)
        # Layer Norm
        self.norm = LayerNorm(layer.size)

    def forward(self, x):
        for layer in self.layers:
            x = layer(x)
        x = self.norm(x)

        return x

class EncoderLayer(nn.Module):
    def __init__(self, size, self_attn, feed_forward, dropout):
        super(EncoderLayer, self).__init__()
        self.self_attn = self_attn
        self.feed_forward = feed_forward

        # The role of SublayerConnection is to connect multi and ffn together
        self.sublayer = clones(SublayerConnection(size, dropout), 2)
        self.size = size

    def forward(self, x):
        x = self.sublayer[0](x, lambda x: self.self_attn(x, x, x))
        return self.sublayer[1](x, self.feed_forward)

class PositionwiseFeedForward(nn.Module):
    def __init__(self, label_size, d_ff, dropout=0.1):
        super(PositionwiseFeedForward, self).__init__()
        self.w_1 = nn.Linear(label_size, d_ff)
        self.w_2 = nn.Linear(d_ff, label_size)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x):
        return self.w_2(self.dropout(F.relu(self.w_1(x))))


def attention(query, key, value, mask=None, dropout=None):
    # Use the last dimensional value of the query matrix as d_k
    d_k = query.size(-1)

    # Swap (transpose) the last two dimensions of the key 
    # in order to multiply it with the query matrix, 
    # and divide it by the equal root of d_k after multiplication

    scores = torch.matmul(query, key.transpose(-2, -1)) / math.sqrt(d_k)

    # If there is content to be masked, replace those parts that are 0 with a large negative number 
    if mask is not None:  # TODO mask
        scores = scores.masked_fill(mask == 0, -1e9)

    # The attention matrix after mask is softmaxed according to the last dimension
    p_attn = F.softmax(scores, dim=-1)

    # If the dropout parameter is set to non-null, then the dropout operation is performed
    if dropout is not None:
        p_attn = dropout(p_attn)

    # Finally, we return the product of the attention matrix and value, and the attention matrix
    return torch.matmul(p_attn, value), p_attn

class MultiHeadedAttention(nn.Module):
    def __init__(self, h, label_size, dropout=0.1):
        super(MultiHeadedAttention, self).__init__()
        # Guaranteed to be divisible
        assert label_size % h == 0
        self.d_k = label_size // h

        # Number of heads 
        self.h = h

        # Define 4 fully connected matrices for query, key, value and attention matrix
        self.linears = clones(nn.Linear(label_size, label_size), 4)
        self.attn = None
        self.dropout = nn.Dropout(p=dropout)

    def forward(self, query, key, value, mask=None):
        if mask is not None:
            mask = mask.unsqueeze(1)

        if query.shape[1] == 1:
            return self.linears[-1](query)
             
        # The first dimension value of the query is batch size
        nbatches = query.size(0)

        # Get the query, key and value matrix
        query, key, value = [l(x).view(nbatches, -1, self.h, self.d_k).transpose(1, 2)
                             for l, x in zip(self.linears, (query, key, value))]

        # Get the attention matrix
        x, self.attn = attention(query, key, value, mask=mask, dropout=self.dropout)

        # Concat up the h multi-headed attention matrix (Change h back to the third dimension first)
        x = x.transpose(1, 2).contiguous().view(nbatches, -1, self.h * self.d_k)
        
        return self.linears[-1](x)

class Transformer(pl.LightningModule):
    def __init__(self, lr=config.lr, n_layers=config.n_layers, n_heads=config.n_heads, 
                d_ff=config.d_ff, label_size=config.label_size, dropout=config.dropout, **kwargs):
        super().__init__()

        # get the module parameters
        self.save_hyperparameters()

        c = copy.deepcopy
        attn = MultiHeadedAttention(self.hparams.n_heads, self.hparams.label_size).to(DEVICE)
        ff = PositionwiseFeedForward(self.hparams.label_size, self.hparams.d_ff, self.hparams.dropout).to(DEVICE)

        self.encoder = Encoder(EncoderLayer(self.hparams.label_size, c(attn), c(ff), self.hparams.dropout).to(DEVICE), self.hparams.n_layers).to(DEVICE)
        self.linear_en = nn.Linear(self.hparams.signal_dim, self.hparams.label_size).to(DEVICE) 
        self.norm = LayerNorm(self.hparams.label_size)

    def encode(self, src):
        return self.encoder(src)

    def forward(self, src):

        src = src/(src[:, :, 0].unsqueeze(2))

        src = self.linear_en(src)
        src = self.norm(src)
        out = self.encode(src)

        return out

    def training_step(self, batch, batch_idx):

        train_input, train_label = batch
        out = self(train_input)

        train_label = train_label*config.mul_label/train_input[:, :, 0].unsqueeze(2)  # TODO

        loss = F.mse_loss(out, train_label)

        self.log("train_loss", loss, on_step=False, on_epoch=True, prog_bar=True)

        return loss
    
    def validation_step(self, batch, batch_idx):

        train_input, train_label = batch
        out = self(train_input)

        train_label = train_label*config.mul_label/train_input[:, :, 0].unsqueeze(2)  # TODO

        loss = F.mse_loss(out, train_label)

        self.log("val_loss", loss, on_step=False, on_epoch=True, prog_bar=True)
    
    def configure_optimizers(self):

        optimizer = torch.optim.Adam(self.parameters(), lr=self.hparams.lr)
        if config.lr_scheduler=='ReduceLR':
            scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'min', patience=3, factor=0.5, verbose=True)  # TODO
        elif config.lr_scheduler=='Cosine':
            scheduler = torch.optim.lr_scheduler.CosineAnnealingWarmRestarts(optimizer, T_0=config.cos_T0, verbose=True, T_mult=2, eta_min=1e-6)

        return {
            'optimizer': optimizer,
            'lr_scheduler': scheduler,
            'monitor':'val_loss'
        }


def make_model():

    model = Transformer(config.lr, config.n_layers, config.n_heads, config.d_ff, config.label_size, config.dropout, 
                        n_fre=config.n_fre, batch_size=config.batch_size, o_fre=config.o_fre, max_fre=config.max_fre,
                        num_D=config.num_D, min_sep=config.min_sep, dB=config.dB, max_b=config.max_b, ratio=config.ratio,
                        sig=config.sig, signal_dim=config.signal_dim, sig_lorz=config.sig_lorz, max_D=config.max_D,
                        Type=config.Type, mul_label=config.mul_label, base_D=config.base_D).to(DEVICE)

    # Initialize the model
    for p in model.parameters():
        if p.dim() > 1:
            nn.init.xavier_uniform_(p)
            
    return model.to(DEVICE)