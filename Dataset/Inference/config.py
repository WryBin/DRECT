

Type='DOSY'
save_Dataset = True

# Transformer
n_layers = 10
n_heads = 7
d_ff = 4096
lr = 1e-3
dropout = 0.1
label_size = 140

# main
batch_size = 48
mul_label = 3
max_epochs=40

# dataset
num_D = 3
n_fre = 300
o_fre = 20
max_fre = 10
num_D = 3
dB =30
min_sep = 0.5
num_samples=30
ratio=0.1
sig=0.04
base_D = 0

sig_lorz=0.06

# set device
gpu = '0'

lr_scheduler='ReduceLR'

if Type == 'DOSY':
    max_b = 0.8
    max_D = 14
    signal_dim=30

elif Type == 'VD':
    max_b = 0.122
    max_D = 2
    signal_dim=12
    base_D = 13
    num_D = 2

elif Type == 'T1T2':
    n_fre = 1
    num_D = 1
    max_b = 20
    max_D = num_D
    signal_dim=10
    min_sep =  0.1
    
    dB = num_D * 30
