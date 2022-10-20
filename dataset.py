import os
import sys

import torch
import numpy as np

from tqdm import trange
from sklearn.model_selection import train_test_split as train_val

import config

import scipy.io as scio

np.random.seed(42)

def gaussian_noise(S, dB):
    """
    Add Gaussian noise to the input signal. The std of the gaussian noise is uniformly chosen between 0 and 1/sqrt(snr).
    """
    snr = np.exp(np.log(10) * float(dB) / 10)
    num_samples, n_fre, signal_dim = np.shape(S)
    noise_S = np.zeros([num_samples, n_fre, signal_dim])
    sigma = np.sqrt(1. / snr)
    weight_b = np.zeros([n_fre, signal_dim])
    weight_f = np.zeros([n_fre, signal_dim])

    snr_nmr = np.zeros([num_samples, n_fre])
    
    for i in trange(num_samples):
        noise = np.random.randn(n_fre, signal_dim)
        mult = sigma * np.linalg.norm(S[i, :, :], 2) / (np.linalg.norm(noise, 2))
        noise = noise * mult
        for j in np.arange(n_fre):
            # weight_b[j, :] = 1 / (np.arange(signal_dim) + 1)
            weight_b[j, :] = -np.linspace(0, 1, config.signal_dim) + 1
        for k in np.arange(signal_dim):
            weight_f[:, k] = S[i, :, 0] / np.max(S[i, :, 0], axis=0)

        if config.Type == 'T1T2':
            noise_S[i, :, :] = S[i, :, :] + noise
        else:
            noise_S[i, :, :] = S[i, :, :] + noise * weight_f
            for j in np.arange(n_fre):
                snr_nmr[i, j] = S[i,j,0]/mult/weight_f[j,0]

    return noise_S


def load_dataloader(batch_size):

    print("Begin to generate simulation signals")

    if config.Type == 'T1T2':
        clean_signals, label = gen_signal_1D()
    else:
        clean_signals, label = gen_signal_2D()
    
    train_input, val_input, train_label, val_label = train_val(clean_signals, label, test_size=config.ratio, random_state=42)

    train_input = torch.from_numpy(train_input).float()
    val_input = torch.from_numpy(val_input).float()
    train_label = torch.from_numpy(train_label).float()
    val_label = torch.from_numpy(val_label).float()

    if config.save_Dataset == True:
        output_dir_dataset = "./Dataset/"
        np.save(os.path.join(output_dir_dataset, "train_input"), train_input)
        np.save(os.path.join(output_dir_dataset, "val_input"), val_input)
        np.save(os.path.join(output_dir_dataset, "train_label"), train_label)
        np.save(os.path.join(output_dir_dataset, "val_label"), val_label)

    train_dataset = torch.utils.data.TensorDataset(train_input, train_label)
    train_loader = torch.utils.data.DataLoader(train_dataset, batch_size=batch_size, num_workers=40, persistent_workers=True, shuffle=True)

    val_dataset = torch.utils.data.TensorDataset(val_input, val_label)
    val_loader = torch.utils.data.DataLoader(val_dataset, batch_size=batch_size, num_workers=40, persistent_workers=True, shuffle=False)

    print('successfully load dataloader')
    
    return train_loader, val_loader


def load_dataloader_exist(batch_size):

    print("Begin to generate simulation signals")

    # train_input = np.load("./Dataset/train_input.npy")
    # train_label = np.load("./Dataset/train_label.npy")
    # val_input = np.load("./Dataset/val_input.npy")
    # val_label = np.load("./Dataset/val_label.npy")

    train_input = np.load("/mnt/DATA1/chenbo/code/DRILT_new/Dataset/train_input.npy")
    train_label = np.load("/mnt/DATA1/chenbo/code/DRILT_new/Dataset/train_label.npy")
    val_input = np.load("/mnt/DATA1/chenbo/code/DRILT_new/Dataset/val_input.npy")
    val_label = np.load("/mnt/DATA1/chenbo/code/DRILT_new/Dataset/val_label.npy")

    train_input = torch.from_numpy(train_input).float()
    val_input = torch.from_numpy(val_input).float()
    train_label = torch.from_numpy(train_label).float()
    val_label = torch.from_numpy(val_label).float()

    train_dataset = torch.utils.data.TensorDataset(train_input, train_label)
    train_loader = torch.utils.data.DataLoader(train_dataset, batch_size=batch_size, num_workers=40, persistent_workers=True, shuffle=True)

    val_dataset = torch.utils.data.TensorDataset(val_input, val_label)
    val_loader = torch.utils.data.DataLoader(val_dataset, batch_size=batch_size, num_workers=40, persistent_workers=True, shuffle=False)

    print('successfully load dataloader')
    
    return train_loader, val_loader

# to generate gaussian distribution
def Gaussian_distribution(max_D, avg, num, sig):

    xgrid = np.linspace(0, max_D, num)
    sqrt_2pi=np.power(2*np.pi,0.5)
    coef=1/(sqrt_2pi*sig)
    powercoef=-1/(2*np.power(sig,2))
    mypow=powercoef*(np.power((xgrid-avg),2))
    result = coef*(np.exp(mypow))

    if config.Type == 'T1T2':
        return result/np.max(result)
        
    return result/np.tile(np.max(result, axis=1).reshape(config.o_fre, 1), [1, config.label_size]) 


def gen_signal_1D():

    num_samples = config.num_samples
    max_D = config.max_D
    label_size = config.label_size
    sig = config.sig
    signal_dim = config.signal_dim
    n_fre = 1
    min_sep = config.min_sep

    dB = config.dB

    b = (pow(2, np.arange(10)) * (config.max_b/(2**9)))

    S = np.zeros([num_samples, signal_dim, n_fre])
    label = np.zeros([num_samples, label_size, n_fre])

    for i in trange(num_samples):

        num_D = np.random.randint(config.num_D) + 1

        D = np.random.randint(0, num_D, [num_D, 1]).astype(float)

        if num_D == 1:
            D[D==0] = np.random.random(num_D) * max_D
            amp = 1
        else:
            while True:  # TODO
                D_value = np.random.random(num_D) * max_D
                D_value = D_value

                D_value = np.sort(D_value)
                D_value_t = np.roll(D_value, 1)

                if np.min(np.abs(D_value-D_value_t)) > min_sep:
                    break
            
            if n_fre == 1:
                for j in np.arange(num_D):
                    D[j] = D_value[j]
            else:
                for j in np.arange(num_D):
                    D[D==j] = D_value[j]

            # amp = (np.random.random() * 0.7) + 0.3
            # amp = np.array([[amp], [1-amp]])
            while True:
                amp = np.random.random([num_D, 1])
                if np.max(amp)/np.min(amp) < 3:
                    amp = amp/np.sum(amp)
                    break
        
        D = D + 2e-2
        signal = np.dot(np.exp(-b/D).T, amp)
        signal = signal/signal[0]

        S[i] = signal

        label[i] = np.sum(Gaussian_distribution(max_D, D, label_size, sig=sig), axis=0).reshape(config.label_size, 1)
        
    S = S.swapaxes(1, 2)
    label = label.swapaxes(1, 2)
    
    print('Add noise')
    noise_S = gaussian_noise(S, dB)

    return noise_S.astype('float32'), label.astype('float32')


def gen_signal_2D():

    num_samples = config.num_samples
    max_D = config.max_D
    label_size = config.label_size
    sig = config.sig
    signal_dim = config.signal_dim
    n_fre = config.n_fre
    o_fre = config.o_fre
    max_fre = config.max_fre
    num_D = config.num_D
    min_sep = config.min_sep
    dB = config.dB
    max_b = config.max_b


    if config.Type == 'VD':
        if os.path.exists("Dataset/" + config.Type + "_net_input.mat"):
            b = scio.loadmat("Dataset/" + config.Type + "_net_input.mat")['b'][0]
        else:
            print('Please generate the input data first, \
                Raw data can be downloaded from https://www.escholar.manchester.ac.uk/uk-ac-man-scw:303820')
            sys.exit(1)

    else:
        b = np.linspace(0, max_b, signal_dim)

    S = np.zeros([num_samples, signal_dim, n_fre])
    label = np.zeros([num_samples, label_size, n_fre])

    for i in trange(num_samples):

        fre = np.tile(np.random.random([o_fre, 1]) * max_fre, [1 ,n_fre])

        D = np.random.randint(0, num_D, [o_fre, 1]).astype(float)

        if num_D == 1:
            D[D==0] = np.random.random(num_D) * max_D
        else:
            while True:  # TODO
                D_value = (np.random.random(num_D) * max_D) + config.base_D

                D_value = np.sort(D_value)
                D_value_t = np.roll(D_value, 1)

                if np.min(np.abs(D_value-D_value_t)) > min_sep:
                    break

            for j in np.arange(num_D):
                D[D==j] = D_value[j]

        x_array = np.tile(np.linspace(0, max_fre, n_fre), [o_fre, 1])

        lorz =config.sig_lorz ** 2/((x_array-fre)**2 + config.sig_lorz**2).reshape(o_fre, n_fre)
        lorz = lorz/(np.max(lorz, axis=1).reshape([o_fre,1]))

        s = np.exp(-D*b).T
        signal = np.dot(s, lorz)
        S[i] = signal

        label[i] = np.dot(Gaussian_distribution(max_D, D - config.base_D, label_size, sig=sig).T, lorz)
        
    S = S.swapaxes(1, 2)
    label = label.swapaxes(1, 2)

    noise_S = gaussian_noise(S, dB)

    return noise_S.astype('float32'), label.astype('float32')

if __name__ == '__main__':
    load_dataloader(config.batch_size)