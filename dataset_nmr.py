def alter(file,old_str,new_str):
    """
    Replace characters in a file
    """
    file_data = ""
    with open(file, "r", encoding="utf-8") as f:
        for line in f:
            if old_str in line:
                line = line.replace(old_str,new_str)
            file_data += line
    with open(file,"w",encoding="utf-8") as f:
        f.write(file_data)

alter("config.py", "='T1T2'", "='"+'DOSY'+"'")
alter("config.py", "='VD'", "='"+'DOSY'+"'")
alter("config.py", "='DOSY'", "='"+'DOSY'+"'")

import os

import torch
import numpy as np

from tqdm import trange
from sklearn.model_selection import train_test_split as train_val

import config

import scipy.io as scio

np.random.seed(42)


def gaussian_noise_nmr(S, snr_ratio):
    """
    Add Gaussian noise to the input signal. The std of the gaussian noise is uniformly chosen between 0 and 1/sqrt(snr).
    """
    num_samples, n_fre, signal_dim = np.shape(S)
    noise_S = np.zeros([num_samples, n_fre, signal_dim])
    weight_f = np.zeros([n_fre, signal_dim])

    for i in range(num_samples):
        noise = np.random.randn(n_fre, signal_dim)
        for k in np.arange(signal_dim):
            weight_f[:, k] = S[i, :, 0] / np.max(S[i, :, 0], axis=0)

        sigma = np.max(S[i, :, 0], axis=0)/snr_ratio
        noise_S[i, :, :] = S[i, :, :] + noise * sigma * weight_f

    return noise_S
    

def load_dataloader():

    print("Begin to generate simulation signals")


    noise_signals, label = gen_signal_2D()
    
    noise_signals = torch.from_numpy(noise_signals).float()
    label = torch.from_numpy(label).float()

    if config.save_Dataset == True:
        output_dir_dataset = "./Dataset/Noise_data_nmr"
        np.save(os.path.join(output_dir_dataset, "noise_signals"), noise_signals)
        np.save(os.path.join(output_dir_dataset, "label"), label)

    print('Save successfully')
    

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


def gen_signal_2D():

    num_samples = 10
    max_D = config.max_D
    label_size = config.label_size
    sig = config.sig
    signal_dim = config.signal_dim
    n_fre = config.n_fre
    o_fre = config.o_fre
    max_fre = config.max_fre
    num_D = config.num_D
    min_sep = config.min_sep
    max_b = config.max_b


    b = np.linspace(0, max_b, signal_dim)

    S = np.zeros([num_samples, signal_dim, n_fre])
    label = np.zeros([num_samples, label_size, n_fre])

    for i in range(num_samples):

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
            
            D_value = np.array([1.1, 4, 8])

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

    Noise_S = np.zeros([4, n_fre, signal_dim])
    snr_ratio = [10, 20, 30, 40]
    for i in range(4):
        noise_S = gaussian_noise_nmr(S, snr_ratio[i])
        Noise_S[i] = noise_S[0]

    return Noise_S.astype('float32'), label.astype('float32')

if __name__ == '__main__':
    load_dataloader()
