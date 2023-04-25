import config
import numpy as np
import scipy.io as scio
import torch
from model import Transformer

def test(scale, S, b, model_path):  
    
    with torch.no_grad():
        module = Transformer.load_from_checkpoint(model_path)
        module.cpu()
        module.eval()
        module.freeze()
    
    test_input = S[np.newaxis, :, :]
        
    test_input = torch.tensor(test_input)
    test_input = test_input.to(torch.float32)
    
    test_out = module(test_input)/config.mul_label
    
    test_out = test_out * test_input[:, :, 0][:, :, np.newaxis]
    test_out = test_out[0].cpu().detach().numpy()
    
    return test_out