# Fast and High-Resolution Reconstruction for Laplace NMR based on Deep Learning
This is the origin Pytorch implementation of AcDpILT in the following paper: Fast and High-Resolution Reconstruction for Laplace NMR based on Deep Learning. 

![image-20221019130652846](https://raw.githubusercontent.com/WryBin/Images/main/image-20221019130652846.png)

## Requirements

- Python 3.9.12
- matplotlib == 3.5.2
- scikit_learn == 1.0.2
- torch == 1.11.0+cu113
- numpy == 1.22.3
- pytorch-lightning == 1.5.10

Dependencies can be installed using the following command:
```bash
pip3 install torch==1.11.0 --extra-index-url https://download.pytorch.org/whl/cu113
pip3 install -r requirements.txt
```

## Reproducibility

#### Test

- Test.ipynb 
(Using MATLAB to plot the last result，the display code for the comparative results of the paper can be found at "DRECT/MATLAB/Result/Results_Graph_Paper.mlx". To facilitate user convenience and verification, an example of how to draw your results graph is available at "DRECT/MATLAB/Result/Results_Graph.m".It is worth noting that you can achieve optimal display results by adjusting the contour_level parameter or by adjusting the parameters in preprocessing for better model output)

#### Train

```shell
# For T1T2
python3 main.py --Type T1T2

# For DOSY
python3 main.py --Type DOSY

# For VD
python3 main.py --Type VD
```

## Pre-trained models
Something wrong with git-lfs, you can get the pre-trained models on [Google Drive](https://drive.google.com/drive/folders/1B-OZLdKW9k4eDrqzUU9UySQUqVvSripd) instead.

## Parameters
The parameters can be configured in config.py. The details about the parameters are described below.

| Parameter name | Description of parameter |
| --- | --- |
| n_layers | Num of encoder layers |
| n_heads | Num of heads                                 |
| d_ff  | Dimension of feed forward modules                            |
| lr    | Learning rate |
| dropout | The probability of dropout |
| label_size | Number of points indicating the value of diffusion coefficient |
| batch_size | The batch size of training input data |
| mul_label | Label magnification |
| max_epochs | Max train epochs                                             |
| n_fre | Number of points of chemical shift |
| o_fre | Number of simulated peaks |
| max_fre | Maximum value of simulated chemical shift |
| num_D | Maximum number of simulated molecular components |
| min_sep | Minimum separation between spikes, normalized by signal_dim |
| dB | Noise level |
| max_b | Maximum value of parameter b                                 |
| num_samples | Number of simulated samples (including training and evaluation) |
| max_D | Maximum value of simulated diffuision coefficient            |
| ratio | Ratio of validation dataset in all data                      |
| sig | Width of the simulated peak of  Gaussian distribution. |
| signal_dim | Dimension of the input signal                                |
| sig_lorz | Width of the simulated peak of  Lorentz distribution. |
| gpu | The gpu no, used for training and test |
| base_D | Minimum value of diffusion coefficient                       |
| lr_scheduler | learning rate scheduler(Cosine or ReduceLR) |

## Others
Raw data for VD available at DOI: 10.15127/1.303820.(https://www.escholar.manchester.ac.uk/uk-ac-man-scw:303820
), It needs to be processed by DOSYToolbox and our preprocessing module(MATLAB/Data/VD/Solve_data.m) after downloading.

## Contact
If you have any questions, please contact us through Email (33320211150312@stu.xmu.edu.cn) or Github issues. Pull requests are highly welcomed!
