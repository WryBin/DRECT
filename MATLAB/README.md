# How to run python functions in Matlab

#### Step-1: install python packages required

```shell
pip3 install torch==1.11.0 --extra-index-url https://download.pytorch.org/whl/cu113
pip3 install -r requirements.txt
```

#### Step-2: Specify the path of python interpreter  in Matlab command window

```matlab
pyversion [your Python interpreter path]
```

for example:

```matlab
pyversion C:\Users\20440\miniconda3\envs\DRILT\python.exe   [windows]
```

