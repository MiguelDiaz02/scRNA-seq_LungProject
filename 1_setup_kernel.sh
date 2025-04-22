#!/bin/bash

# Env Name
ENV_NAME="theLungPr"

# Verify conda availability
if ! command -v conda &> /dev/null
then
    echo "Conda is not install or is not in PATH, please install Miniconda or Anaconda first."
    exit 1
fi

# Create conda env with Python 3.10 
echo "Creating Conda: $ENV_NAME"
conda create -y -n $ENV_NAME python=3.10

# Activate the env (we use 'conda run' instead of conda init)
echo "Installing ipykernel and other packages"
conda run -n $ENV_NAME pip install ipykernel

pip install \
scanpy \
notebook \
numpy \
scipy \
matplotlib \
pandas \
scvi-tools \
seaborn \
sc_toolbox \
rpy2 \
anndata2ri \
loompy \
celltypist \
mygene

# Add the kernel to Jupyter
echo "Installing the environment as kernel to Jupyter"
conda run -n $ENV_NAME python -m ipykernel install --user --name $ENV_NAME --display-name "Python ($ENV_NAME)"

echo "Done. '$ENV_NAME' kernel is now available in Jupyter."
