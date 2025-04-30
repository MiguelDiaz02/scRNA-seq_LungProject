#!/bin/bash

# Env Name
ENV_NAME="theLungPr"

# Verify conda availability
if ! command -v conda &> /dev/null; then
    echo "Conda is not installed or is not in PATH. Please install Miniconda or Anaconda first."
    exit 1
fi

# Create conda env with Python 3.10 
echo "Creating Conda environment: $ENV_NAME"
conda create -y -n $ENV_NAME python=3.10

# Install packages INSIDE the environment using conda run
echo "Installing packages in $ENV_NAME"
conda run -n $ENV_NAME pip install \
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
