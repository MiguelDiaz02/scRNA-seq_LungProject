# scRNA-seq_LungProject
We are developing a scRNA-seq pipeline for testing lung fibrosis data and refer them with the HCA data

# Requiered packages (they can be all installed with pip) 

pip install \
  scanpy==1.11.1 \
  numpy==2.2.5 \
  scipy==1.15.2 \
  matplotlib==3.10.1 \
  pandas==2.2.3 \
  scvi-tools==1.3.0 \
  seaborn \
  sc-toolbox==0.12.3 \
  rpy2==3.5.17 \
  anndata2ri==1.3.2 \
  loompy==3.0.8 \
  celltypist==1.6.3 \
  mygene==3.2.2 \
  ipykernel==6.29.5

For an easier installation we advice the user to execute the 'setup_kernel.sh' and start working, this file creates a conda environment with this packages and converts it into a Jupyter compatible kernel.
