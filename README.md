# scRNA-seq_LungProject
We are developing a scRNA-seq pipeline for testing lung fibrosis data and refer them with the HCA data

## Python:

- scanpy (v. 1.11.1)  
- numpy (v. 2.2.5)  
- scipy (v. 1.15.2)  
- matplotlib (v. 3.10.1)  
- pandas (v. 2.2.3)  
- scvi-tools (v. 1.3.0)  
- seaborn (v. 0.13.2)  
- sc-toolbox (v. 0.12.3)  
- rpy2 (v. 3.5.17)  
- anndata2ri (v. 1.3.2)  
- loompy (v. 3.0.8)  
- celltypist (v. 1.6.3)  
- mygene (v. 3.2.2)  
- ipykernel (v. 6.29.5)  

You can install these dependencies using pip:

```bash
pip install scanpy==1.11.1 numpy==2.2.5 scipy==1.15.2 matplotlib==3.10.1 pandas==2.2.3 \
  scvi-tools==1.3.0 seaborn==0.13.2 sc-toolbox==0.12.3 rpy2==3.5.17 anndata2ri==1.3.2 \
  loompy==3.0.8 celltypist==1.6.3 mygene==3.2.2 ipykernel==6.29.5


For an easier installation we advice the user to execute the 'setup_kernel.sh' and start working, this file creates a conda environment with this packages and converts it into a Jupyter compatible kernel.
