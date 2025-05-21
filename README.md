# scRNA-seq_LungProject
We are developing a scRNA-seq pipeline for testing lung fibrosis data and refer them with the HCA data

## How to run the pipeline?

We advice the user to run as follow the scripts:

1. TheLungProject.ipynb
2. CHOIR_part0(pre).ipynb
3. run_CHOIR.R
4. Training_model_celltypist-nativeIPFmodel.ipynb
5. run_pseudotime_M3.R

NOTE: Before start playing with data it will be necessary to run ".sh" files to set the kernel and download the files (3_run_cellbender.sh won be necessary to run till authors already removed ambient RNA)

## References and Acknowledgments

This project is based on the best practices and tutorials from the scientific community. We especially thank the following works and publications for their contribution:

1. M.D. Luecken, F.J. Theis, "Current best practices in single-cell RNA-seq analysis: a tutorial", Molecular Systems Biology 15(6) (2019): e8746
2. Heumos, L., Schaar, A.C., Lance, C. et al. Best practices for single-cell analysis across modalities. Nat Rev Genet (2023). https://doi.org/10.1038/s41576-023-00586-w
3. Sanborn, Mark (2022). Sanbomics scripts. https://github.com/mousepixels/sanbomics_scripts

## Python:

You can install independently of the kernel file these dependencies using pip as:

(EXAMPLE)

```bash
pip install scanpy==1.11.1 numpy==2.2.5 scipy==1.15.2 matplotlib==3.10.1 pandas==2.2.3 \
  scvi-tools==1.3.0 seaborn==0.13.2 sc-toolbox==0.12.3 rpy2==3.5.17 anndata2ri==1.3.2 \
  loompy==3.0.8 celltypist==1.6.3 mygene==3.2.2 ipykernel==6.29.5



