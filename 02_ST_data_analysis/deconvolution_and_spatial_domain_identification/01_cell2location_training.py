import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import subprocess
import cell2location
import os
import scvi
from matplotlib import rcParams
import seaborn as sns
from scipy.sparse import csr_matrix
from cell2location.utils.filtering import filter_genes

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.serif'] = ['Arial']
sc.settings.set_figure_params(dpi=150,dpi_save=200,facecolor='white',fontsize=10,vector_friendly=True,figsize=(2,2))



ref_file_path = './data/RNA.merged.analysed.h5ad'
celltype_key='refined_celltype'
batch_key='sample'
reference_model_path= './output/reference_signatures'

### Load dataset from scRNA
adata_ref = sc.read_h5ad(ref_file_path)
adata_ref.var_names_make_unique()
adata_ref.obs_names_make_unique()
adata_ref.X = csr_matrix(adata_ref.X)

# remove cells and genes with 0 counts everywhere
sc.pp.filter_genes(adata_ref, min_cells=1)
sc.pp.filter_cells(adata_ref, min_genes=1)
adata_ref.var['SYMBOL'] = adata_ref.var_names
adata_ref.var.set_index('gene_ids', drop=True, inplace=True)

from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
# filter the object
adata_ref = adata_ref[:, selected].copy()


cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                        # 10X reaction / sample / batch
                        batch_key=batch_key,
                        # cell type, covariate used for constructing signatures
                        labels_key=celltype_key,
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        # categorical_covariate_keys=[]
                       )


# create the regression model
from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref)
mod.train(max_epochs=250, batch_size=2500, train_size=1, lr=0.002, use_gpu=True) 


# summary of the posterior distribution
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)

# Save model
mod.save(reference_model_path, overwrite=True)

# Save anndata object with results
adata_file = f"{reference_model_path}/sc.h5ad"
adata_ref.write(adata_file)
adata_file

print(reference_model_path)
