import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import subprocess
import os
import cell2location
import scvi
import torch
torch.cuda.is_available()
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 
import seaborn as sns
from scipy.sparse import csr_matrix
from cell2location.utils.filtering import filter_genes

spatial_file_path = 'Visium_data.h5ad'
celltype_key='refined_celltype'
batch_key='sample'
N_cells_per_location = 20
reference_model_path = './output/reference_signatures' #scRNA model
output_file_path = './output/deconvolution'
# gpu_use ='1'

# np.random.seed()
def get_freer_gpu():
    out = subprocess.getoutput('nvidia-smi -q -d Memory |grep -A4 GPU|grep Free')
    memory_available = [int(x.split()[2]) for x in out.split('\n')]
    max_idx = np.where(memory_available == np.max(memory_available))[0]
    return np.random.permutation(max_idx)[0]


### Load dataset from visium
adata_vis = sc.read_h5ad(spatial_file_path)
adata_vis.var_names_make_unique()
adata_vis.obs_names_make_unique()
adata_vis.X = csr_matrix(adata_vis.X)
adata_vis.var['SYMBOL'] = adata_vis.var_names
adata_vis.var.set_index('gene_ids', drop=True, inplace=True)

# find mitochondria-encoded (MT) genes
adata_vis.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_vis.var['SYMBOL']]
# remove MT genes for spatial mapping (keeping their counts in the object)
adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]

# load model
adata_file = f"{reference_model_path}/sc.h5ad"
adata_ref = sc.read_h5ad(adata_file)
print('scRNA data:',adata_ref.shape)

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:5]


# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
sample_names = adata_vis.obs['sample'].unique()
for sample in sample_names:
    adata_vis_subsets = adata_vis[adata_vis.obs['sample'] == sample, :]
    cell2location.models.Cell2location.setup_anndata(adata=adata_vis_subsets)
    mod = cell2location.models.Cell2location(
    adata_vis_subsets, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location = N_cells_per_location,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=20 
)
    mod.train(max_epochs=30000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
          use_gpu=True,
         )
    adata_vis_subsets = mod.export_posterior(
    adata_vis_subsets, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
)
    # Save model
    new_path = output_file_path+'/'+sample
    os.makedirs(new_path, exist_ok=True)
    mod.save(new_path, overwrite=True)
    adata_file = f"{new_path}/sp_{sample}.h5ad"
    adata_vis_subsets.write(adata_file)
    adata_file
    adata_vis_subsets.obsm['q05_cell_abundance_w_sf'].to_csv(f"{new_path}/q05_cell_abundance_w_sf_{sample}.csv")


