import sys
from SPACEL.setting import set_environ_seed
set_environ_seed(42)
from SPACEL import Splane
import numpy as np
import os
import scanpy as sc
import matplotlib
import pandas as pd
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.serif'] = ['Arial']
sc.settings.set_figure_params(dpi=100,dpi_save=300,facecolor='white',fontsize=10,vector_friendly=True,figsize=(3,3))

import argparse
parser = argparse.ArgumentParser(description='SPACEL parameter')
parser.add_argument('--k', type=int, default=1)
parser.add_argument('--c', type=int, default=7)
parser.add_argument('--d_l', type=float, default=0.67)
parser.add_argument('--gnn_dropout', type=float)
args = parser.parse_args()

k = args.k
c = args.c
d_l = args.d_l
gnn_dropout = args.gnn_dropout
print(k, c, d_l, gnn_dropout)

st_ad_list = []
for i, sample in enumerate(samples):
    print(sample)

    adata = sc.read_h5ad(f'./data/ST/{sample}.h5ad')  
    cell_abundance = pd.read_csv(f'./deconvo/{sample}.csv', index_col=0, header=0)  
    cell_abundance = cell_abundance.T
    cell_abundance = cell_abundance.apply(lambda x: x/np.sum(x))
    cell_abundance = cell_abundance.T
    cell_abundance.columns = [i.split('_')[4] for i in cell_abundance.columns]
    adata.uns['celltypes'] = cell_abundance.columns
    adata.obs = pd.concat([adata.obs, cell_abundance],axis=1)
    st_ad_list.append(adata)

splane = Splane.init_model(st_ad_list, n_clusters=c, k=k, gnn_dropout=gnn_dropout, use_gpu=False)
splane.train(d_l=d_l, simi_l=0)
splane.identify_spatial_domain()

domain_df = None
for ad in st_ad_list:
    if domain_df is None:
        domain_df = ad.obs[['spatial_domain']]
    else:
        domain_df = pd.concat([domain_df, ad.obs[['spatial_domain']]], axis=0)      
        
domain_df.to_csv(f'./output/domain_df.csv')