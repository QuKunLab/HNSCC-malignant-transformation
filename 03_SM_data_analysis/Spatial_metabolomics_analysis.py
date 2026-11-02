import scanpy as sc
import scanpy.external as sce
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rc_context
import warnings
from collections import Counter
from typing import List
import os
warnings.filterwarnings('ignore')
from  tqdm import tqdm
import anndata as ad
import random
import math

from scipy.sparse import csr_matrix, coo_matrix

import PyComplexHeatmap
from PyComplexHeatmap import *

import scipy

sns.set_style('ticks')

sc.settings.verbosity = 3 
sc.settings.set_figure_params(dpi=100, dpi_save=300,
                              facecolor='white', figsize=(5,4))

merged_adata = sc.read_h5ad('./data/merged_SM.h5ad')

metabolites = ['Glutathione', 'Glutathione disulfide', 'Adrenic acid','Arachidonate',  'PE(38:4)', 'PE(40:4)',
'PE(18:0/20:4(5Z,8Z,11Z,14Z;15OOH))']

merged_adata.obs['domain_refined'] = merged_adata.obs['spatial_domain'].astype('object').copy() 
merged_adata.obs['domain_refined'] = (merged_adata.obs['domain_refined'].isin(['D2'])).map({True:'D2', False:'Other domains'})

expr_df = merged_adata[:, metabolites].to_df()
expr_df['GSH/GSSG'] = expr_df['Glutathione']/expr_df['Glutathione disulfide']
expr_df['section'] = merged_adata.obs['section']
expr_df['domain_refined'] = merged_adata.obs['domain_refined']
expr_df = expr_df.groupby(['section', 'domain_refined'], observed=True).mean().reset_index()
expr_df.columns.name = None

df = expr_df.copy()
meta_cols = ['Glutathione', 'GSH/GSSG', 'Adrenic acid','Arachidonate',  'PE(38:4)', 'PE(40:4)', 'PE(18:0/20:4(5Z,8Z,11Z,14Z;15OOH))',]
wide = df.pivot_table(
    index='section',
    columns='domain_refined',
    values=meta_cols,
    aggfunc='first'  
)

eps = 1e-12
ratio_df = wide.xs('D2', axis=1, level=1) / (wide.xs('Other domains', axis=1, level=1) + eps)
ratio_df = ratio_df.astype('float')

top_anno = pd.DataFrame(dict(group=ratio_df.index.values))
top_anno.index = top_anno.group.values

plt.figure(figsize=(12, 3))
cm = ClusterMapPlotter(data=ratio_df[['Glutathione', 'GSH/GSSG', 'Adrenic acid','Arachidonate',  'PE(38:4)', 'PE(40:4)','PE(18:0/20:4(5Z,8Z,11Z,14Z;15OOH))']].T, 
                       vmax=2,vmin=0, 
                       cmap='RdYlBu_r', show_rownames=True, show_colnames=True, 
                       row_cluster=False, 
                       col_cluster=False, 
                       legend_width=2, linecolor='white',linewidths=0.4, legend_gap=8, label='D2/Non-D2',
                       col_split_gap=0.8,
                       xticklabels=True,
                      )
plt.savefig('./SM_MALDI/figures/Heatmap_metabolites_signal_FC_between_D2_and_Non-D2.pdf', dpi=300)

plt.show()


metabolites = [
'Glutathione', 'Glutathione disulfide', 'Arachidonate', 'Adrenic acid', 
    'PE(38:4)', 'PE(40:4)','PE(18:0/20:4(5Z,8Z,11Z,14Z;15OOH))',  'L-Glutamate',  'L-Glutamine','Lactate'
]
meta_df = merged_adata.to_df()
meta_df['GSH/GSSG'] = meta_df['Glutathione']/meta_df['Glutathione disulfide']
meta_df['spatial_domain'] = merged_adata.obs['spatial_domain'].values
meta_df['section'] = merged_adata.obs['section'].values
meta_df = meta_df.groupby(['section', 'spatial_domain'])[metabolites+['GSH/GSSG']].mean().reset_index()

plot_df = meta_df[metabolites+['GSH/GSSG']].copy()
plot_df['spatial_domain'] = meta_df['spatial_domain']
plot_df = plot_df.groupby('spatial_domain').mean()

plt.figure(figsize=(6.5, 4))
cm = ClusterMapPlotter(data=plot_df[['Glutathione','GSH/GSSG','Arachidonate', 'Adrenic acid', 
    'PE(38:4)', 'PE(40:4)','PE(18:0/20:4(5Z,8Z,11Z,14Z;15OOH))',  'L-Glutamate',  'L-Glutamine','Lactate']], 
                       z_score=1,
                       vmax=1,vmin=-1, 
                       cmap='RdYlBu_r', show_rownames=True, show_colnames=True, 
                       row_cluster=False, 
                       col_cluster=False, 
                       linecolor='white',
                       linewidths=1.1, 
                       label='intensity z-score',
                      )

plt.show()


def plot_spatial(anndata, sample, key='section', **kwargs):
    sc.pl.spatial( anndata[anndata.obs[key]==sample], **kwargs)

merged_adata_smoothed = sc.read_h5ad('./data/merged_SM_smoothed.h5ad')

plot_spatial(merged_adata_smoothed, 'ST36P+T', color=['Glutathione','Glutathione disulfide',  'Arachidonate', 'Adrenic acid', 'PE(38:4)',  'PE(40:4)',  
                'L-Glutamine', 'Lactate','L-Glutamate',  'PE(18:0/20:4(5Z,8Z,11Z,14Z;15OOH))'], 
             cmap='magma',   frameon=False, ncols=10, wspace=0.01,
             img_key=None, spot_size=280,
 vmax=[18, 2.6, 4.4, 0.7, 25, 6.2, 3.1, 1.1, 12, 1.4],
            vmin = [0, 1, 2.4 , 0.1,0.1,0.1,0.1,0.1,1, 0.1 ], 

            )

plot_spatial(merged_adata_smoothed, 'ST30P+T', color=[ 'Glutathione','Glutathione disulfide',  'Arachidonate', 'Adrenic acid', 'PE(38:4)',  'PE(40:4)',   'L-Glutamine', 'Lactate','L-Glutamate',
                'PE(18:0/20:4(5Z,8Z,11Z,14Z;15OOH))'],
             cmap='magma', frameon=False, ncols=10, wspace=0.01,
              img_key=None, spot_size=250,
                 vmax=[40, 6, 2, 0.7, 11, 3, 1.8, 2.6, 7, 14],
            vmin = [0.1, 0.1, 0.8 , 0.3, 2.5, 0.8, 0.5,0.5, 1,1.5]
            )
