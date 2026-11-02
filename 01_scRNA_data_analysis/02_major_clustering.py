import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams["font.sans-serif"] = "Arial"
import os
from os.path import dirname, join
import scanpy.external as sce
from matplotlib.pyplot import rc_context
import warnings
warnings.filterwarnings('ignore')

sc.settings.verbosity = 3 

sc.settings.set_figure_params(dpi=100, dpi_save=300,
                              facecolor='white', figsize=(4,4))

merged_adata_raw = sc.read_h5ad('./raw_data/RNA_merged_raw.h5ad')

merged_adata_raw.var['mt'] = merged_adata_raw.var_names.str.startswith('MT-') 
sc.pp.calculate_qc_metrics(merged_adata_raw, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

merged_adata = merged_adata_raw.copy()

sc.pp.filter_cells(merged_adata, min_genes=500)   
sc.pp.filter_genes(merged_adata, min_cells=50)
merged_adata.var['mt'] = merged_adata.var_names.str.startswith('MT-')  
merged_adata.var['RP'] = merged_adata.var_names.str.contains("^RP[SL]|^RPLP|^RPSA",regex=True)
sc.pp.calculate_qc_metrics(merged_adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pp.calculate_qc_metrics(merged_adata, qc_vars=['RP'], percent_top=None, log1p=False, inplace=True)

with rc_context({'figure.figsize': (20, 4)}):
    sc.pl.violin(merged_adata, ['n_genes_by_counts'],groupby='sample' , 
             jitter=0.4, multi_panel=True, stripplot=False, rotation=45)

with rc_context({'figure.figsize': (20, 4)}):
    sc.pl.violin(merged_adata, ['total_counts'],groupby='sample' , 
             jitter=0.4, multi_panel=True, stripplot=False, rotation=45)

with rc_context({'figure.figsize': (20, 4)}):
    sc.pl.violin(merged_adata, ['pct_counts_mt'],groupby='sample' , 
             jitter=0.4, multi_panel=True, stripplot=False, rotation=45)

DoubletDetection =  pd.read_csv('./remove_doublet/douletdetection/outputs/Merged_doublets.DoubletDetection.result.csv', index_col=0)

#remove dobulet
adata = merged_adata[merged_adata.obs_names.isin(DoubletDetection.index[DoubletDetection.labels == 'Singlet'])]

#remove genes
RPL = adata.var_names[adata.var_names.str.contains("^RP[SL]|^RPLP|^RPSA",regex=True)].tolist()  
MT = adata.var_names[adata.var_names.str.startswith("MT-")].tolist()
adata = adata[:,~adata.var_names.isin(RPL+MT)] 

#remove cells
adata = adata[adata.obs.n_genes_by_counts < 6000, :]
adata = adata[adata.obs.pct_counts_mt < 20, :] 
adata = adata[adata.obs.total_counts <40000 , :]

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, n_top_genes=3000)  
adata.raw = adata
adata = adata[:, adata.var.highly_variable]    
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])  
sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, svd_solver='arpack')

sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50)

sce.pp.harmony_integrate(adata, key='sample', sigma= 0.07, lamb=1.0, theta=0.5)   
sc.pp.neighbors(adata, use_rep='X_pca_harmony', n_pcs=40,  n_neighbors=15)   
sc.tl.umap(adata)   

adata.uns['coarse_condition_colors'] = ['#7ed4ad', '#d76c82', '#3d0301', '#fff4b7']

with rc_context({'figure.figsize': (7, 7)}):
    sc.pl.umap(adata, 
               color=[ 'coarse_condition','tissue'],size=2)  

for i in [0.01,0.05, 0.07, 0.09, 0.1,0.2,0.3,0.4,0.5,0.6]:     
    sc.tl.leiden(adata, resolution=i, key_added='leiden_'+str(i)) 

with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(adata, 
                legend_loc='on data',
               color=[ 'leiden_'+str(i) for i in [0.6]], size=3) 

sc.tl.rank_genes_groups(adata, 'leiden_0.6', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, ncols=5)

print(adata.shape)

celltype_map={
'0': 'CD8 T',   
'1':  'CD4 T',  
'2':   'Endothelium',  
'3':  'Epithelium',
'4': 'Neutrophil',  
'5':   'Fibroblast', 
'6' :    'Macrophage',   
'7':  'NK',   
'8':  'Monocyte',
'9':  'SMC',
'10':   'Treg', 
'11':  'DC',  
'12':  'Mast cell', 
'13':  'B cell'  ,
'14': 'Epithelium',  
'15':   'Plasma cell', 
'16':  'SG Epithelium',  
'17':   'Myoblast', 
'18':  'Cycling T',     
'19':  'Endothelium',
'20':  'Platelet',  
}

celltype_order = ['CD8 T', 'CD4 T', 'Treg', 'Cycling T', 'B cell', 'Plasma cell','NK',
                        'Monocyte', 'Macrophage', 'Neutrophil','DC', 'Mast cell', 'Platelet',
                        'Fibroblast', 'SMC', 'Myoblast', 'Epithelium', 'SG Epithelium', 'Endothelium']

adata.obs['celltype'] = adata.obs['leiden_0.6'].map(celltype_map)
adata.obs['celltype'] = adata.obs['celltype'].astype('category')
adata.obs['celltype'] = adata.obs.celltype.cat.reorder_categories(celltype_order)

with rc_context({'figure.figsize': (7, 7)}):
    sc.pl.umap(adata, color=['celltype'], 
               legend_loc = 'on data', 
               legend_fontsize=8, size=4
              )

adata.obs['major_celltype'] = adata.obs['celltype'].copy()
adata.obs['major_celltype']  = adata.obs['major_celltype'].astype('object')

adata.obs.loc[adata.obs.major_celltype.isin(['CD8 T', 'CD4 T', 'Treg', 'Cycling T','NK']), 'major_celltype'] = 'NK/T'
adata.obs.loc[adata.obs.major_celltype.isin(['Monocyte', 'Macrophage',  'Neutrophil', 'DC']), 'major_celltype'] = 'Myeloid cell'

with rc_context({'figure.figsize': (7, 7)}):
    sc.pl.umap(adata, color=['major_celltype'], 
               legend_loc = 'on data',
               legend_fontsize=8, size=4
              )

adata.write_h5ad('./data/RNA.merged.analysed.h5ad')
