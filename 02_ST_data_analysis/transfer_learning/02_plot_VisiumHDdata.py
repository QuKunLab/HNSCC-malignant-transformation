import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
import seaborn as sns
import os
import scipy

# import warnings
# warnings.filterwarnings('ignore')

plt.rcParams['pdf.fonttype'] = 42  
plt.rcParams["font.sans-serif"] ='Arial' 


from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all"

sc.settings.verbosity = 1
sc.set_figure_params(figsize=(4, 4), dpi=100, dpi_save=300)
sns.set_style("ticks")
sample = '...'
adata = sc.read_h5ad(f'./VisiumHD/data/sample_data/cell_{sample}_adata_analysed.h5ad')
groups = [
        'GPX2+ Epi','MKI67+ Epi','SBSN+ Epi','COL17A1+ Epi'
    ]
fig, ax = plt.subplots(figsize=(4, 4), dpi=100)
sc.pl.embedding(adata, basis='spatial', color=['celltype'], size=2,
                groups=groups,
                ax=ax,  show=False )
ax.invert_yaxis()


x=12000;y=4300
width=4000 
height=4000
rect = plt.Rectangle((x, y), width, height, fill=False, edgecolor ='black',linewidth=1)
ax.add_patch(rect)
ax.set_aspect('equal', adjustable='box')

x=19200;y=12400
width=1500 
height=1500
rect = plt.Rectangle((x, y), width, height, fill=False, edgecolor ='black',linewidth=1)
ax.add_patch(rect)
ax.set_aspect('equal', adjustable='box')

# fig.savefig('/path/to/save/figure.pdf', bbox_inches='tight', dpi=300)
plt.show()
x=12000;y=4300
width=4000 
height=4000

tmp = adata[(adata.obs.X > x)&(adata.obs.X < x+width)&(adata.obs.Y >y)&(adata.obs.Y <y+height)].copy()
_ = sc.pl.embedding(tmp, basis='spatial', color=[
    'celltype',
     'SBSN', 'GPX2', 'MKI67', 'COL17A1'], 
                    size=40,
                groups=groups,
                    cmap='RdPu', 
                    vmin=[None, 0, 0.02, 0, 0],
                    vmax=[None, 'p96.5', 'p97', 'p99', 'p98'], 
                    title=[None, 'SBSN', 'GPX2', 'MKI67', 'COL17A1'], ncols=5, legend_loc=None,
                show=False)
for i in _:
    i.invert_yaxis()

# plt.savefig('/path/to/save/figure.pdf', bbox_inches='tight', dpi=300)
plt.show()
x=19200;y=12400
width=1500 
height=1500

tmp = adata[(adata.obs.X > x)&(adata.obs.X < x+width)&(adata.obs.Y >y)&(adata.obs.Y <y+height)].copy()
_ = sc.pl.embedding(tmp, basis='spatial', color=[
    'celltype',
     'SBSN', 'GPX2', 'MKI67', 'COL17A1'], 
                    size=100,
                groups=groups,
                    cmap='RdPu', 
                    vmin=[None, 0, 0, 0, 0.04],
                    vmax=[None, 'p97', 'p97', 'p97', 0.1], 
                    title=[None, 'SBSN', 'GPX2', 'MKI67', 'COL17A1'], ncols=5, legend_loc=None,
                show=False)
for i in _:
    i.invert_yaxis()

# plt.savefig('/path/to/save/figure.pdf', bbox_inches='tight', dpi=300)
plt.show()
sample = '...'
adata = sc.read_h5ad(f'./VisiumHD/data/sample_data/cell_{sample}_adata_analysed.h5ad')
groups = [
        'GPX2+ Epi','MKI67+ Epi','SBSN+ Epi','COL17A1+ Epi'
    ]
fig, ax = plt.subplots(figsize=(4, 4), dpi=100)
sc.pl.embedding(adata, basis='spatial', color=['celltype'], size=2,
                groups=groups,
                ax=ax,  show=False )
ax.invert_yaxis()


x=2500;y=10800
width=3000 
height=3000
rect = plt.Rectangle((x, y), width, height, fill=False, edgecolor ='black',linewidth=1)
ax.add_patch(rect)
ax.set_aspect('equal', adjustable='box')

x=13800;y=13400
width=2000 
height=2000
rect = plt.Rectangle((x, y), width, height, fill=False, edgecolor ='black',linewidth=1)
ax.add_patch(rect)
ax.set_aspect('equal', adjustable='box')

# fig.savefig('/path/to/save/figure.pdf', bbox_inches='tight', dpi=300)
plt.show()
x=2500;y=10800
width=3000 
height=3000

tmp = adata[(adata.obs.X > x)&(adata.obs.X < x+width)&(adata.obs.Y >y)&(adata.obs.Y <y+height)].copy()

_ = sc.pl.embedding(tmp, basis='spatial', color=[
    'celltype',
    'SBSN', 'GPX2', 'MKI67', 'COL17A1'], 
                    size=40,
                groups=groups,
                    cmap='RdPu', 
                    vmin=[None, 0, 0.3, 0, 0.1],
                    vmax=[None, 'p97', 'p97', 'p96.5', 'p95'], 
                    title=[None, 'SBSN', 'GPX2', 'MKI67', 'COL17A1'], ncols=5, legend_loc=None,
                show=False)
for i in _:
    i.invert_yaxis()

# plt.savefig('/path/to/save/figure.pdf', bbox_inches='tight', dpi=300)
plt.show()
x=13800;y=13400
width=2000 
height=2000


tmp = adata[(adata.obs.X > x)&(adata.obs.X < x+width)&(adata.obs.Y >y)&(adata.obs.Y <y+height)].copy()

_ = sc.pl.embedding(tmp, basis='spatial', color=[
    'celltype',
    'SBSN', 'GPX2', 'MKI67', 'COL17A1'], 
                    size=80,
                groups=groups,
                    cmap='RdPu', 
                    vmin=[None, 0, 0, 0, 0],
                    vmax=[None, 'p97', 'p97', 'p97', 'p97'], 
                    title=[None, 'SBSN', 'GPX2', 'MKI67', 'COL17A1'], ncols=5, legend_loc=None,
                show=False)
for i in _:
    i.invert_yaxis()

# plt.savefig('/path/to/save/figure.pdf', bbox_inches='tight', dpi=300)
plt.show()
