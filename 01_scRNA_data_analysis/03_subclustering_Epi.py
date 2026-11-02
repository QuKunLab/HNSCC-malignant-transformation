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
# import scirpy as ir
from matplotlib.pyplot import rc_context
import warnings
warnings.filterwarnings('ignore')

sc.settings.set_figure_params(dpi=100, dpi_save=300,
                              facecolor='white', figsize=(4,4))

import PyComplexHeatmap
from PyComplexHeatmap import *

adata = sc.read_h5ad('./data/RNA.merged.analysed.h5ad')

adata_raw = adata.raw.to_adata()

cell_cycle_genes = pd.read_csv('./external_data/cell_cycle_genes.csv', header=0)

def func(x):
    if isinstance(x, float):
        return ''
    return x.strip()
cell_cycle_genes = cell_cycle_genes.applymap(func)

def preprocessing_step1(adata_tmp, n_top_genes=2000):
    sc.pp.highly_variable_genes(adata_tmp, n_top_genes=n_top_genes)
    adata_tmp.raw = adata_tmp
    adata_tmp = adata_tmp[:, adata_tmp.var.highly_variable]
    sc.pp.scale(adata_tmp, max_value=10)
    sc.tl.score_genes_cell_cycle(adata_tmp, s_genes=G1_S , g2m_genes=G2_M)
    sc.pp.regress_out(adata_tmp, ['S_score','G2M_score','total_counts', 'pct_counts_mt'])  

    sc.pp.scale(adata_tmp, max_value=10)
    sc.tl.pca(adata_tmp, svd_solver='arpack')
    sc.pl.pca_variance_ratio(adata_tmp, log=True, n_pcs=50)
    return adata_tmp


def preprocessing_step2(adata_tmp, key, n_pcs, n_neighbors, sigma=0.1, lamb=1.0, theta=1.0):
    sce.pp.harmony_integrate(adata_tmp, key='sample', sigma=sigma, lamb=lamb, theta=theta)
    sc.pp.neighbors(adata_tmp, use_rep='X_pca_harmony', n_pcs=n_pcs, n_neighbors=n_neighbors)
    sc.tl.umap(adata_tmp)
    for i in [0.2, 0.3, 0.4, 0.5]:  
        sc.tl.leiden(adata_tmp, resolution=i, key_added=key+'_leiden_'+str(i)) 
    return adata_tmp

cellcycle_geneset = pd.read_excel('./external_data/cell_cycle_geneset.xlsx', sheet_name='Sheet1', skiprows=6)
G1_S = [ i.strip() for i in cellcycle_geneset['G1/S'][~cellcycle_geneset['G1/S'].isna()].tolist()]
G2_M = [ i.strip() for i in cellcycle_geneset['G2/M'][~cellcycle_geneset['G2/M'].isna()].tolist()]

adata_Epi = adata_raw[adata_raw.obs.celltype.isin(['Epithelium'])].copy() 

print(adata_Epi.shape)

adata_Epi = preprocessing_step1(adata_Epi, n_top_genes=3000)

adata_Epi = preprocessing_step2(adata_Epi,  key='Epi', n_pcs=50, n_neighbors=15,sigma=0.1, lamb=1.0, theta=0.3) 

sc.pl.umap(adata_Epi, color=['Epi_leiden_0.2'],legend_loc='on data')  

Epi_map = {
'0': 'SBSN+ Epi',
'1': 'COL17A1+ Epi',
'2': 'GPX2+ Epi',
'3': 'MKI67+ Epi',
'4': 'COL17A1+ Epi'
}

order = ['SBSN+ Epi',  'COL17A1+ Epi', 'MKI67+ Epi',  'GPX2+ Epi']

adata_Epi.obs['Epi_celltype'] = adata_Epi.obs['Epi_leiden_0.2'].map(Epi_map)
adata_Epi.obs['Epi_celltype'] = adata_Epi.obs['Epi_celltype'].astype('category')
adata_Epi.obs['Epi_celltype'] = adata_Epi.obs.Epi_celltype.cat.reorder_categories(order)

sc.pl.dotplot(adata_Epi, var_names = ['SBSN',  'KRT6B', 
'COL17A1', 'IGFBP7',
'MKI67','TOP2A',
'GPX2', 'NTS', 'GSTA1'], groupby='Epi_celltype', standard_scale='var')

adata_Epi.uns['Epi_celltype_colors'] = ['#3EC0AA','#FFC600','#3E6DB6','#EE1133']

sc.pl.umap(adata_Epi, color=['Epi_celltype'], size=5)

adata_Epi.write_h5ad('./data/adata_Epi.h5ad')


sns.set_style('ticks')

adata_Epi = adata_Epi.raw.to_adata()

gene_list = [
    "IL1B", "IL13", "HLA-DPA1", "HLA-DQA1", "HLA-DQA2", "HLA-DPB1",
    "CXCL9", "CXCL10", "HLA-B", "HLA-C", "IL2", "CD40LG",
    "VEGFB", "HLA-A", "CD70", "IL10", "IL1A", "IFNG",
    "HLA-DQB2", "HLA-DQB1", "CCL5", "HLA-DRA", "HLA-DRB5", "HLA-DRB1"
]

expr_df = adata_Epi[:,gene_list].to_df()
expr_df['Epi_celltype'] = adata_Epi.obs['Epi_celltype'].copy()
expr_df = expr_df.groupby('Epi_celltype').mean()

top_anno = pd.DataFrame({'Epi_celltype':expr_df.index})
top_anno.index = top_anno.Epi_celltype.values

celltyep_map = dict(zip(adata_Epi.obs['Epi_celltype'].cat.categories, adata_Epi.uns['Epi_celltype_colors']))

plt.figure(figsize=(2, 10))
top_annotation = HeatmapAnnotation(
    Epi_celltype = anno_simple(top_anno['Epi_celltype'], add_text=False, legend_kws={'frameon':False}, colors= celltyep_map, legend=True),  axis=1)
cm2 = ClusterMapPlotter(data=expr_df.T, 
                       z_score=0,
                       vmax=1, vmin=-1, 
                        top_annotation =top_annotation,
                       cmap='viridis',   #BuPu  PRGn_r RdYlBu_r
                       show_rownames=True, show_colnames=True,
                       row_cluster=False, col_cluster=False,  legend=True,
                        legend_width=2, linecolor='white',linewidths=0.8, label='Normalized expression')
plt.savefig('./figures/Heatmap_Epi_feature.pdf', dpi=300, bbox_inches='tight')

plt.show()


keratinocyte_pathway = pd.read_excel('./pathway/keratinocyte_genesets.xlsx', header=None, index_col=0)
keratinocyte_pathway = keratinocyte_pathway.T
keratinocyte_pathway = keratinocyte_pathway.iloc[1:,]

print(keratinocyte_pathway)

for item in keratinocyte_pathway.columns:
    print(item)
    genes = keratinocyte_pathway[item][~keratinocyte_pathway[item].isna()].tolist()
    print(len(genes))
    sc.tl.score_genes(adata_Epi, gene_list = genes, score_name = item)
    print('---------------------------------------------------------------------')


stem_pathway = pd.read_excel('./18 干性相关功能通路.xlsx', header=None) #GOBP  HALLMARK   KEGG\
stem_pathway.index = stem_pathway[0]
stem_pathway = stem_pathway.iloc[:,1:]
stem_pathway = stem_pathway.T
stem_pathway.shape
stem_pathway = stem_pathway.loc[:, [len(i) > 10 for i in stem_pathway.columns]]
print(stem_pathway.shape)

for item in ['GOBP_REGULATION_OF_STEM_CELL_POPULATION_MAINTENANCE']:#stem_pathway.columns:
    print(item)
    genes = stem_pathway[item][~stem_pathway[item].isna()].tolist()
    print(len(genes))
    sc.tl.score_genes(adata_Epi, gene_list = genes, score_name = item)
    print('---------------------------------------------------------------------')

pathways = pd.read_excel('./14 MSigDB HALLMARK KEGG.xlsx', header=None) #GOBP  HALLMARK   KEGG\

pathways.index = pathways[0]
pathways = pathways.iloc[:,1:]
pathways = pathways.T
pathways.shape
pathways = pathways.loc[:, [len(i) > 10 for i in pathways.columns]]
print(pathways.shape)

for item in pathways.columns:
    print(item)
    genes = pathways[item][~pathways[item].isna()].tolist()
    print(len(genes))
    sc.tl.score_genes(adata_Epi, gene_list = genes, score_name = item)
    print('---------------------------------------------------------------------')

Ferroptosis = pd.read_csv('./WP_FERROPTOSIS.v2023.1.Hs.tsv', sep='\t')
NRF2 = pd.read_csv('./WP_NRF2_PATHWAY.v2023.1.Hs.tsv', sep='\t')

Ferroptosis = Ferroptosis['WP_FERROPTOSIS'][16].strip(',').split(',')
NRF2 = NRF2['WP_NRF2_PATHWAY'][16].strip(',').split(',')

sc.tl.score_genes(adata_Epi, gene_list=Ferroptosis, score_name='FERROPTOSIS')  
sc.tl.score_genes(adata_Epi, gene_list=NRF2, score_name='NRF2_PATHWAY')  

# https://panglaodb.se/markers.html?cell_type=%27Basal%20cells%27
Basal_cells = pd.read_csv('./Basal_cells_PanglaoDB_markers_27_Mar_2020.tsv', sep='\t')
Basal_cells = Basal_cells.loc[Basal_cells['cell type']=='Basal cells',]

sc.tl.score_genes(adata_Epi, gene_list = Basal_cells['official gene symbol'].tolist(), score_name = 'Basal cell signature')

pEMT = pd.read_excel('./12 2022 HNSCC EMT 相关marker（PMID 36076232).xlsx')

sc.tl.score_genes(adata_Epi, gene_list = pEMT['pEMT (Puram et al.)'][~pEMT['pEMT (Puram et al.)'].isna()].tolist(), score_name = 'pEMT signature')

selected_pathways_list = [

'GOBP_KERATINOCYTE_DIFFERENTIATION',
'HALLMARK_TGF_BETA_SIGNALING',
'pEMT signature',
'Basal cell signature',
'HALLMARK_HYPOXIA',
'HALLMARK_TNFA_SIGNALING_VIA_NFKB',
'HALLMARK_ANGIOGENESIS',
'HALLMARK_INFLAMMATORY_RESPONSE',
'KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION',
'HALLMARK_G2M_CHECKPOINT',
'HALLMARK_E2F_TARGETS',
'HALLMARK_FATTY_ACID_METABOLISM',
'HALLMARK_OXIDATIVE_PHOSPHORYLATION',
'KEGG_GLUTATHIONE_METABOLISM',
'FERROPTOSIS',    
'NRF2_PATHWAY',
'GOBP_REGULATION_OF_STEM_CELL_POPULATION_MAINTENANCE'
]

def prepare_heatmap_data(pathways, key=['Epi_celltype']):  
    df = adata_Epi.obs[pathways+key] 
    df = df.groupby(key).mean().reset_index()  
    df_heatmap = df.iloc[:,len(key):].T
    df_color = df[key] 
    for i in key:
        df_color[i] = df_color[i] .astype('object')
    return df_heatmap, df_color

df_heatmap, df_color = prepare_heatmap_data(pathways=selected_pathways_list, key=['Epi_celltype'])

Epi_celltype_map = {'SBSN+ Epi': '#3EC0AA', 'COL17A1+ Epi': '#FFC600', 'MKI67+ Epi': '#3E6DB6', 'GPX2+ Epi': '#EE1133'}   

plt.figure(figsize=(2.5, 9))
top_annotation = HeatmapAnnotation(Epi_celltype = anno_simple(df_color['Epi_celltype'], add_text=False,legend_kws={'frameon':False}, colors=Epi_celltype_map), axis=1)
cm = ClusterMapPlotter(data=df_heatmap, z_score=0,
                       vmax=1.5,vmin=-1.5, 
                       top_annotation=top_annotation,   
                       cmap='RdBu_r',show_rownames=True, 
                       row_cluster=False, col_cluster=False, 
                       linecolor='white',linewidths=2, label='pathway score')
plt.savefig('./figures/Heatmap_of_Epi_celltype_pathway_score.pdf', dpi=300, bbox_inches='tight')

plt.show()
