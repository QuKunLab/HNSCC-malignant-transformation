import tifffile
from tifffile import TiffWriter, TiffReader

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
import os
import scipy

import warnings
warnings.filterwarnings('ignore')

plt.rcParams['pdf.fonttype'] = 42   
plt.rcParams["font.sans-serif"] ='Arial' 


import xml.etree.ElementTree as ET

import sys

sys.path.append('/mIF_TMA_analysis/')

from src.plot import  plot_roi
from src.plot.plot_roi import plot_roi_cell, plot_roi_region

sc.set_figure_params(dpi=100, figsize=(4, 4), dpi_save=500)

def add_rectangle(
    ax,
    crop_para
):
    '''xxx'''
    x, y, width, height = crop_para
    
    rect = plt.Rectangle((x, y), width, height, fill=False, edgecolor='red', linewidth=1)
    ax.add_patch(rect)

def plot_spatial(
    adata,
    sample,
    **kwargs
):
    '''xxx'''
    _ = sc.pl.embedding(adata[adata.obs['sample'].isin([sample])], basis='spatial', show=False, **kwargs)

    if not isinstance(_, list):
        _ = [_]
    
    for i in _:
        i.invert_yaxis()

ct2col = plot_roi.default_cell_type_color()

celltype_colors = {
    "CD8T":         "#0077B6",
    "PDCD1+ CD8T":  "#48CAE4",

    "CD4T":         "#C68A00",
    "CXCL13+ CD4T": "#FFD166",

    "B":             "#006D77",  
    "CXCR5+ B": "#00BFA6",  
    "IgG PC":       "#FF6D00",

    "CXCL12+ Fibro":"#C44536",
    "Macrophage":   "#6FAF9F",

    "Epithelium": "#aec7e8",
    
    'Myofibroblast': "#D46AA5"
}

ct2col.update(celltype_colors)

sample = '...'
adata = sc.read_h5ad(f'/data/{sample}/adata_{sample}_analysed.h5ad')
cell_mask = tifffile.imread(f'/data/cell_mask.tif')
boundary_mask = tifffile.imread(f'/data/boundary_mask.tif')
tiff = TiffReader(f'/ome_tiff/protein_stain_{sample}.ome.tif')

img = tiff.asarray()

metadata = tiff.ome_metadata  ## 
img.shape
root = ET.fromstring(metadata)
namespace = {'ome':'http://www.openmicroscopy.org/Schemas/OME/2016-06'}
channel_names = []
for channel in root.findall('.//ome:Channel', namespace):
    name = channel.get('Name')
    channel_names.append(name)

crop_para = [10900, 2220, 1400, 1400]

### zoom-in
colors = {
     "DAPI":   "#0000FF",
    "COL1A1":  "#2CA25F",  
    "CK":      "#AFAFAF",  
    "CD8":     "#00C2FF",  
    "CD4":     "#F4E04D",  
    "CD20":    "#008C95",
    "IgG":     "#E76F00",  
    "CXCL12":  "#E31A1C",  
    "FAP":     "#E15AAE",  
}

vmax = {
    'DAPI': 220,
    'COL1A1': 90,
    'CK': 100,
    'CD8': 55,
    'CD4': 70,
    'CD20': 180,
    'IgG': 180,
    'CXCL12': 70,
     'FAP': 30
}

vmin = {
    'DAPI': 8,
    'COL1A1': 8,
    'CK': 8,
    'CD8': 6,
    'CD4': 12,
    'CD20': 15,
    'IgG': 15,
    'CXCL12': 8,
     'FAP': 15
}

raw_channel_order = channel_names

plot_roi.plot_roi_signal_with_legend(
    img, 
    cell_mask=cell_mask,
    crop_para=crop_para,
    channels_color=colors,
    channel_order=raw_channel_order,
    vmin=vmin,
    vmax=vmax,
    bbox_to_anchor=(1.3, 1),
    thre_c2=[1],
    um_per_pixel=0.5,
    scale_bar_um=100,
    clip_circle=True
)

plt.show()

cellid2celltype_dict = plot_roi.get_cellid_to_celltype(adata, crop_para=crop_para, celltype_key='celltype')

plot_roi.plot_roi_mask_with_legend(cell_mask, boundary_mask, crop_para=crop_para, cellid2celltype=cellid2celltype_dict, 
                          groups=[ 'CD8T', 'PDCD1+ CD8T', 'CD4T', 'CXCL13+ CD4T', 'B', 'CXCR5+ B', 
                                   'IgG PC', 'CXCL12+ Fibro','Macrophage', 
                                   'Myofibroblast'
                                 ], 
                                   ct2col=ct2col,
                                   um_per_pixel=0.5,
                                   scale_bar_um=100,
                                   show_legend=True,
                                   clip_circle=True,
                                   bbox_to_anchor=(1.4, 1)
                                  )

plt.show()
