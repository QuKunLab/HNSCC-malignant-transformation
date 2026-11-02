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

adata_all = sc.read_h5ad('./mIF_TMA_analysis/data/adata_analyased_array.h5ad')

sample = '...'

cell_mask = tifffile.imread(f'/h5ad_out/{sample}/cell_mask.tif')
boundary_mask = tifffile.imread(f"/h5ad_out/{sample}/boundary_mask.tif")
tiff = TiffReader(f'/core_ome/{sample}.ome.tif')
img = tiff.asarray()

Epi_groups = ['GPX2+ Epi', 'SBSN+ Epi', 'COL17A1+ Epi', 'MKI67+ Epi']

crop_para=[  int(round(1594/0.325)),
             int(round(1107/0.325)),
             int(round(700/0.325)),
             int(round(700/0.325))
          ]

raw_channel_order = adata_all.var_names.tolist()

colors = {'DAPI':'#0000ff',  
          'PAN-CK':'#999999',  
          "KI67":        "#00ffff",
    "TGM1":        "#3EC0AA",
    "NRF2(abcam)": "#ff0000",
    "COL17A1":     "#FFC600",
          "ITGB1":     "#000000",
                                 }

vmax = {
    'DAPI': 70,
    'PAN-CK': 140,
    'KI67': 185,
    'TGM1': 140,
    'NRF2(abcam)': 35,
    'COL17A1': 45,
    'ITGB1': 40
}

vmin = {
    'DAPI': 3,
    'PAN-CK': 18,
    'KI67': 70,
    'TGM1': 20,
    'NRF2(abcam)': 30,
    'COL17A1': 8,
    'ITGB1': 8
}

plot_roi.plot_roi_signal_with_legend(
    img, 
    cell_mask=cell_mask,
    crop_para=crop_para,
    channels_color=colors,
    channel_order=raw_channel_order,
    vmin=vmin,
    vmax=vmax,
    bbox_to_anchor=(1.4, 1),
    um_per_pixel=0.325,
    scale_bar_um=100,
    clip_circle=True,
)

plt.show()

roi_adata = adata_all[adata_all.obs['sample'].isin([sample])]
roi = roi_adata[(roi_adata.obs.pixel_X > crop_para[0]) & (roi_adata.obs.pixel_X < crop_para[0] + crop_para[2]) &
            (roi_adata.obs.pixel_Y > crop_para[1]) & (roi_adata.obs.pixel_Y < crop_para[1] + crop_para[3])]

cellid2celltype_dict = dict(zip(roi.obs["cell_id_label"].astype(int), roi.obs["celltype"]))

ct2col = plot_roi.default_cell_type_color()

plot_roi.plot_roi_mask_with_legend(cell_mask, boundary_mask, crop_para=crop_para, cellid2celltype=cellid2celltype_dict, 
                          groups=Epi_groups, 
                                   ct2col=ct2col,
                                   um_per_pixel=0.325,
                                   scale_bar_um=100,
                                   show_legend=True,
                                   clip_circle=True,
                                   bbox_to_anchor=(1.4, 1)
                                  )

plt.show()
