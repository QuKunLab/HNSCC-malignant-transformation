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
import scirpy as ir
from matplotlib.pyplot import rc_context

samples = [...] 

def build_anndata(sample):
    """
    loading scRNA and scTCR data at locaton: `./raw_data/%s' % sample`
    """
    _adata = sc.read_10x_mtx('./raw_data/{}/{}_filtered_feature_bc_matrix'.format(sample,sample) )
    _adata.obs['sample'] = sample
    print(_adata.shape)
    _adata_ir = ir.io.read_10x_vdj('./raw_data/{}/match_contigs.csv'.format(sample) )
    ir.pp.merge_with_ir(_adata, _adata_ir)
    ir.tl.chain_qc(_adata)
    print('-----%s-----' % sample)
    print(_adata.obs['chain_pairing'].value_counts())
    _adata = _adata[_adata.obs["multi_chain"] != "True", :]
    print('----------------------')
    return _adata

def build_anndatas(samples):
    adata_list = []
    for sample in samples:
        _adata = build_anndata(sample)
        adata_list.append(_adata)

    if len(adata_list)==1:
        return adata_list[0]
    else:
        return adata_list[0].concatenate(adata_list[1:], batch_key='batch', batch_categories=samples)

merged_adata_raw = build_anndatas(samples)  

merged_adata_raw.write('./raw_data/RNA_merged_raw.h5ad')
