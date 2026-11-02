### Import packages
import os
import numpy as np
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import classification_report
from scarches.models.scpoli import scPoli
import warnings
warnings.filterwarnings('ignore')
import argparse
import torch
import random

### Settings

sc.settings.set_figure_params(dpi=100, frameon=False)
sc.set_figure_params(dpi=100)
sc.set_figure_params(figsize=(3, 3))

def set_seed(seed: int = 42):
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

set_seed(42)

parser = argparse.ArgumentParser(description="Train and predict with sample data.")
parser.add_argument('sample', type=str, help="Sample identifier")
args = parser.parse_args()


####--- training settings ---###
ref_scRNA_counts_path = '/RNA.merged.analysed_raw_count.h5ad'
ref_scRNA_norm_path = '/RNA.merged.analysed.h5ad'

### ----- VisiumHD ---- ###
sample = args.sample 
visiumHD_path = f'/cell_{sample}_adata_raw.h5ad'


celltype_key = 'celltype'   
condition_key = 'sample'

early_stopping_kwargs = {
    "early_stopping_metric": "val_prototype_loss",
    "mode": "min",
    "threshold": 0,
    "patience": 20,
    "reduce_lr": True,
    "lr_patience": 13,
    "lr_factor": 0.1,
}

n_HVG = [5500]  
n_epochs =50
pretraining_epochs =40
eta =10
latent_dim_ls = [40]  
embedding_dim_list = [30]  

### Read in datas
ref_scRNA_counts = sc.read_h5ad(ref_scRNA_counts_path)
ref_scRNA_norm = sc.read_h5ad(ref_scRNA_norm_path)
query_HD_counts = sc.read_h5ad(visiumHD_path)


ref_norm_base = ref_scRNA_norm.raw.to_adata() if ref_scRNA_norm.raw is not None else ref_scRNA_norm

for n_hvg in n_HVG:
    ref_norm = ref_norm_base.copy()

    # 1) HVG
    sc.pp.highly_variable_genes(ref_norm, n_top_genes=n_hvg, batch_key='sample')
    hvg_names = ref_norm.var_names[ref_norm.var.highly_variable].tolist()

    ref_scRNA_counts_sub = ref_scRNA_counts[:, hvg_names].copy()

    # 2) overlap
    overlap_hvgs = list(np.intersect1d(ref_scRNA_counts_sub.var['gene_ids'], query_HD_counts.var['gene_ids']))
    ref_scRNA_counts_sub = ref_scRNA_counts_sub[:, ref_scRNA_counts_sub.var['gene_ids'].isin(overlap_hvgs)].copy()
    query_HD_counts_sub = query_HD_counts[:, query_HD_counts.var['gene_ids'].isin(overlap_hvgs)].copy()

    print(f"[n_hvg={n_hvg}] ref shape={ref_scRNA_counts_sub.shape}, query shape={query_HD_counts_sub.shape}, overlap={len(overlap_hvgs)}")

    # --- train ref model ---
    for latent_dim in latent_dim_ls:
        for embedding_dim in embedding_dim_list:
            ref_scRNA_model_dir = f"./level1_output/nHVG_{n_hvg}/{sample}_level1_ref_model/latent_dim_{latent_dim}_embedding_dim_{embedding_dim}"
            os.makedirs(ref_scRNA_model_dir, exist_ok=True)

            ref_scpoli_model = scPoli(
                adata=ref_scRNA_counts_sub,
                condition_keys=condition_key,
                cell_type_keys=celltype_key,
                embedding_dims=embedding_dim,
                latent_dim=latent_dim,
                recon_loss='zinb',
            )
            ref_scpoli_model.train(
                n_epochs=n_epochs,
                pretraining_epochs=pretraining_epochs,
                early_stopping_kwargs=early_stopping_kwargs,
                eta=eta,
            )
            ref_scpoli_model.save(ref_scRNA_model_dir, overwrite=True, save_anndata=True)

            # eval once
            ref_results_dict = ref_scpoli_model.classify(ref_scRNA_counts_sub, scale_uncertainties=True)
            preds = ref_results_dict[celltype_key]["preds"]
            classification_df = pd.DataFrame(
                classification_report(ref_scRNA_counts_sub.obs[celltype_key], preds, output_dict=True)
            ).transpose()
            classification_df.to_csv(f"{ref_scRNA_model_dir}/level1_ref_model_evaluation.csv")

            # --- query transfer ---
            save_path = f"./level1_output/nHVG_{n_hvg}/{sample}_pred_results/latent_dim_{latent_dim}_embedding_dim_{embedding_dim}"
            os.makedirs(save_path, exist_ok=True)

            query = query_HD_counts_sub.copy()
            query.obs[celltype_key] = "unlabeled"

            scpoli_query = scPoli.load_query_data(adata=query, reference_model=ref_scRNA_model_dir, labeled_indices=[])
            scpoli_query.train(n_epochs=n_epochs, pretraining_epochs=pretraining_epochs, eta=eta)

            results_dict = scpoli_query.classify(query, scale_uncertainties=True)

            data_latent = scpoli_query.get_latent(query, mean=True)
            adata_latent = sc.AnnData(data_latent)
            adata_latent.obs = query.obs.copy()
            adata_latent.obs["cell_type_pred"] = results_dict[celltype_key]["preds"].tolist()
            adata_latent.obs["cell_type_uncert"] = results_dict[celltype_key]["uncert"].tolist()

            adata_latent.obs[["cell_type_pred","cell_type_uncert"]].to_csv(f"{save_path}/prediction_meta.csv")
            adata_latent.write_h5ad(f"{save_path}/adata_latent.h5ad")
