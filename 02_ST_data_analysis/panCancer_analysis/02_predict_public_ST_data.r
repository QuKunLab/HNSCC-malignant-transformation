### Load packages
library(scPred)
library(Seurat)
library(magrittr)
library(ggplot2)
library(patchwork)
library(tidyverse)


### Settings
setwd('~/workspace/pan-cancer_ST/')
###--- Global Functions ---###
run_scPred_predict <- function(query, reference, is.processed=FALSE, threshold=0.55, recompute_alignment=TRUE){
    # using the same normalization method for both the reference and the query datasets.
    if (!is.processed){
        query <- query %>%
            NormalizeData()
    }

    # predict (Aligned data using harmony)
    query <- scPredict(new = query, reference = reference, threshold = threshold, recompute_alignment = recompute_alignment)

    return(query)
}

# models: avNNet
ref_avNNet <- readRDS('./Cell_classifiers/scPred/results/ref_all_10X_trained_avNNet.rds')
model_avNNet <- get_scpred(ref_avNNet)

tissue_type <- 'Tumor'
out_base_dir <- './Cell_classifiers/scPred/results/prediction_avNNet/'  

cancers_to_predict <- c('BLCA', 'CM', 'PDAC', 'PRAD', 'COAD', 'BRCA', 'HCC', 'KIRC', 'GBM', 'NPC', 'LUSC', 'LUAD', 'CSCC', 'ESCC')
input_base_dir <- './Cell_classifiers/scPred/data/other_cancer_type_data/'

for (cancer_type in cancers_to_predict) {
    print(cancer_type)
    
    samples_ls <- list.files(paste0(input_base_dir, cancer_type))
    for (sp in samples_ls) {
        print(sp)
        

        out_path <- file.path(out_base_dir, cancer_type)
        if (!dir.exists(out_path)) {
            dir.create(out_path, recursive = TRUE)   
        }

        fout_path <- file.path(out_path, tissue_type, sp)
        if (!dir.exists(fout_path)) {
            dir.create(fout_path, recursive = TRUE)  
        }
        
        data_dir <- file.path(input_base_dir, cancer_type, sp)
        query <- CreateSeuratObject(counts = Read10X(data.dir = data_dir), min.cells = 0, min.features = 0, assay = 'RNA')
        query <- run_scPred_predict(query = query, is.processed = FALSE, reference = model_avNNet)   
        
        query@meta.data$barcodes <- colnames(query)
        write_csv(query@meta.data, file = file.path(fout_path, 'niche_prediction.csv'))
    }
}


