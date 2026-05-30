library(scPred)
library(Seurat)
library(magrittr)
library(ggplot2)
library(patchwork)
library(tidyverse)


setwd('~/workspace/pan-cancer_ST/')


###--- Global Functions ---###
run_scPred_train <- function(seu_obj, is.processed=FALSE, classify_key='Niche', reduction='pca', model='svmRadial', reclassify = NULL){
    #preprocess raw counts data
    if (!is.processed){
        seu_obj <- seu_obj %>%
            NormalizeData() %>%
            FindVariableFeatures() %>%
            ScaleData() %>%
            RunPCA() %>%
            RunUMAP(dims = 1:30)
    }

    # Training classifiers with scPred
    ## 1. FeatureSpace
    seu_obj <- getFeatureSpace(object = seu_obj, pvar = classify_key, reduction = reduction)

    ## 2.training
    seu_obj <- trainModel(object = seu_obj, model = model, reclassify = reclassify)

    return(seu_obj)
}

    
process_our_data <- function(data.dir, niche.meta.file){
    # seurat object
    mtx <- Read10X(data.dir = data.dir)
    seurat_obj <- CreateSeuratObject(counts = mtx, min.cells = 0, min.features = 0, assay = 'RNA')

    # niche metadata
    niche_metadata <- read_csv(file = niche.meta.file, show_col_types = FALSE)
    
    niche_metadata <- niche_metadata |>
        column_to_rownames(var = colnames(niche_metadata)[1])

    seurat_obj <- AddMetaData(seurat_obj, metadata = niche_metadata, col.name = 'Niche')

    return(seurat_obj)
}


reference <- process_our_data(
    data.dir = './Cell_classifiers/scPred/data/all_10X_data/', 
    niche.meta.file = './Cell_classifiers/scPred/data/all_10X_data/_niche_metadata.csv'
)

#1.avNNet
### ---training--- ###
ref_trained_avNNet <- run_scPred_train(reference, is.processed=FALSE, classify_key = 'Niche', reduction = 'pca', model = 'avNNet')

saveRDS(object = ref_trained_avNNet, file = './Cell_classifiers/scPred/results/ref_all_10X_trained_avNNet.rds')

png('./Cell_classifiers/scPred/figures/avNNet_model_probabilities.png')
plot_probabilities(ref_trained_avNNet)
dev.off()

