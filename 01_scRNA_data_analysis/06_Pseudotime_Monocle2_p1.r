library(monocle)  
library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyr) 
library(ggridges)
library(RColorBrewer)
library(scales)
library(ggsci)
library(patchwork)

library(ggrastr)

warnings('off')

options(list(ggplot2.discrete.fill = c('#a6cee3','#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#ffff99','#b15928')))

New_theme = theme_classic() + theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1,size = 12, colour = '#000000'), 
      axis.text.y = element_text(,size = 15, colour = '#000000'), 
      axis.title = element_text(size = 8),
    plot.title = element_text(size=8 ,  hjust = 0.5),
     panel.grid=element_blank())

theme_set(New_theme)

counts = read.csv('./input/raw_counts.csv', row.names=1, ,check.names =F )
data <- as(as.matrix(counts), 'sparseMatrix') 
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
meta = read.csv('./input/meta.csv', row.names=1, ,check.names =F )
pd <- new('AnnotatedDataFrame', data = meta)

cds <- newCellDataSet(data,    
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.1, 
                      expressionFamily = negbinomial.size())  

cds <- estimateSizeFactors(cds)  
cds <- estimateDispersions(cds) 

cds <- detectGenes(cds, min_expr = 0.1 )   
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))
cds <- cds[expressed_genes,]  

## Step1: choosing genes that define progress
 ordering_genes = read.csv('./input/Fibro_subtypes_DEG.csv', header = F) ## using firboblast subtypes top DEG as ordering genes
ordering_genes = ordering_genes[['V1']]
print(length(ordering_genes))
monocle_cds <- setOrderingFilter(cds, ordering_genes)

## step2
monocle_cds <- reduceDimension(monocle_cds, max_components = 2, 
                                 reduction_method = 'DDRTree')  

## step3
monocle_cds <- orderCells(monocle_cds) 

options(repr.plot.height=6, repr.plot.width=6)
celltype_colour=c('#00C9A7','#F3C5FF','#845EC2'  )
plot_cell_trajectory(monocle_cds, color_by = "Fibro_celltype_refined",show_tree=T, show_branch_points=F, size=1) + 
scale_color_manual(values = celltype_colour)  + 
theme(legend.text = element_text(size = 8)) + guides(color = guide_legend(override.aes = list(size = 5)))

options(repr.plot.height=6, repr.plot.width=6)
plot_cell_trajectory(monocle_cds, color_by = "Pseudotime",show_tree=T, show_branch_points=F, size=2) + 
scale_fill_gradientn(colors = brewer.pal(9, "Purples"))
