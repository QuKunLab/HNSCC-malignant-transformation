library(infercnv)
raw_counts_matrix = './inputs/raw_counts.txt'
annotations_file = './inputs/meta.txt'  
out_dir_1 = './outputs/'
out_dir_2 = './outputs/figures'

ref_celltypes = c('CD8 T', 'CD4 T', 'Treg', 'Cycling T', 'B cell', 'Plasma cell', 'NK', 'Monocyte', 'Neutrophil', 'Macrophage', 'DC', 'Mast cell')

print('********************* CreateInfercnvObject start ***********************')

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=raw_counts_matrix,  
                                    annotations_file=annotations_file,     
                                    delim="\t",
                                    gene_order_file='./hg38_gen_pos.deduplicates.ordered.txt',
                                    ref_group_names= ref_celltypes)  
print('********************* CreateInfercnvObject finished ***********************')

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,   
                             out_dir=out_dir_1,
                             cluster_by_groups=TRUE,  
                             plot_steps=FALSE,
                             denoise=TRUE,
                              HMM=FALSE,
                             num_threads=25,
                             cluster_references=FALSE,
                             output_format = 'pdf',
                              png_res = 300)
save.image(file = paste0(out_dir_1,'/inferCNV.result.RData'))

print('********************* infercnv::run finished ***********************')
infercnv::plot_cnv(infercnv_obj, 
                   out_dir=out_dir_2, 
                   output_filename='infercnv_heatmap',   
                   x.range="auto",
                   x.center=1,
                   title = '', 
                   cluster_by_groups = T,      
                   color_safe_pal = FALSE,
                   output_format='pdf',
                   png_res = 300,
                  cluster_references=FALSE)

print('********************* infercnv::plot_cnv finished ***********************')
print('********************* over! ***********************')