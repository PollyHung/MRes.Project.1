library(Seurat)
library(ggplot2)
library(magrittr)
library(dplyr)
library(patchwork)
library(Signac)
library(SingleR)
library(celldex)
library(stringr)
library(infercnv)
library(Matrix)
library(SeuratObject)
library(SeuratDisk)
library(SeuratData)
library(infercnv)


T_cells <- readRDS("/rds/general/user/ph323/ephemeral/GSE180661/cancerTcells.rds")


## load in seurat object for count matrix
sub_matrix <- GetAssayData(T_cells, layer="data")
metadata <- T_cells@meta.data
metadata$cell_subtype[which(is.na(metadata$cell_subtype))] <- "cancer"

sub_metadata <- metadata %>% dplyr::select(cell, cell_subtype) %>% 
  dplyr::rename(barcode=cell, newCelliD=cell_subtype) %>% na.omit(cell_subtype)

## remove single occurence cells 
cell_id_table <- table(sub_metadata$newCelliD)
single_occurrences <- cell_id_table[cell_id_table == 1]
sub_metadata <- sub_metadata %>% dplyr::filter(!newCelliD %in% names(single_occurrences))
newCelliDs <- sub_metadata$newCelliD
sub_metadata <- sub_metadata %>% dplyr::filter(barcode %in% colnames(sub_matrix))
sub_matrix <- sub_matrix[, sub_metadata$barcode]

head(sub_metadata)
print(sub_matrix[1:5, 1:5])

## write out file 
write.table(sub_metadata, "~/MRes_project_1/docs/SC_SEQ/ova/GSE180661/inferCNV/T_cell/T_cell_annotation_file.txt", sep="\t", 
            col.names = FALSE, row.names = FALSE, quote=FALSE)
saveRDS(sub_matrix, "/rds/general/user/ph323/ephemeral/GSE180661/sub_matrix.rds")

## reference group names 
cancers_id = sub_metadata$newCelliD[grep(pattern="dysfunc", sub_metadata$newCelliD)] %>% unique
ref_group_names <- setdiff(unique(sub_metadata$newCelliD), cancers_id)

head(cancers_id)
head(ref_group_names)

#load("/rds/general/user/ph323/ephemeral/GSE131907/inferCNV.RData")
## file 5: define output directory for this sample 
output_dir <- "~/MRes_project_1/docs/SC_SEQ/ova/GSE180661/inferCNV/T_cell/"
gene_order_file <- "~/MRes_project_1/codes/6_single_cell/reference/gof.txt"
meta.data <- "~/MRes_project_1/docs/SC_SEQ/ova/GSE180661/inferCNV/T_cell/T_cell_annotation_file.txt"

if(FALSE){  
  ## inferCNV object 
  infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = sub_matrix, 
                                       annotations_file = meta.data, 
                                       gene_order_file = gene_order_file, 
                                       ref_group_names = ref_group_names) ## set to various normal cell types 
  
  infercnv_obj = infercnv::run(infercnv_obj, 
                               cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                               out_dir=output_dir,  # dir is auto-created for storing outputs
                               cluster_by_groups=TRUE,   # cluster
                               HMM_report_by = "cell", 
                               denoise=TRUE,
                               sd_amplifier = 3, 
                               noise_logistic = TRUE, 
                               BayesMaxPNormal=0.5, 
                               HMM=TRUE, 
                               num_threads = 20, 
                               write_expr_matrix = TRUE)
}