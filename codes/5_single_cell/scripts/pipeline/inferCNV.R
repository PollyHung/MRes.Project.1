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


run_inferCNV <- function(seurat_obj, ## path to seurat object
                         output_dir, ## path to result directory 
                         gene_order_file ## path to gene order file 
                         ){
  ## change directory to output_dir 
  setwd(output_dir)
  
  ## read in 
  BSC <- readRDS(seurat_obj)

  ## generating annotations_file for each run uniquely 
  sub_metadata <- BSC@meta.data
  sub_metadata$barcode <- row.names(sub_metadata)
  sub_metadata$newCelliD <-  paste0(sub_metadata$sample, "_", sub_metadata$CellType)
  BSC@meta.data <- sub_metadata
  cell_id_table <- table(sub_metadata$newCelliD)
  single_occurrences <- cell_id_table[cell_id_table == 1]
  sub_metadata2 <- sub_metadata[-which(sub_metadata$newCelliD %in% names(single_occurrences)), ]
  newCelliDs <- sub_metadata2$newCelliD
  BSC <- subset(BSC, subset = newCelliD %in% newCelliDs) ## resulting in a matrix with single occurances removed 
  write.table(sub_metadata2[, c("barcode", "newCelliD")], file = "total_cellAnnotations.txt", 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  ## generating raw_counts_matrix
  sub_matrix <- as.matrix(GetAssayData(BSC, assay = "RNA"))
  
  ## gene_order_file 
  ## this is the same across different runs, already generated elsewhere for global use 
  
  ## getting ref_group_names 
  cancers_id <-  metadata$newCelliD[grepl(pattern="*_Cancer", metadata$newCelliD)] %>% unique
  ref_group_names <- setdiff(unique(metadata$newCelliD), cancers_id)
  
  ## create infercnv_obj
  infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = sub_matrix, 
                                       annotations_file = "total_cellAnnotations.txt", 
                                       gene_order_file = gene_order_file, 
                                       ref_group_names = ref_group_names) 
  
  ## run infercnv
  infercnv_obj = infercnv::run(infercnv_obj, 
                               cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                               out_dir=output_dir,  # dir is auto-created for storing outputs
                               cluster_by_groups=TRUE,   # cluster
                               denoise=TRUE,
                               sd_amplifier = 3, 
                               noise_logistic = TRUE, 
                               BayesMaxPNormal=0.5, 
                               HMM=TRUE, 
                               num_threads = 20, 
                               write_expr_matrix = TRUE)
  
}



