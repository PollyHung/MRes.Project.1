library(Seurat)
library(SeuratObject)
library(ggplot2)
library(magrittr)
library(dplyr)
library(patchwork)
library(tidyverse)
library(Matrix)


add_CNV <- function(cnv_file, ## path to cnv_file, example: "~/MRes_project_1/docs/SC_SEQ/ova/lambrechtslab/inferCNV/72_hr_runtime/broad_values_by_arm.txt"
                    seurat_obj, ## path to seurat_obj, example: "~/MRes_project_1/docs/SC_SEQ/ova/lambrechtslab/results/umapped_obj.rds"
                    chr_arm ## the chromosome arm you wish to get, starts with chr, example: chr2q
                    ){
  
  ## read in 
  cnv <- read.table(cnv_file, sep="\t", header=TRUE, row.names = 1)
  cnv$sample_id <- rownames(cnv)
  cnv_sub <- cnv %>% dplyr::select(sample_id , all_of(chr_arm)) 
  BSC <- readRDS(seurat_obj)
  
  ## update metadata
  metadata <- BSC@meta.data
  metadata <- left_join(metadata, cnv_sub, by="sample_id")
  metadata[paste0("round_", chr_arm)]<- round(metadata[chr_arm])
  rownames(metadata) <- metadata$cell
  metadata <- metadata[rownames(BSC@meta.data), ]
  
  BSC@meta.data <- metadata
}


