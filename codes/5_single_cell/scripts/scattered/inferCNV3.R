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

save.rds <- "~/MRes_project_1/docs/SC_SEQ/lung/raw/lambrechtslab/lung_singleOccurRemoved.rds"
meta.data <- "~/MRes_project_1/docs/SC_SEQ/lung/raw/lambrechtslab/inferCNV/total_cellAnnotations.txt"
gene_order_file <- "~/MRes_project_1/codes/6_single_cell/reference/gof.txt"

## load in seurat object for count matrix
BSC <- readRDS(save.rds) 
sub_matrix <- as.matrix(GetAssayData(BSC, assay = "RNA"))
metadata <- read.table(meta.data, sep="\t") 
colnames(metadata) <- c("barcodes", "newCelliD")

## load in metadata 
#metadata <- read.csv(meta.data, row.names = "Cell") 
#metadata <- metadata[rownames(BSC@meta.data), ] 
#BSC <- AddMetaData(BSC, metadata = metadata, col.name = colnames(metadata))

## annotation file
#sub_metadata <- BSC@meta.data
#sub_metadata$barcode <- row.names(sub_metadata)
#sub_metadata$newCelliD = paste0(sub_metadata$sample, "_", sub_metadata$CellType)
#BSC@meta.data <- sub_metadata

## remove single occurence cells 
#cell_id_table <- table(sub_metadata$newCelliD)
#single_occurrences <- cell_id_table[cell_id_table == 1]
#sub_metadata2 <- sub_metadata[-which(sub_metadata$newCelliD %in% names(single_occurrences)), ]
#newCelliDs <- sub_metadata2$newCelliD
## re subset the seurat object 
#BSC <- subset(BSC, subset = newCelliD %in% newCelliDs)

#saveRDS(BSC, "~/MRes_project_1/docs/SC_SEQ/lung/raw/lambrechtslab/lung_singleOccurRemoved.rds")

#write.table(sub_metadata2[, c("barcode", "newCelliD")], 
#            file = paste0(inferCNV.dir, "total_cellAnnotations.txt"), 
#            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

## reference group names 
cancers_id = metadata$newCelliD[grepl(pattern="*_Cancer", metadata$newCelliD)] %>% unique
ref_group_names <- setdiff(unique(metadata$newCelliD), cancers_id)

## file 5: define output directory for this sample 
output_dir <- "~/MRes_project_1/docs/SC_SEQ/lung/raw/lambrechtslab/inferCNV/72_hr_runtime/"

## inferCNV object 
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = sub_matrix, 
                                     annotations_file = meta.data, 
                                     gene_order_file = gene_order_file, 
                                     #max_cells_per_group = 20, 
                                     ref_group_names = ref_group_names) ## set to various normal cell types 

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




