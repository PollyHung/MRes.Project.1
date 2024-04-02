## This script is part of the series of script that performs sample level single cell sequencing 
## The data is downloaded online and is a collection of different samples 
## 5 ovarian cancer patents (45,114 cells), 3' prime scRNA-seq
## https://lambrechtslab.sites.vib.be/en/pan-cancer-blueprint-tumour-microenvironment-0
## Sample and Data Relationship Format (SDRF)
## Pipeline: 
## Step 1: 
rm(list = ls())

## Section 0: Environment  ============================================================
## Packages 
library(infercnv)
library(dplyr)
library(Matrix)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(SeuratData)
library(ggplot2)
library(magrittr)
library(patchwork)
library(Signac)
library(SingleR)
library(celldex)
library(stringr)

## Controls 
count.matrix <- "~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/count_matrix/"
meta.data <- "~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/metadata.csv.gz"
data.dir <- "~/MRes_project_1/docs/SC_SEQ/lung/raw/lambrechtslab/"
inferCNV.dir <- "~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/InferCNV"
gene.order.file <- "~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/gene_order_file_exp.txt"
output.dir <- "~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/InferCNV/result/"
rds.dir <- "~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/ova.rds"

## Section 2: Perpare object for inferCNV ==========================================================
## load the general dataset and intialize the seurat object 
#BSC <- Read10X(data.dir = count.matrix)
#BSC <- CreateSeuratObject(counts = BSC, min.cells = 3, min.features = 200)
BSC <- readRDS(rds.dir)

## Read in metadata 
metadata <- read.csv(meta.data)
rownames(metadata) <- metadata$Cell
obj_names <- c("nGene", "nUMI", "CellFromTumor", "PatientNumber", "TumorType", "TumorSite", "CellType")
metadata <- metadata[obj_names] ## select out interested groups 

## read in sdrf
sdrf <- read.delim("~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/E-MTAB-8107.sdrf.txt")

## append 
metadata <- metadata[rownames(BSC@meta.data), ] ## reorder by BSC 
for(i in 1:7){
  named_vector <- metadata[[i]]
  names(named_vector) <- rownames(BSC@meta.data)
  BSC <- AddMetaData(object=BSC, 
                     metadata=named_vector, 
                     col.name=obj_names[[i]])
}



## Split Seurat objects based on Samples 
obj.list <- SplitObject(BSC, split.by = "orig.ident")


## perform InferCNV ==============================================================================
## example with scrSOL007
gene_order <- read.delim("~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/gene_order_file_exp.txt", header = FALSE)
((gene_order$V1) %in% (obj.list$scrSOL007 %>% Features())) %>% sum

sub_matrix <- obj.list$scrSOL004@assays$RNA$counts %>% as.matrix()
sub_metadata <- obj.list$scrSOL004@meta.data
sub_metadata$barcode <- row.names(sub_metadata)
sub_metadata <- sub_metadata[, c("barcode", "CellType")]
sub_metadata$CellType <- ifelse(sub_metadata$CellType == "Cancer", "cancer", "non_cancer")
rownames(sub_metadata) <- NULL
write.table(sub_metadata, paste0(inferCNV.dir, "cellAnnotations_alt.txt"), 
            sep = '\t', quote = F, row.names = F, col.names = F)

infercnv_obj <-  CreateInfercnvObject(raw_counts_matrix=sub_matrix, 
                                      annotations_file=paste0(inferCNV.dir, "cellAnnotations_alt.txt"),
                                      delim="\t",
                                      gene_order_file=gene.order.file,
                                      ref_group_names=c("non_cancer")) ## set to various normal cell types 

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir=paste0(output.dir, "scrSOL004_alt/"),  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T, 
                             num_threads = 10)




















