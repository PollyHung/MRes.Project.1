## Library packages 
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(SeuratObject)


list.files("~/MRes_project_1/docs/SC_SEQ/ova/raw/Olbrecht")

rds <- readRDS("~/MRes_project_1/docs/SC_SEQ/ova/raw/Olbrecht/SOL_counts_matrix.rds")
metadata <- read.csv("~/MRes_project_1/docs/SC_SEQ/ova/raw/Olbrecht/SOL_metadata.csv")

BSC <- CreateSeuratObject(rds, assay = "RNA", meta.data = metadata)
