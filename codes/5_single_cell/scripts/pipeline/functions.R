library(Seurat)
library(SeuratObject)
library(SeuratData)
library(SeuratWrappers)
library(scCustomize)
library(Nebulosa)
library(EnhancedVolcano)
library(magrittr)
library(dplyr)
library(monocle3)
library(pheatmap)
library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)

setwd("~/MRes_project_1/codes/6_single_cell/RNotebook/subcluster_analysis/cancer/")
cancer <- readRDS("cancer.rds")
Idents(cancer) <- cancer$chr2qStatus

DE.markers <- FindMarkers(cancer, ident.1 = "amp", ident.2 = "wt", test.use = "MAST")

write.csv(DE.markers, "DE_markers.csv")