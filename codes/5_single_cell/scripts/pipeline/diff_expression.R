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

## small functions 
quick_seurat <- function(cell.type){
  object <- subset(BSC, subset=cell_type %in% cell.type)
  object <- RunPCA(object, features=VariableFeatures(object=object))
  object <- FindNeighbors(object, dims = 1:10)
  object <- FindClusters(object, resolution = 1)
  object <- RunUMAP(object, dims = 1:10, reduction="pca")
  return(object)
}

## setting up 
work.dir = "~/MRes_project_1/codes/6_single_cell/RNotebook/subcluster_analysis/T_cells/"

## DE function 
differential_expression <- function(work.dir, 
                                    cell_type, #"T.cell"
                                    ){
  setwd(work.dir)
  BSC <- readRDS("~/MRes_project_1/docs/SC_SEQ/ova/GSE180661/results/harmony_obj.rds")
  Idents(BSC) <- BSC@meta.data$cell_type
  
  sub.obj <- quick_seurat(cell_type)
  sub.obj@meta.data$chr2qStatus <- ifelse(sub.obj@meta.data$chr2q > 0, "amp", "wt")
  sub.obj <- subset(sub.obj, subset = chr2qStatus %in% c("wt", "amp"))
  
  p1 <- DimPlot(sub.obj, label=TRUE, pt.size = 0.3) + NoLegend() + ggtitle(paste0(cell_type, " subclusters"))
  ggsave(filename = "cell_subtype.png", plot = p1, width = 3.5, height = 3.6, units = "in") 
  
  
  
}
BSC <- readRDS("~/MRes_project_1/docs/SC_SEQ/ova/GSE180661/results/harmony_obj.rds") 
























