library(SeuratDisk)
library(RColorBrewer)
library(magrittr)
library(CellChat)
library(patchwork)
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(dplyr)
library(celldex)
library(harmony)
library(SingleR)
library(scCustomize)
library(clustree)
library(tidyverse)
library(ggrepel)

setwd("~/MRes_project_1/codes/6_single_cell/RNotebook/subcluster_analysis/myeloid")
seurat_object <- readRDS("myeloid_2.rds")
DE.markers <- FindAllMarkers(seurat_object, test.use = "MAST")
write.csv(DE.markers, "DE_All_markers_MAST.csv")

#setwd("~/MRes_project_1/docs/SC_SEQ/ova/GSE180661/cellChat/cancer_myeloid/")

## interferon signalling? 
#cellchat_amp <- readRDS("cellchat_amp.rds")
#cellchat_wt <- readRDS("cellchat_wt.rds")

#reticulate::use_python("/rds/general/user/ph323/home/anaconda3/bin/python3")
#reticulate::py_install(packages = 'numpy', )
#reticulate::py_install(packages = 'umap-learn')

## Identify signaling groups based on their functional similarity
#signalling_group <- function(cellchat){
#  cellchat <- computeNetSimilarity(cellchat, type = "functional")
#  cellchat <- netEmbedding(cellchat, type = "functional")
#  cellchat <- netClustering(cellchat, type = "functional")
#  return(cellchat)
#}
#cellchat_amp <- signalling_group(cellchat_amp)
#cellchat_wt <- signalling_group(cellchat_wt)

#saveRDS(cellchat_amp, "cellchat_amp.rds")
#saveRDS(cellchat_wt, "cellchat_wt.rds")

#netVisual_embedding(cellchat_amp, type = "functional", label.size = 3.5)
#netVisual_embedding(cellchat_wt, type = "functional", label.size = 3.5)