library(Seurat)
library(ggalluvial)
library(patchwork)
library(CellChat)
library(dplyr)
library(magrittr)
library(SeuratObject)
library(SeuratWrappers)
library(NMF)
library(monocle3)



## cancer to M5 cells ========================================================================================
setwd("~/MRes_project_1/docs/SC_SEQ/ova/GSE180661/cellChat/cancer_myeloid/")

#load cell chat object 
cellchat_amp <- readRDS("cellchat_amp.rds")
cellchat_wt <- readRDS("cellchat_wt.rds")
cellchat <- readRDS("cellchat_combined.rds")

# pathways 
cellchat_amp.pathways <- cellchat_amp@netP$pathways
cellchat_wt.pathways <- cellchat_wt@netP$pathways
#setdiff(cellchat_amp.pathways, cellchat_wt.pathways)
common <- intersect(cellchat_amp.pathways, cellchat_wt.pathways)

# from cancer to M5? 
pdf("plots.pdf")
for(i in common){
  #netVisual_chord_gene(cellchat_amp, sources.use = 1, targets.use = c(2:12), lab.cex = 0.5, signaling = i, 
  #                     legend.pos.y = 5, title.name = paste0("cancer to M ", i, " in chr2q amp"))
  #netVisual_chord_gene(cellchat_wt, sources.use = 1, targets.use = 8, lab.cex = 0.5, signalling = i, 
  #                     legend.pos.y = 5, title.name = paste0("cancer to M ", i, " in chr2q wt"))
  netVisual_aggregate(cellchat_amp, signaling = i, layout = "chord", sources.use = c(1, 8), targets.use = c(1:12))
  netVisual_aggregate(cellchat_wt, signaling = i, layout = "chord", sources.use = c(1, 8), targets.use = c(1:12))
}
dev.off()

# select K
amp_selectK <- selectK(cellchat_amp, pattern = "outgoing") # 7
wt_selectK <- selectK(cellchat_wt, pattern = "outgoing") # 

# cellchat_communication patterns  
cellchat_amp <- identifyCommunicationPatterns(cellchat_amp, pattern = "outgoing", k = 7)
cellchat_wt <- identifyCommunicationPatterns(cellchat_wt, pattern = "outgoing", k = 8)

# river plot 
netAnalysis_river(cellchat_amp, pattern = "outgoing")
netAnalysis_river(cellchat_wt, pattern = "outgoing")

# compute centrality 
cellchat_amp <- netAnalysis_computeCentrality(cellchat_amp, slot.name = "netP") 
cellchat_wt <- netAnalysis_computeCentrality(cellchat_wt, slot.name = "netP") 
netAnalysis_signalingRole_network(cellchat_amp, signaling = "CCL", width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat_wt, signaling = "CCL", width = 8, height = 2.5, font.size = 10)

# dot plot 
netVisual_bubble(cellchat_amp, sources.use = 1, targets.use = c(8), remove.isolate = FALSE)
netVisual_bubble(cellchat_wt, sources.use = 1, targets.use = c(8), remove.isolate = FALSE)

# upregulate and downregulate 
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(3, 8),  comparison = c(1, 2), angle.x = 45)

object.list <- 
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat_amp <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cellchat_wt <-  netAnalysis_computeCentrality(cellchat_wt, slot.name = "netP")
gg1 <- netAnalysis_signalingChanges_scatter(cellchat_amp, idents.use = "cancer", signaling.exclude = "MIF")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "cDC1", signaling.exclude = c("MIF"))


## M5 to T cells ========================================================================================
setwd("~/MRes_project_1/docs/SC_SEQ/ova/GSE180661/cellChat/myeloid_tcells/")

#load cell chat object 
cellchat_amp <- readRDS("cellchat_amp.rds")
cellchat_wt <- readRDS("cellchat_wt.rds")

# pathways 
cellchat_amp.pathways <- cellchat_amp@netP$pathways
cellchat_wt.pathways <- cellchat_wt@netP$pathways
#setdiff(cellchat_amp.pathways, cellchat_wt.pathways)

# cells 
cells <- cellchat_wt@netP$prob %>% rownames()

# signalling M5 to T cells? 
netVisual_chord_gene(cellchat_amp, sources.use = 22, targets.use = c(1:11, 15), title.name = "M5 to T cells, chr2q amp", 
                     legend.pos.x = 1, legend.pos.y = 6, slot.name = "netP", lab.cex = 0.4)
netVisual_chord_gene(cellchat_wt, sources.use = 22, targets.use = c(1:11, 15), title.name = "M5 to T cells, chr2q wt", 
                     legend.pos.x = 1, legend.pos.y = 6, slot.name = "netP", lab.cex = 0.4)

# heatmap 
pathway <- "CXCL"
netVisual_heatmap(cellchat_amp, signaling = pathway, color.heatmap = "Reds", 
                  row.show = c(1:11, 15, 16:26, 32), col.show = c(1:11, 15, 16:26, 32), 
                  title.name = paste0(pathway, " signalling in chr2q amp"))
netVisual_heatmap(cellchat_wt, signaling = pathway, color.heatmap = "Reds", 
                  row.show = c(1:11, 15, 16:26, 32), col.show = c(1:11, 15, 16:26, 32), 
                  title.name = paste0(pathway, " signalling in chr2q wt"))



## monocle3 ==============================================================================
setwd("~/MRes_project_1/codes/6_single_cell/RNotebook/subcluster_analysis/myeloid")
seurat_object <- readRDS("myeloid_2.rds")

monocle <- as.cell_data_set(seurat_object)
monocle <- cluster_cells(cds = monocle, reduction_method = "UMAP")
monocle <- learn_graph(monocle, use_partition = TRUE)
monocle <- order_cells(monocle, reduction_method = "UMAP")

plot_cells(
  cds = monocle,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)








