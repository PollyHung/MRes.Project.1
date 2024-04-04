library(CellChat)
library(patchwork)
library(Seurat)


do_cellchat <- function(seurat_object ## an seurat object with Ident = cell type 
                        ){
  cellChat <- createCellChat(object = seurat_object, group.by = "ident", assay = "RNA") 
  groupSize <- as.numeric(table(cellChat@idents))
  CellChatDB <- CellChatDB.human 
  CellChatDB.use <- subsetDB(CellChatDB) 
  cellChat@DB <- CellChatDB.use
  cellchat <- subsetData(cellChat) 
  future::plan("multisession", workers = 4) 
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  options(future.globals.maxSize = 1024 * 1024^2)
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  df.net.ligand.receotor <- subsetCommunication(cellchat)
  df.net.signalling <- subsetCommunication(cellchat, slot.name = "netP") 
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  groupSize <- as.numeric(table(cellchat@idents))
  
  return(cellchat)
}