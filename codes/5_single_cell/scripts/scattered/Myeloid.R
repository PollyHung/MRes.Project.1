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
SPLIT <- FALSE


## prepare ---------------------------------------------------------------------------------------------------
if(!file.exists("~/MRes_project_1/codes/6_single_cell/RNotebook/subcluster_analysis/myeloid/myeloid_2.rds")){
  # read in file 
  setwd("~/MRes_project_1/codes/6_single_cell/RNotebook/subcluster_analysis/myeloid/")
  myeloid <- readRDS("myeloid.rds")
  myeloid <- subset(myeloid, subset = cell_type == "Myeloid.cell")
  
  # rerun seurat pipeline 
  myeloid <- RunPCA(myeloid, features=VariableFeatures(myeloid))
  myeloid <- FindNeighbors(myeloid, dims = 1:15)
  myeloid <- FindClusters(myeloid, resolution = 0.4)
  myeloid <- RunUMAP(myeloid, dims = 1:15, reduction="pca")
  myeloid@meta.data$chr2qStatus <- ifelse(myeloid@meta.data$chr2q > 0, "amp", "wt")
  myeloid <- subset(myeloid, subset = chr2qStatus %in% c("wt", "amp"))
  
  p1 <- DimPlot(myeloid, label = TRUE, pt.size = 0.7) + NoLegend() 
  p2 <- DimPlot(myeloid, label = TRUE, pt.size = 0.7, split.by = "chr2qStatus")  + NoLegend() 
  p3 <- clustree(myeloid)
  
  myeloid$seurat_clusters[which(myeloid$seurat_clusters == 10)] <- 0
  myeloid$seurat_clusters[which(myeloid$seurat_clusters == 11)] <- 10
  
  color <- brewer.pal(length(unique(Idents(myeloid))), "Paired")
  Idents(myeloid) <- myeloid$seurat_clusters
  p2 <- DimPlot(myeloid, label = TRUE, split.by = "chr2qStatus", cols = color)  + NoLegend() 
  
  # differential expression 
  DE.TcellAllMarkers <- FindAllMarkers(myeloid)
  DE.TcellAllMarkers <- Add_Pct_Diff(DE.TcellAllMarkers)
  DE.TcellAllMarkers <- DE.TcellAllMarkers %>% dplyr::filter(p_val_adj < 0.05) %>% 
    dplyr::arrange(cluster, desc(avg_log2FC), desc(pct_diff))
  write.csv(DE.TcellAllMarkers, "DE.TcellAllMarkers.csv")
  
  markers <- c("CD1C", "FCER1A", "CLEC10A", "HLA-DQA1", "HLA-DQB1", "FITM1", "IRF8", "CLEC9A", "XCR1", "CLNK",
               "HLA-DPB1", "HLA-DPA1", "LAMP3", "IDO1", "IRF4", "LGALS2", "S100A12", "LYZ", "VCAN", "S100A9",
               "ITGB2", "ITGAM", "SERPINA1", "CX3CR1", "LILRA5", "LILRB1", "EREG", "CD68", "CCL20", "AREG",
               "CXCL8", "IL1RN", "VEGFA", "THBS1", "APOBEC3A", "CXCL10", "FCAR", "FN1", "PLIN2", "TIMP1",
               "FCGR3A", "CD163", "FABP5", "C1QA", "RGS2", "IL4I1", "MAFB", "RNASE2", "FCGR3B", "CCL4",
               "TNF", "DUSP2", "C3", "MMP9", "SPP1", "TGFBI", "PLTP", "MARCO", "HAMP", "NUPR1",
               "FOLR2", "GPNMB", "CD14", "MKI67") 
  
  p6 <- DotPlot(myeloid, features = markers, dot.scale = 4) + theme_bw() + 
    scale_colour_gradientn(colours = c("blue", "white", "red"), values = scales::rescale(c(-2, 0, 2))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          legend.position = "top") + 
    ggtitle("Myeloid subtype by chr 2q status for DE features")
  
  cell_type <- c(1:length(unique(myeloid@meta.data$seurat_clusters))) %>% as.data.frame
  cell_type$cell_type <- c("M0", "M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10")
  colnames(cell_type) <- c("clusters", "cell_type")
  
  Idents(myeloid) <- myeloid@meta.data$seurat_clusters
  new.cluster.ids <- cell_type[, 2]
  names(new.cluster.ids) <- levels(myeloid)
  myeloid <- RenameIdents(myeloid, new.cluster.ids)
}
if(!file.exists("/rds/general/user/ph323/ephemeral/GSE180661/cellChat/cancer_myeloid.rds")){
  cancer <- readRDS("~/MRes_project_1/codes/6_single_cell/RNotebook/subcluster_analysis/cancer/cancer.rds")
  myeloid <- readRDS("~/MRes_project_1/codes/6_single_cell/RNotebook/subcluster_analysis/myeloid/myeloid_2.rds")
  
  cancer$cell_subtype <- "cancer"
  myeloid$cell_subtype <- Idents(myeloid)
  
  seurat_object <-  merge(cancer, y = myeloid, project = "GSE180661")
  seurat_object <- JoinLayers(seurat_object)
  Idents(seurat_object) <- seurat_object$cell_subtype
  
  saveRDS(seurat_object, "/rds/general/user/ph323/ephemeral/GSE180661/cellChat/cancer_myeloid.rds")
  
  Idents(myeloid) <- myeloid$chr2qStatus
  DE.markers <- FindMarkers(myeloid, ident.1 = "amp", ident.2 = "wt")
  #DE.markers <- FindAllMarkers(myeloid)
  #DE.markers <- Add_Pct_Diff(DE.markers)
  #DE.markers$pct_diff <- round(DE.markers$pct_diff, digits = 1)
  DE.markers <- DE.markers %>% #dplyr::filter(p_val_adj < 0.05) %>% 
    #dplyr::arrange(desc(pct_diff), desc(avg_log2FC)) %>% 
    dplyr::arrange(desc(avg_log2FC)) %>% 
    dplyr::mutate(gene = rownames(DE.markers)) %>% 
    dplyr::select(gene, avg_log2FC)
  write.table(DE.markers, "DE.markers_new.txt", row.names = F, col.names = F, quote = F, sep = "\t")
}
if(!file.exists("~/MRes_project_1/codes/6_single_cell/RNotebook/subcluster_analysis/myeloid/13_Clustered_dotPlot_1.png")){
## gene markers 
markers_2 <- DE.markers %>% group_by(cluster) %>% top_n(5) %>% pull(gene)
markers <- c(markers, markers_2) %>% unique
markers <- markers[order(markers)]

## david clustering 
DAVID <- readxl::read_xlsx("~/MRes_project_1/codes/6_single_cell/RNotebook/subcluster_analysis/myeloid/DAVID.xlsx", sheet = 3)
DAVID <- DAVID %>% dplyr::filter(Category != "DRUGBANK") %>% 
  dplyr::mutate(FDR = as.numeric(FDR), 
                annotation_updated = sapply(str_split(Term, pattern = "[:~]"), function(x) if (length(x) > 1) x[2] else x[1])) %>% 
  dplyr::filter(FDR<0.10)
genes_expanded <- DAVID %>% mutate(Genes = strsplit(as.character(Genes), ", ")) %>% 
  unnest(Genes) %>% distinct(annotation_updated, Genes) %>% 
  dplyr::rename(gene = Genes, annotation = annotation_updated)
genes <- DAVID$Genes %>% strsplit(split = ", ") %>% unlist() %>% trimws() %>% unique() ## filtered genes 
immune_genes <- genes_expanded %>% dplyr::filter(annotation == "Immune System") %>% pull(gene)
genes_annotations <- genes_expanded %>% group_by(gene) %>% summarise(annotations = paste(annotation, collapse = ", ")) %>% ungroup()
genes_annotations_2 <- genes_annotations %>% group_by(annotations) %>% summarise(genes = paste(gene, collapse = ", ")) %>% ungroup()

# coloring 
start_color <- "#1B3C73" 
middle_color <- "#FFFFFF"
end_color <- "#FF407D"
palette_start_to_mid <- colorRampPalette(colors = c(start_color, middle_color))(10)
palette_mid_to_end <- colorRampPalette(colors = c(middle_color, end_color))(10)
full_palette <- c(palette_start_to_mid, palette_mid_to_end[-1])

# plotting
p7 <- Clustered_DotPlot(myeloid, genes, cluster_feature = TRUE, cluster_ident = TRUE, 
                        colors_use_exp = full_palette, flip = TRUE, x_lab_rotate = 90, 
                        colors_use_idents = color) 
png("13_Clustered_dotPlot_1.png", width = 15, height = 3, units = "in", res = 600)
p7
dev.off()
}
if(!file.exists("~/MRes_project_1/docs/SC_SEQ/ova/GSE180661/cellChat/cancer_myeloid/cellchat_amp.rds")){
  # load 
  setwd("~/MRes_project_1/docs/SC_SEQ/ova/GSE180661/cellChat/cancer_myeloid/")
  seurat_object <- readRDS("/rds/general/user/ph323/ephemeral/GSE180661/cellChat/cancer_myeloid.rds")
  
  # split 
  amp <- subset(seurat_object, subset = chr2qStatus == "amp")
  wt <- subset(seurat_object, subset = chr2qStatus == "wt")
  
  # set up function 
  prepare_cellChat <- function(seurat_object){
    data.input <- seurat_object[["RNA"]]$data
    labels <- seurat_object$cell_subtype
    meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
    
    # Create a CellChat object ---------------------------------------------------------------------------------------
    cellChat <- createCellChat(object = seurat_object, group.by = "cell_subtype", assay = "RNA") # ident = cell subtypes 
    head(cellChat@meta)
    levels(cellChat@idents) # show factor levels of the cell labels
    groupSize <- as.numeric(table(cellChat@idents)) # number of cells in each cell group
    
    # Set the ligand-receptor interaction database --------------------------------------------------------------------
    CellChatDB <- CellChatDB.human 
    CellChatDB.use <- subsetDB(CellChatDB) # use all except for "Non-protein Signaling" (i.e., metabolic and synaptic signaling)
    cellChat@DB <- CellChatDB.use
    
    # Preprocessing the expression data for cell-cell communication analysis -------------------------------------------
    cellchat <- subsetData(cellChat) # This step is necessary even if using the whole database
    future::plan("multisession", workers = 4) # do parallel
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    
    # Compute the communication probability and infer cellular communication network ----------------------------------------
    options(future.globals.maxSize = 2 * 1024^3)
    cellchat <- computeCommunProb(cellchat, type = "triMean")
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    
    ## Extract the inferred cellular communication network as a data frame --------------------------------------------------
    df.net.ligand.receotor <- subsetCommunication(cellchat) #all the inferred cell-cell communications at the level of ligands/receptors
    df.net.signalling <- subsetCommunication(cellchat, slot.name = "netP") #inferred communications at the level of signaling pathways
    
    # Infer the cell-cell communication at a signaling pathway level --------------------------------------------------
    cellchat <- computeCommunProbPathway(cellchat)
    
    # Calculate the aggregated cell-cell communication network --------------------------------------------------
    cellchat <- aggregateNet(cellchat)
    groupSize <- as.numeric(table(cellchat@idents))
    
    # return object 
    return(cellchat)
  }
  
  # apply function to amp, create cellchat object 
  cellchat_amp <- prepare_cellChat(amp)
  saveRDS(cellchat_amp, "cellchat_amp.rds")
  
  # apply function to wt, create cellchat object 
  cellchat_wt <- prepare_cellChat(wt)
  saveRDS(cellchat_wt, "cellchat_wt.rds")
  
  # stack bar plot 
  amp2q <- table(Idents(subset(seurat_object, subset=chr2qStatus=="amp")))
  wt2q <- table(Idents(subset(seurat_object, subset=chr2qStatus=="wt")))
  
  skBarPlot <- data.frame(amp2q, wt2q) %>% dplyr::select(Var1, Freq, Freq.1) %>% 
    dplyr::rename(Cell_Type = Var1, amp = Freq, wt = Freq.1)
  colsum <- colSums(skBarPlot[2:3])
  skBarPlot[2:3] <- sweep(skBarPlot[2:3], 2, colsum, FUN = "/")
  skBarPlot[2:3] <- skBarPlot[2:3]*100
  long_data <- gather(skBarPlot, condition, count, -Cell_Type)
  color <- viridis::mako(n = nrow(skBarPlot))
  color <- sample(color)
  p5 <- ggplot(long_data, aes(x = condition, y = count, fill = Cell_Type)) + 
    geom_bar(stat = "identity") + theme_bw() +
    labs(x = NULL, y = "percent of cells") +
    scale_fill_manual(values = color) + theme(legend.position = "right") + guides(fill = guide_legend(ncol = 1))
}

# cancer and myeloid -------------------------------------------------------------------------------------------------------
setwd("~/MRes_project_1/docs/SC_SEQ/ova/GSE180661/cellChat/cancer_myeloid/plots/wildType")
cellchat_wt_cancer <- readRDS("~/MRes_project_1/docs/SC_SEQ/ova/GSE180661/cellChat/cancer_myeloid/cellchat_wt.rds")
cellchat_amp_cancer <- readRDS("~/MRes_project_1/docs/SC_SEQ/ova/GSE180661/cellChat/cancer_myeloid/cellchat_amp.rds")
cellchat_merged_cancer <- mergeCellChat(list(chr2q_amp = cellchat_amp_cancer, chr2q_wt = cellchat_wt_cancer), 
                                 add.names = c("chr2q_amp", "chr2q_wt"))

cellchat_wt <- netAnalysis_computeCentrality(cellchat_wt, slot.name = "netP")
cellchat_amp <- netAnalysis_computeCentrality(cellchat_amp, slot.name = "netP")

groupSize_wt <- as.numeric(table(cellchat_wt@idents))
groupSize_amp <- as.numeric(table(cellchat_amp@idents))

mat_wt <- cellchat_wt@net$weight
mat_amp <- cellchat_amp@net$weight

sig_pathway_amp <- cellchat_amp@netP$pathways
sig_pathway_wt <- cellchat_wt@netP$pathways

netVisual_chord_gene(cellchat_amp, sources.use = c(2, 5:8), targets.use = c(1), slot.name = "netP", legend.pos.x = 10)



# myeloid and t cells -------------------------------------------------------------------
setwd("~/MRes_project_1/docs/SC_SEQ/ova/GSE180661/cellChat/myeloid_tcells/")

cellchat_wt <- readRDS("~/MRes_project_1/docs/SC_SEQ/ova/GSE180661/cellChat/myeloid_tcells/cellchat_wt.rds")
cellchat_amp <- readRDS("~/MRes_project_1/docs/SC_SEQ/ova/GSE180661/cellChat/myeloid_tcells/cellchat_amp.rds")

cellchat_merged <- mergeCellChat(list(chr2q_amp = cellchat_amp, chr2q_wt = cellchat_wt), 
                                 add.names = c("chr2q_amp", "chr2q_wt"))

cellchat_wt <- netAnalysis_computeCentrality(cellchat_wt, slot.name = "netP")
cellchat_amp <- netAnalysis_computeCentrality(cellchat_amp, slot.name = "netP")

groupSize_wt <- as.numeric(table(cellchat_wt@idents))
groupSize_amp <- as.numeric(table(cellchat_amp@idents))

mat_wt <- cellchat_wt@net$weight
mat_amp <- cellchat_amp@net$weight

sig_pathway_amp <- cellchat_amp@netP$pathways
sig_pathway_wt <- cellchat_wt@netP$pathways





## plotting ----------------------------------------------------------------------------
group.cellType <- c("CD4.dysfunc", rep("CD4.mem", 2), "CD4.helper", rep("CD4.reg", 2), 
                    "CD8.cytotoxic", rep("CD8.mem", 2), "CD8.dysfunc", "CD8.ISG", "CD8.cycling", 
                    rep("NK.cycling", 2), "gd.T.cells", "M0", "M1", "M10", "M2", "M3", "M4", 
                    "M5", "M6", "M7", "M8", "M9", "NK.cytotoxic", rep("NK.reg", 4), "CD4.reg") 
names(group.cellType) <- levels(cellchat_wt@idents)
group.cellType <- as.factor(group.cellType)


pathways.show <- c("IFN-II") 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat_wt, signaling = pathways.show, layout = "circle", #vertex.receiver = seq(1,10), 
                    group = group.cellType)


netVisual_chord_cell(cellchat_amp, signaling = pathways.show, group = group.cellType)


cellchat_merged <- mergeCellChat(list(chr2q_amp = cellchat_amp, chr2q_wt = cellchat_wt), 
                                 add.names = c("chr2q_amp", "chr2q_wt"))
gg1 <- rankNet(cellchat_merged, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat_merged, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
gg1 + gg2


netVisual_bubble(cellchat_merged, targets.use = 7, sources.use = c(16:17, 19:25), 
                 comparison = c(1, 2), angle.x = 90, thresh = 0.20, 
                 title.name = "M0/M1, M2:M8 to CD8 cytotoxic T cells") 
ggsave("~/MRes_project_1/docs/SC_SEQ/ova/GSE180661/cellChat/myeloid2tcells.png", gg, width = 22, height = 7, units = "in")

netVisual_bubble(cellchat_merged_cancer, sources.use = 1, 
                 comparison = c(1, 2), angle.x = 45, thresh = 0.20, 
                 title.name = "Cancer to Macrophages") + coord_flip()



#c(7, 10, 11, 15)


pathway.union <- union(cellchat_merged[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))





