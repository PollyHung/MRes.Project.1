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


## function 1
quick_seurat <- function(seurat_object, 
                         cell.type){ ## select cell type 
  ## set idents 
  Idents(seurat_object) <- seurat_object@meta.data$cell_type
  
  ## subset cells of that cell type 
  object <- subset(seurat_object, subset=cell_type %in% cell.type) 
  
  ## redo seurat standard pipeline 
  object <- RunPCA(object, features=VariableFeatures(object=object))
  object <- FindNeighbors(object, dims = 1:10)
  object <- FindClusters(object, resolution = 1)
  object <- RunUMAP(object, dims = 1:10, reduction="pca")
  
  ## add metadata on chromosome 2q status 
  object@meta.data$chr2qStatus <- ifelse(object@meta.data$chr2q > 0, "amp", "wt")
  
  ## quality control to remove samples with NA in chr2q status 
  object <- subset(object, subset = chr2qStatus %in% c("wt", "amp"))
  return(object)
}



## function 2 
differential_expression <- function(seurat_object, 
                                    result.dir, ## "~/MRes_project_1/docs/SC_SEQ/ova/GSE180661/results/differential_exp/"
                                    name){ ## cancer
  ## environment -------------------------------------------------------------------
  if(!dir.exists(result.dir)){
    dir.create(result.dir)
  }
  setwd(result.dir)
  
  ## create DE between clusters ----------------------------------------------------
  ## rename idents 
  Idents(seurat_object) <- seurat_object$seurat_clusters
  
  ## DE between clusters 
  DE.markers.by.cluster <- FindAllMarkers(seurat_object) ## perform differential expression on seurat clusters 
  DE.markers.by.cluster <- Add_Pct_Diff(DE.markers.by.cluster) ## add pct_diff 
  DE.markers.by.cluster <- DE.markers.by.cluster %>% group_by(cluster) %>% 
    dplyr::filter(p_val_adj<0.05) %>% 
    dplyr::arrange(cluster, desc(avg_log2FC), desc(pct_diff)) %>% 
    dplyr::select(cluster, gene, avg_log2FC, pct.1, pct.2, pct_diff, p_val_adj) 
  write.csv(DE.markers.by.cluster, paste0("00_markers_by_cluster_",name,".csv"))
  
  ## top 50 
  DE.markers.by.cluster.top <- DE.markers.by.cluster %>% group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group = TRUE) %>%
    slice_head(n = 50) %>% ungroup()
  DE.markers.by.cluster.top <- DE.markers.by.cluster.top %>% dplyr::mutate(cluster = paste0("cluster_", cluster)) %>%
    dplyr::select(cluster, gene) %>% group_by(cluster) %>% mutate(row = row_number()) %>% ungroup() %>%
    tidyr::pivot_wider(names_from = cluster, values_from = gene) %>% dplyr::select(-row)
  write.csv(DE.markers.by.cluster.top, paste0("01_top50_markers_by_cluster_wider_", name, ".csv"))

  ## create DE between chr2q status for gsea ------------------------------------------ 
  ## rename idents 
  Idents(seurat_object) <- seurat_object$chr2qStatus
  
  ## prepare DE for gsea 
  DE.gsea <- FindMarkers(seurat_object, ident.1 = "amp", ident.2 = "wt")
  DE.gsea <- DE.gsea %>% 
    dplyr::mutate(gene = rownames(DE.gsea)) %>% 
    dplyr::arrange(desc(avg_log2FC)) %>% 
    dplyr::select(gene, avg_log2FC)
  write.table(x = DE.gsea, file =  paste0("02_DE_gsea_rankedList_", name,".txt"), sep = "\t", quote = F, row.names = F, col.names = F)
}



make_plots <- function(seurat_object, 
                       stacked_bar_plot = TRUE){
  # plot functions ---------------------------------------------------------------------
  stacked_bar_plot <- function(seurat_object){
    Idents(seurat_object) <- seurat_object$cell_type
    amp <- table(Idents(subset(seurat_object, subset=chr2q > 0)))
    wt <- table(Idents(subset(seurat_object, subset=chr2q <= 0)))
    
    skBarPlot <- data.frame(amp, wt) 
    skBarPlot <- skBarPlot %>% dplyr::select(Var1, Freq, Freq.1) %>% 
      dplyr::rename(Cell_Type = Var1, amp2q = Freq, wt2q = Freq.1) 
    colsum <- colSums(skBarPlot[, 2:3])
    skBarPlot[, 2:3] <- sweep(skBarPlot[, 2:3], 2, colsum, FUN="/")
    skBarPlot[, 2:3] <- skBarPlot[, 2:3]*100
    
    long_data <- gather(skBarPlot, condition, count, -Cell_Type)
    color <- viridis::mako(n = nrow(skBarPlot)) %>% sample
    p1 <- ggplot(long_data, aes(x = condition, y = count, fill = Cell_Type)) + 
      geom_bar(stat = "identity") + theme_bw() +
      labs(x = NULL, y = "percent of cells") +
      scale_fill_manual(values = color) + theme(legend.position = "right") + guides(fill = guide_legend(ncol = 1))
    ggsave("01_StackedBarPlot.png", p1, width = 4, height = 6, units = "in")
  }
  
  # plot 1 -----------------------------------------------------------------------------
  if(stacked_bar_plot){
    stacked_bar_plot(seurat_object)
    }
}










