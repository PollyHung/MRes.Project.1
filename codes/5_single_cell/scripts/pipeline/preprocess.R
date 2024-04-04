library(Seurat)
library(SeuratObject)
library(ggplot2)
library(magrittr)
library(dplyr)
library(patchwork)
library(Signac)
library(SingleR)
library(celldex)
library(harmony)
library(multtest)
library(AnnotationHub)
library(tidyverse)
library(Matrix)
library(RCurl)
library(scales)
library(cowplot)
library(scCustomize)

## this is a function to perform preprocessing of seurat objects  
## function 
preprocess <- function(count_mtx,  ## a path to the dataframe of count matrix 
                       result_dir,  ## a path to the result directory 
                       metadata_external, ## a path to the external metadata 
                       names_field, ## names.field = c(2, 3) for lung GSE, 3 for ovarian GSE, not needed for lambrechelab
                       doHarmony = FALSE, ## whether or not to perform Harmony normalisation? 
                       var_to_regress ## other variables to regress according to metadata table 
                       ){ 

  ## create directory if it does not exist  
  if(!dir.exists(result_dir)){
    print("create directory")
    dir.create(file.path(result_dir)) 
  }
  
  ## set new working directory 
  setwd(result_dir)
  getwd()
  
  ## check if we have done umapped_obj before? If not, we will create this rds 
  if(!file.exists("umapped_obj.rds")){ 
    ## Read in 
    if(file_test("-f", count_mtx)){ ## check is this a file? 
      print("a count matrix file provided")
      count_mtx <- readRDS(count_mtx)
      BSC <- CreateSeuratObject(counts = count_mtx, min.cells = 10, min.features = 200, assay="RNA", names.field=names_field)
    } else {
      print("a count matrix folder including matrix and cell barcode provided")
      BSC <- Read10X(data.dir = count_mtx)
      BSC <- CreateSeuratObject(counts = BSC, min.cells = 10, min.features = 200, assay="RNA")
    }
    BSC@meta.data$cell <- BSC@meta.data %>% row.names() 
    metadata <- BSC@meta.data
    head(metadata)
    
    ## QC metrics 
    BSC <- Add_Cell_QC_Metrics(seurat_object = BSC, species = "Human") ## function from package scCustomize
    BSC[["sample"]] <- BSC@meta.data$orig.ident
    head(BSC@meta.data)
    
    ## Get QC plots 
    p1 <- BSC@meta.data %>% ggplot(aes(x=sample, fill=sample)) + geom_bar() + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1), legend.position = "none") + ggtitle("Cell Count per Sample for combined dataset") 
    p2 <- BSC@meta.data %>% ggplot(aes(color=sample, x=nCount_RNA, fill= sample)) + geom_density(alpha = 0) + theme_bw() + theme(legend.position = "none") + scale_x_log10() + geom_vline(xintercept = 200) + ggtitle("UMI per cells")  
    p3 <- BSC@meta.data %>% ggplot(aes(color=sample, x=nFeature_RNA, fill= sample)) + geom_density(alpha = 0) + theme_bw() + theme(legend.position = "none") + scale_x_log10() + geom_vline(xintercept = 2500) + ggtitle("Gene Per Cells") 
    p4 <- BSC@meta.data %>% ggplot(aes(x=sample, y=log10GenesPerUMI, fill=sample)) + geom_violin() + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1), legend.position = "none") + ggtitle("NCells vs NGenes") 
    if("percent_mito" %in% colnames(BSC@meta.data)){
      p5 <- QC_Plot_UMIvsGene(seurat_object = BSC, meta_gradient_name = "percent_mito", low_cutoff_gene = 600, high_cutoff_gene = 6000, high_cutoff_UMI = 50000, meta_gradient_low_cutoff = 15, combination = TRUE)
    } else {
      p5 <- QC_Plot_UMIvsGene(seurat_object = BSC, low_cutoff_gene = 600, high_cutoff_gene = 6000, high_cutoff_UMI = 50000, low_cutoff_UMI = 300)
    }
    
    pdf("QualityControl.pdf", width=10, height=5)
    print(p1)
    print(p2)
    print(p3)
    print(p4)
    print(p5)
    dev.off()
    
    ## Standard work flow 1 
    if("percent_mito" %in% colnames(BSC@meta.data)){
      BSC <- subset(BSC, subset= nCount_RNA>200 & nFeature_RNA<6000 & percent_mito <15 & log10GenesPerUMI>0.75)
    } else {
      BSC <- subset(BSC, subset= nCount_RNA>200 & nFeature_RNA<6000 & log10GenesPerUMI>0.75)}
    
    BSC <- NormalizeData(BSC, normalization.method = "LogNormalize", scale.factor = 10000)
    BSC <- FindVariableFeatures(BSC, selection.method = "vst", nfeatures = 3000)
    
    if("percent_mito" %in% colnames(BSC@meta.data)){
      BSC <- ScaleData(BSC, vars.to.regress = c("percent_mito", "nCount_RNA", "S.Score", "G2M.Score"))
    } else {
      BSC <- ScaleData(BSC, vars.to.regress = c("nCount_RNA", "S.Score", "G2M.Score"))}
    
    BSC <- JoinLayers(BSC)
    BSC <- RunPCA(BSC, features = VariableFeatures(object = BSC))
    pca_stdev <- BSC@reductions$pca@stdev
    
    p7 <- ElbowPlot(BSC, ndims = length(pca_stdev))
    ggsave("SelectPCA.pdf", p7, width=7, height=5, units="in")
    
    pca_above_2 <- sum(pca_stdev > 2)
    BSC <- FindNeighbors(BSC, dims = 1:pca_above_2, reduction="pca")
    BSC <- FindClusters(BSC, resolution = 0.5, reduction="pca")
    BSC <- RunUMAP(BSC, dims = 1:pca_above_2, reduction="pca")
    saveRDS(BSC, "interim.rds")
    print("find cell clusters and run UMAP")
    
    table(BSC$orig.ident, BSC@meta.data$seurat_clusters)
    Idents(object = BSC) <- BSC@meta.data$seurat_clusters
    saveRDS(BSC, "interim.rds")
    
    ## update the metadata for BSC with external inputs 
    metadata <- BSC@meta.data
    metadata_ext <- read.table(metadata_external, sep="\t", header=TRUE)
    metadata <- left_join(metadata, metadata_ext, by="cell")
    rownames(metadata) <- metadata$cell
    metadata <- metadata[rownames(BSC@meta.data), ]
    BSC@meta.data <- metadata
    
    ## save the plot 
    pdf("ClusterResult_noHarmony.pdf", width = 8, height = 7)
    p5 <- DimPlot(BSC, reduction="umap", label=TRUE, pt.size = 0.1) + ggtitle("DimPlot labelled by seurat sluster")
    p6 <- DimPlot(BSC, reduction="umap", group.by = "cell_type") + ggtitle("DimPlot labelled by annotated cells") 
    p7 <- DimPlot(BSC, reduction="umap", group.by = "sample_id") + ggtitle("DimPlot labelled by sample id") + NoLegend()
    print(p5)
    print(p6)
    print(p7)
    dev.off()
    
    ## save the RDS object 
    saveRDS(BSC, "umapped_obj.rds")
  }
  
  
  ## Do Harmony?
  if(doHarmony){
    setwd(result_dir)
    getwd()
    
    ## read in rds 
    BSC <- readRDS("umapped_obj.rds")
    
    ## run Harmony
    BSC <- RunHarmony(BSC, 
                      group.by.vars = var_to_regress, 
                      reduction = "pca", 
                      assay.use = "RNA", 
                      reduction.save = "harmony", 
                      lambda = length(var_to_regress))
    harmony_stdev <- BSC@reductions$harmony@stdev
    harmony_above_2 <- sum(harmony_stdev > 2)
    
    BSC <- FindNeighbors(BSC, dims = 1:harmony_above_2, reduction="harmony")
    BSC <- FindClusters(BSC, resolution = 0.5, reduction="harmony")
    BSC <- RunUMAP(BSC, dims = 1:harmony_above_2, reduction="harmony")
    table(BSC$orig.ident, BSC@meta.data$seurat_clusters)
    Idents(object = BSC) <- BSC@meta.data$seurat_clusters
    
    ## save plots 
    pdf("ClusterResult_Harmony.pdf", width = 8, height = 7)
    p5 <- DimPlot(BSC, reduction="umap", label=TRUE, pt.size = 0.1) + ggtitle("DimPlot labelled by seurat sluster")
    p6 <- DimPlot(BSC, reduction="umap", group.by = "cell_type") + ggtitle("DimPlot labelled by annotated cells") 
    p7 <- DimPlot(BSC, reduction="umap", group.by = "sample_id") + ggtitle("DimPlot labelled by sample id") 
    print(p5)
    print(p6)
    print(p7)
    dev.off()
    
    ## save RDS
    saveRDS(BSC, "harmony_obj.rds")
  }
}

