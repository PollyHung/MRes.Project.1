## Section 1: Environment: ======================================================================================
## Data info 
## GSE180661
## Ovarian cancer mutational processes drive site-specific immune evasion

## Set control 
LIBRARY <- TRUE
SET_DIRECTORY <- TRUE
SAVE_IMAGE <- TRUE
TYPE <- "lung"
LOAD <- TRUE
SAVE_RDS <- TRUE
PREPROCESS <- TRUE
PROCESS <- TRUE
INFERCNV <- TRUE

## Library 
library(Seurat)
library(ggplot2)
library(magrittr)
library(dplyr)
library(patchwork)
library(Signac)
library(SingleR)
library(celldex)
library(stringr)
library(infercnv)
library(Matrix)
library(SeuratObject)
library(SeuratDisk)
library(SeuratData)



if(SET_DIRECTORY){
  ## for seurat 
  save.dir <-  paste0("~/MRes_project_1/docs/SC_SEQ/", TYPE,"/raw/bischoff/", TYPE, ".RData")
  save.rds <- paste0("~/MRes_project_1/docs/SC_SEQ/", TYPE,"/raw/lambrechtslab/", TYPE, ".rds") ## TYPE = blood or skin 
  result.dir <- paste0("~/MRes_project_1/docs/SC_SEQ/", TYPE,"/processed/lambrechtslab/")
  
  ## for infer cnv 
  count.matrix <- paste0("~/MRes_project_1/docs/SC_SEQ/", TYPE,"/raw/lambrechtslab/count_matrix/")
  meta.data <- paste0("~/MRes_project_1/docs/SC_SEQ/", TYPE,"/raw/lambrechtslab/metadata.csv.gz")
  data.dir <- "~/MRes_project_1/docs/SC_SEQ/lung/raw/bischoff/"
  inferCNV.dir <- "~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/InferCNV"
  gene.order.file <- "~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/gene_order_file_exp.txt"
  output.dir <- "~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/InferCNV/result/"
}

print("setting up environment done!")




list.files("/rds/general/user/ph323/")

## Section 2: Preprocess ======================================================================================
setwd("~/MRes_project_1/docs/SC_SEQ/lung/raw/bischoff/") 

if(PREPROCESS){
  
  file_names <- list.files("cellranger")
  mylist <- list()
  
  for(i in 1:length(file_names)){
    BSC_data <- Read10X(paste0("cellranger/", file_names[[i]], "/filtered_feature_bc_matrix"))
    mylist[[i]] <- CreateSeuratObject(counts = BSC_data, min.cells = 3, min.features = 200)
  }
  
  BSC.big <- merge(mylist[[1]], y = c(mylist[[2]],mylist[[3]],mylist[[4]],mylist[[5]],mylist[[6]],mylist[[7]],mylist[[8]],mylist[[9]],mylist[[10]],mylist[[11]],mylist[[12]],mylist[[13]],mylist[[14]],mylist[[15]],mylist[[16]],mylist[[17]],mylist[[18]],mylist[[19]],mylist[[20]]), 
                   add.cell.ids = c(file_names[1],file_names[2],file_names[3],file_names[4],file_names[5],file_names[6],file_names[7],file_names[8],file_names[9],file_names[10],file_names[11],file_names[12],file_names[13],file_names[14],file_names[15],file_names[16],file_names[17],file_names[18],file_names[19],file_names[20]), 
                   project = "BSC")
  
  ## load the PBMC dataset 
  BSC_data <- Read10X(data.dir = "~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/count_matrix/")
  ## Initialize the Seurat oject with the raw (non-normalized data)
  BSC <-  CreateSeuratObject(counts = BSC_data, min.cells = 3, min.features = 200)
  ## build rownames 
  BSC@meta.data$sample_ids <- BSC@meta.data %>% row.names() %>% str_extract("^[^_]*")   
  print("section 2 set up data done!")
  
  
  
  ## Section 2.2: quality control and normalization  
  ## calculate mitochondrial QC metrics (percentage of counts originating from a set of features)
  BSC[["percent.mt"]] <- PercentageFeatureSet(BSC, pattern = "^MT-")
  
  ## visualize QC metrics and use these to filter cells 
  ## visualisation 
  pdf(file = paste0(result.dir, "QC_violin_plot.pdf"), width = 15, height = 5)   # The directory you want to save the file in
  VlnPlot(BSC, features = "nFeature_RNA", ncol = 1, pt.size = 0) 
  VlnPlot(BSC, features = "nCount_RNA", ncol = 1, pt.size = 0)
  VlnPlot(BSC, features ="percent.mt", ncol = 1, pt.size = 0)
  dev.off() 
  
  ## filter cells 
  BSC <- subset(BSC, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15) ## nFeature_originalexp
  
  ## normalizing the data 
  ## normalization based on assumption that each cell originally contains same number of RNA molecules 
  BSC <- NormalizeData(BSC, normalization.method = "LogNormalize", scale.factor = 10000)
  
  print("QC and normalisation done!")
  if(SAVE_IMAGE){save.image(file = save.dir)}
  
  
  
  ## Section 2.3: identification of highly variable features (feature selection)  
  ## calculate subset of features that exhibit high cell-to-cell variation in the dataset 
  ## by default we return 2000 features per dataset for downstream analysis like PCA, here we identify 3000 
  BSC <- FindVariableFeatures(BSC, selection.method = "vst", nfeatures = 3000)
  
  ## identify 10 most highly variable genes 
  top10 <- head(VariableFeatures(BSC), 10)
  
  ## plot variable features with and without labels 
  plot1 <- VariableFeaturePlot(BSC)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  
  ## save the file 
  pdf(file = paste0(result.dir, "VariableFeaturePlot.pdf"), width = 8, height = 4) 
  print(plot1) 
  print(plot2)
  dev.off() 
  
  print("feature selection done!")
  if(SAVE_IMAGE){save.image(file = save.dir)}
  
  
  
  ## Section 2.4: scaling the data  
  ## apply linear transformation prior to dimensional reduction techniques 
  ## shift the expression of each gene so mean expression is 0 and variabce across cells is 1 
  ## result stored in BSC[["RNA"]]$scale.data
  ## Scaling the data 
  BSC <- ScaleData(BSC, vars.to.regress = c("nCount_RNA",
                                            "percent.mt"),
                   do.par = TRUE,
                   num.cores = 3)
  
  ## Cell cycle adjust
  regev_lab_cell_cycle_genes <- read.delim("~/MRes_project_1/codes/6_single_cell/cell_cycle_adjust/regev_lab_cell_cycle_genes.txt", 
                                           header=F)
  s.genes <- regev_lab_cell_cycle_genes[1:43,1]
  g2m.genes <- regev_lab_cell_cycle_genes[44:97,1]
  BSC_join <- JoinLayers(BSC)
  BSC_join <- CellCycleScoring(object = BSC_join,
                               s.features = s.genes,
                               g2m.features = g2m.genes,
                               set.ident = TRUE)
  
  if(SAVE_IMAGE){save.image(file = save.dir)}
  print("data scaling done!")
  
  
  
  ## Section 2.5: perform linear dimensional reduction  
  BSC_join <- RunPCA(BSC_join, features = c(s.genes, g2m.genes))
  
  ## print plots 
  pdf(file = paste0(result.dir,"dimplot_beforeccadj.pdf"), width = 8, height = 8)
  VizDimLoadings(BSC_join, dims = 1:2, reduction = "pca") ## visualise features 
  DimPlot(BSC_join) ## visualise cells 
  DimHeatmap(BSC_join, dims = 1:9, cells = 500, balanced = TRUE) ## explore if feature is included in downstream analysis 
  ElbowPlot(BSC_join) ## determine dimensionality of the dataset
  dev.off()
  
  if(SAVE_IMAGE){save.image(file = save.dir)}
  print("linear dimensional reduction done!")
  
  
  
  ## Section 2.6: scaling again against set features  
  BSC_join <- ScaleData(BSC_join,vars.to.regress = c("nCount_RNA",
                                                     "percent.mt",
                                                     "S.Score",
                                                     "G2M.Score"),
                        do.par = TRUE,
                        num.cores = 3)
  ## Run pca again 
  BSC_join <- RunPCA(BSC_join, features = c(s.genes, g2m.genes))
  
  ## print plots 
  pdf(file = paste0(result.dir,"dimplot_var.pdf"), width = 8, height = 4)
  VizDimLoadings(BSC_join, dims = 1:2, reduction = "pca") ## visualise features 
  DimPlot(BSC_join, reduction = "pca") ## visualise cells 
  DimHeatmap(BSC_join, dims = 1:15, cells = 500, balanced = TRUE) ## explore if feature is included in downstream analysis 
  ElbowPlot(BSC_join) ## determine dimensionality of the dataset
  dev.off()
  
  ## pca against automatically selected features 
  BSC_join <- RunPCA(BSC_join, features = VariableFeatures(object = BSC_join))
  
  ## print plots 
  pdf(file = paste0(result.dir,"dimplot_var_2.pdf"), width = 8, height = 4)
  VizDimLoadings(BSC_join, dims = 1:2, reduction = "pca") ## visualise features 
  DimPlot(BSC_join, reduction = "pca") ## visualise cells 
  DimHeatmap(BSC_join, dims = 1:15, cells = 500, balanced = TRUE) ## explore if feature is included in downstream analysis 
  ElbowPlot(BSC_join) ## determine dimensionality of the dataset
  dev.off()
  
  print("re-scaling and pca done!")
  if(SAVE_IMAGE){save.image(file = save.dir)}
  
  
  
  ## Section 2.7: cell clustering  
  BSC_join <- FindNeighbors(BSC_join, dims = 1:15)
  BSC_join <- FindClusters(BSC_join, resolution = 0.5)
  head(Idents(BSC_join), 5) ## view the cluster IDS of the first 5 cells 
  ## we know that there is 25 levels 
  
  print("cell clustering done!")
  if(SAVE_IMAGE){save.image(file = save.dir)} 
  
  
  
  ## Section 2.8: run non-linear dimensional reduction (UMAP/tSNE)  
  BSC_join <- RunUMAP(BSC_join, dims = 1:15)
  
  ## plot 
  pdf(file = paste0(result.dir, "umap_BSC.pdf"))
  DimPlot(BSC_join, reduction = "umap") %>% print
  dev.off()
  
  ## table the clusters 
  table(BSC_join$old.ident, BSC_join@meta.data$seurat_clusters)
  
  ## re-enter the identities 
  Idents(object = BSC_join) <- BSC_join@meta.data$seurat_clusters
  pdf(file = paste0(result.dir, "umap_BSC_sampleid.pdf"))
  DimPlot(BSC_join, reduction = "umap", label = T, pt.size = 0.2) %>% print
  dev.off()
  
  print("non-linear dimensional reduction done!")
  if(SAVE_IMAGE){save.image(file = save.dir)} 
  if(SAVE_RDS){saveRDS(BSC_join, file = save.rds)} ## save the object 
}


