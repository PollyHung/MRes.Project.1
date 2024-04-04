## Section 1: Environment: ======================================================================================
## Data info: Lambrechtslab 
## data downlaoded from https://cellxgene.cziscience.com/collections/4796c91c-9d8f-4692-be43-347b1727f9d8

## Set control 
LIBRARY <- TRUE
SET_DIRECTORY <- TRUE
SAVE_IMAGE <- TRUE
TYPE <- "ova"
LOAD <- TRUE
SAVE_RDS <- TRUE
PREPROCESS <- FALSE
PROCESS <- FALSE
INFERCNV <- TRUE
RUN_INFERCNV <- FALSE
BATCH_EFFECT <- TRUE

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
  save.dir <-  paste0("~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/", TYPE, ".RData")
  save.rds <- paste0("~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/", TYPE, ".rds") ## TYPE = blood or skin 
  result.dir <- "~/MRes_project_1/docs/SC_SEQ/ova/processed/lambrechtslab/"
  
  ## for infer cnv 
  count.matrix <- "~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/count_matrix/"
  meta.data <- "~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/metadata.csv.gz"
  data.dir <- "~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/"
  inferCNV.dir <- "~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/InferCNV"
  gene.order.file <- "~/MRes_project_1/codes/6_single_cell/reference/gene_order_file.txt"
  output.dir <- "~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/InferCNV/result/"
}

setwd("~/MRes_project_1/docs/SC_SEQ/ova/raw/") ## set directory 
list.files()

print("setting up environment done!")






## Section 2: Preprocess ======================================================================================
if(PREPROCESS){
  ## Section 2.1: Setup the Seurat Object 
  
  ## load the PBMC dataset 
  BSC_data <- Read10X(data.dir = "~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/count_matrix/")
  ## Initialize the Seurat oject with the raw (non-normalized data)
  BSC <-  CreateSeuratObject(counts = BSC_data, min.cells = 3, min.features = 200)
  ## add metadata 
  BSC@meta.data$sample_ids <- BSC@meta.data %>% row.names() %>% str_extract("^[^_]*")   
  metadata <- read.csv(meta.data) ## read metadata 
  rownames(metadata) <- metadata$Cell 
  obj_names <- c("nGene", "nUMI", "CellFromTumor", "PatientNumber", "TumorType", "TumorSite", "CellType")
  metadata <- metadata[rownames(BSC@meta.data), obj_names] ## select columns and rows of interest 
  for(col_name in colnames(metadata)){ ## add metadata to BSC object 
    BSC <- AddMetaData(BSC, metadata = metadata[, col_name, drop = FALSE], col.name = col_name)
  }
  
  print("section 2 set up data done!")
  
  
  
  ## Section 2.2: quality control and normalization  
  ## calculate mitochondrial QC metrics (percentage of counts originating from a set of features)
  BSC[["percent.mt"]] <- PercentageFeatureSet(BSC, pattern = "^MT-")
  
  ## visualize QC metrics and use these to filter cells 
  ## visualisation 
  pdf(file = paste0(result.dir, "preprocess/QC_violin_plot.pdf"), width = 15, height = 5)   # The directory you want to save the file in
  print(VlnPlot(BSC, features = "nFeature_RNA", ncol = 1, pt.size = 0))
  print(VlnPlot(BSC, features = "nCount_RNA", ncol = 1, pt.size = 0))
  print(VlnPlot(BSC, features ="percent.mt", ncol = 1, pt.size = 0))
  dev.off() 
  
  ## filter cells 
  BSC <- subset(BSC, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15) ## nFeature_originalexp
  
  ## normalizing the data 
  ## normalization based on assumption that each cell originally contains same number of RNA molecules 
  BSC <- NormalizeData(BSC, normalization.method = "LogNormalize", scale.factor = 10000, verbose = TRUE)
  
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
  pdf(file = paste0(result.dir, "preprocess/VariableFeaturePlot.pdf"), width = 8, height = 4) 
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
  pdf(file = paste0(result.dir,"preprocess/dimplot_beforeccadj.pdf"), width = 8, height = 8)
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
  pdf(file = paste0(result.dir,"preprocess/dimplot_var.pdf"), width = 8, height = 4)
  VizDimLoadings(BSC_join, dims = 1:2, reduction = "pca") ## visualise features 
  DimPlot(BSC_join, reduction = "pca") ## visualise cells 
  DimHeatmap(BSC_join, dims = 1:15, cells = 500, balanced = TRUE) ## explore if feature is included in downstream analysis 
  ElbowPlot(BSC_join) ## determine dimensionality of the dataset
  dev.off()
  
  ## pca against automatically selected features 
  BSC_join <- RunPCA(BSC_join, features = VariableFeatures(object = BSC_join))
  
  ## print plots 
  pdf(file = paste0(result.dir,"preprocess/dimplot_var_2.pdf"), width = 8, height = 4)
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
  pdf(file = paste0(result.dir, "preprocess/umap_BSC.pdf"))
  DimPlot(BSC_join, reduction = "umap") %>% print
  dev.off()
  
  ## table the clusters 
  table(BSC_join$old.ident, BSC_join@meta.data$seurat_clusters)
  
  ## re-enter the identities 
  Idents(object = BSC_join) <- BSC_join@meta.data$seurat_clusters
  pdf(file = paste0(result.dir, "preprocess/umap_BSC_sampleid.pdf"))
  DimPlot(BSC_join, reduction = "umap", label = T, pt.size = 0.2) %>% print
  dev.off()
  
  print("non-linear dimensional reduction done!")
  if(SAVE_IMAGE){save.image(file = save.dir)} 
  if(SAVE_RDS){saveRDS(BSC_join, file = save.rds)} ## save the object 
}




## Section 3: Batch effects ======================================================================================
if(BATCH_EFFECT){
  BSC <- readRDS(save.rds)
  metadata <- read.csv(meta.data) ## read metadata 
  rownames(metadata) <- metadata$Cell 
  
  obj_names <- c("nGene", "nUMI", "CellFromTumor", "PatientNumber", "TumorType", "TumorSite", "CellType")
  metadata <- metadata[rownames(BSC@meta.data), obj_names] ## select columns and rows of interest 
  
  ## add metadata to BSC object 
  for(col_name in colnames(metadata)){
    BSC <- AddMetaData(BSC, metadata = metadata[, col_name, drop = FALSE], col.name = col_name)
  }
  
  ## examine batch effects 
  batch_by_sample <- DimPlot(BSC, reduction = "umap", group.by = "sample_ids", alpha = 0.1)
  batch_by_patient <- DimPlot(BSC, reduction = "umap", group.by = "PatientNumber", alpha = 0.1)
  elbow_plot <- ElbowPlot(BSC)
  
  ## plot 
  pdf(file = paste0(result.dir, "batch_effect/batch_effects.pdf"), width = 8, height = 8)
  print(batch_by_sample)
  print(batch_by_patient)
  dev.off()
  
  ## remove batch effects 
  ## regress out by patient number? 
  BSC_batch.cor2 <- ScaleData(BSC, vars.to.regress = c("PatientNumber", "sample_ids"))
  BSC_batch.cor2 <- FindNeighbors(BSC_batch.cor2, reduction = "pca", dims = 1:20)
  BSC_batch.cor2 <- FindClusters(BSC_batch.cor2, resolution = 1)
  BSC_batch.cor2 <- RunUMAP(BSC_batch.cor2, dims = 1:20)
  
}













## Section 3: Process ======================================================================================
if(PROCESS){
  ## a small function to remove all non-vector objects from enviroment 
  for(i in ls()){
    if(!is.vector(get(i))){
      print(paste(i, "is not vector"))
      rm(list=i)
    }
  }
  print("Environment cleared! Start from Section 3!")
  
  ## Section 3.0 load data 
  BSC <- readRDS(save.rds)
  gc()
  print("finish loading BSC_join into environment")
  
  #BSC@meta.data$sample_ids <- sub("_.*", "", rownames(BSC@meta.data)) 
  
  ## Section 3.2 automated cluster annotation 
  ## Get reference datasets from celldex package. Note that there are two cell type assignments, label.main and label.fine. 
  ## We’re only going to run the annotation against the Monaco Immune Database
  MI <- celldex::MonacoImmuneData(ensembl = FALSE, cell.ont = "all") ## monaco immune 
  HPCA <- celldex::HumanPrimaryCellAtlasData(ensembl = FALSE, cell.ont = "all") ## human primary cell atlas 
  DICE <- celldex::DatabaseImmuneCellExpressionData(ensembl = FALSE, cell.ont = "all") ## database immune cell expression 
  BE <- celldex::BlueprintEncodeData(ensembl = FALSE, cell.ont = "all") ## blueprint encode data 
  
  MI$label.main <- paste0("MI.", MI$label.main) ## monaco immune data 
  HPCA$label.main <- paste0("HPCA.", HPCA$label.main) ## human primary cell atlas data 
  DICE$label.main <- paste0("DICE.", DICE$label.main) ## database immune cell expression 
  BE$label.main <- paste0("BE.", BE$label.main) ## blueprint encode data 
  
  shared <- rownames(MI) %>% intersect(rownames(DICE)) %>% intersect(rownames(BE))  %>% intersect(rownames(HPCA))
  combined <- cbind(MI[shared, ], DICE[shared, ], HPCA[shared, ], BE[shared, ])  
  
  ## Let’s convert our Seurat object to single cell experiment (SCE) for convenience. After this, using SingleR becomes very easy:
  #BSC_sce <- as.SingleCellExperiment(DietSeurat(BSC))
  
  ## Get parameters for experiment 
  testdata <- GetAssayData(BSC, layer="data")
  clusters <- BSC@meta.data$seurat_clusters
  
  ## perform SingleR
  cellpred <- SingleR(test =  testdata,  ## BSC_sce,  ## cellpred
                      ref = combined, 
                      labels = combined$label.fine, 
                      method = "cluster", 
                      clusters = clusters, 
                      assay.type.ref = "logcounts", 
                      assay.type.test = "logcounts")
  
  if(SAVE_IMAGE){save.image(file = save.dir)}
  print("finish finding markers")
  
  
  
  ## Section 3.3 view automated clustering result 
  celltype <- data.frame(ClusterID=rownames(cellpred), 
                         celltype=cellpred$labels, 
                         stringAsFactors=FALSE)
  ## save dataframe 
  write.csv(celltype, paste0(result.dir, "celltype_BSC_singleR.csv"), row.names=F)
  
  ## code to see individual features 
  #FeaturePlot(BSC, "CD38") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
  if(SAVE_IMAGE){save.image(file = save.dir)}
  print("finish viewing automated clustering result")
  
  ## Section 3.4 cell labelling 
  Idents(object = BSC) <- BSC@meta.data$seurat_clusters
  celltype <- read.csv(file = paste0(result.dir, "celltype_BSC_singleR.csv"), row.names = 1)
  
  new.cluster.ids <- celltype$celltype 
  names(new.cluster.ids) <- levels(BSC)
  #names(new.cluster.ids) <- BSC@meta.data$seurat_clusters
  
  BSC <- RenameIdents(BSC, new.cluster.ids)
  
  pdf(file = paste0(result.dir, "umap_BSC_celllabel.pdf"), width = 10, height = 10)
  print(DimPlot(BSC, reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend())
  dev.off()
  
  if(SAVE_IMAGE){save.image(file = save.dir)}
  if(SAVE_RDS){saveRDS(BSC, file = save.rds)}
}




## Section 4: Infer CNV ====================================================================================================
if(INFERCNV){
  
  BSC <- readRDS(save.rds) ## load in data 
  
  metadata <- read.csv(meta.data) ## read metadata 
  rownames(metadata) <- metadata$Cell 
  
  obj_names <- c("nGene", "nUMI", "CellFromTumor", "PatientNumber", "TumorType", "TumorSite", "CellType")
  metadata <- metadata[rownames(BSC@meta.data), obj_names] ## select columns and rows of interest 
  
  ## add metadata to BSC object 
  for(col_name in colnames(metadata)){
    BSC <- AddMetaData(BSC, metadata = metadata[, col_name, drop = FALSE], col.name = col_name)
  }
  
  ## Split Seurat objects based on orig.ident
  obj.list <- SplitObject(BSC, split.by = "orig.ident")
  
  ## Initialize list of store updated objects 
  updated.obj.list <- list()
  
  ## file 1: gene order file 
  gene_order_file <- "~/MRes_project_1/codes/6_single_cell/reference/gene_order_file.txt"
  sdrf <- read.table("~/MRes_project_1/codes/6_single_cell/reference/E-MTAB-8107.sdrf.txt", sep = "\t", header = TRUE)
  
  for(sample_name in names(obj.list)){ ## 
    sample_obj <- obj.list[[sample_name]] ## get the object 
    print(sample_obj) 
    
    ## file 2: count matrix 
    sub_matrix <- as.matrix(GetAssayData(sample_obj, assay = "RNA"))
    print(dim(sub_matrix)) ## dimension? 
    
    ## file 3: annotation file 
    sub_metadata <- sample_obj@meta.data
    sub_metadata$barcode <- row.names(sub_metadata)
    write.table(sub_metadata[, c("barcode", "CellType")], 
                file = paste0(inferCNV.dir, sample_name, "_cellAnnotations.txt"), 
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    ## file 4: reference group names 
    ref_group_names <- setdiff(unique(sub_metadata$CellType), "Cancer") ## remove the cancer cells, because the reference has to be normal cell types 
    
    ## file 5: define output directory for this sample 
    output_dir <- paste0(output.dir, sample_name)
    
    ## Create infer CNV object 
    if(RUN_INFERCNV){
      infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = sub_matrix, 
                                           annotations_file = paste0(inferCNV.dir, sample_name, "_cellAnnotations.txt"), 
                                           gene_order_file = gene_order_file, 
                                           ref_group_names = ref_group_names) ## set to various normal cell types 
      
      ## Run infer CNV analysis 
      infercnv_obj = infercnv::run(infercnv_obj,
                                   cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                   out_dir=output_dir,  # dir is auto-created for storing outputs
                                   cluster_by_groups=T,   # cluster
                                   denoise=T,
                                   HMM=T, 
                                   num_threads = 10, 
                                   write_expr_matrix = TRUE)
    }
    ## Add inferCNV results to Seurat object 
    updated_sample_obj <- add_to_seurat(seurat_obj = sample_obj, infercnv_output_path = output_dir)
    
    ## Store updated object 
    updated.obj.list[[sample_name]] <- updated_sample_obj
  }
  
  BSC <- merge(updated.obj.list[[1]],y = updated.obj.list[2:10], 
               add.cell.ids = c("BT1303", "BT1304", "BT1305", "BT1306", "BT1307", "scrSOL001", "scrSOL003", "scrSOL004", "scrSOL006", "scrSOL007"), 
               project = "SeuratProject")
  #saveRDS(BSC, "cnv_ova_lamb.rds")
}











