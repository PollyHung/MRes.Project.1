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
library(fuzzyjoin)

TYPE <- "ova"
CHR_AMR <- "chr16p"

cyto.band <- "~/MRes_project_1/codes/6_single_cell/reference/cytoBand.txt"
cnv.file <- "HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat"
output.dir <- "~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/InferCNV/result"

save.dir <-  "~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/ova.RData"
save.rds <- "~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/feb_15_ova.rds" 
result.dir <- "~/MRes_project_1/docs/SC_SEQ/ova/processed/lambrechtslab/"
rdata <- "~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/RData/"
count.matrix <- "~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/count_matrix/"
meta.data <- "~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/metadata.csv.gz"
data.dir <- "~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/"
inferCNV.dir <- "~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/InferCNV/"
gene.order.file <- "~/MRes_project_1/codes/6_single_cell/reference/gene_order_file.txt"



summary_table <- read.table("~/MRes_project_1/docs/SC_SEQ/ova/processed/lambrechtslab/InferCNV/chr2q10p.txt", 
                            header = TRUE)
chr2q <- summary_table %>% dplyr::filter(chr_arm == "chr2q")
chr10p <- summary_table %>% dplyr::filter(chr_arm == "chr10p")

chr2q$status <- ifelse(chr2q$aneu_score>mean(chr2q$aneu_score), "amplification", "neutral")
chr10p$status <- ifelse(chr10p$aneu_score>mean(chr2q$aneu_score), "amplification", "neutral")

summary_table <- chr2q

amplification <- summary_table %>% dplyr::filter(status == "amplification") %>% dplyr::select(sample) %>% unlist()
neutral <- summary_table %>% dplyr::filter(status == "neutral") %>% dplyr::select(sample) %>% unlist()


BSC_data <- Read10X(data.dir="~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/count_matrix/")
BSC_orig <- CreateSeuratObject(counts=BSC_data, 
                               min.cells=3, 
                               min.features=200)
nFeature_RNA <- rowSums(BSC_orig@assays$RNA$counts > 0)
nCount_RNA <- rowSums(BSC_orig@assays$RNA$counts)
BSC_orig@meta.data$nFeature_RNA <- nFeature_RNA
BSC_orig@meta.data$nCount_RNA <- nCount_RNA


BSC_mut <- subset(x = BSC_orig, subset = orig.ident %in% amplification) 
BSC_wild <- subset(x = BSC_orig, subset = orig.ident %in% neutral) 


for(dataset in c("BSC_mut", "BSC_wild")){
## auto annotation -----------------------------------------------------------------------------------

  print(paste("This run is for", dataset))
  BSC <- get(dataset, envir = globalenv())
  
  BSC[["percent.mt"]] <- PercentageFeatureSet(BSC, pattern = "^MT-")
  plot.1.1 <- VlnPlot(BSC, features = "percent.mt", ncol = 1, pt.size = 0)
  
  BSC <- subset(BSC, subset=nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)
  BSC <- NormalizeData(BSC, normalization.method = "LogNormalize", scale.factor = 10000, verbose = TRUE)
  
  BSC <- FindVariableFeatures(BSC, selection.method = "vst", nfeatures = 3000)
  top10 <- head(VariableFeatures(BSC), 10)
  plot.2.1 <- VariableFeaturePlot(BSC)
  plot.2.2 <- LabelPoints(plot = plot.2.1, points = top10, repel = TRUE)
  
  BSC <- ScaleData(BSC, vars.to.regress = c("nCount_RNA", "percent.mt"), 
                   do.par = TRUE, num.cores = 6)
  
  cell_cycles <- read.delim("~/MRes_project_1/codes/6_single_cell/cell_cycle_adjust/regev_lab_cell_cycle_genes.txt", 
                            header = FALSE)
  s.genes <- cell_cycles[1:43, 1]
  g2m.genes <- cell_cycles[44:97, 1]
  
  BSC <- JoinLayers(BSC)
  BSC <- CellCycleScoring(object = BSC, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  head.1 <- head(BSC@meta.data)
  
  BSC <- RunPCA(BSC, features = c(s.genes, g2m.genes))
  plot.3.1 <- VizDimLoadings(BSC, dims = 1:2, reduction = "pca")
  plot.3.2 <- DimPlot(BSC)
  plot.3.3 <- DimHeatmap(BSC, dims = 1:9, cells = 500, balanced = TRUE)
  plot.3.4 <- ElbowPlot(BSC)
  
  BSC <- ScaleData(BSC, vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"), 
                   do.par = TRUE, num.cores = 6)
  BSC <- RunPCA(BSC, features = c(s.genes, g2m.genes))
  plot.4.1 <- VizDimLoadings(BSC, dims = 1:2, reduction = "pca")
  plot.4.2 <- DimPlot(BSC, reduction = "pca")
  plot.4.3 <- DimHeatmap(BSC, dims = 1:15, cells = 500, balanced = TRUE)
  plot.4.4 <- ElbowPlot(BSC)
  
  BSC <- RunPCA(BSC, features = VariableFeatures(BSC))
  plot.5.1 <- VizDimLoadings(BSC, dims = 1:2, reduction = "pca")
  plot.5.2 <- DimPlot(BSC, reduction = "pca")
  plot.5.3 <- DimHeatmap(BSC, dims = 1:15, cells = 500, balanced = TRUE)
  plot.5.4 <- ElbowPlot(BSC)
  
  BSC <- FindNeighbors(BSC, dims = 1:30)
  BSC <- FindClusters(BSC, resolution = 0.5)
  head(Idents(BSC), 5) ## view the cluster ids for the first 5 cells, there are 25 levels 
  
  BSC <- RunUMAP(BSC, dims = 1:30)
  plot.6.1 <- DimPlot(BSC, reduction = "umap")
  
  table(BSC$old.ident, BSC@meta.data$seurat_clusters)
  Idents(BSC) <- BSC@meta.data$seurat_clusters
  plot.6.2 <- DimPlot(BSC, reduction = "umap", label = TRUE, pt.size = 0.2)
  plot.6.3 <- DimPlot(BSC, reduction = "umap", group.by = "orig.ident", alpha = 0.1)
  
  saveRDS(BSC, paste0(paste0(rdata, dataset, "_", CHR_AMR, ".rds"))) ## 存档！！
  
  BSC.markers <- FindAllMarkers(BSC, only.pos = TRUE) 
  BSC.markers.filt <- BSC.markers %>% group_by(cluster) %>%
    dplyr::filter(p_val_adj < 0.005) %>% ungroup() %>% 
    arrange(desc(cluster), desc(avg_log2FC)) 
  
  metadata <- read.csv(meta.data, row.names = 1)
  metadata <- metadata[rownames(BSC@meta.data), ]
  BSC <- AddMetaData(BSC, metadata = metadata, col.name = colnames(metadata))
  plot.7.1 <- DimPlot(BSC, reduction = "umap", group.by = c("PatientNumber", "orig.ident", 
                                                            "seurat_clusters", "CellType"), alpha = 0.1)
  
  MI <- celldex::MonacoImmuneData(ensembl = FALSE, cell.ont = "all") 
  HPCA <- celldex::HumanPrimaryCellAtlasData(ensembl = FALSE, cell.ont = "all") 
  DICE <- celldex::DatabaseImmuneCellExpressionData(ensembl = FALSE, cell.ont = "all") 
  BE <- celldex::BlueprintEncodeData(ensembl = FALSE, cell.ont = "all") 
  MI$label.main <- paste0("MI.", MI$label.main) 
  HPCA$label.main <- paste0("HPCA.", HPCA$label.main) 
  DICE$label.main <- paste0("DICE.", DICE$label.main) 
  BE$label.main <- paste0("BE.", BE$label.main) 
  shared <- rownames(MI) %>% intersect(rownames(DICE)) %>% intersect(rownames(BE)) %>% intersect(rownames(HPCA))
  combined <- cbind(MI[shared, ], DICE[shared, ], HPCA[shared, ], BE[shared, ])  
  
  testdata <- GetAssayData(BSC, layer="data")
  clusters <- BSC@meta.data$seurat_clusters
  
  cellpred <- SingleR(test = testdata, 
                      ref = combined, 
                      labels = combined$label.fine, 
                      method = "cluster", 
                      clusters = clusters, 
                      assay.type.ref = "logcounts", 
                      assay.type.test = "logcounts")
  
  celltype <- data.frame(ClusterID=rownames(cellpred), 
                         celltype=cellpred$labels, 
                         stringAsFactors=FALSE)
  
  Idents(object = BSC) <- BSC@meta.data$seurat_clusters
  new.cluster.ids <- celltype$celltype 
  names(new.cluster.ids) <- levels(BSC)
  
  BSC <- RenameIdents(BSC, new.cluster.ids)
  plot.7.2 <- DimPlot(BSC, reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend()
  
  plot_names <- ls(pattern = "^plot.*", envir = .GlobalEnv)
  pdf("~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/InferCNV/examples/chr16p_mut.pdf")
  for (obj_name in plot_names) {
    obj <- get(obj_name)
    print(obj)
  }
  dev.off()
  
  
  #save.image(paste0(rdata, dataset, "_", CHR_AMR, ".RData"))
}


load("~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/RData/BSC_mut_chr16p.RData")





