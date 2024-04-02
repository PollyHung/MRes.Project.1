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


## Perform normalisation separately and then perform analysis with integration 
setwd("/rds/general/user/ph323/ephemeral/GSE131907/")
list.files()

count.mtx <- readRDS("GSE131907_Lung_Cancer_raw_UMI_matrix.rds")
metadata <- read.table("GSE131907_Lung_Cancer_cell_annotation.txt", sep="\t", header=TRUE, row.names = 1) 
feature.summary <- readxl::read_xlsx("GSE131907_Lung_Cancer_Feature_Summary.xlsx", 
                                     skip = 1, n_max = 59, col_names=TRUE) %>% as.data.frame
BSC <- CreateSeuratObject(counts = count.mtx, assay = "RNA", project = "GSE131907", 
                          min.cells = 3, min.features = 200)

setwd("~/MRes_project_1/docs/SC_SEQ/lung/raw/GSE131907/")
list.files()


## Perform normalisation separately and then perform analysis with integration 
BSC[["log10GenesPerUMI"]] <- log10(BSC$nFeature_RNA)/log10(BSC$nCount_RNA)
BSC[["percent.mt"]] <- PercentageFeatureSet(BSC, pattern = "^MT-")
BSC[["mitoRatio"]] <- BSC@meta.data$percent.mt/100 
BSC[["sample"]] <- BSC@meta.data$orig.ident
BSC@meta.data <- BSC@meta.data %>% dplyr::rename(nUMI=nCount_RNA, nGene=nFeature_RNA)

plot <- BSC@meta.data %>% ggplot(aes(x=sample, fill=sample)) + geom_bar() + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) + 
  ggtitle("Cell Count per Sample for combined dataset") + NoLegend()
ggsave("plots/cellCounts.png", plot = plot, width = 20, height = 5, units = "in") 

plot <- BSC@meta.data %>% ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0) + scale_x_log10() + theme_bw() + ylab("Cell density") +
  geom_vline(xintercept = 200)
ggsave("plots/UMIperCell.png", plot = plot, width = 15, height = 5, units = "in")

plot <- BSC@meta.data %>% ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0) + theme_bw() + scale_x_log10() + geom_vline(xintercept = 2500)
ggsave("plots/GenePerCells.png", plot = plot, width = 15, height = 5, units = "in")

plot <- BSC@meta.data %>% ggplot(aes(x=sample, y=log10GenesPerUMI, fill=sample)) + 
  geom_boxplot() + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) + ggtitle("NCells vs NGenes")
ggsave("plots/boxGenePerCells.png", plot = plot, width = 15, height = 5, units = "in")

BSC <- subset(BSC, subset=nUMI>300 & nGene<3000 & percent.mt<15 & log10GenesPerUMI>0.8)
BSC <- NormalizeData(BSC, normalization.method = "LogNormalize", scale.factor = 10000)
BSC <- FindVariableFeatures(BSC, selection.method = "vst", nfeatures = 3000)
BSC <- ScaleData(BSC, vars.to.regress = c("percent.mt", "nUMI"), 
                 do.par = TRUE, num.cores = 3)

regev_lab_cell_cycle_genes <- read.delim("~/MRes_project_1/codes/6_single_cell/reference/regev_lab_cell_cycle_genes.txt",header=F)
s.genes <- regev_lab_cell_cycle_genes[1:43,1]
g2m.genes <- regev_lab_cell_cycle_genes[44:97,1]
BSC <- JoinLayers(BSC)
BSC <- CellCycleScoring(object = BSC, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

## perform first pca -----------------------------------------------------------------------------------
BSC <- RunPCA(BSC, features = c(s.genes, g2m.genes))

## scale data again for cell cycle factors and regress out these and variable features using pca 
BSC <- ScaleData(BSC, vars.to.regress = c("S.Score", "G2M.Score"),
                 do.par = TRUE, num.cores = 3)
BSC <- RunPCA(BSC, features = c(s.genes, g2m.genes))
BSC <- RunPCA(BSC, features = VariableFeatures(object = BSC))

## run umap      
BSC <- RunUMAP(BSC, dims = 1:15)
plot.BeforeHarmony_sample <- DimPlot(BSC, reduction = "umap", group.by = "sample") + ggtitle('DimPlot Prior to Harmony') 
ggsave("plots/BeforeHarmony_sample.png", plot = plot.BeforeHarmony_sample, width = 10, height = 5, units = "in")

saveRDS(BSC, "/rds/general/user/ph323/ephemeral/GSE180661/updated_lung.rds")

## add labels 
#table(BSC$old.ident, BSC@meta.data$seurat_clusters)
#Idents(object = BSC) <- BSC@meta.data$seurat_clusters
#plot.BeforeHarmony_cluster <- DimPlot(BSC, reduction = "umap", label=TRUE) + ggtitle('DimPlot Prior to Harmony')
#ggsave("plots/BeforeHarmony_cluster.png", plot=plot.BeforeHarmony_cluster, width = 10, height = 5, units = "in")

## run harmony to remove batch effect from sequencing company and sample type -----------------------------
BSC <- RunHarmony(object = BSC,
                  group.by.vars = c("sample"),
                  assay.use = "RNA", 
                  plot_convergence = TRUE)

## redo standard workflow with reduction changed to harmony 
BSC <- FindNeighbors(BSC, dims = 1:15, reduction="harmony")
BSC <- FindClusters(BSC, resolution = 0.5, reduction="harmony")
BSC <- RunUMAP(BSC, dims = 1:15, reduction="harmony")
plot.AfterHarmony_sample <- DimPlot(BSC, reduction="umap", group.by = "sample") + ggtitle('DimPlot after harmony')
ggsave("plots/AfterHarmony_sample.png", plot=plot.AfterHarmony_sample, width = 10, height = 5, units = "in")

saveRDS(BSC, "/rds/general/user/ph323/ephemeral/GSE180661/updated_lung.rds")

## add labels 
table(BSC$old.ident, BSC@meta.data$seurat_clusters)
Idents(object = BSC) <- BSC@meta.data$seurat_clusters
plot.AfterHarmony_cluster <- DimPlot(BSC, reduction="umap", group.by = "sample") + ggtitle('DimPlot after harmony')
ggsave("AfterHarmony_cluster.png", plot=plot.AfterHarmony_cluster, width = 10, height = 5, units = "in")

saveRDS(BSC, "/rds/general/user/ph323/ephemeral/GSE180661/updated_lung.rds")



















