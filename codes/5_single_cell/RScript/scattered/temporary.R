library(Seurat)
library(SeuratObject)
library(ggplot2)
library(magrittr)
library(dplyr)
library(patchwork)
library(Signac)
library(SingleR)
library(celldex)
library(presto)
library(harmony)
library(multtest)
library(AnnotationHub)
library(ensembldb)
library(glmGamPoi)
library(tidyverse)
library(Matrix)
library(RCurl)
library(scales)
library(cowplot)
library(harmony)
library(knitr)

setwd("~/MRes_project_1/docs/SC_SEQ/ova/processed/combined/plot")

## in total of 5 samples 
GSE <- readRDS("/rds/general/user/ph323/ephemeral/GSE180661/ova_GSE180661.rds")
count_mtx <- GetAssayData(GSE)
metadata_GSE <- GSE@meta.data
GSE <- CreateSeuratObject(counts=count_mtx, assay="RNA", project="GSE180661")

LAM <- readRDS("/rds/general/user/ph323/ephemeral/skin_rds/ov_sc2.rds")
count_mtx <- GetAssayData(LAM)
metadata_LAM <- LAM@meta.data
LAM <- CreateSeuratObject(counts=count_mtx, assay="RNA", project="lambrechtslab")

LAM@meta.data$orig.ident %>% unique %>% print
GSE@meta.data$orig.ident %>% unique %>% print

## Merge datasets 
BSC <- merge(GSE, y=LAM, project="ovarian")

saveRDS(BSC, "~/MRes_project_1/docs/SC_SEQ/ova/processed/combined/ovarian.rds")

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
ggsave("cellCounts.png", plot = plot, width = 20, height = 5, units = "in") 

plot <- BSC@meta.data %>% ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0) + scale_x_log10() + theme_bw() + ylab("Cell density") +
  geom_vline(xintercept = 200)
ggsave("UMIperCell.png", plot = plot, width = 15, height = 5, units = "in")

plot <- BSC@meta.data %>% ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0) + theme_bw() + scale_x_log10() + geom_vline(xintercept = 2500)
ggsave("GenePerCells.png", plot = plot, width = 15, height = 5, units = "in")

plot <- BSC@meta.data %>% ggplot(aes(x=sample, y=log10GenesPerUMI, fill=sample)) + 
  geom_boxplot() + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) + ggtitle("NCells vs NGenes")
ggsave("boxGenePerCells.png", plot = plot, width = 15, height = 5, units = "in")

BSC <- subset(BSC, subset=nUMI>200 & nGene<2500 & percent.mt<15 & log10GenesPerUMI>0.8)
BSC <- NormalizeData(BSC, normalization.method = "LogNormalize", scale.factor = 10000)
BSC <- FindVariableFeatures(BSC, selection.method = "vst", nfeatures = 3000)
BSC <- ScaleData(BSC, vars.to.regress = c("nUMI", "percent.mt", "log10GenesPerUMI"), 
                   do.par = TRUE, num.cores = )

regev_lab_cell_cycle_genes <- read.delim("~/scRNA/files/regev_lab_cell_cycle_genes.txt",header=F)
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
ggsave("BeforeHarmony_sample.png", plot = plot.BeforeHarmony_sample, width = 10, height = 5, units = "in")

## add labels 
table(BSC$old.ident, BSC@meta.data$seurat_clusters)
Idents(object = BSC) <- BSC@meta.data$seurat_clusters
plot.BeforeHarmony_cluster <- DimPlot(BSC, reduction = "umap", label=TRUE) + ggtitle('DimPlot Prior to Harmony')
ggsave("BeforeHarmony_cluster.png", plot=plot.BeforeHarmony_cluster, width = 10, height = 5, units = "in")

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
ggsave("AfterHarmony_sample.png", plot=plot.AfterHarmony_sample, width = 10, height = 5, units = "in")

## add labels 
table(BSC$old.ident, BSC@meta.data$seurat_clusters)
Idents(object = BSC) <- BSC@meta.data$seurat_clusters
plot.AfterHarmony_cluster <- DimPlot(BSC, reduction="umap", group.by = "sample") + ggtitle('DimPlot after harmony')
ggsave("AfterHarmony_cluster.png", plot=plot.AfterHarmony_cluster, width = 10, height = 5, units = "in")

saveRDS(BSC, "~/MRes_project_1/docs/SC_SEQ/ova/processed/combined/ovarian.rds")



















