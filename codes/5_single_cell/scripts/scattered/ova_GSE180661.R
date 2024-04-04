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


if(FALSE){
  setwd("/rds/general/user/ph323/ephemeral/GSE180661/")
  
  BSC <- readRDS("ova_GSE180661.rds")
  count.mtx <- GetAssayData(BSC)
  metadata <- BSC@meta.data
  rm(BSC)
  
  BSC <- CreateSeuratObject(counts = count.mtx, assay = "RNA", project = "GSE180661", 
                            min.cells = 3, min.features = 200)
  
  setwd("~/MRes_project_1/docs/SC_SEQ/ova/raw/GSE180661/")
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
  
  BSC <- RunPCA(BSC, features = c(s.genes, g2m.genes))
  
  ## scale data again for cell cycle factors and regress out these and variable features using pca 
  BSC <- ScaleData(BSC, vars.to.regress = c("S.Score", "G2M.Score"),
                   do.par = TRUE, num.cores = 3)
  BSC <- RunPCA(BSC, features = c(s.genes, g2m.genes))
  BSC <- RunPCA(BSC, features = VariableFeatures(object = BSC))
  
  #saveRDS(BSC, "/rds/general/user/ph323/ephemeral/GSE180661/updated_ova.rds")
  BSC <- readRDS("/rds/general/user/ph323/ephemeral/GSE180661/updated_ova.rds")
  
  
  ## run umap      
  BSC <- FindNeighbors(BSC, dims = 1:15, reduction="pca")
  BSC <- FindClusters(BSC, resolution = 0.5, reduction="pca")
  BSC <- RunUMAP(BSC, dims = 1:15, reduction="pca")
  
  table(BSC$old.ident, BSC@meta.data$seurat_clusters)
  Idents(object = BSC) <- BSC@meta.data$seurat_clusters
  
  plot.AfterHarmony_sample <- DimPlot(BSC, reduction="umap", label = TRUE) + ggtitle('DimPlot after harmony')
  ggsave("plots/AfterHarmony_sample.png", plot=plot.AfterHarmony_sample, width = 10, height = 5, units = "in")
  
  saveRDS(BSC, "/rds/general/user/ph323/ephemeral/GSE180661/updated_ova.rds")
  
  ## step 2: --------------------------------------------------------------------------------------------------
  BSC <- readRDS("/rds/general/user/ph323/ephemeral/GSE180661/ova_GSE180661.rds")
  metadata <- BSC@meta.data
  
  rm(BSC)
  gc()
  
  BSC <- readRDS("/rds/general/user/ph323/ephemeral/GSE180661/updated_ova.rds")
  
  
  ## step 3: ---------------------------------------------------------------------------------------------------
  metadata_2 <- BSC@meta.data
  metadata_2$cell_id <- rownames(metadata_2)
  
  rownames(metadata) <- metadata$cell_id
  metadata <- metadata[metadata_2$cell_id, ] ## realign 
  
  pdf("~/MRes_project_1/docs/SC_SEQ/ova/raw/GSE180661/plots/compare.pdf", width=8, height=8)
  plot(metadata$G2M.Score, metadata_2$G2M.Score, type="p", pch='.')
  plot(metadata$percent.mt, metadata_2$percent.mt, type="p", pch='.')
  plot(metadata$S.Score, metadata_2$S.Score, type="p", pch='.')
  dev.off()
  
  metadata.1 <- metadata %>% dplyr::select("cell_id", "cell_type", "cell_type_super", 
                                           "sort_parameters", "percent.rb", "percent.mt", 
                                           "CC.Diff", "isabl_experiment_system_id") 
  
  metadata_2 <- metadata_2 %>% dplyr::select("cell_id", "orig.ident", "nUMI", "nGene", "log10GenesPerUMI", 
                                             "S.Score", "G2M.Score", "Phase", "RNA_snn_res.0.5", "seurat_clusters") %>% 
    dplyr::rename(new_clusters = seurat_clusters)
  
  metadata.1 <- merge(metadata.1, metadata_2, by="cell_id")
  metadata.1$sort_parameters <- gsub("singlet, live, ", "", metadata.1$sort_parameters)
  rownames(metadata.1) <- metadata.1$cell_id
  metadata.1 <- metadata.1[metadata_2$cell_id, ]
  
  cnv <- read.table("~/MRes_project_1/docs/SC_SEQ/ova/raw/GSE180661/CNVbyGISTIC/broad_values_by_arm.txt", sep="\t", header=TRUE)
  cnv <- cnv %>% t %>% as.data.frame 
  colnames(cnv) <- cnv[1, ]
  cnv <- cnv[2:nrow(cnv), ]
  cnv_arm2 <- cnv[, c("2p", "2q")]
  cnv_arm2$sample_id <- gsub("\\.", "-", rownames(cnv_arm2))
  
  clinical.mtx <- read_tsv("~/MRes_project_1/docs/SC_SEQ/ova/raw/GSE180661/clinical_data.tsv") %>% as.data.frame
  colnames(clinical.mtx) <- gsub(" ", "_", colnames(clinical.mtx)) %>% tolower
  clinical.mtx <- clinical.mtx %>% dplyr::select("sample_id", "patient_display_name", 
                                                 "patient_age_at_diagnosis", "fraction_genome_altered", 
                                                 "tmb_(nonsynonymous)", "tumor_purity")
  
  cnv_arm2 <- merge(clinical.mtx, cnv_arm2, by="sample_id")
  cnv_arm2 <- cnv_arm2 %>% dplyr::rename(orig.ident=patient_display_name)
  cnv_arm2 <- cnv_arm2[which(cnv_arm2$orig.ident %in% unique(metadata.1$orig.ident)), ] 
  cnv_arm2.WGS <- cnv_arm2[which(grepl("^SHAH", cnv_arm2$sample_id)), ]
  cnv_arm2.WGS$tumor_purity <- NULL
  
  metadata.2 <- left_join(metadata.1, cnv_arm2.WGS, by="orig.ident")
  rownames(metadata.2) <- metadata.2$cell_id
  metadata.2 <- metadata.2[rownames(metadata_2), ]
  
  BSC@meta.data <- metadata.2
  
  plot.Labelled <- DimPlot(BSC, group.by = "cell_type") + ggtitle('DimPlot after harmony')
  ggsave("plots/labelled.png", plot=plot.Labelled, width = 8, height = 5, units = "in")
  
  write.table(metadata.2, "~/MRes_project_1/codes/6_single_cell/reference/ova_GSE180661_cnv_metadata.txt", 
              sep="\t", quote=FALSE, row.names = TRUE, col.names = TRUE)
  saveRDS(BSC, "/rds/general/user/ph323/ephemeral/GSE180661/cnv_appended_ova.rds")
  
  
  ## Step 5: Add the cnv level, this is ovarian so we should focus on: 2q! -----------------------------------
  BSC <- readRDS("/rds/general/user/ph323/ephemeral/GSE180661/cnv_appended_ova.rds")
  metadata = read.table("~/MRes_project_1/codes/6_single_cell/reference/ova_GSE180661_cnv_metadata.txt", header=TRUE, row.names=1)
  
  
  BSC@meta.data$`2q` <- as.numeric(BSC@meta.data$`2q`)
  amp_2q <- subset(BSC, subset = `2q` > 0)
  wt_2q <- subset(BSC, subset=`2q` <= 0)
  
  saveRDS(amp_2q, "/rds/general/user/ph323/ephemeral/GSE180661/amp2q.rds")
  saveRDS(wt_2q, "/rds/general/user/ph323/ephemeral/GSE180661/wt2q.rds")
}


## step 6: subsetting -------------------------------------------------------------------------------------------
amp_2q <- readRDS("/rds/general/user/ph323/ephemeral/GSE180661/amp2q.rds")
wt_2q <- readRDS("/rds/general/user/ph323/ephemeral/GSE180661/wt2q.rds")
print("read in complete")

amp_2q.Tcell <- subset(x=amp_2q, subset=cell_type %in% c("T.cell", "Dendritic.cell"))
amp_2q.Bcell <- subset(x=amp_2q, subset=cell_type %in% c("B.cell", "Plasma.cell"))
amp_2q.Myeloid <- subset(x=amp_2q, subset=cell_type %in% c("Mast.cell", "Myeloid.cell"))
amp_2q.Cancer <- subset(x=amp_2q, subset=cell_type %in% c("Ovarian.cancer.cell"))
print("amp_2q subset complete")

wt_2q.Tcell <- subset(x=wt_2q, subset=cell_type %in% c("T.cell", "Dendritic.cell"))
wt_2q.Bcell <- subset(x=wt_2q, subset=cell_type %in% c("B.cell", "Plasma.cell"))
wt_2q.Myeloid <- subset(x=wt_2q, subset=cell_type %in% c("Mast.cell", "Myeloid.cell"))
wt_2q.Cancer <- subset(x=wt_2q, subset=cell_type %in% c("Ovarian.cancer.cell"))
print("wt_2q subset complete")

save(list=c("wt_2q.Tcell", "wt_2q.Bcell", "wt_2q.Myeloid", "wt_2q.Cancer"), 
     file = "/rds/general/user/ph323/ephemeral/GSE180661/wt_2q_immune.RData")
save(list=c("amp_2q.Tcell", "amp_2q.Bcell", "amp_2q.Myeloid", "amp_2q.Cancer"), 
     file = "/rds/general/user/ph323/ephemeral/GSE180661/mut_2q_immune.RData")
print("saving the RData")


## Step 7: re-normalise T cell ----------------------------------------------------------------------------------------
if(FALSE){
  amp_2q.Tcell <- NormalizeData(object = amp_2q.Tcell, normalization.method = "LogNormalize", scale.factor = 10000)
  amp_2q.Tcell <- FindVariableFeatures(object = amp_2q.Tcell, selection.method = "vst", nfeatures = 3000)
  amp_2q.Tcell <- ScaleData(object = amp_2q.Tcell, vars.to.regress = c("percent.mt", "nUMI"))
  amp_2q.Tcell <- RunPCA(object = amp_2q.Tcell, features = VariableFeatures(object = amp_2q.Tcell))
  amp_2q.Tcell <- FindNeighbors(object = amp_2q.Tcell, dims = 1:15, reduction="pca")
  amp_2q.Tcell <- FindClusters(object = amp_2q.Tcell, resolution = 0.5, reduction="pca")
  amp_2q.Tcell <- RunUMAP(amp_2q.Tcell, dims = 1:15, reduction="pca")
  print("re-do seurat on amp_2q_Tcells")
  
  wt_2q.Tcell <- NormalizeData(object = wt_2q.Tcell, normalization.method = "LogNormalize", scale.factor = 10000)
  wt_2q.Tcell <- FindVariableFeatures(object = wt_2q.Tcell, selection.method = "vst", nfeatures = 3000)
  wt_2q.Tcell <- ScaleData(object = wt_2q.Tcell, vars.to.regress = c("percent.mt", "nUMI"))
  wt_2q.Tcell <- RunPCA(object = wt_2q.Tcell, features = VariableFeatures(object = wt_2q.Tcell))
  wt_2q.Tcell <- FindNeighbors(object = wt_2q.Tcell, dims = 1:15, reduction="pca")
  wt_2q.Tcell <- FindClusters(object = wt_2q.Tcell, resolution = 0.5, reduction="pca")
  wt_2q.Tcell <- RunUMAP(wt_2q.Tcell, dims = 1:15, reduction="pca")
  print("re-do seurat on wt_2q_Tcells")
  
  saveRDS(amp_2q.Tcell, "/rds/general/user/ph323/ephemeral/GSE180661/amp_2q_Tcell.rds")
  saveRDS(wt_2q.Tcell, "/rds/general/user/ph323/ephemeral/GSE180661/wt_2q_Tcell.rds")
  
}



## Step 8: process amp_2q.Tcell ----------------------------------------------------------------------------
amp_2q.Tcell <- readRDS("/rds/general/user/ph323/ephemeral/GSE180661/amp_2q_Tcell.rds")
wt_2q.Tcell <- readRDS("/rds/general/user/ph323/ephemeral/GSE180661/wt_2q_Tcell.rds")

plot.amp_2q_Tcell_clusters <- DimPlot(amp_2q.Tcell, label = TRUE)
plot.wt_2q_Tcell_clusters <- DimPlot(wt_2q.Tcell, label = TRUE)
ggsave("~/MRes_project_1/docs/SC_SEQ/ova/raw/GSE180661/plots/amp_2q_Tcell_clusters.png", plot.amp_2q_Tcell_clusters, 
       width = 6, height = 5, units = "in")
ggsave("~/MRes_project_1/docs/SC_SEQ/ova/raw/GSE180661/plots/wt_2q_Tcell_clusters.png", plot.wt_2q_Tcell_clusters, 
       width = 6, height = 5, units = "in")

amp_2q.Tcell.markers <- FindAllMarkers(object = amp_2q.Tcell, 
                                       slot = "data")
amp_2q.Tcell.markers$diff_pct <- abs(amp_2q.Tcell.markers$pct.1 - amp_2q.Tcell.markers$pct.2)
amp_2q.Tcell.markers.filt <- amp_2q.Tcell.markers %>% 
  dplyr::filter(p_val_adj<0.01, avg_log2FC>1) %>% 
  dplyr::arrange(cluster, desc(pct.1), desc(avg_log2FC)) 
write.csv(amp_2q.Tcell.markers.filt, "~/MRes_project_1/docs/SC_SEQ/ova/raw/GSE180661/amp_2q_Tcell_markers_filt.csv")
amp_2q.Tcell.markers.filt <- amp_2q.Tcell.markers.filt %>% group_by(cluster) %>% slice_head(n = 10)

wt_2q.Tcell.markers <- FindAllMarkers(object = wt_2q.Tcell, 
                                      slot = "data")
wt_2q.Tcell.markers$diff_pct <- abs(wt_2q.Tcell.markers$pct.1 - wt_2q.Tcell.markers$pct.2)
wt_2q.Tcell.markers.filt <- wt_2q.Tcell.markers %>% 
  dplyr::filter(p_val_adj<0.01, avg_log2FC>1) %>% 
  dplyr::arrange(cluster, desc(pct.1), desc(avg_log2FC)) 
write.csv(wt_2q.Tcell.markers.filt, "~/MRes_project_1/docs/SC_SEQ/ova/raw/GSE180661/wt_2q_Tcell_markers_filt.csv")
wt_2q.Tcell.markers.filt <- wt_2q.Tcell.markers.filt %>% group_by(cluster) %>% slice_head(n = 10)




## Step 9: SingleR -------------------------------------------------------------------------------------
if(FALSE){
  MI <- celldex::MonacoImmuneData(ensembl = FALSE, cell.ont = "all") 
  testdata <- GetAssayData(amp_2q.Tcell, layer="data")
  clusters <- amp_2q.Tcell@meta.data$seurat_clusters
  
  cellpred <- SingleR(test = testdata, 
                      ref = MI, 
                      labels = MI$label.fine, 
                      method = "cluster", 
                      clusters = clusters, 
                      assay.type.ref = "logcounts", 
                      assay.type.test = "logcounts")
  
  celltype.MI <- data.frame(ClusterID=rownames(cellpred), 
                            celltype=cellpred$labels, 
                            stringAsFactors=FALSE)
  
  Idents(object = amp_2q.Tcell) <- amp_2q.Tcell@meta.data$seurat_clusters
  new.cluster.ids <- celltype.MI$celltype 
  names(new.cluster.ids) <- levels(amp_2q.Tcell)
  
  amp_2q.Tcell <- RenameIdents(amp_2q.Tcell, new.cluster.ids)
  DimPlot(amp_2q.Tcell, reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend()
}
## Step 10 Manual Cluster deductions --------------------------------------------------------------------
## T cell cell types: 
## Naïve T, central memory T, effector memory T, GZMK CD4 T, IFN-act T, Cytotoxic T, GZMK CD8 T, Exhausted T

#"GZMK", "LEF1", "IL7R", "TRBC2", "CD3D", "CD3G", "CD3E" -------- haonan given markers 
#"CD3D", "CD79A", "CD14", "FCER1A", "FCGR3A", "IL3RA", "NKG7", "PPBP", "MZB1" ------- canonical markers 
#"IGHD", "CD27", "IGHM", "IGHG1", "IGHG3", "IGHA1", "MZB1" -------- B cell expression 
#"CCR7", "SELL" ---------- naive T cells 
#"PCNA", "MCM5", "MKI67", "STMN1", "CDK1"----------- mitotic TEM clusters 

## generage dot plots 
markers.filt <- markers %>% dplyr::filter(pct.1>0.60) 

top_markers <- markers.filt %>% group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) %>%
  pull(gene)
top_markers <- unique(top_markers)

p <- DotPlot(Tcell, features = top_markers) +
  scale_colour_gradient2(low = "#1D2B53", mid = "#fcfcfc", high = "#FF004D") +
  RotatedAxis() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Top 5 Genes per Cluster by Average log2(FC)") + 
  ylab("Seurat Cluster ID") + 
  ggtitle("Gene Expression in Ovarian Cancer T-cell Subcluster with Wild Type 2q")
ggsave("~/MRes_project_1/docs/SC_SEQ/ova/raw/GSE180661/plots/Tcell_wt2q_dotplot_long.png", 
       p, width=15, height=5, units = "in")



#amplification <- c("amp_2q.Bcell", "amp_2q.Myeloid", "wt_2q.Bcell", "wt_2q.Myeloid")
#"/rds/general/user/ph323/ephemeral/GSE180661/",i,".rds"


















