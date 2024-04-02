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

source("~/MRes_project_1/codes/6_single_cell/RScript/pipeline/cell_chat.R")
## example application
#> cellchat_wt <- do_cellchat(wildtype)
#[1] "Create a CellChat object from a Seurat object"
#The `meta.data` slot in the Seurat object is used as cell meta information 
#Set cell identities for the new CellChat object 
#The cell groups used for CellChat analysis are  Cancer, Endothelial, Fibroblast, T cell, Myeloid, B cell 
#The number of highly variable ligand-receptor pairs used for signaling inference is 707 
#triMean is used for calculating the average gene expression per cell group. 
#[1] ">>> Run CellChat on sc/snRNA-seq data <<< [2024-03-21 20:37:43]"
#|==========================================================================================================================| 100%
#[1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2024-03-21 20:48:42]"
#Warning message:
#  In createCellChat(object = seurat_object, group.by = "ident", assay = "RNA") :
#  The 'meta' data does not have a column named `samples`. We now add this column and all cells are assumed to belong to `sample1`!

## read in data 
amplified <- readRDS("~/MRes_project_1/codes/6_single_cell/RNotebook/subcluster_analysis/wholeSample/amplified.rds")
wildtype <- readRDS("~/MRes_project_1/codes/6_single_cell/RNotebook/subcluster_analysis/wholeSample/wildtype.rds")
cellchat_amp <- readRDS("~/MRes_project_1/codes/6_single_cell/RNotebook/subcluster_analysis/wholeSample/cellchat_amp.rds")
cellchat_wt <- readRDS("~/MRes_project_1/codes/6_single_cell/RNotebook/subcluster_analysis/wholeSample/cellchat_wt.rds")

groupSize_amp = as.numeric(table(cellchat_amp@idents))
groupSize_wt = as.numeric(table(cellchat_wt@idents))

## cell chat analysis 
netVisual_circle(cellchat_amp@net$count, vertex.weight = groupSize_amp, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_wt@net$count, vertex.weight = groupSize_wt, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")



netVisual_circle(cellchat_amp@net$weight, vertex.weight = groupSize_amp, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
netVisual_circle(cellchat_wt@net$weight, vertex.weight = groupSize_wt, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")


common_pathways <- intersect(cellchat_amp@netP$pathways, cellchat_wt@netP$pathways)

pdf("~/MRes_project_1/codes/5_plots/plots/plot_3/wholesample.pdf", width = 20, height = 10)
for(i in common_pathways){
  par(mfrow=c(1,1))
  netVisual_aggregate(cellchat_amp, signaling = i, layout = "chord", color.use = color) 
  netVisual_aggregate(cellchat_wt, signaling = i, layout = "chord", color.use = color) 
}
dev.off()

#> cellchat_amp@netP$pathways
#[1] "MHC-II"     "MHC-I"      "COLLAGEN"   "MK"         "PECAM1"     "LAMININ"    "VEGF"       "FN1"        "ANNEXIN"   
#[10] "CLDN"       "ADGRE"      "CCL"        "ADGRG"      "CDH5"       "ApoE"       "COMPLEMENT" "SELL"       "ESAM"      
#[19] "CLEC"       "BAFF"       "GRN"        "SPP1"       "NOTCH"      "PARs"       "NECTIN"     "HSPG"       "GAP"       
#[28] "GAS"        "CXCL"       "ApoA"       "THBS"       "LAIR1"      "PERIOSTIN"  "JAM"        "GALECTIN"   "KLK"       
#[37] "ICAM"       "AGRN"       "LCK"        "PTN"        "PROS"       "ANGPTL"     "CDH"        "CD86"       "PTPRM"     
#[46] "SEMA3"      "DESMOSOME"  "SN"         "ANGPT"      "PDGF"       "SEMA6"      "IGF"       
#> cellchat_wt@netP$pathways
#[1] "MHC-I"      "MHC-II"     "COLLAGEN"   "MK"         "PECAM1"     "ANNEXIN"    "FN1"        "ApoE"       "LAMININ"   
#[10] "VEGF"       "CLDN"       "CCL"        "COMPLEMENT" "ADGRE"      "CDH5"       "ESAM"       "CLEC"       "ADGRG"     
#[19] "GRN"        "PARs"       "NOTCH"      "SPP1"       "HSPG"       "GAP"        "GAS"        "LCK"        "NECTIN"    
#[28] "THBS"       "PERIOSTIN"  "PLAU"       "LAIR1"      "ApoA"       "JAM"        "PTN"        "CXCL"       "GALECTIN"  
#[37] "ICAM"       "ANGPTL"     "AGRN"       "PTPRM"      "PROS"       "CD86"       "SN"         "SEMA3"      "DESMOSOME" 
#[46] "ANGPT"      "KLK"        "SEMA6"      "CADM"       "CDH"        "PDGF"       "IGF"       
#> setdiff(cellchat_amp@netP$pathways, cellchat_wt@netP$pathways)
#[1] "SELL" "BAFF"
#> setdiff(cellchat_wt@netP$pathways, cellchat_amp@netP$pathways)
#[1] "PLAU" "CADM"

color = c("#F5786D", "#07BA35", "#B4A100", "#5E9EFE", "#00BFC3", "#F464E2")

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_chord_gene(cellchat_amp, #sources.use = c(1,4,5,6), targets.use = c(1,4,5,6), 
                     signaling = c("CCL", "CXCL", "SPP1", "COMPLEMENT"), title.name = "chr2q amp", color.use = color, 
                     legend.pos.x = 10)
netVisual_chord_gene(cellchat_wt, #sources.use = c(1,4,5,6), targets.use = c(1,4,5,6), 
                     signaling = c("CCL", "CXCL", "SPP1", "COMPLEMENT"), title.name = "chr2q wt", color.use = color, 
                     legend.pos.x = 10)



table(Idents(amplified)) 
table(Idents(wildtype))

table <- data.frame(table(Idents(amplified)), table(Idents(wildtype))) %>% 
  dplyr::select(Var1, Freq, Freq.1) %>% 
  dplyr::rename(cell_type = Var1, amp = Freq, wt = Freq.1)

sums <- colSums(table[, 2:3])

table_2 <- sweep(table[, 2:3], MARGIN = 2, STATS = sums, FUN = "/")
table_2 <- table_2*100
table_2$cell_type <- table$cell_type
table_3 <- table_2
#table_2[2, ] <- table_2[3, ]
#table_2[3, ] <- table_3[2, ]

long_data <- gather(table_2, condition, count, -cell_type)

p2 <- ggplot(long_data, aes(x = condition, y = count, fill = cell_type)) + 
  geom_bar(stat = "identity") + theme_bw() + #coord_flip() + 
  labs(x = "chr2q Status", y = "Percent of Cells") + 
  scale_fill_manual(values = color) + 
  theme(legend.position = "right", legend.title = element_blank()) + 
  #ggtitle("Cell Composition (pct)")
  #guides(fill = guide_legend(nrow = 2)) + 
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 6, margin(t = 0, unit = "pt")),
        axis.ticks.length = unit(1, "mm"), 
        title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(2.5, 'mm'), 
        legend.margin = margin(t = 6, unit = "pt"), 
        strip.text.x = element_text(size = 6, face = "bold"), 
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "mm"))

ggsave(plot = p2, filename = "~/MRes_project_1/codes/5_plots/plots/plot_3/stackbarPlot_wholesample.png", width = 1.6, height = 1.6, units = "in", dpi = 600)



## Plot 2 ======================================================================


color = c("#F5786D", "#07BA35", "#B4A100", "#5E9EFE", "#00BFC3", "#F464E2")

# #F5786D = red = cancer
# #07BA35 = green = endothelial 
# #B4A100 = yellow = fibroblast 
# #5E9EFE = dark blue = T cell 
# #00BFC3 = light blue = myeloid 
# #F464E2 = purple = b cell 














