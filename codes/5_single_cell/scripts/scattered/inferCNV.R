## This script is part of the series of script that performs sample level single cell sequencing 
## The data is downloaded online and is a collection of different samples 
## 5 ovarian cancer patents (45,114 cells), 3' prime scRNA-seq
## https://lambrechtslab.sites.vib.be/en/pan-cancer-blueprint-tumour-microenvironment-0
## Sample and Data Relationship Format (SDRF)
## Pipeline: 
## Step 1: 
rm(list = ls())

## Section 0: Environment  ============================================================
## Packages 
library(infercnv)
library(dplyr)
library(Matrix)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(SeuratData)
library(ggplot2)
library(magrittr)
library(patchwork)
library(Signac)
library(SingleR)
library(celldex)
library(stringr)

## Controls 
count.matrix <- "~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/count_matrix/"
meta.data <- "~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/metadata.csv.gz"
data.dir <- "~/MRes_project_1/docs/SC_SEQ/lung/raw/lambrechtslab/"
inferCNV.dir <- "~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/InferCNV"
gene.order.file <- "~/MRes_project_1/codes/6_single_cell/reference/gene_order_file.txt"
output.dir <- "~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/InferCNV/result/"
rds.dir <- "~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/ova.rds"

## Section 1: Download Data ========================================================================
"https://www.ebi.ac.uk/biostudies/files/E-MTAB-8107/BT1303.counts.csv"
"https://www.ebi.ac.uk/biostudies/files/E-MTAB-8107/BT1304.counts.csv"
"https://www.ebi.ac.uk/biostudies/files/E-MTAB-8107/BT1305.counts.csv"
"https://www.ebi.ac.uk/biostudies/files/E-MTAB-8107/BT1306.counts.csv"
"https://www.ebi.ac.uk/biostudies/files/E-MTAB-8107/BT1307.counts.csv"
"https://www.ebi.ac.uk/biostudies/files/E-MTAB-8107/scrSOL001.counts.csv"
"https://www.ebi.ac.uk/biostudies/files/E-MTAB-8107/scrSOL003.counts.csv"
"https://www.ebi.ac.uk/biostudies/files/E-MTAB-8107/scrSOL004.counts.csv"
"https://www.ebi.ac.uk/biostudies/files/E-MTAB-8107/scrSOL006.counts.csv"
"https://www.ebi.ac.uk/biostudies/files/E-MTAB-8107/scrSOL007.counts.csv"


## Section 2: Perpare object for inferCNV ==========================================================
## load the general dataset and intialize the seurat object 
#BSC <- Read10X(data.dir = count.matrix)
#BSC <- CreateSeuratObject(counts = BSC, min.cells = 3, min.features = 200)
BSC <- readRDS(rds.dir)

## Read in metadata 
metadata <- read.csv(meta.data)
rownames(metadata) <- metadata$Cell
obj_names <- c("nGene", "nUMI", "CellFromTumor", "PatientNumber", "TumorType", "TumorSite", "CellType")
metadata <- metadata[obj_names] ## select out interested groups 

## read in sdrf
sdrf <- read.delim("~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/E-MTAB-8107.sdrf.txt")

## append 
metadata <- metadata[rownames(BSC@meta.data), ] ## reorder by BSC 
for(i in 1:7){
  named_vector <- metadata[[i]]
  names(named_vector) <- rownames(BSC@meta.data)
  BSC <- AddMetaData(object=BSC, 
                             metadata=named_vector, 
                             col.name=obj_names[[i]])
}



## Split Seurat objects based on Samples 
obj.list <- SplitObject(BSC, split.by = "orig.ident")


## perform InferCNV ==============================================================================
## example with scrSOL007
gene_order <- read.delim("~/MRes_project_1/codes/6_single_cell/reference/gene_order_file.txt", header = FALSE)

sub_matrix <- obj.list$scrSOL004@assays$RNA$counts %>% as.matrix()
sub_metadata <- obj.list$scrSOL004@meta.data
sub_metadata$barcode <- row.names(sub_metadata)
sub_metadata <- sub_metadata[, c("barcode", "CellType")]

rownames(sub_metadata) <- NULL
write.table(sub_metadata, paste0(inferCNV.dir, "cellAnnotations.txt"), 
            sep = '\t', quote = F, row.names = F, col.names = F)

infercnv_obj <-  CreateInfercnvObject(raw_counts_matrix=sub_matrix, 
                                    annotations_file=paste0(inferCNV.dir, "cellAnnotations.txt"),
                                    delim="\t",
                                    gene_order_file=gene.order.file,
                                    ref_group_names=c("T_cell", "Myeloid", "B_cell", "Fibroblast", "EC")) ## set to various normal cell types 

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir=paste0(output.dir, "scrSOL004/"),  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T, 
                             num_threads = 10)

#infercnv_obj_run <-  infercnv::run(infercnv_obj, 
                                   #cutoff=0.1, ## becase we uses 10x Genomics 
                                   #out_dir=paste0(output.dir, "scrSOL004/"), ## output directory 
                                   #cluster_by_groups=TRUE, ## basically set in all examples 
                                   #denoise=TRUE, ## same as above 
                                   #plot_steps = FALSE, 
                                   #sd_amplifier = 3, ## sigmoidal function 
                                   #noise_logistic = TRUE, ## turns gradient filtering on 
                                   #HMM=TRUE, ## run HMM based CNV prediction 
                                   #BayesMaxPNormal=0.4, ## posterior probabilies for cnv states as per bayesian network 
                                   #analysis_mode = "subclusters", ## partition cells into groups having consistent patterns of CNV
                                   #tumor_subcluster_partition_method = "random_trees", ## if a tree height is found statistically sig according to tumor_subclister_pval = 0.05, the tree is bifurcated 
                                   #tumor_subcluster_pval = 0.05, ## setting pvalue for random_trees, 
                                   #num_threads = 8, ## to speed up the process because random_trees is very slow. 
#)


#scrSOL004 <-  add_to_seurat(seurat_obj = obj.list$scrSOL004,
                            #infercnv_output_path = paste0(output.dir, "scrSOL004/"))

#temp1 <- read.table("~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/InferCNV/result/scrSOL004/HMM_CNV_predictions.HMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.4.pred_cnv_genes.dat", header = TRUE)
#temp1$cell_group_name %>% unique

#temp1$cell_name <- gsub("\\..*","", temp1$cell_group_name)
#summarise_temp1 <- temp1 %>% group_by(cell_name, chr) %>% summarise(temp1_mean = mean(state),
                                                                    #temp1_IQR = IQR(state), 
                                                                    #temp1_max = max(state), 
                                                                    #temp1_min = min(state))



temp1 <- read.table("~/MRes_project_1/codes/6_single_cell/reference/gene_order_file.txt")
















