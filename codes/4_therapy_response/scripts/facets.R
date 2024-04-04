## environment ============================================================================
## controls whether to run 
INSTALL <- FALSE
WRITE <- FALSE
UNZIP <- FALSE
EXECUTE <- FALSE

## paths 
DATA <- "/rds/general/user/ph323/ephemeral"
RESULT <- "~/MRes_project_1/docs/HH_ova/facet_input"
EXAMPLE <- "~/MRes_project_1/Codes/3_therapy_resp/example"
GATK <- "~/MRes_project_1/otherCodes/gatk-4.2.5.0/gatk" ## this is also a file exec
SNP_PILEUP <- "~/MRes_project_1/Codes/3_therapy_resp/pileup_exe/snp-pileup" ## this is a file exec
VCF_FILE <- "/MRes_project_1/Codes/3_therapy_resp/vcf_file"

## library packages 
library(ABSOLUTE)
library(argparse)
library(ASCAT)
library(caroline)
library(data.table)
library(data.tree)
library(dplyr)
library(DoAbsolute)
library(egg)
library(facets)
library(facetsSuite)
library(ggplot2)
library(httr)
library(igvR)
library(magrittr)
library(nnet)
library(purrr)
library(readr)
library(readxl)
library(Rsamtools)
library(stringr)
library(survival)
library(tibble)


setwd("~/MRes_project_1/Codes/3_therapy_response/HH_lung/facets")
rm(list = ls())
load("facets.RData")

## List files ============================================================================
## list the bam files 
file_list <- list.files(paste0(DATA, "/HH_lung/Alignment_lung"), pattern = "_sorted_nodup.bam$")
sample_ids <- gsub("_sorted_nodup.bam", "", file_list)
sample_ids <- sample_ids[grepl("^MDL", sample_ids)]

## tumor and normal sample lists 
tumour <- file_list[grepl("^MDL", file_list)]  ## 52 lung cancer samples in total 
normal <- rep("X991N_sorted_nodup.bam", length(tumour))

## write the output to the folder 
if(WRITE){
  write.table(tumour, "tumour_order.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(normal, "normal_order.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(sample_ids, "sample_ids.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
}


## Quality Control ============================================================================




## SNP-pileup ============================================================================
## list the chromosome order files comparing the lung and vcf files 
## Method: choose 2 files and make it to the chromosome order, compare with the vcf chr order file we currently have 
## Location: /rds/general/user/ph323/ephemeral/HH_lung/Alignment_lung 
## Code: 
## produced by samtools idxstats MDL18-2281_sorted_nodup.bam | cut -f 1 | uniq > MDL18-2281.txt
## produced by samtools idxstats MDL16-4473_sorted_nodup.bam | cut -f 1 | uniq > MDL16-4473.txt
chr_order <- read.table("/rds/general/user/ph323/ephemeral/HH_lung/Alignment_lung/MDL18-2281.txt")
chr_order$sample2 <- read.table("/rds/general/user/ph323/ephemeral/HH_lung/Alignment_lung/MDL16-4473.txt") %>% unlist 
chr_order$X991 <- read.table("X991N_order.txt") %>% unlist
identical(chr_order$V1, chr_order$sample2)
identical(chr_order$X991, chr_order$sample2) 
## Result: confirmed that samples have the same chromosome order 
vcf_order <- read.table("vcf_new_order.txt", header = TRUE)
chr_order[1:25, ]
## about the same, except for an extra M between 9 and X 


## Execute SNP-pileup 
if(EXECUTE){
  system("chmod +x snp_pileup.sh") ## make the snp_pileup.sh executable 
  system("./snp_pileup.sh") ## apply 
}

## reorder the files 
for(cancer in sample_ids){
  ## output file: 
  setwd("~/MRes_project_1/docs/HH_lung/facets/facet_input/")
  output <- paste0(cancer, "_ordered.csv.gz")
  
  if (!file.exists(output)) {
    ## reorder the snp_pileup files to 1, 2, 3... 22, X, Y format 
    temp <- read.csv(paste0(cancer, ".csv")) ## read in files 
    chrom_order <- c(as.character(1:22), "X", "Y")  ## set up the order 
    temp <- temp[order(match(temp$Chromosome, chrom_order), temp$Position), ] ## reorder 
    
    ## save the ordered csv 
    write.csv(temp, paste0(cancer, "_ordered.csv"), 
              quote = FALSE, col.names = TRUE, row.names = FALSE)  ## save the files 
    
    ## write the file 
    write.csv(temp, gzfile(paste0(cancer, "_ordered.csv.gz")), 
              row.names = FALSE)
    
    ## echo
    print(paste0(cancer, " has been reordered and gunzipped"))
  } else {
    ## if it has already been reorganized and gunzipped, echo skip 
    print(paste0(cancer, " has already been reordered and gunzipped, no action performed"))
  }
}
list.files(".", pattern = "ordered.csv$") %>% length()

## examine the reordered files.    
facet_inputs_summary <- matrix(NA, nrow = 1, ncol = 13) %>% as.data.frame()
colnames(facet_inputs_summary) <- c("sample", "ncol", "nrow", 
                                    "Ref_null", "Alt_null", 
                                    "File1R", "File1A", "File1E", "File1D", 
                                    "File2R", "File2A", "File2E", "File2D")

# Function to summarize each file
summarize_facet_input <- function(sample_id) {
  file_path <- paste0(sample_id, "_ordered.csv")
  temp <- read.csv(file_path)
  
  # Create summary for this sample
  tibble(
    sample = sample_id,
    ncol = ncol(temp),
    nrow = nrow(temp),
    Ref_null = sum(temp$Ref == "."),
    Alt_null = sum(temp$Alt == "."),
    File1R = median(temp$File1R, na.rm = TRUE),
    File1A = median(temp$File1A, na.rm = TRUE),
    File1E = median(temp$File1E, na.rm = TRUE),
    File1D = median(temp$File1D, na.rm = TRUE),
    File2R = median(temp$File2R, na.rm = TRUE),
    File2A = median(temp$File2A, na.rm = TRUE),
    File2E = median(temp$File2E, na.rm = TRUE),
    File2D = median(temp$File2D, na.rm = TRUE)
  )
}

# Apply the function to each sample_id and combine the results
facet_inputs_summary <- lapply(sample_ids, summarize_facet_input) %>% bind_rows()


## FACETS ============================================================================
## Execute facets 
if(EXECUTE){
  setwd("~/MRes_project_1/Codes/3_therapy_response/HH_lung/facets")
  system("chmod +x facets.sh") ## make the snp_pileup.sh executable 
  system("./facets.sh") ## apply 
}

## quality control of facet runs 
## set working directory 
if(EXECUTE){
  setwd("~/MRes_project_1/docs/HH_lung/facets/facet_output")
  quality_control <- read.table(paste0(sample_ids[1], "/", sample_ids[1], ".qc.txt"), header = TRUE)
  summary_table <- read.table(paste0(sample_ids[1], "/", sample_ids[1], ".txt"), header = TRUE, sep = "\t")     
  arm_level <- read.table(paste0(sample_ids[1], "/", sample_ids[1], ".arm_level.txt"), header = TRUE, sep = "\t") 
  gene_level <- read.table(paste0(sample_ids[1], "/", sample_ids[1], ".gene_level.txt"), header = TRUE, sep = "\t") 
  ## loop thorugh all the files 
  for(i in 2:length(sample_ids)) {
    ## make file path 
    file_path_qc <- paste0(sample_ids[i], "/", sample_ids[i], ".qc.txt")
    file_path_s <- paste0(sample_ids[i], "/", sample_ids[i], ".txt")
    file_path_a <- paste0(sample_ids[i], "/", sample_ids[i], ".arm_level.txt")
    file_path_g <- paste0(sample_ids[i], "/", sample_ids[i], ".gene_level.txt")
    
    if(file.exists(file_path_qc)) {
      ## quality_control 
      temp_qc <- read.table(file_path_qc, header = TRUE)
      quality_control <- rbind(quality_control, temp_qc)
    }
    if(file.exists(file_path_s)) {
      ## summary_table
      temp_s <- read.table(file_path_s, header = TRUE, sep = "\t")
      summary_table <- rbind(summary_table, temp_s)
    }
    if(file.exists(file_path_a)) {
      ## arm level 
      temp_a <- read.table(file_path_a, header = TRUE, sep = "\t")
      arm_level <- rbind(arm_level, temp_a)
    }
    if(file.exists(file_path_g)) {
      ## gene level 
      temp_g <- read.table(file_path_g, header = TRUE, sep = "\t")
      gene_level <- rbind(gene_level, temp_g)
    }
  }
  summary_table$flags <- NULL
}
## write output 
if(WRITE){
  setwd("~/MRes_project_1/docs/HH_ova/facets/quality_control")
  write.csv(quality_control, "quality_control.csv", row.names = FALSE, quote = FALSE)
  write.table(summary_table, "summary_table.txt", row.names = FALSE, quote = FALSE)
  write.csv(arm_level, "arm_level.csv", row.names = FALSE, quote = FALSE)
  write.csv(gene_level, "gene_level.csv", row.names = FALSE, quote = FALSE)
}


## Segmentation file 
facet_process <- function(sample_ids, seg_file_type, facet_dir, seg_file_dir){
  ## sample_id from sample_ids 
  ## seg_file_type: _hisens_diplogR.adjusted.seg, _hisens_diplogR.unadjusted.seg, 
  ##                _purity_diplogR.adjusted.seg, _purity_diplogR.unadjusted.seg
  ## facet_directory: facet_output
  ## seg_file_dir: seg_cval_50, seg_cval_25
  
  ## build an empty dataframe 
  temp <- matrix(NA, nrow = 1, ncol = 6) %>% as.data.frame
  colnames(temp) <- c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")
  
  
  ## loop through each directory and get the dataframe 
  for(sample_id in sample_ids){
    ## set directory and switch to directory 
    directory <- file.path(paste0("~/MRes_project_1/docs/HH_lung/facets/", facet_dir,"/", sample_id))
    setwd(directory)
    
    ## import hisens_adj segmented file 
    segmented_file <- read.table(paste0(sample_id, seg_file_type), header = TRUE)
    temp <- rbind(temp, segmented_file)
  }
  ## remove the first NA row 
  temp <- temp[2:nrow(temp), ]
  
  ## print the result 
  print(gsub("_|\\.", " ", seg_file_type))
  print(summary(temp))
  
  ## write the file 
  output_dir <- file.path(paste0("~/MRes_project_1/docs/HH_lung/facets/", seg_file_dir))
  setwd(output_dir)
  write.table(temp, seg_file_type, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  ## return the file 
  return(temp)
}

if(WRITE){
  hisens_adj <- facet_process(sample_ids, "_hisens_diplogR.adjusted.seg", "facet_output", "seg_cval_50")
  hisens_unadj <- facet_process(sample_ids, "_hisens_diplogR.unadjusted.seg", "facet_output", "seg_cval_50")
  purity_adj <- facet_process(sample_ids, "_purity_diplogR.adjusted.seg", "facet_output", "seg_cval_50")
  purity_unadj <- facet_process(sample_ids, "_purity_diplogR.unadjusted.seg", "facet_output", "seg_cval_50")
}


## DoABSOLUTE ============================================================================ 
## find ~/MRes_project_1/docs/HH_lung/facets/facet_output -name "*_hisens_diplogR.adjusted.seg" -exec cp {} ~/MRes_project_1/docs/HH_lung/absolute/absolute_input/ \;
## the above command copies the set of segmentation file into a new folder for following process 

# segmentation file
if(EXECUTE){
  setwd("~/MRes_project_1/docs/HH_lung/absolute/absolute_input/")
  files = list.files(".", pattern = "seg", full.names = TRUE)
  files
  Seg = rbindlist(lapply(files, fread))
  
  Seg
  Seg2 = Seg[, list(Sample = ID, Chromosome = chrom, Start = loc.start, End = loc.end, 
                    Num_Probes = num.mark, Segment_Mean = seg.mean)]
  
  # test function
  DoAbsolute(Seg = Seg2, platform = "Illumina_WES", copy.num.type = "total",
             results.dir = "~/MRes_project_1/docs/HH_lung/absolute/absolute_output/", 
             keepAllResult = TRUE, verbose = TRUE)
  ## this is not working 
}

## GISTIC2 ============================================================================ 
setwd("~/MRes_project_1/docs/HH_lung/gistic/gistic_output/")
if(UNZIP){
  system("unzip hisens_adj.zip") 
  system("mv 558627 hisens_adj") 
  
  system("unzip hisens_unadj.zip")
  system("mv 558628 hisens_unadj") 
  
  system("unzip purity_adj.zip")
  system("mv 558629 purity_adj") 
  
  system("unzip purity_unadj.zip")
  system("mv 558630 purity_unadj") 
}

## which segmentation resulted gistic output should we use? 
distribution_plot <- function(file_directory, name){
  ## file_directory of the file 
  ## name = focal_tcga, broad_tcga 
  
  temp <- read.table(file_directory, sep = "\t", header = TRUE)
  
  pdf(paste0("~/MRes_project_1/docs/HH_lung/gistic/distribution_plots/", name, ".pdf"))
  temp[, 11:(ncol(temp)-1)] %>% colMeans %>% hist(main = paste0(name, "_colMean"))
  temp[, 11:(ncol(temp)-1)] %>% rowMeans %>% hist(main = paste0(name, "_rowMean"))
  temp[, 11:(ncol(temp)-1)] %>% as.matrix %>% colVars(useNames = FALSE) %>% hist(main = paste0(name, "_colVar"))
  temp[, 11:(ncol(temp)-1)] %>% as.matrix %>% rowVars(useNames = FALSE) %>% hist(main = paste0(name, "_rowVar"))
  dev.off()
}

if(WRITE){
  distribution_plot("~/MRes_project_1/GISTIC2_by_genePattern/tcga/lung/all_lesions.conf_95.txt", "focal_tcga")
  distribution_plot("~/MRes_project_1/GISTIC2_by_genePattern/tcga/lung/broad_values_by_arm.txt", "broad_tcga")
  distribution_plot("hisens_adj/all_lesions.conf_95.txt", "focal_his_a")
  distribution_plot("hisens_adj/broad_values_by_arm.txt", "broad_his_a")
  distribution_plot("purity_adj/broad_values_by_arm.txt", "broad_pur_a")
}
## After examing the distribution we decided that hisens_adj is the dataset we decide to proceed with.    


## read in clinical data 
clinical_doc <- read.csv("~/MRes_project_1/docs/HH_lung/clinical_data/HH_lung_clinical.csv")
#radiomics <- read.csv("~/MRes_project_1/docs/HH_lung/clinical_data/Clinical_data_HH_Immunoradiomics_full_2.csv")
#tumour_cellularity <- read_xlsx("~/MRes_project_1/docs/HH_lung/clinical_data/HH_tumourcellularity.xlsx")
colnames(clinical_doc)

## read in focal and broad data 
focal <- read.table("~/MRes_project_1/docs/HH_lung/gistic/gistic_output/hisens_adj/all_lesions.conf_95.txt",
                    sep = "\t", header = TRUE)
focal <- focal[grepl("CN values$", focal$Unique.Name), ]
focal$Unique.Name <- gsub("lification Peak|- CN values|etion Peak|\\s+", "", focal$Unique.Name) ## clean name 
focal$Unique.Name <- paste0(focal$Unique.Name, "_", focal$Descriptor) %>% tolower ## create unique identifier 
rownames(focal) <- focal$Unique.Name ## rename 

broad <- read.table("~/MRes_project_1/docs/HH_lung/gistic/gistic_output/hisens_adj/broad_values_by_arm.txt", 
                    sep = "\t", header = TRUE)
rownames(broad) <- broad$Chromosome.Arm
## focal analysis     
clinical_cleaned <- clinical_doc[c("Sample.Name", "sex", "AgeDx", "PD.L1.expression", "M_sites", "OVERALL.SURVIVAL", 
                                   "Best.Response", "RESPONSE.RATE", "PROGRESSION.RATE", "Progression", 
                                   "PROGRESSION.FREE.SURVIVAL", "OS_days", "PFS_days", "OS_event")]
clinical_cleaned$Sample.Name <- gsub("-", ".", clinical_cleaned$Sample.Name)
sample_ids <- intersect(clinical_cleaned$Sample.Name, colnames(focal))
rownames(clinical_cleaned) <- clinical_cleaned$Sample.Name
clinical_cleaned <- clinical_cleaned[sample_ids, ]
immune <- clinical_cleaned[sample_ids, "PD.L1.expression"] %>% as.data.frame()
broad <- broad[, sample_ids]
focal <- focal[, sample_ids]

gistic <- focal %>% as.matrix()

for(i in 1:ncol(immune)){
  ## build an empty dataframe 
  output <- array(NA, c(nrow(gistic), 5))
  rownames(output) <- rownames(gistic)
  colnames(output) <- c("estimate", "std_error", "t_value", "p_values", "r_sq")
  output <- output %>% as.data.frame
  
  ## multivariate analysis accounting for age and stage 
  for(j in 1:nrow(gistic)){
    coxphmodel <- lm(immune[, i]~gistic[j, ]+clinical_cleaned$AgeDx+clinical_cleaned$sex, data = as.data.frame(gistic))
    temp <- summary(coxphmodel)
    output$estimate[j] <- temp$coefficients[2, 1]
    output$std_error[j] <- temp$coefficients[2, 2]
    output$t_value[j] <- temp$coefficients[2, 3]
    output$p_values[j] <- temp$coefficients[2, 4]
    output$r_sq[j] <- temp$adj.r.squared
  }
  
  output <- as.data.frame(output)
  output$FDR <- p.adjust(output$p_values, method = "fdr")
  output <- output[order(output$FDR, output$p_values, decreasing = FALSE), ] 
  
  os_survival_name <- make.names("focal_PDL1")
  assign(os_survival_name, output, envir = .GlobalEnv)
}









