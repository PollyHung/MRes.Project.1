## We have 7 datasets, for each dataset we shall extract the estimate, p-value, and FDR 
library(openxlsx)


## HH_ova===========================================================================================================
setwd("~/MRes_project_1/GISTIC2/output/tcga_lung/result_imm_lm_b_lung.rds")
read_list <- list.files(pattern = "^broad_m")

chromosomes <- c()
for(i in 1:22) {
  chromosomes <- c(chromosomes, paste(i, "p", sep=""), paste(i, "q", sep=""))}
chromosomes <- c(chromosomes, "Xp", "Xq")

p_value <- chromosomes %>% as.data.frame()
colnames(p_value) <- "chr_arm"

estimate <- chromosomes %>% as.data.frame()
colnames(estimate) <- "chr_arm"

fdr <- chromosomes %>% as.data.frame()
colnames(fdr) <- "chr_arm"

for(i in read_list){
  setwd("~/MRes_project_1/GISTIC2/output/cbioportal_mel/")
  temp <- read.table(i)
  temp$chr_arm <- rownames(temp)
  name <- gsub("broad_m_|_1.txt", "", i)
  
  ## p-values 
  p_values <- temp[c("chr_arm", "p_values")] %>% as.data.frame()
  colnames(p_values) <- c("chr_arm", name)
  p_value <- merge(p_value, p_values, by="chr_arm", all=TRUE)
  
  ## estimate 
  Estimates <- temp[c(ncol(temp), 1)] %>% as.data.frame()
  colnames(Estimates) <- c("chr_arm", name)
  estimate <- merge(estimate, Estimates, by="chr_arm", all=TRUE)
  
  ## FDR
  FDR <- temp[c("chr_arm", "FDR")] %>% as.data.frame()
  colnames(FDR) <- c("chr_arm", name)
  fdr <- merge(fdr, FDR, by="chr_arm", all=TRUE)
}

## TCGA系列===================================================================================================================
tcga <- readRDS("~/MRes_project_1/GISTIC2/output/tcga_ova/result_imm_lm_b_ova.rds")
tcga_os <- read.table("~/MRes_project_1/GISTIC2/output/tcga_ova/ova_os.txt")
tcga_pfs <- read.table("~/MRes_project_1/GISTIC2/output/tcga_ova/ova_pfi.txt")
tcga_age <- read.table("~/MRes_project_1/GISTIC2/output/tcga_ova/ova_age.txt")
tcga_sex <- read.table("~/MRes_project_1/GISTIC2/output/tcga_ova/ova_sex.txt")
tcga_stage <- read.table("~/MRes_project_1/GISTIC2/output/tcga_ova/ova_stage.txt")
tcga_therapy_response <- read.table("~/MRes_project_1/GISTIC2/output/tcga_ova/ova_therapy_resp.txt")

## reorder dataframe 
for(i in names(tcga)){
  temp <- tcga[[i]] 
  temp <- temp[chromosomes, ]
  rownames(temp) <- chromosomes
  tcga[[i]] <- temp
}

# Initialize an empty list to store the estimate vectors
estimate_list <- list()
p_value_list <- list()
fdr_list <- list()

# Loop through the list 'tcga'
for (name in names(tcga)) {
  # Extract the estimate component from each sublist
  estimate_list[[name]] <- tcga[[name]]$estimate
  p_value_list[[name]] <- tcga[[name]]$p_values
  fdr_list[[name]] <- tcga[[name]]$FDR
}

# Combine the estimates into a single dataframe
# Each sublist's estimates become a column in the new dataframe
estimate <- do.call(cbind, estimate_list) %>% as.data.frame()
p_value <- do.call(cbind, p_value_list) %>% as.data.frame()
fdr <- do.call(cbind, fdr_list) %>% as.data.frame()

estimate$chr_arm <- chromosomes
p_value$chr_arm <- chromosomes
fdr$chr_arm <- chromosomes

## clinical features 
list_files <- ls(pattern = "^tcga_")
for(i in list_files){
  files <- get(i, envir = .GlobalEnv)
  files$chr_arm <- rownames(files)
  name <- gsub("tcga_", "", i)
  
  ## p-values 
  p_values <- files[c("chr_arm", "p_values")] %>% as.data.frame()
  colnames(p_values) <- c("chr_arm", name)
  p_value <- merge(p_value, p_values, by="chr_arm", all=TRUE)
  
  ## estimate 
  Estimates <- files[c(ncol(files), 1)] %>% as.data.frame()
  colnames(Estimates) <- c("chr_arm", name)
  estimate <- merge(estimate, Estimates, by="chr_arm", all=TRUE)
  
  ## FDR
  FDR <- files[c("chr_arm", "FDR")] %>% as.data.frame()
  colnames(FDR) <- c("chr_arm", name)
  fdr <- merge(fdr, FDR, by="chr_arm", all=TRUE)
}

estimate <- estimate[chromosomes, ]


## save ===================================================================================================================
# Create a new workbook
wb <- createWorkbook()

addWorksheet(wb, "estimate")
addWorksheet(wb, "p_value")
addWorksheet(wb, "fdr")

# Write the dataframes to their respective sheets
writeData(wb, sheet = "estimate", estimate)
writeData(wb, sheet = "p_value", p_value)
writeData(wb, sheet = "fdr", fdr)

# Save the workbook to a file
saveWorkbook(wb, "~/MRes_project_1/0_RESULTS/xlsx/tcga_ova.xlsx", overwrite = TRUE)








