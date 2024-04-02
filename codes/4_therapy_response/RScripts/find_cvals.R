## This is a script from a series of scripts to test all kinds of possibilities for facets package 
## This script focuses on finding the optimal cval 
library(dplyr)
library(magrittr)
library(httr)
library(Rsamtools)
library(facets)
library(ASCAT)
library(data.tree)
library(ABSOLUTE)
library(caroline)
library(DoAbsolute)
library(facetsSuite)
library(argparse)
library(ggplot2)
library(egg)
library(purrr)
library(tibble)
library(igvR)
library(readr)
library(readxl)
library(stringr)
library(survival)
library(nnet)



## change to current working directories 
setwd("/rds/general/user/ph323/home/MRes_project_1/Codes/3_therapy_response/example")

## read in sample ids 
sample_ids <- read.table("sample_ids.txt", header = FALSE) %>% unlist
samples <- sample(sample_ids, 10, replace = FALSE)

## build a dataframe to contain 
cvals <- c(50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100)
seg_length <- matrix(NA, nrow = length(cvals), ncol = 10) %>% as.data.frame
colnames(seg_length) <- samples
rownames(seg_length) <- cvals

## find optimal cvals 
for(i in 1:length(samples)){
  cancer <- samples[i]
  
  ## read in the snp-piled up file 
  for(j in 1:length(cvals)){
    CVAL <- cvals[j]
    
    datafile <- paste0("/rds/general/user/ph323/home/MRes_project_1/docs/HH_ova/facets/facet_input/", cancer, "_ordered.csv")
    
    ## preprocessing 
    set.seed(1234)
    rcmat <-  readSnpMatrix(datafile)
    xx <- preProcSample(rcmat, gbuild = "hg38")
    oo <- procSample(xx, cval = CVAL, dipLogR = xx$dipLogR)
    
    ## fitting 
    fit <- emcncf(oo)
    seg_length[j, i] <- nrow(fit$cncf)
  }
    print("finish with ", cancer)
} 

write.table(seg_length, "find_cval.txt", quote = FALSE, col.names = TRUE, row.names = TRUE, 
            sep = "\t")












