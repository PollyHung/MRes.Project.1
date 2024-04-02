## This is a script that converts facets result to absolute result 
## The Facets result is from original FACETS package 

## I wouldn't use Method 2, if you want to adjust for the diploid log ratio, use cnlr.median - diplogR.
## A third method is to use the log tcn values, which corresponds to also adjusting for the purity.

## Section 1: environment ======================================================
## Control ---------------------------------------------------------------------
OVA <- "~/MRes_project_1/docs/HH_ova/facets_original/dataframe/"
LUNG <- "~/MRes_project_1/docs/HH_lung/facets_original/dataframe/"
ABS_OVA <- "~/MRes_project_1/docs/ABSOLUTE/HH_ova/"
ABS_LUNG <- "~/MRes_project_1/docs/ABSOLUTE/HH_lung/"
COLUMNS <- c("chrom", "start", "end", "num.mark", "cnlr.median")

## Library ---------------------------------------------------------------------
library(ABSOLUTE)
library(dplyr)
library(magrittr)
library(httr)
library(facets)


## Section 2: ABSOLUTE =========================================================
## list files    
facets_output <- list.files(OVA)
facets_output <- gsub(".csv", "", facets_output)
facets_output <- c("X103", "X1032", "X1040", "X2362")

for(i in facets_output) {
  ## set working directory -----------------------------------------------------   
  setwd(ABS_OVA)
  output <- paste0(ABS_OVA, i,".ABSOLUTE.RData") ## this is where your output will be stored 
  x <- read.csv(paste0(OVA, i, ".csv"))  ## read in table 
  x <- x[, COLUMNS] ## edit the table by selecting desirable columns 
  names(x) <- c("Chromosome","Start","End","Num_Probes","Segment_Mean") ## rename dataframe 
  x <- x[x[, 1] %in% c(1:22),] ## only select chromosome 1 to 22 
  
  ## write the table 
  write.table(x, "x", sep = "\t", row.names = FALSE)
  
  ## run absolute --------------------------------------------------------------
  RunAbsolute("x", 
              sigma.p=0, 
              max.sigma.h=0.015, min.ploidy=0.95, max.ploidy=10, 
              primary.disease="ov", platform="Illumina_WES", 
              sample.name= i, 
              results.dir=ABS_OVA, 
              max.as.seg.count=1500, max.non.clonal=0.05, max.neg.genome=0.005, 
              copy_num_type="total", 
              maf.fn=NULL, min.mut.af=NULL, 
              output.fn.base=NULL, verbose=FALSE)
  
  ## set working directory 
  setwd(ABS_OVA)
  load(paste0(i,".ABSOLUTE.RData"))
  write.csv(seg.dat$segtab, file = paste0(i, ".csv"))

  ## Create the purity and ploidy scores for the absolute ----------------------
  CreateReviewObject(obj.name = i, 
                     absolute.files = paste0(i, ".ABSOLUTE.RData"), 
                     indv.results.dir = paste0(ABS_OVA, i), 
                     copy_num_type = "total", 
                     plot.modes = TRUE)
  
  print(paste("finished RunAbsolute for", i))
}