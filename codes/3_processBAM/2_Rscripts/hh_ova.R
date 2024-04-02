## Section 0: Environment ========================================================================
setwd("~/MRes_project_1/docs/absolute")

## library things 
library(ABSOLUTE)
library(dplyr)


## Section 1: read in and set up files ===========================================================
hh_ova.dir <- "~/MRes_project_1/docs/HH_ova/facets/seg_cval_50/_hisens_diplogR.adjusted.seg"


## Section 2: runABSOLUTE =========================================================================
## Section 2.1: HH_ova ----------------------------------------------------------------------------
segmentation <- read.table(hh_ova.dir)
colnames(segmentation) <- c("ID", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")
segmentation_list <- split(segmentation, segmentation$ID)

for(i in names(segmentation_list)){ 
  setwd("~/MRes_project_1/docs/absolute/")
  x <-  segmentation_list[[i]]
  x <-  x[, c(2:6)]
  names(x) <- c("Chromosome","Start","End","Num_Probes","Segment_Mean")
  x <- x[x[, 1] %in% c(1:22), ]
  write.table(x, "x", sep = "\t", row.names = FALSE)
  
  RunAbsolute("x", sigma.p = 0, max.sigma.h = 0.015, min.ploidy = 0.95, 
              max.ploidy = 10, primary.disease = "ov", 
              platform = "Illumina_WES", sample.name = i, 
              results.dir = "~/MRes_project_1/docs/absolute/HH_ova", 
              max.as.seg.count = 1500, max.non.clonal = 0.05, max.neg.genome = 0.005, 
              copy_num_type = "total", maf.fn = NULL, min.mut.af = NULL, output.fn.base = NULL, 
              verbose = FALSE)
  
  setwd("~/MRes_project_1/docs/absolute/HH_ova")
  load(paste0(i, ".ABSOLUTE.RData"))
  write.csv(seg.dat$segtab, file = paste0(i, ".csv"))
  
  print(paste0("finish running ABSOLUTE algorithm on sample ", i))
}
