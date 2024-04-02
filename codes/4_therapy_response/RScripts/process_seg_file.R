library(ABSOLUTE)

## read in segmentation file 
segmentation <- read.table("~/MRes_project_1/docs/cBioPortal_lung/msk.seg")
colnames(segmentation) <- c("Sample", "Chromosome","Start","End","Num_Probes","Segment_Mean")
segmentation$Sample <- gsub("-", ".", segmentation$Sample)

## make unique sample ids 
sample_ids <- unique(segmentation$Sample)

## split the dataframe 
list_df <- split(segmentation, as.factor(segmentation$Sample))
list2env(list_df, envir = .GlobalEnv)

## list files 
filelist <- ls(pattern="P.")

## loop through absolute 
for(i in filelist) {
  ## preprocess the file 
  setwd("~/MRes_project_1/ABSOLUTE/temp")
  seg <- get(i, envir = .GlobalEnv)
  seg <- seg[, c(2:6)]
  colnames(seg) <- c("Chromosome","Start","End","Num_Probes","Segment_Mean")
  seg <- seg[seg[, 1] %in% c(1:22), ]
  write.table(seg, "seg", sep = "\t", row.names = FALSE)

  ## run absolute 
  RunAbsolute("seg", sigma.p=0, max.sigma.h=0.015, min.ploidy=0.95, max.ploidy=10, 
              primary.disease="lung", platform="Illumina_WES", sample.name= i, 
              results.dir="~/MRes_project_1/docs/cBioPortal_lung/absolute/", 
              max.as.seg.count=1500, max.non.clonal=0.05, max.neg.genome=0.005, 
              copy_num_type="total", maf.fn=NULL, min.mut.af=NULL, 
              output.fn.base=NULL, verbose=FALSE)
  setwd("~/MRes_project_1/docs/cBioPortal_lung/absolute/")
  load(paste0(i, ".ABSOLUTE.RData"))
  write.csv(seg.dat$segtab, file=paste0(i,".csv"))

  CreateReviewObject( 
    obj.name = i, 
    absolute.files = paste0(i, ".ABSOLUTE.RData"), 
    indv.results.dir = "~/MRes_project_1/docs/cBioPortal_lung/absolute/purity", 
    copy_num_type = "total", 
    plot.modes = FALSE, verbose = FALSE)
  
  print(paste("Finish Doing ABSOLUTE for", i))
}



