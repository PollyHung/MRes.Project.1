library(data.table)
library(ABSOLUTE)
library(DoAbsolute)
# Load Test Data ----------------------------------------------------------

# segmentation file
Seg <- file.path("~/MRes_project_1/docs/HH_ova/facets/segmentation_file/_hisens_diplogR.adjusted.seg")
Seg <- fread(Seg)
colnames(Seg) <- c("Sample" ,"Chromosome","Start","End","Num_Probes","Segment_Mean")

# test function
DoAbsolute(Seg = Seg, 
           sigma.p = 0, 
           max.sigma.h = 0.015, 
           min.ploidy = 0.95, 
           max.ploidy = 10, 
           primary.disease = "ov", 
           platform = "Illumina_WES", 
           results.dir = "~/MRes_project_1/docs/HH_ova/absolute/do_absolute", 
           max.as.seg.count = 5000, 
           max.non.clonal = 0.5, 
           max.neg.genome = 0.9, 
           copy.num.type = "total", 
           keepAllResult = TRUE, verbose = TRUE)