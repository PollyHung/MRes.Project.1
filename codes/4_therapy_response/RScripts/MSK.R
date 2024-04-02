library(DoAbsolute)

msk_seg <- read.table("~/MRes_project_1/docs/MSK/msk.seg")
colnames(msk_seg) <- c("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")

DoAbsolute(Seg = msk_seg, platform = "Illumina_WES", copy.num.type = "total",
           results.dir = "~/MRes_project_1/docs/MSK/absolute", 
           keepAllResult = TRUE, verbose = TRUE)








