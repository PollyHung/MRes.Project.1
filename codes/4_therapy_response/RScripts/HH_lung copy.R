library(DoAbsolute)

hh_lung_seg <- read.table("~/MRes_project_1/docs/HH_lung/facets/seg_cval_50/_hisens_diplogR.adjusted.seg")
colnames(hh_lung_seg) <- c("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")

DoAbsolute(Seg = hh_lung_seg, platform = "Illumina_WES", copy.num.type = "total",
           results.dir = "~/MRes_project_1/docs/HH_lung/absolute/jan_30", 
           keepAllResult = TRUE, verbose = TRUE)