## segmentation files 
seg_file <- "~/MRes_project_1/docs/HH_ova/facets/seg_cval_50/hhova_hisens_diplogR.adjusted.seg"
segfile <- read.table(seg_file)
colnames(segfile) <- c("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")
write.table(segfile, "~/MRes_project_1/codes/7_biscut/BISCUT-py3/docs/hhova_hisens_diplogR.adjusted.seg", 
            quote = F, col.names = T, row.names = F, sep = "\t")

## command:
cmd <- "python BISCUT_preprocessing.py --input hhova_hisens_diplogR.adjusted.seg"
system(cmd)
## produced a folder called breakpoint_files_2024_03_17/hhova

