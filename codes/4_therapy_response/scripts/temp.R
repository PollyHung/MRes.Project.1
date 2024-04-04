format_igv_seg = function(facets_output, sample_id, normalize = TRUE) {
  
  if (!all(c('snps', 'segs', 'dipLogR') %in% names(facets_output))) {
    stop(paste0('Input is missing segs, snps or dipLogR ojbect.'), call. = FALSE)
  }
  
  seg = facets_output$segs
  
  if (normalize) { seg = mutate(seg, adj.seg.mean = cnlr.median - facets_output$dipLogR) }
  data.frame(seg)
}

hh_ova <- list.files(path = "~/MRes_project_1/docs/HH_lung/facets/facet_cval_50/")
segmentation <- matrix(NA, nrow = 1, ncol = 19) %>% as.data.frame
colnames(segmentation) <- c(column_names)
hh_ova <- setdiff(hh_ova, "README.md")

for(i in hh_ova){
  temp <- readRDS(paste0("~/MRes_project_1/docs/HH_lung/facets/facet_cval_50/", i,"/", i, "_hisens.rds"))
  seg <- format_igv_seg(temp, i)
  seg$sample <- i
  segmentation <- rbind(segmentation, seg)
}

segmentation <- segmentation[2:nrow(segmentation), ]

write.table(segmentation, file = "~/hh_lung_comp.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)







