segmentation <- read.table("~/MRes_project_1/GISTIC2/cbioportal/lung/segmentation.seg", header = TRUE)
segmentation$ID <- gsub("-", ".", segmentation$ID)
segmentation_list <- split(segmentation, segmentation$ID)

for (i in names(segmentation_list)) {
  temp_df <- segmentation_list[[i]]
  dup_ids <- temp_df %>% duplicated
  segmentation_list[[i]] <- temp_df[!dup_ids, ]
}

segmentation_cleaned <- do.call(rbind, segmentation_list)
write.table(segmentation_cleaned, "~/MRes_project_1/GISTIC2/cbioportal/lung/segmentation.seg", 
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
