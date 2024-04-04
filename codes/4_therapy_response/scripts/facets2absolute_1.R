## This is a script that converts facets result to absolute result 
## The Facets result is from original FACETS package 

## Section 1: environment ======================================================
OVA <- "~/MRes_project_1/docs/HH_ova/facets_original/dataframe/"
LUNG <- "~/MRes_project_1/docs/HH_lung/facets_original/dataframe/"
COLUMNS <- c("chrom", "start", "end", "num.mark", "cnlr.median")

## Section 2: list files and read in ===========================================
hh_ova <- list()
facets_ova <- list.files(OVA)
facets_lung <- list.files(LUNG)

setwd(LUNG)
for(i in 1:length(facets_lung)){
  dataframe <- read.csv(facets_lung[i])
  dataframe$sample <- gsub(".csv", "", facets_lung[i])
  hh_ova[[i]] <- dataframe 
}

hh_ova_df <- do.call(rbind, hh_ova)

hh_ova_df <- hh_ova_df[, c("sample", COLUMNS)]
colnames(hh_ova_df) <- c("Sample", "Chromosome", "Start", "End", "NumMarkers", "Seg.CN")
write.table(hh_ova_df, "~/hh_lung_example.seg", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")


hh_ova_df$NumMarkers %>% summary


















