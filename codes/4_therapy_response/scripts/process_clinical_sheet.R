## read in and combine files 
sample_ids <- read.table("~/MRes_project_1/Codes/3_therapy_resp/sample_ids.txt") %>% unlist
radio_omics <- read.csv("~/MRes_project_1/docs/HH_ova/clinical_data/Clinical_data_HH_Immunoradiomics_full_2.csv", 
                        header = TRUE)
radio_omics$sample <- radio_omics$Row.names
# Step 1: Filter out duplicates where Unique_code is not 1
radio_omics_filtered <- radio_omics %>%
  group_by(sample) %>%
  filter(if(any(Unique_code == 1)) Unique_code == 1 else TRUE) %>%
  ungroup()
# step 2 remove unwanted columns 
radio_omics_filtered <- subset(radio_omics_filtered, select = -c(X.1, Row.names, ID1, Grade_di, ID2, CD163, PDL1, Alpha.SMA, Ki67, CD14, Cohort, 
                                                                 Sample.Origin_2020, Sample.ID.short, Unique_code, Serous, X, S..no.))
colnames(radio_omics_filtered) <- c("RPV", "age", "stage", "surgery_outcome", "post_opt_tumour_free", "os_event", "os_time", 
                                    "pfs_time", "pfs_event", "ca125", "surgery_type", "molecular_subtype", "primary_chemo_outcome", 
                                    "thickness", "storage", "CD20", "CD4", "CD8", "PD1", "CD4_FOXP3", "CD8_PD1", "sample")
radio_omics_filtered <- radio_omics_filtered %>% as.data.frame
rownames(radio_omics_filtered) <- radio_omics_filtered$sample
radio_omics_filtered <- radio_omics_filtered[sample_ids, ]

## cellularity
cellularity <- readxl::read_xlsx("~/MRes_project_1/docs/HH_ova/clinical_data/HH_tumourcellularity.xlsx")
colnames(cellularity) <- c("sample_name", "sample", "QC") ## rename columns 
cellularity <- filter(cellularity, cellularity$sample %in% sample_ids) ## select exist samples 
cellularity <- na.omit(cellularity) ## remove nas 
cellularity <- distinct(cellularity, sample, .keep_all = TRUE) %>% as.data.frame ## remove duplicates 
rownames(cellularity) <- cellularity$sample
cellularity_2 <- cellularity[sample_ids, "QC"] %>% as.data.frame()
rownames(cellularity_2) <- sample_ids
colnames(cellularity_2) <- "QC"
## combine 
hh_ova <- merge(radio_omics_filtered, cellularity_2, by = "row.names")
row.names(hh_ova) <- hh_ova$Row.names
hh_ova_immune <- hh_ova[, c("ca125", "CD20", "CD4", "CD8", "PD1", "CD4_FOXP3", "CD8_PD1")]

write.table(hh_ova, "~/MRes_project_1/docs/HH_ova/clinical_data/hh_ova.txt", quote = FALSE, 
            sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(hh_ova_immune, "~/MRes_project_1/docs/HH_ova/clinical_data/hh_ova_immune.txt", quote = FALSE, 
            sep = "\t", row.names = TRUE, col.names = TRUE)

temp <- read.table("~/MRes_project_1/docs/HH_ova/clinical_data/hh_ova.txt", header = TRUE)





