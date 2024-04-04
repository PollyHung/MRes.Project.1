## environment =================================================================
library(dplyr)
library(magrittr)
library(tibble)
library(survival)

## HH_ova ======================================================================
## process the radiomics dataframe 
radio_omics <- read.csv("~/MRes_project_1/docs/HH_ova/clinical_data/Clinical_data_HH_Immunoradiomics_full_2.csv", 
                        header = TRUE)
radio_omics <- radio_omics %>% group_by(Row.names) %>% filter(if(any(Unique_code == 1)) Unique_code == 1 else TRUE) %>% 
  ungroup() ## remove rows with unique code = 1 
radio_omics$sample <- radio_omics$Row.names

## process the cellularity dataframe 
cellularity <- readxl::read_xlsx("~/MRes_project_1/docs/HH_ova/clinical_data/HH_tumourcellularity.xlsx")
cellularity$sample <- cellularity$...2
cellularity <- cellularity[c("sample", "QC", "Sample name")]
cellularity <- filter(cellularity, cellularity$sample %in% radio_omics$sample)
cellularity <- na.omit(cellularity)
cellularity <- distinct(cellularity, sample, .keep_all = TRUE) %>% as.data.frame
rownames(cellularity) <- cellularity$sample

## combine the two dataframes 
combined <- merge(radio_omics, cellularity, by="sample", all.x=TRUE)

## data cleaning 
hh_ova_new <- subset(combined, select = -c(X.1, Row.names, ID1, Grade_di, ID2, CD163, PDL1, Alpha.SMA, Ki67, CD14, Cohort, 
                                           Sample.Origin_2020, Sample.ID.short, Unique_code, Serous, X, S..no., `Sample name`))
colnames(hh_ova_new) <- c("sample", "RPV", "age", "stage", "surgery_outcome", "post_opt_tumour_free", "os_event", "os_time", 
                          "pfs_time", "pfs_event", "ca125", "surgery_type", "molecular_subtype", "primary_chemo_outcome", 
                          "thickness", "storage", "CD20", "CD4", "CD8", "PD1", "CD4_FOXP3", "CD8_PD1", "cellularity")
hh_ova_new$di_CD20 <- as.integer(hh_ova_new$CD20 > (hh_ova_new$CD20 %>% na.omit %>% median))
hh_ova_new$di_CD8 <- as.integer(hh_ova_new$CD8 > (hh_ova_new$CD8 %>% na.omit %>% median))
hh_ova_new$di_CD4_FOXP3 <- as.integer(hh_ova_new$CD4_FOXP3 > (hh_ova_new$CD4_FOXP3 %>% na.omit %>% median))

hh_ova_new$storage <- gsub("Box ", "", hh_ova_new$storage) %>% as.numeric
hh_ova_new$primary_chemo_outcome[which(hh_ova_new$primary_chemo_outcome == "complete.response")] <- 2
hh_ova_new$primary_chemo_outcome[which(hh_ova_new$primary_chemo_outcome == "partial.response")] <- 1
hh_ova_new$primary_chemo_outcome[which(hh_ova_new$primary_chemo_outcome == "partial response ")] <- 1
hh_ova_new$primary_chemo_outcome[which(hh_ova_new$primary_chemo_outcome == "stable.disease")] <- 0
hh_ova_new$primary_chemo_outcome[which(hh_ova_new$primary_chemo_outcome == "progressive.disease")] <- -1
hh_ova_new$primary_chemo_outcome[which(hh_ova_new$primary_chemo_outcome == "")] <- NA
hh_ova_new$molecular_subtype <- gsub("C", "", hh_ova_new$molecular_subtype) %>% as.numeric

RowNames <- hh_ova_new$sample
hh_ova_new$sample <- NULL
hh_ova_new <- apply(hh_ova_new, 2, as.numeric) %>% as.data.frame()
rownames(hh_ova_new) <- RowNames

clinical <- hh_ova_new[c("ca125", "CD20", "CD4", "CD8", "PD1", "CD4_FOXP3", "CD8_PD1", "cellularity", "RPV", "age", 
                       "stage", "surgery_outcome", "post_opt_tumour_free", "surgery_type", "molecular_subtype", 
                       "primary_chemo_outcome", "thickness", "storage", "di_CD20", "di_CD8", "di_CD4_FOXP3")]
survival <- hh_ova_new[c("os_event", "os_time", "pfs_time", "pfs_event")]



write.table(clinical, "~/MRes_project_1/docs/HH_ova/clinical_data/ho_clinical.txt", sep = "\t", 
            row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(survival, "~/MRes_project_1/docs/HH_ova/clinical_data/ho_survival.txt", sep = "\t", 
            row.names = TRUE, col.names = TRUE, quote = FALSE)






## HH_lung ======================================================================
lung_clinical <- read.csv("~/MRes_project_1/docs/HH_ova/clinical_data/HH_lung_clinical.csv")
lung_clinical <- subset(lung_clinical, select = -c(X.10, Seg.ID, Sample.Number, Tube.No., 
                                                   Volume..L., Total.Mass..g., Library.Type, DEAD, 
                                                   Test.Result, Remark, X, Id, hospNr, DOB, DORx_start, 
                                                   DODx, DODth.DOLFU, IMM_Type, Current_Imm, 
                                                   ProgressionDate, PROGRESSION.RATE, CONTROL.FUNCTION, 
                                                   CaType, Histopath.N, Molec.path.N, M_Other_Location, 
                                                   X1L, X2L, X3L, Chemo_Type, X.1, X.2, X.3, X.4, 
                                                   X.5, X.6, X.7, X.8, X.9, FHC_code, Others_Grade, 
                                                   Others_Type, Others))
colnames(lung_clinical) <- c("sample", "sample_conc", "sample_integrity", "FHC", "sex", "age", "immunotherapy", 
                             "line_of_therapy", "subtype", "PS", "PD_L1", "stage", "M_Lungs", "M_Bone", "M_LN", 
                             "M_Adrenal", "M_Liver", "M_CNS", "M_Other", "surgery", "chemotherapy", 
                             "os", "best_response", "response_rate", "pfs_event", "pfs", 
                             "rash", "diarrhoea", "fatigue", "liver_tox", "thyroid_tox", "pituitary_tox", 
                             "nausea", "pneumo_tox", "max_side_effect_grade", "os_time", "pfs_time", "os_event")
## 
lung_clinical$immunotherapy[which(lung_clinical$immunotherapy == 6)] <- 1 #"atezolizumab"
lung_clinical$immunotherapy[which(lung_clinical$immunotherapy == 4)] <- 1 #"durvalumab"
lung_clinical$immunotherapy[which(lung_clinical$immunotherapy == 2)] <- 1 #"nivolumab"
lung_clinical$immunotherapy[which(lung_clinical$immunotherapy == 1)] <- 1 #"pembrolizumab"

lung_clinical$subtype[which(lung_clinical$subtype == "ADK")] <- 5 #"adenocarcinoma"
lung_clinical$subtype[which(lung_clinical$subtype == "adk")] <- 4 #"adenocarcinoma"
lung_clinical$subtype[which(lung_clinical$subtype == "SCC")] <- 3 #"small_cell"
lung_clinical$subtype[which(lung_clinical$subtype == "LAGE CELL")] <- 2 #"large_cell"
lung_clinical$subtype[which(lung_clinical$subtype == "pleomorph")] <- 1 #"pleomorphic"
lung_clinical$subtype[which(lung_clinical$subtype == "NOS-poor diff")] <- 0 #"poor_diff"

lung_clinical$sample_integrity[which(lung_clinical$sample_integrity == "DNA fragments > 500bp")] <- 5
lung_clinical$sample_integrity[which(lung_clinical$sample_integrity == "500bp > DNA fragments > 250bp")] <- 4
lung_clinical$sample_integrity[which(lung_clinical$sample_integrity == "DNA fragments < 250bp")] <- 3
lung_clinical$sample_integrity[which(lung_clinical$sample_integrity == "Degraded slightly")] <- 2
lung_clinical$sample_integrity[which(lung_clinical$sample_integrity == "Degraded Moderate")] <- 1

lung_clinical$FHC[which(lung_clinical$FHC == "negative")] <- -1
lung_clinical$FHC[which(lung_clinical$FHC == "high")] <- 2
lung_clinical$FHC[which(lung_clinical$FHC == "low")] <- 1
lung_clinical$FHC[which(lung_clinical$FHC == "")] <- NA

RowNames <- gsub("-", ".", lung_clinical$sample)
lung_clinical$sample <- NULL
lung_clinical <- apply(lung_clinical, 2, as.numeric) %>% as.data.frame()
row.names(lung_clinical) <- RowNames

clinical <- subset(lung_clinical, select = -c(os_time, os_event, pfs_time, pfs_event))
survival <- lung_clinical[c("os_time", "os_event", "pfs_time", "pfs_event")]


write.table(clinical, "~/MRes_project_1/docs/HH_ova/clinical_data/hl_clinical.txt", sep = "\t", 
            row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(survival, "~/MRes_project_1/docs/HH_ova/clinical_data/hl_survival.txt", sep = "\t", 
            row.names = TRUE, col.names = TRUE, quote = FALSE)

sample_ids <- read.table("~/MRes_project_1/Codes/3_therapy_response/HH_lung/facets/sample_ids.txt") %>% unlist
sample_ids <- gsub("-", ".", sample_ids)

write.table(sample_ids, "~/MRes_project_1/Codes/3_therapy_response/HH_lung/facets/sample_ids_2.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)





