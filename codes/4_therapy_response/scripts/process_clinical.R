################################################################################
clinical <- read.delim("~/MRes_project_1/GISTIC2/clinical_docs/combined_study_clinical_data (3).tsv")
clinical$Patient.ID <- NULL
clinical$CIBERSORT.Absolute <- NULL
clinical$CIBERSORT.Correlation <- NULL
clinical$CIBERSORT.P.value <- NULL
clinical$CIBERSORT.RMSE <- NULL
clinical$Treatment <- NULL
clinical$age <- c(na.omit(clinical$Age.at.Diagnosis.1), na.omit(clinical$Age.at.Diagnosis))
clinical$Age.at.Diagnosis <- NULL
clinical$Age.at.Diagnosis.1 <- NULL
clinical$ESTIMATE.Immune.Score <- NULL
clinical$ESTIMATE.Score <- NULL
clinical$ESTIMATE.Stromal.Score <- NULL

columns <- c("sample", "os_event", "mut_load", 
             "neo_antigen_load", "CD8", "M0", 
             "M1", "M2", "CD4_mem_act", 
             "dendritic_act", "clinical_benefit", "stage", 
             "sex", "ploidy", "purity", 
             "response_duration", "TMB", "treatment_response", 
             "B_cell_memory", "B_cell_naive", "dendritic_rest", 
             "esoinophils", "mast_act", "mast_rest", 
             "monocyte", "neutrophils", "NK_act", 
             "NK_rest", "plasma_cells", "CD4_mem_rest", 
             "CD4_naive", "follicular_helper_T", "gamma_delta_T", "Tregs", "age")

colnames(clinical) <- columns 

rownames(clinical) <- clinical$sample
clinical$sample <- NULL
clinical$response_duration <- NULL
clinical$treatment_response <- NULL
clinical$os_event <- NULL

clinical$stage <- gsub("M0", 0, clinical$stage)
clinical$stage <- gsub("M1a|M1b|M1c", 1, clinical$stage)

clinical$sex <- gsub("Male", 0, clinical$sex)
clinical$sex <- gsub("Female", 1, clinical$sex)

clinical$clinical_benefit <- gsub("PD", -1, clinical$clinical_benefit)
clinical$clinical_benefit <- gsub("CR|LB", 2, clinical$clinical_benefit)
clinical$clinical_benefit <- gsub("PR", 1, clinical$clinical_benefit)
clinical$clinical_benefit <- gsub("SD|NB", 0, clinical$clinical_benefit)
clinical$clinical_benefit <- as.numeric(clinical$clinical_benefit)

#clinical <- apply(clinical, 2, as.numeric) %>% as.data.frame()

write.table(clinical, "~/MRes_project_1/GISTIC2/clinical_docs/cbioportal_mel_pheno.txt", 
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

################################################################################
survival <- read.delim("~/MRes_project_1/GISTIC2/clinical_docs/combined_study_clinical_data (4).tsv")
survival <- survival[c('Sample.ID', 'Overall.Survival.Status', 'Overall.Survival..Months.',
                       'Disease.Free.status', 'Disease.Free.Survival..Months.', 'Treatment.Response', 
                       'Response.Duration..Weeks.')]
rownames(survival) <- survival$Sample.ID
survival$Sample.ID <- NULL
colnames(survival) <- c("os_event", "os_time", "pfs_event", "pfs_time", "resp", "resp_months")

survival$pfs_event[is.na(survival$pfs_event)] <- survival$resp[is.na(survival$pfs_event)]
survival$resp_months <- (survival$resp_months)/4.33
survival$pfs_time[is.na(survival$pfs_time)] <- survival$resp_months[is.na(survival$pfs_time)]
survival$pfs_event <- gsub("^response|^Responder|^0:DiseaseFree", "0", survival$pfs_event) ## pfs 
survival$pfs_event <- gsub("^nonresponse|^Non-responder|^1:Recurred/Progressed", "1", survival$pfs_event) ## pfs 
survival$os_event <- gsub(":LIVING", "", survival$os_event) ## 0 = living 
survival$os_event <- gsub(":DECEASED", "", survival$os_event) ## 1 = deceased 

survival <- survival[c("os_event", "os_time", "pfs_event", "pfs_time")]

#survival <- apply(survival, 2, as.numeric) %>% as.data.frame()

write.table(survival, "~/MRes_project_1/GISTIC2/clinical_docs/cbioportal_mel_surv.txt", 
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)


sample_ids <- rownames(survival)
write.table(sample_ids, "~/MRes_project_1/GISTIC2/sample_ids/cbioportal_mel_ids.txt", sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)


clinical <- read.table("~/MRes_project_1/docs/HH_ova/clinical_data/ho_clinical.txt", sep = "\t", 
                       row.names = 1, header = TRUE)
clinical$sex <- 1
write.table(clinical, "~/MRes_project_1/docs/HH_ova/clinical_data/ho_clinical.txt", sep = "\t", 
            row.names = TRUE, col.names = TRUE, quote = FALSE)
