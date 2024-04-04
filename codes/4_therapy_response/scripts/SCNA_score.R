## associating SNCA scores with various 
setwd("~/MRes_project_1/GISTIC2/SCNA_scores")

folder_name <- list.files()

for(i in folder_name){
  cancer <- list.files(i)
  for(j in cancer){
    file_path <- paste0(i, "/", j, "/normalised_SCNA_score.csv")
    if(file.exists(file_path)){
      scna <- read.csv(file_path, row.names = 1, header = TRUE)
      
      scna_name <- paste0(i, "_", j) %>% make.names
      assign(scna_name, scna, envir = .GlobalEnv)
    }

  }
}

## HH_lung
clinical <- read.table("~/MRes_project_1/GISTIC2/clinical_docs/hh_lung_pheno.txt", 
                       sep = "\t", row.names = 1, header = TRUE)
survival <- read.table("~/MRes_project_1/GISTIC2/clinical_docs/hh_lung_surv.txt", 
                       sep = "\t", row.names = 1, header = TRUE)
temp <- merge(clinical, survival, by="row.names")
rownames(temp) <- temp$Row.names
temp$Row.names <- NULL
rownames(hh_lung) <- gsub("-", ".", rownames(hh_lung))
hh_lung <- merge(temp, hh_lung, by="row.names")
rownames(hh_lung) <- hh_lung$Row.names
hh_lung$Row.names <- NULL
write.table(hh_lung, "~/MRes_project_1/GISTIC2/clinical_docs/hh_lung.txt", 
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

## HH_ova
clinical <- read.table("~/MRes_project_1/GISTIC2/clinical_docs/hh_ova_pheno.txt", 
                       sep = "\t", row.names = 1, header = TRUE)
survival <- read.table("~/MRes_project_1/GISTIC2/clinical_docs/hh_ova_surv.txt", 
                       sep = "\t", row.names = 1, header = TRUE)
temp <- merge(clinical, survival, by="row.names")
rownames(temp) <- temp$Row.names
temp$Row.names <- NULL
hh_ova <- merge(hh_ova, temp, by = "row.names")
rownames(hh_ova) <- hh_ova$Row.names
hh_ova$Row.names <- NULL
write.table(hh_ova, "~/MRes_project_1/GISTIC2/clinical_docs/hh_ova.txt", 
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)


## cBioPortal melanoma 
clinical <- read.table("~/MRes_project_1/GISTIC2/clinical_docs/cbioportal_mel_pheno.txt", 
                       sep = "\t", row.names = 1, header = TRUE)
survival <- read.table("~/MRes_project_1/GISTIC2/clinical_docs/cbioportal_mel_surv.txt", 
                       sep = "\t", row.names = 1, header = TRUE)
temp <- merge(clinical, survival, by="row.names")
rownames(temp) <- temp$Row.names
temp$Row.names <- NULL
cbioportal_mel <- merge(cbioportal_mel, temp, by = "row.names")
rownames(cbioportal_mel) <- cbioportal_mel$Row.names
cbioportal_mel$Row.names <- NULL
write.table(cbioportal_mel, "~/MRes_project_1/GISTIC2/clinical_docs/cbioportal_mel.txt", 
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

## cBioPortal lung 
clinical <- read.table("~/MRes_project_1/GISTIC2/clinical_docs/nsclc_pd1_msk_2018_clinical_data.tsv", 
                       sep = "\t", header = TRUE)
survival <- read.table("~/MRes_project_1/GISTIC2/clinical_docs/KM_Plot__Progression_Free_(months).txt", 
                       sep = "\t", header = TRUE)
clinical <- merge(clinical, survival, by="Patient.ID")

clinical <- select(clinical, -c(Patient.ID, Progression.Free.Status, Somatic.Status, Study.ID))
rownames(clinical) <- clinical$Sample.ID
clinical$Sample.ID <- NULL

colnames(clinical) <- c("mut_count", "FGA", "sex", "clinical_benefit", "lines_of_treatment", 
                        "mutation_rate", "PD_L1", "smoker", "TMB", "treatment_type", 
                        "age", "pfs_event", "pfs_time")
clinical$pfs_event <- ifelse(clinical$pfs_event == "0:Not Progressed", 0, 1)
clinical$sex <- ifelse(clinical$sex == "Female", 1, 0) ## female = 1, male = 0
clinical$clinical_benefit <- ifelse(clinical$clinical_benefit == "YES", 1, 0) ## yes = 1, no = 0
clinical$smoker <- ifelse(clinical$smoker == "Ever", 1, 0) ## ever = 1, never = 0
clinical$treatment_type <- ifelse(clinical$treatment_type == "Monotherapy", 0, 1) ## monotherapy = 0, combinationa = 1

rowNames <- rownames(clinical)
clinical <- apply(clinical, 2, as.numeric) %>% as.data.frame()
rownames(clinical) <- rowNames

write.table(clinical, "~/MRes_project_1/GISTIC2/clinical_docs/cbioportal_lung.txt", 
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)


## corrleation 
cbioportal_mel <- read.table("~/MRes_project_1/GISTIC2/clinical_docs/cbioportal_mel.txt", 
                             sep = "\t", row.names = 1, header = TRUE)
hh_ova <- read.table("~/MRes_project_1/GISTIC2/clinical_docs/hh_ova.txt", 
                     sep = "\t", row.names = 1, header = TRUE)
hh_lung <- read.table("~/MRes_project_1/GISTIC2/clinical_docs/hh_lung.txt", 
                     sep = "\t", row.names = 1, header = TRUE)
save.image("~/MRes_project_1/GISTIC2/SCNA_scores/scna_scores.RData")


lm(hh_lung$rank_norm_chrom ~ hh_lung$best_response) %>% summary
lm(hh_lung$rank_norm_chrom ~ hh_lung$os_time) %>% summary

hh_lung$os_surv <- Surv(hh_lung$os_time, hh_lung$os_event)
hh_lung$pfs_surv <- Surv(hh_lung$pfs_time, hh_lung$pfs_event)

coxph(hh_lung$os_surv~hh_lung$sum_score)
coxph(hh_lung$pfs_surv~hh_lung$rank_norm_chrom)






