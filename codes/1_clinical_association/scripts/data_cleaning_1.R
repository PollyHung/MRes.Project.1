## dataset 
tcga_ova_cna <- read.table("~/MRes_project_1/docs/GISTIC2/tcga/ova/broad_values_by_arm.txt", sep = "\t", header = TRUE, row.names = 1)
tcga_lung_cna <- read.table("~/MRes_project_1/docs/GISTIC2/tcga/lung/broad_values_by_arm.txt", sep = "\t", header = TRUE, row.names = 1)
tcga_ova_mtx <- read.table("~/MRes_project_1/docs/00_mitéra/clinical/raw/tcgaOvaPheno.txt", sep = "\t", header = TRUE)
rownames(tcga_ova_mtx) <- gsub("-", ".", tcga_ova_mtx$sampleID)
tcga_lung_mtx <- read.table("~/MRes_project_1/docs/00_mitéra/clinical/raw/tcgaLungPheno.txt", sep = "\t", header = TRUE)
rownames(tcga_lung_mtx) <- gsub("-", ".", tcga_lung_mtx$sampleID)

## sample id 
sample_ova <- intersect(rownames(tcga_ova_mtx), colnames(tcga_ova_cna))
sample_lung <- intersect(rownames(tcga_lung_mtx), colnames(tcga_lung_cna))

## organize 
tcga_ova_cna <- tcga_ova_cna[, sample_ova]
tcga_lung_cna <- tcga_lung_cna[, sample_lung]
tcga_ova_mtx <- tcga_ova_mtx[sample_ova, ]
tcga_lung_mtx <- tcga_lung_mtx[sample_lung, ]

## tmb 
tcga_ova_tmb <- read.table("~/MRes_project_1/docs/GISTIC2/tcga/ova/ova_tmb.tsv", sep = "\t", header = TRUE)
rownames(tcga_ova_tmb) <- gsub("-", ".", tcga_ova_tmb$Sample.ID)
tcga_ova_tmb <- tcga_ova_tmb[sample_ova, ] 
tcga_lung_tmb <- read.table("~/MRes_project_1/docs/GISTIC2/tcga/lung/lung_tmb.tsv", sep = "\t", header = TRUE)
rownames(tcga_lung_tmb) <- gsub("-", ".", tcga_lung_tmb$Sample.ID)
tcga_lung_tmb <- tcga_lung_tmb[sample_lung, ] 

## 
tcga_ova_mtx <- merge(tcga_ova_mtx, tcga_ova_tmb, by = "row.names")
tcga_lung_mtx <- merge(tcga_lung_mtx, tcga_lung_tmb, by = "row.names")



write.csv(tcga_ova_mtx, "~/MRes_project_1/docs/00_mitéra/clinical/raw/tcga_ova_mtx.csv")
write.csv(tcga_lung_mtx, "~/MRes_project_1/docs/00_mitéra/clinical/raw/tcga_lung_mtx.csv")





tcga_ova_mtx <- read.table("~/MRes_project_1/docs/00_mitéra/clinical/raw/tcga_ova_mtx.csv", sep = ",", header = TRUE, row.names = "Row.names")
tcga_lung_mtx <- read.table("~/MRes_project_1/docs/00_mitéra/clinical/raw/tcga_lung_mtx.csv", sep = ",", header = TRUE, row.names = "Row.names")

tcga_lung_mtx <- tcga_lung_mtx %>% dplyr::select("age_at_initial_pathologic_diagnosis", "Stage", "gender", 
                                                 "primary_therapy_outcome_success", "Person.Cigarette.Smoking.History.Pack.Year.Value", 
                                                 "Fraction.Genome.Altered", 
                                                 "Overall.Survival..Months.", "Overall.Survival.Status", "TMB..nonsynonymous.")
tcga_lung_mtx <- tcga_lung_mtx %>% dplyr::rename(age = age_at_initial_pathologic_diagnosis, 
                                                 stage = Stage, 
                                                 sex = gender, 
                                                 therapy.resp = primary_therapy_outcome_success, 
                                                 pack.yr = Person.Cigarette.Smoking.History.Pack.Year.Value, 
                                                 FGA = Fraction.Genome.Altered, 
                                                 TMB = TMB..nonsynonymous., 
                                                 os.time = Overall.Survival..Months., 
                                                 os.event = Overall.Survival.Status)
tcga_lung_mtx$stage <- roman2int(gsub("A|B|C", "", tcga_lung_mtx$stage)) 
tcga_lung_mtx$sex <- ifelse(tcga_lung_mtx$sex == "FEMALE", 0, 1)
tcga_lung_mtx$therapy.resp <- case_when(tcga_lung_mtx$therapy.resp == "Complete Remission/Response" ~ 2, 
                                        tcga_lung_mtx$therapy.resp == "Partial Remission/Response" ~ 1, 
                                        tcga_lung_mtx$therapy.resp == "Stable Disease" ~ 0, 
                                        tcga_lung_mtx$therapy.resp == "Progressive Disease" ~ -1, 
                                        TRUE ~ NA_integer_)
tcga_lung_mtx$resp.rate <- ifelse(tcga_lung_mtx$therapy.resp > 0, 1, 0)
tcga_lung_mtx$dfi.event <- ifelse(tcga_lung_mtx$dfi.event == "0:ALIVE OR DEAD TUMOR FREE", 0, 1) ## 0 = tumour free, 1 = tumour 
tcga_lung_mtx$os.event <- ifelse(tcga_lung_mtx$os.event == "0:LIVING", 0, 1) # 0 = living, 1 = dead 
tcga_lung_mtx$pfs.event <- ifelse(tcga_lung_mtx$pfs.event == "0:CENSORED", 0, 1) # 0 = censored, 1 = progression 

write.csv(tcga_lung_mtx, file = "~/MRes_project_1/docs/00_mitéra/clinical/processed/tcga_ova_mtx.csv")









tcga_ova_mtx <- read.table("~/MRes_project_1/docs/00_mitéra/clinical/raw/tcga_ova_mtx.csv", sep = ",", header = TRUE, row.names = "Row.names")
tcga_lung_mtx <- read.table("~/MRes_project_1/docs/00_mitéra/clinical/raw/tcga_lung_mtx.csv", sep = ",", header = TRUE, row.names = "Row.names")

tcga_ova_mtx <- tcga_ova_mtx %>% dplyr::select("age_at_initial_pathologic_diagnosis", "clinical_stage", "gender", 
                                               "primary_therapy_outcome_success", "Aneuploidy.Score", 
                                               "Months.of.disease.specific.survival",  "Disease.specific.Survival.status", 
                                               "Fraction.Genome.Altered", "MSI.MANTIS.Score", 
                                               "Overall.Survival..Months.", "Overall.Survival.Status", 
                                               "Progress.Free.Survival..Months.", "Progression.Free.Status", "TMB..nonsynonymous.")
tcga_ova_mtx <- tcga_ova_mtx %>% dplyr::rename(age = age_at_initial_pathologic_diagnosis, 
                                               stage = clinical_stage, 
                                               sex = gender, 
                                               therapy.resp = primary_therapy_outcome_success, 
                                               aneuploidy.score = Aneuploidy.Score, 
                                               dfi.time = Months.of.disease.specific.survival, 
                                               dfi.event = Disease.specific.Survival.status, 
                                               FGA = Fraction.Genome.Altered, 
                                               MSI = MSI.MANTIS.Score, 
                                               TMB = TMB..nonsynonymous., 
                                               os.time = Overall.Survival..Months., 
                                               os.event = Overall.Survival.Status, 
                                               pfs.time = Progress.Free.Survival..Months., 
                                               pfs.event = Progression.Free.Status)
tcga_ova_mtx$stage <- roman2int(gsub("Stage |A|B|C", "", tcga_ova_mtx$stage)) 
tcga_ova_mtx$sex <- ifelse(tcga_ova_mtx$sex == "FEMALE", 0, NA)
tcga_ova_mtx$therapy.resp <- case_when(tcga_ova_mtx$therapy.resp == "Complete Remission/Response" ~ 2, 
                                       tcga_ova_mtx$therapy.resp == "Partial Remission/Response" ~ 1, 
                                       tcga_ova_mtx$therapy.resp == "Stable Disease" ~ 0, 
                                       tcga_ova_mtx$therapy.resp == "Progressive Disease" ~ -1, 
                                       TRUE ~ NA_integer_)
tcga_ova_mtx$resp.rate <- ifelse(tcga_ova_mtx$therapy.resp > 0, 1, 0)
tcga_ova_mtx$dfi.event <- ifelse(tcga_ova_mtx$dfi.event == "0:ALIVE OR DEAD TUMOR FREE", 0, 1) ## 0 = tumour free, 1 = tumour 
tcga_ova_mtx$os.event <- ifelse(tcga_ova_mtx$os.event == "0:LIVING", 0, 1) # 0 = living, 1 = dead 
tcga_ova_mtx$pfs.event <- ifelse(tcga_ova_mtx$pfs.event == "0:CENSORED", 0, 1) # 0 = censored, 1 = progression 

write.csv(tcga_ova_mtx, file = "~/MRes_project_1/docs/00_mitéra/clinical/processed/tcga_ova_mtx.csv")



tcga_ova_cna <- read.table("~/MRes_project_1/docs/GISTIC2/tcga/ova/broad_values_by_arm.txt", sep = "\t", header = TRUE, row.names = 1)
tcga_lung_cna <- read.table("~/MRes_project_1/docs/GISTIC2/tcga/lung/broad_values_by_arm.txt", sep = "\t", header = TRUE, row.names = 1)

## sample id 
sample_ova <- intersect(rownames(tcga_ova_mtx), colnames(tcga_ova_cna))
sample_lung <- intersect(rownames(tcga_lung_mtx), colnames(tcga_lung_cna))

## organize 
tcga_ova_cna <- tcga_ova_cna[, sample_ova]
tcga_lung_cna <- tcga_lung_cna[, sample_lung]



tcga_ova_cna <- t(tcga_ova_cna) %>% as.data.frame()
tcga_lung_cna <- t(tcga_lung_cna) %>% as.data.frame()




tcga_ova <- merge(tcga_ova_mtx, tcga_ova_cna, by="row.names")
rownames(tcga_ova) <- tcga_ova$Row.names
tcga_ova$Row.names <- NULL
tcga_lung <- merge(tcga_lung_mtx, tcga_lung_cna, by="row.names")
rownames(tcga_lung) <- tcga_lung$Row.names
tcga_lung$Row.names <- NULL



colnames(tcga_ova) <- gsub("X", "chr", colnames(tcga_ova))
colnames(tcga_lung) <- gsub("X", "chr", colnames(tcga_lung))












tcga_ova_mtx <- read.table("~/MRes_project_1/docs/00_mitéra/clinical/raw/tcga_ova_mtx.csv", sep = ",", header = TRUE, row.names = "Row.names")
tcga_ova_mtx <- tcga_ova_mtx %>% dplyr::select("age_at_initial_pathologic_diagnosis", "clinical_stage", "gender", 
                                               "primary_therapy_outcome_success", "Aneuploidy.Score", "sample_type", 
                                               "Months.of.disease.specific.survival",  "Disease.specific.Survival.status", 
                                               "Fraction.Genome.Altered", "MSI.MANTIS.Score", "X_PANCAN_miRNA_PANCAN", 
                                               "Overall.Survival..Months.", "Overall.Survival.Status", 
                                               "Progress.Free.Survival..Months.", "Progression.Free.Status", "TMB..nonsynonymous.")
tcga_ova_mtx <- tcga_ova_mtx %>% dplyr::rename(age = age_at_initial_pathologic_diagnosis, 
                                               stage = clinical_stage, 
                                               sex = gender, 
                                               therapy.resp = primary_therapy_outcome_success, 
                                               aneuploidy.score = Aneuploidy.Score, 
                                               dfi.time = Months.of.disease.specific.survival, 
                                               dfi.event = Disease.specific.Survival.status, 
                                               miRNA = X_PANCAN_miRNA_PANCAN, 
                                               FGA = Fraction.Genome.Altered, 
                                               MSI = MSI.MANTIS.Score, 
                                               TMB = TMB..nonsynonymous., 
                                               os.time = Overall.Survival..Months., 
                                               os.event = Overall.Survival.Status, 
                                               pfs.time = Progress.Free.Survival..Months., 
                                               pfs.event = Progression.Free.Status)
tcga_ova_mtx$stage <- roman2int(gsub("Stage |A|B|C", "", tcga_ova_mtx$stage)) 
tcga_ova_mtx$sex <- ifelse(tcga_ova_mtx$sex == "FEMALE", 0, NA)
tcga_ova_mtx$therapy.resp <- case_when(tcga_ova_mtx$therapy.resp == "Complete Remission/Response" ~ 2, 
                                       tcga_ova_mtx$therapy.resp == "Partial Remission/Response" ~ 1, 
                                       tcga_ova_mtx$therapy.resp == "Stable Disease" ~ 0, 
                                       tcga_ova_mtx$therapy.resp == "Progressive Disease" ~ -1, 
                                       TRUE ~ NA_integer_)
tcga_ova_mtx$resp.rate <- ifelse(tcga_ova_mtx$therapy.resp > 0, 1, 0)
tcga_ova_mtx$dfi.event <- ifelse(tcga_ova_mtx$dfi.event == "0:ALIVE OR DEAD TUMOR FREE", 0, 1) ## 0 = tumour free, 1 = tumour 
tcga_ova_mtx$os.event <- ifelse(tcga_ova_mtx$os.event == "0:LIVING", 0, 1) # 0 = living, 1 = dead 
tcga_ova_mtx$pfs.event <- ifelse(tcga_ova_mtx$pfs.event == "0:CENSORED", 0, 1) # 0 = censored, 1 = progression 
tcga_ova_mtx$sample_type <- ifelse(tcga_ova_mtx$sample_type == "Primary Tumor", 0, 1) # 0 = primary, 1 = recurrent 

write.csv(tcga_ova_mtx, file = "~/MRes_project_1/docs/00_mitéra/clinical/processed/tcga_ova_mtx.csv")



tcga_ova_cna <- read.table("~/MRes_project_1/docs/GISTIC2/tcga/ova/broad_values_by_arm.txt", sep = "\t", header = TRUE, row.names = 1)
sample_ova <- intersect(rownames(tcga_ova_mtx), colnames(tcga_ova_cna))
tcga_ova_cna <- tcga_ova_cna[, sample_ova]
tcga_ova_cna <- t(tcga_ova_cna) %>% as.data.frame()
tcga_ova <- merge(tcga_ova_mtx, tcga_ova_cna, by="row.names")
rownames(tcga_ova) <- tcga_ova$Row.names
tcga_ova$Row.names <- NULL
write.table(tcga_ova, "~/MRes_project_1/docs/00_mitéra/clinical/processed/tcga_ova.txt", quote = F, sep = "\t", row.names = T, col.names = T)
colnames(tcga_ova) <- gsub("X", "chr", colnames(tcga_ova))

tcga_ova$miRNA <- as.numeric(as.factor(tcga_ova$miRNA))


tcga_ova <- tcga_ova %>% drop_na(c("age", "stage", "TMB", "FGA", "sample_type", "miRNA", "os.time", "os.event"))


tcga_ova$age_rounded <- round(tcga_ova$age, digit = -1)

tcga_ova_transformed <- tcga_ova %>%
  mutate(across(starts_with("chr"), ~ifelse(. > 0, "amp", "wt")))


colnames(tcga_ova)[colnames(tcga_ova) == "chrp"] <- "chrXp"
colnames(tcga_ova)[colnames(tcga_ova) == "chrq"] <- "chrXq"
write.table(tcga_ova, "~/MRes_project_1/docs/00_mitéra/clinical/processed/tcga_ova_dichotomized.txt", quote = F, sep = "\t", row.names = T, col.names = T)


tcga_ova <- read.table("~/MRes_project_1/docs/00_mitéra/clinical/processed/tcga_ova_dichotomized.txt", sep = "\t", header = TRUE, row.names = 1)
chromosome_cols <- colnames(tcga_ova)[grep(pattern = "^chr", colnames(tcga_ova))]

os_result <- matrix(data = NA, nrow = length(chromosome_cols), ncol = 5) %>% as.data.frame()
rownames(os_result) <- chromosome_cols
colnames(os_result) <- c("HR", "upper_CI", "lower_CI", "p_values", "FDR")

for(i in chromosome_cols){
  test <- coxph(Surv(os.time, os.event)~tcga_ova[, i]+age_rounded+stage+TMB+FGA+sample_type+miRNA, data = tcga_ova)
  test <- summary(test)
  
  os_result[i, "HR"] <- test$coefficients["tcga_ova[, i]wt", "exp(coef)"]
  os_result[i, "p_values"] <- test$coefficients["tcga_ova[, i]wt", "Pr(>|z|)"]
  os_result[i, "upper_CI"] <- test$conf.int["tcga_ova[, i]wt", "upper .95"]
  os_result[i, "lower_CI"] <- test$conf.int["tcga_ova[, i]wt", "lower .95"]
}

os_result$FDR <- p.adjust(os_result$p_values)
rownames(os_ova)[rownames(os_ova) == "chrq"] <- "chrXq"
rownames(os_ova)[rownames(os_ova) == "chrp"] <- "chrXp"
write.table(os_result, "~/MRes_project_1/codes/5_plots/files/tcga_ova_os.txt", sep = "\t", quote = F, row.names = T, col.names = T)

pfs_result <- matrix(data = NA, nrow = length(chromosome_cols), ncol = 5) %>% as.data.frame()
rownames(pfs_result) <- chromosome_cols
colnames(pfs_result) <- c("HR", "upper_CI", "lower_CI", "p_values", "FDR")

for(i in chromosome_cols){
  test <- coxph(Surv(pfs.time, pfs.event)~tcga_ova[, i]+age_rounded+stage+TMB+FGA+sample_type+miRNA, data = tcga_ova)
  test <- summary(test)
  
  pfs_result[i, "HR"] <- test$coefficients["tcga_ova[, i]wt", "exp(coef)"]
  pfs_result[i, "p_values"] <- test$coefficients["tcga_ova[, i]wt", "Pr(>|z|)"]
  pfs_result[i, "upper_CI"] <- test$conf.int["tcga_ova[, i]wt", "upper .95"]
  pfs_result[i, "lower_CI"] <- test$conf.int["tcga_ova[, i]wt", "lower .95"]
}

pfs_result$FDR <- p.adjust(pfs_result$p_values)
write.table(pfs_result, "~/MRes_project_1/codes/5_plots/files/tcga_ova_pfs.txt", sep = "\t", quote = F, row.names = T, col.names = T)


tcga_ova <- read.table("~/MRes_project_1/docs/00_mitéra/clinical/processed/tcga_ova.txt", sep = "\t", header = TRUE, row.names = 1)
chromosome_cols <- colnames(tcga_ova)[grep(pattern = "^chr", colnames(tcga_ova))]
therapy.resp_result <- matrix(data = NA, nrow = length(chromosome_cols), ncol = 5) %>% as.data.frame()
rownames(therapy.resp_result) <- chromosome_cols
colnames(therapy.resp_result) <- c("estimate", "Std.Error", "t_values", "p_values", "FDR")

for(i in chromosome_cols){
  test <- lm(therapy.resp~tcga_ova[, i]+age+stage+TMB+FGA+sample_type+miRNA, data = tcga_ova)
  test <- summary(test)
  
  therapy.resp_result[i, "estimate"] <- test$coefficients["tcga_ova[, i]", "Estimate"]
  therapy.resp_result[i, "p_values"] <- test$coefficients["tcga_ova[, i]", "Pr(>|t|)"]
  therapy.resp_result[i, "Std.Error"] <- test$coefficients["tcga_ova[, i]", "Std. Error"]
  therapy.resp_result[i, "t_values"] <- test$coefficients["tcga_ova[, i]", "t value"]
}

therapy.resp_result$FDR <- p.adjust(therapy.resp_result$p_values)
write.table(therapy.resp_result, "~/MRes_project_1/codes/5_plots/files/tcga_ova_therapy.txt", sep = "\t", quote = F, row.names = T, col.names = T)






tcga_ova <- tcga_ova %>%
  mutate(across(starts_with("chr"), ~ifelse(. > 0, "amp", ifelse(. < 0, "del", "wt"))))

tcga_ova <- read.table("~/MRes_project_1/docs/00_mitéra/clinical/processed/tcga_ova.txt", sep = "\t", header = TRUE, row.names = 1)
chromosome_cols <- colnames(tcga_ova)[grep(pattern = "^chr", colnames(tcga_ova))]

os_result <- matrix(data = NA, nrow = length(chromosome_cols), ncol = 5) %>% as.data.frame()
rownames(os_result) <- chromosome_cols
colnames(os_result) <- c("HR", "upper_CI", "lower_CI", "p_values", "FDR")

for(i in chromosome_cols){
  test <- coxph(Surv(os.time, os.event)~tcga_ova[, i]+age_rounded+stage+TMB+FGA+sample_type+miRNA, data = tcga_ova)
  test <- summary(test)
  
  os_result[i, "HR"] <- test$conf.int["tcga_ova[, i]wt", "exp(-coef)"]
  os_result[i, "p_values"] <- test$coefficients["tcga_ova[, i]wt", "Pr(>|z|)"]
  os_result[i, "upper_CI"] <- test$conf.int["tcga_ova[, i]wt", "upper .95"]
  os_result[i, "lower_CI"] <- test$conf.int["tcga_ova[, i]wt", "lower .95"]
}

os_result$FDR <- p.adjust(os_result$p_values)
rownames(os_ova)[rownames(os_ova) == "chrq"] <- "chrXq"
rownames(os_ova)[rownames(os_ova) == "chrp"] <- "chrXp"
write.table(os_result, "~/MRes_project_1/codes/5_plots/files/tcga_ova_os.txt", sep = "\t", quote = F, row.names = T, col.names = T)







