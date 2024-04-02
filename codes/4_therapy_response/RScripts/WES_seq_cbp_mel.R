library(AnnotationDbi)
library(org.Hs.eg.db) 
library(tidyr)
library(dplyr)
rpkm2tpm <- function(x) {
  x <- as.matrix(x)
  return(t(t(x)/colSums(x))*10^6)
}



## section 1: combine mRNA seq from 3 studies and convert them from rkpm to tpm -----------------------------------------
setwd("/rds/general/user/ph323/home/MRes_project_1/docs/WES")

dfci <- read.table("dfci/data_mrna_seq_rpkm.txt", header = TRUE)
msk <- read.table("msk/data_mrna_seq_rpkm.txt", header = TRUE)
ucla <- read.table("ucla/data_mrna_seq_rpkm.txt", header = TRUE)

dataframe <- merge(dfci, msk, by="Entrez_Gene_Id")
dataframe <- merge(dataframe, ucla, by="Entrez_Gene_Id")
dataframe$Entrez_Gene_Id <- as.character(dataframe$Entrez_Gene_Id)

dataframe$Entrez_Gene_Id <- mapIds(org.Hs.eg.db,
                                   keys = dataframe$Entrez_Gene_Id,
                                   column = "SYMBOL",
                                   keytype = "ENTREZID",
                                   multiVals = "first")

dataframe <- dataframe[-which(duplicated(dataframe$Entrez_Gene_Id)), ]
dataframe <- dataframe %>% drop_na(Entrez_Gene_Id)
rownames(dataframe) <- dataframe$Entrez_Gene_Id
dataframe$Entrez_Gene_Id <- NULL

dataframe <- rpkm2tpm(dataframe) ## transform to TPM 

write.table(dataframe, "~/MRes_project_1/docs/WES/cbp_mel/mRNAseq.txt", 
            row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)




## section 2: cleaning up and adding additional columns for clinical.mtx --------------------------------
## read in clinical matrix 
clinical.mtx <- read_tsv("~/MRes_project_1/docs/WES/cbp_mel/clinical_mtx.tsv") 
colnames(clinical.mtx) <- gsub(" ", "_", colnames(clinical.mtx)) %>% tolower
rownames(clinical.mtx) <- clinical.mtx$patient_id
clinical.mtx <- clinical.mtx[colnames(dataframe), ]
clinical.mtx <- clinical.mtx %>% 
  dplyr::select(c("patient_id", "mutation_load", "neo-antigen_load", "cytolytic_score", 
                  "estimate_immune_score", "estimate_stromal_score", 
                  "durable_clinical_benefit", "m_stage", "sex", "tmb_(nonsynonymous)", 
                  "age_at_diagnosis...50", "age_at_diagnosis...30"), starts_with("cbset"))

## create tls signature 
tls_signature <- c("CD79B", "CD1D", "CCR6", "LAT", "SKAP1", "CETP", "EIF1AY", "RBP5", "PTGDS", "CCL19", 
                   "CCL21", "CXCL13", "CCR7", "CXCR5", "SELL", "LAMP3")
tls <- dataframe[tls_signature, ]
tls <- tls %>% t %>% as.matrix %>% scale(center=F) %>% as.data.frame
tls$tls <- apply(tls, 1, mean)
tls$patient_id <- rownames(tls)
tls <- tls[, 17:18]
clinical.mtx <- merge(clinical.mtx, tls, by="patient_id")

## add other genes and MHCI and MHCII
genes <- c("PRF1", "GZMB", "GNLY", "CD274", 
           "HLA-A", "HLA-B", "HLA-C", ## MHCI
           "HLA-DRB1", "HLA-DQB1", "HLA-DPB1", "HLA-DRB5", "HLA-DOB", "HLA-DQB2") ## MHCII
genes <- dataframe[genes, ] %>% t %>% as.data.frame
genes <- genes %>%
  mutate(MHCI = rowMeans(.[, c("HLA-A", "HLA-B", "HLA-C")], na.rm = TRUE),
         MHCII = rowMeans(.[, c("HLA-DRB1", "HLA-DQB1", "HLA-DPB1", "HLA-DRB5", "HLA-DOB", "HLA-DQB2")], na.rm = TRUE))
genes <- genes %>% dplyr::select(c("PRF1", "GZMB", "GNLY", "CD274", "MHCI", "MHCII"))
genes$patient_id <- rownames(genes)
clinical.mtx <- merge(clinical.mtx, genes, by="patient_id")

clinical.mtx$age_at_diagnosis...30 <- ifelse(is.na(clinical.mtx$age_at_diagnosis...30), 
                                             clinical.mtx$age_at_diagnosis...50, clinical.mtx$age_at_diagnosis...30)
clinical.mtx$age_at_diagnosis...50 <- NULL

## rename column names 
colnames(clinical.mtx) <- c("patient_id", "mut.load", "neo.antigen", 
                            "cytolytic.score", "immune.score", "stromal.score", 
                            "treatment.resp", "stage", "sex", 
                            "tmb", "age", "cd8.t", 
                            "M0", "M1", "M2", 
                            "cd4.t.mem.act", "dendritic.act", "b.mem", 
                            "b.naive", "dendritic.rest", "eosinophils", 
                            "mast.act", "mast.rest", "monocytes", 
                            "neutrophils", "nk.act", "nk.rest", 
                            "plasma.cell", "cd4.t.mem.rest", "cd4.t.naive", 
                            "t.help.follicular", "t.γδ", "t.regs", 
                            "tls", "PRF1", "GZMB", "GNLY", "CD274", "MHCI", "MHCII")

survival <- read.table("~/MRes_project_1/docs/WES/cbp_mel/survival.txt", sep="\t", header=TRUE)
survival <- survival[, 2:4]
colnames(survival) <- c("patient_id", "os_event", "os_time")

clinical.mtx <- merge(clinical.mtx, survival, by="patient_id")
clinical.mtx$os_event <- ifelse(clinical.mtx$os_event == "0:LIVING", 0, 1) ## 1 survived, 0 deceased 

## if LB (long term benefit), PR (partial response), or CR (complete response) then its response (1) else no response (0)
clinical.mtx$treatment.resp <- ifelse(clinical.mtx$treatment.resp %in% c("LB", "PR", "CR"), 1, 0)
clinical.mtx$sex <- ifelse(clinical.mtx$sex == "Male", 1, 0) ## male = 1, female = 0
clinical.mtx$stage <- ifelse(clinical.mtx$stage == "M0", 3, 4) ## M0 = stage 3, all others with M1 is stage 4 

rownames(clinical.mtx) <- clinical.mtx$patient_id


write.table(clinical.mtx, "~/MRes_project_1/docs/GISTIC2/clinical_docs/cbp_mel_87obs.txt", sep="\t", quote=FALSE, 
            row.names=FALSE, col.names=TRUE)


## section 3: linear model analysis of selected immune features against cna --------------------------------
# read in CNV object 
cnv <- read.table("~/MRes_project_1/docs/GISTIC2/cbioportal/mel/broad_values_by_arm.txt", 
                  header = TRUE, sep="\t", row.names = 1)

# align object patient ids 
cnv <- cnv[, which(colnames(cnv) %in% clinical.mtx$patient_id)]
clinical.mtx <- clinical.mtx[-78, ]

# create survival object 
surv_obj <- Surv(clinical.mtx$os_time, clinical.mtx$os_event)

# univariate analysis to find confounding variables 
result <- matrix(NA, nrow=ncol(clinical.mtx[, 2:40]), ncol=3) %>% as.data.frame
rownames(result) <- colnames(clinical.mtx[, 2:40])
colnames(result) <- c("HR", "p_val", "FDR")

for (i in colnames(clinical.mtx[, 2:40])) {
  clinical_obj <- as.numeric(clinical.mtx[, i])
  test <- coxph(surv_obj~clinical_obj) %>% summary
  result[i, "HR"] <- test$coefficients[1]
  result[i, "p_val"] <- test$coefficients[5]
}

result$FDR <- p.adjust(result$p_val, method = "BH")
result <- result %>% arrange(FDR)

# multivariate analysis 
survival_df <- matrix(NA, nrow=nrow(cnv), ncol=5) %>% as.data.frame
rownames(survival_df) <- rownames(cnv)
colnames(survival_df) <- c("HR", "lower.95", "upper.95", "p_val", "FDR")

for(i in rownames(cnv)){
  cnv_obj <- as.numeric(cnv[i, ])
  test <- coxph(surv_obj~cnv_obj+clinical.mtx$tmb+clinical.mtx$sex+clinical.mtx$stage+clinical.mtx$age) %>% summary

  survival_df[i, "HR"] <- test$coefficients["cnv_obj", 2]
  survival_df[i, "p_val"] <- test$coefficients["cnv_obj", 5]
  survival_df[i, "lower.95"] <- test$conf.int["cnv_obj", 3]
  survival_df[i, "upper.95"] <- test$conf.int["cnv_obj", 4]
}

survival_df$FDR <- p.adjust(survival_df$p_val, method = "BH")
survival_df <- survival_df %>% arrange(FDR)





  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  