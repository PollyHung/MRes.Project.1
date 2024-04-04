library(survival)
library(survminer)
library(ggplot2)
library(gridExtra)


## cnv 
cbioportal_lung <- read.table("~/MRes_project_1/docs/GISTIC2/cbioportal/lung/broad_values_by_arm.txt", 
                              header = TRUE, sep="\t", row.names = 1)
hh_lung <- read.table("~/MRes_project_1/docs/GISTIC2/hh/lung/broad_values_by_arm.txt", 
                      header = TRUE, sep="\t", row.names = 1)
#lung <- merge(cbioportal_lung, hh_lung, by="row.names")
#rownames(lung) <- rownames(cbioportal_lung)
#lung <- lung[, 2:ncol(lung)]

cbioportal_lung_surv2 <- read.table("~/MRes_project_1/docs/GISTIC2/clinical_docs/cbioportal_lung_surv.txt")
cbioportal_lung_surv2 <- cbioportal_lung_surv2[colnames(cbioportal_lung), ]
hh_lung_surv2 <- read.table("~/MRes_project_1/docs/GISTIC2/clinical_docs/hh_lung_surv.txt")
hh_lung_surv2 <- hh_lung_surv2[colnames(hh_lung), ]

hh_ova <- read.table("~/MRes_project_1/docs/GISTIC2/hh/ova/broad_values_by_arm.txt", 
                     header = TRUE, sep="\t", row.names = 1)
cbioportal_mel <- read.table("~/MRes_project_1/docs/GISTIC2/cbioportal/mel/broad_values_by_arm.txt", 
                             header = TRUE, sep="\t", row.names = 1)

## kaplan meier plot 
dataset <- cbioportal_lung
survival <- cbioportal_lung_surv
hh_ova_surv <- hh_ova_surv[colnames(hh_ova), ]

hh_ova_surv$CD8 %>% summary

#plot.list <- list()
chrosomes <- c("2q", "10p")
for(i in chrosomes){
  chr <- dataset["15q", ] %>% as.numeric
  names(chr) <- colnames(dataset)
  chr <- cut(chr, c(-10,-0.0001,+0.0001, 10))
  survival$chr <- chr
  surv_obj <- Surv(time=survival[, "pfs_time"], event=survival[, "pfs_event"])
  fit <- survfit(surv_obj ~ chr, data = survival)
  
  ggsurvplot(fit, data = survival, pval = TRUE, conf.int = TRUE,
                  palette = c("blue", "red", "green"),
                  xlab = "Time", ylab = "Survival probability",
                  title = paste("Kaplan-Meier Survival Curve Based on 11q Copy Number Level"))
}

plot.list


surv_obj





## clinical docs 
cbioportal_lung_surv <- read.table("~/MRes_project_1/docs/GISTIC2/clinical_docs/cbioportal_lung.txt")
rownames(cbioportal_lung_surv) <- gsub("-", ".", rownames(cbioportal_lung_surv))
cbioportal_lung_surv <- cbioportal_lung_surv[colnames(cbioportal_lung), ]
hh_lung_surv <- read.table("~/MRes_project_1/docs/GISTIC2/clinical_docs/hh_lung.txt") %>% 
  dplyr::rename(clinical_benefit=response_rate)
lung_surv <- rbind(cbioportal_lung_surv[c("pfs_event", "pfs_time", "clinical_benefit")], 
                   hh_lung_surv[c("pfs_event", "pfs_time", "clinical_benefit")])
rownames(lung_surv) <- gsub("-", ".", rownames(lung_surv))
lung_surv <- lung_surv[colnames(lung), ]

cbioportal_mel_surv <- read.table("~/MRes_project_1/docs/GISTIC2/clinical_docs/cbioportal_mel.txt")
cbioportal_mel_surv <- cbioportal_mel_surv[c("os_event", "os_time", "clinical_benefit")]
cbioportal_mel_surv <- cbioportal_mel_surv[colnames(cbioportal_mel), ]

hh_ova_surv <- read.table("~/MRes_project_1/docs/GISTIC2/clinical_docs/hh_ova.txt")
hh_ova_surv <- hh_ova_surv[c("os_event", "os_time", "surgery_outcome")] %>% 
  dplyr::rename(clinical_benefit=surgery_outcome)
hh_ova_surv <- hh_ova_surv[colnames(hh_ova), ]


hh_ova_surv$patient_id <- rownames(hh_ova_surv)
cbioportal_mel_surv$patient_id <- rownames(cbioportal_mel_surv)
lung_surv$patient_id <- rownames(lung_surv)

hh_ova_surv$cancer_type <- "ovarian"
cbioportal_mel_surv$cancer_type <- "melanoma"
lung_surv$cancer_type <- "lung"

lung_surv$chr2q <- as.numeric(lung["2q", ]) 
lung_surv$chr10p <- as.numeric(lung["10p", ])

cbioportal_mel_surv$chr2q <- as.numeric(cbioportal_mel["2q", ]) 
cbioportal_mel_surv$chr10p <- as.numeric(cbioportal_mel["10p", ]) 

hh_ova_surv$chr2q <- as.numeric(hh_ova["2q", ]) 
hh_ova_surv$chr10p <- as.numeric(hh_ova["10p", ]) 


## barplots 

barplots <- rbind(lung_surv[c("patient_id", "cancer_type", "chr2q", "chr10p", "clinical_benefit")], 
                  hh_ova_surv[c("patient_id", "cancer_type", "chr2q", "chr10p", "clinical_benefit")], 
                  cbioportal_mel_surv[c("patient_id", "cancer_type", "chr2q", "chr10p", "clinical_benefit")])
melted_barplots <- reshape2::melt(barplots, 
                  id.vars = c('patient_id', 'cancer_type', 'clinical_benefit'), 
                  measure.vars = c('chr2q', 'chr10p'), 
                  variable.name = 'chr_arm', 
                  value.name = 'value')
melted_barplots$clinical_benefit <- ifelse(melted_barplots$clinical_benefit == 1, "response", "no response")
melted_barplots <- melted_barplots %>% drop_na(clinical_benefit)
p <- ggplot(melted_barplots, aes(x = clinical_benefit, y = value, fill = clinical_benefit)) +
  geom_boxplot() + facet_grid(chr_arm ~ cancer_type) + 
  labs(x = "Treatment Response", y = "Chromosome Arm-Level Aberration", 
       title = "Comparison of Chromosome Arm-Level Aberration by Treatment Response") +
  theme_bw()

qc <- read.table("~/MRes_project_1/docs/HH_ova/facets/quality_control/summary_table.txt", header=TRUE)
qc_hisens <- qc %>% dplyr::filter(run_type == "hisens")



hist(qc_hisens$purity, main = "purity", xlab="purity") 



ggplot(qc_hisens, aes(x = as.factor(clinical_benefit), y = chr15q, group = clinical_benefit)) +
  geom_boxplot() + theme_bw() + scale_x_discrete(labels = c("0" = "No Response", "1" = "Response")) +
  labs(x = "Treatment Response", y = "Chromosome Arm-Level Aberration", 
       title = "Treatment Response by 15q") 













