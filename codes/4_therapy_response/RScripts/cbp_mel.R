library(survival)
library(survminer)
library(ggplot2)
library(gridExtra)

#dataset <- read.table("~/MRes_project_1/docs/GISTIC2/cbioportal/lung/broad_values_by_arm.txt", header = TRUE, sep="\t", row.names = 1)
dataset <- read.table("~/MRes_project_1/docs/GISTIC2/cbioportal/mel_2/broad_values_by_arm.txt", header = TRUE, sep="\t", row.names = 1)
survival <- read.table("~/MRes_project_1/docs/GISTIC2/clinical_docs/cbioportal_mel.txt")
survival <- survival[colnames(dataset), ] ## align 


## cox regression 
cbp_mel <- matrix(data=NA, nrow=nrow(dataset), ncol=5) %>% as.data.frame
colnames(cbp_mel) <- c("exp.Coef", "lower.95", "upper.95", "p_val", "FDR")
rownames(cbp_mel) <- rownames(dataset)

surv_obj <- Surv(time=survival[, "os_time"], event=survival[, "os_event"])
for(i in 1:nrow(dataset)){
  chr <- dataset[i, ] %>% as.numeric
  model <- coxph(surv_obj~chr+survival$age+survival$sex+survival$stage+survival$TMB) %>% summary
  cbp_mel[i, "exp.Coef"] <- model$coefficients[1, "exp(coef)"]
  cbp_mel[i, "p_val"] <- model$coefficients[1, "Pr(>|z|)"]
  cbp_mel[i, "lower.95"] <- model$conf.int[1, "lower .95"]
  cbp_mel[i, "upper.95"] <- model$conf.int[1, "upper .95"]
}

cbp_mel$FDR<-p.adjust(cbp_mel$p_val,method="fdr")

write.table(cbp_mel, file="~/MRes_project_1/docs/GISTIC2/statical_analysis/cbp_mel/coxph(pfs~chr)_new.txt", 
            sep="\t", row.names = TRUE, col.names = TRUE, quote=FALSE)


## kaplan meier curve 
chr <- dataset["3q", ] %>% as.numeric
names(chr) <- colnames(dataset)
chr <- cut(chr, c(-10,-0.0001,+0.0001, 10))
survival$chr <- chr
surv_obj <- Surv(time=survival[, "os_time"], event=survival[, "os_event"])
fit <- survfit(surv_obj ~ chr, data = survival)

#png("~/MRes_project_1/results/plot/cbp_mel/chr15survival.png", width=5, height=4, units = "in")
ggsurvplot(fit, data = survival, pval = TRUE, conf.int = FALSE,
           palette = "npg",
           xlab = "Time", ylab = "Survival probability",
           title = paste("Survival Curve on 3q Copy Number Level"), 
           legend=("top"), legend.title="3q", 
           legend.labs = c("deletion", "neutral", "amplification"))
#dev.off()


## box plot 
survival$response <- ifelse(survival$clinical_benefit < 1, 0, 1) %>% as.factor
response <- survival$response

index <- which(is.na(survival$response))
survival <- survival[-index, ]
dataset <- dataset[, -index]


chr3q <-dataset["3q", ] %>% as.numeric

ggplot(survival, aes(x = as.factor(response), y = chr3q, group = response)) +
  geom_boxplot() + theme_bw() + scale_x_discrete(labels = c("0" = "No Response", "1" = "Response")) +
  labs(x = "Treatment Response", y = "Chromosome Arm-Level Aberration", 
       title = "Treatment Response by 3q") 
