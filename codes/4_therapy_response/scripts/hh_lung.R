library(survival)
library(survminer)
library(ggplot2)
library(gridExtra)

dataset <- read.table("~/MRes_project_1/docs/GISTIC2/hh/lung_2/broad_values_by_arm.txt", header = TRUE, sep="\t", row.names = 1)
survival <- read.table("~/MRes_project_1/docs/GISTIC2/clinical_docs/hh_lung.txt")
survival <- survival[colnames(dataset), ] 

hh_lung <- matrix(data=NA, nrow=nrow(dataset), ncol=5) %>% as.data.frame
colnames(hh_lung) <- c("exp.Coef", "lower.95", "upper.95", "p_val", "FDR")
rownames(hh_lung) <- rownames(dataset)

surv_obj <- Surv(time=survival[, "os_time"], event=survival[, "os_event"])
for(i in 1:nrow(dataset)){
  chr <- dataset[i, ] %>% as.numeric
  model <- coxph(surv_obj~chr+survival$sex+survival$age+survival$PD_L1+survival$stage) %>% summary
  hh_lung[i, "exp.Coef"] <- model$coefficients[1, "exp(coef)"]
  hh_lung[i, "p_val"] <- model$coefficients[1, "Pr(>|z|)"]
  hh_lung[i, "lower.95"] <- model$conf.int[1, "lower .95"]
  hh_lung[i, "upper.95"] <- model$conf.int[1, "upper .95"]
}

hh_lung$FDR<-p.adjust(hh_lung$p_val,method="fdr")

write.table(hh_lung, file="~/MRes_project_1/docs/GISTIC2/statical_analysis/hh_lung/coxph(os~chr).txt", 
            sep="\t", row.names = TRUE, col.names = TRUE, quote=FALSE)


## survival 
chr <- dataset["15q", ] %>% as.numeric
names(chr) <- colnames(dataset)
chr <- cut(chr, c(-10, +0.0001, 10))
survival$chr <- chr
surv_obj <- Surv(time=survival[, "pfs_time"], event=survival[, "pfs_event"])
fit <- survfit(surv_obj ~ chr, data = survival)

ggsurvplot(fit, data = survival, pval = TRUE, conf.int = FALSE,
           palette = "npg",
           xlab = "Time", ylab = "Survival probability",
           title = paste("Survival Curve on 15q Copy Number Level"), 
           legend=("top"), legend.title="15q", 
           legend.labs = c("deletion", "amplification"))





## volcano plot 
volcano <- matrix(data=NA, nrow=nrow(dataset), ncol=5) %>% as.data.frame
colnames(volcano) <- c("mean_0", "mean_1", "log2FC", "p_val", "FDR")
rownames(volcano) <- rownames(dataset)

clinical_benefit <- survival$response_rate %>% as.factor

for(i in 1:nrow(dataset)){
  chr <- dataset[i, ] %>% as.numeric
  model <- t.test(chr~clinical_benefit) 
  volcano[i, "mean_0"] <- model$estimate[1] ## no response 
  volcano[i, "mean_1"] <- model$estimate[2] ## response 
  FC <- model$estimate[2] / model$estimate[1]
  volcano[i, "log2FC"] <- log2(abs(FC)) * sign(FC)
  volcano[i, "p_val"] <- model$p.value
}

volcano$FDR <- p.adjust(volcano$p_val,method="fdr")

volcano$neg_log10_p_val <- -log10(volcano$p_val)
volcano$significant <- volcano$p_val < 0.05
label_candidates <- volcano[volcano$significant & abs(volcano$log2FC) > 1,]


p <- ggplot(volcano, aes(x = log2FC, y = neg_log10_p_val)) +
  geom_point(aes(color = significant), alpha = 0.7) +
  scale_color_manual(values = c('FALSE' = 'grey', 'TRUE' = 'red')) +
  labs(x = "log10(FC)", y = "-log10(p_val)", 
       title = "Chromosome Arm Aneuploidy Status 
       in Hammersmith Lung Dataset") +
  theme_bw() +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "blue") + 
  geom_text_repel(label=rownames(volcano)) 

ggsave("~/MRes_project_1/results/plot/hh_lung/volcanoPlot.png", p, width=5, height=4, units = "in") 



## box plot 
chr15q <-dataset["15q", ] %>% as.numeric

index <- which(is.na(survival$response_rate))
survival <- survival[-index, ]
dataset <- dataset[, -index]

ggplot(survival, aes(x = as.factor(response_rate), y = chr15q, group = response_rate)) +
  geom_boxplot() + theme_bw() + scale_x_discrete(labels = c("0" = "No Response", "1" = "Response")) +
  labs(x = "Treatment Response", y = "Chromosome Arm-Level Aberration", 
       title = "Treatment Response by 15q") 


