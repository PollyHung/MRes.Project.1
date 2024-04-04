library(survival)
library(survminer)
library(ggplot2)
library(gridExtra)

#dataset <- read.table("~/MRes_project_1/docs/GISTIC2/cbioportal/lung/broad_values_by_arm.txt", header = TRUE, sep="\t", row.names = 1)
dataset <- read.table("~/MRes_project_1/docs/GISTIC2/hh/ova_2/broad_values_by_arm.txt", header = TRUE, sep="\t", row.names = 1)
survival <- read.table("~/MRes_project_1/docs/GISTIC2/clinical_docs/hh_ova.txt")
survival <- survival[colnames(dataset), ] ## align 


## cox regression 
hh_ova <- matrix(data=NA, nrow=nrow(dataset), ncol=5) %>% as.data.frame
colnames(hh_ova) <- c("exp.Coef", "lower.95", "upper.95", "p_val", "FDR")
rownames(hh_ova) <- rownames(dataset)

surv_obj <- Surv(time=survival[, "os_time"], event=survival[, "os_event"])
for(i in 1:nrow(dataset)){
  chr <- dataset[i, ] %>% as.numeric
  model <- coxph(surv_obj~chr+survival$age+survival$stage+survival$CD4+survival$CD8+survival$CD20) %>% summary
  hh_ova[i, "exp.Coef"] <- model$coefficients[1, "exp(coef)"]
  hh_ova[i, "p_val"] <- model$coefficients[1, "Pr(>|z|)"]
  hh_ova[i, "lower.95"] <- model$conf.int[1, "lower .95"]
  hh_ova[i, "upper.95"] <- model$conf.int[1, "upper .95"]
}

hh_ova$FDR<-p.adjust(hh_ova$p_val,method="fdr")

write.table(hh_ova, file="~/MRes_project_1/docs/GISTIC2/statical_analysis/hh_ova/coxph(pfs~chr)_new.txt", 
            sep="\t", row.names = TRUE, col.names = TRUE, quote=FALSE)


## kaplan meier curve 
chr <- dataset["13q", ] %>% as.numeric
names(chr) <- colnames(dataset)
chr <- cut(chr, c(-10,-0.0001,+0.0001, 10))
survival$chr <- chr
surv_obj <- Surv(time=survival[, "os_time"], event=survival[, "os_event"])
fit <- survfit(surv_obj ~ chr, data = survival)

#png("~/MRes_project_1/results/plot/hh_ova/chr15survival.png", width=5, height=4, units = "in")
ggsurvplot(fit, data = survival, pval = TRUE, conf.int = FALSE,
           palette = "npg",
           xlab = "Time", ylab = "Survival probability",
           title = paste("Survival Curve on 2q Copy Number Level"), 
           legend=("top"), legend.title="2q", 
           legend.labs = c("deletion", "neutral", "amplification"))
#dev.off()


## volcano plot 
volcano <- matrix(data=NA, nrow=nrow(dataset), ncol=5) %>% as.data.frame
colnames(volcano) <- c("mean_0", "mean_1", "log2FC", "p_val", "FDR")
rownames(volcano) <- rownames(dataset)

clinical_benefit <- survival$primary_chemo_outcome %>% as.factor

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
volcano$significant <- volcano$p_val < 0.20
label_candidates <- volcano[volcano$significant & abs(volcano$log2FC) > 1,]


p <- ggplot(volcano, aes(x = log2FC, y = neg_log10_p_val)) +
  geom_point(aes(color = significant), alpha = 0.7) +
  scale_color_manual(values = c('FALSE' = 'grey', 'TRUE' = 'red')) +
  labs(x = "log10(FC)", y = "p-value", 
       title = "Chromosome Arm Aneuploidy Status 
       in cBioPortal Lung Dataset") +
  theme_bw() +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "blue") + 
  geom_text_repel(label=rownames(volcano)) 

ggsave("~/MRes_project_1/results/plot/hh_ova/volcanoPlot.png", p, width=5, height=4, units = "in") 

# Display the plot
print(volcano_plot)




## box plot 


survival$response <- ifelse(survival$primary_chemo_outcome < 1, 0, 1) %>% as.factor

index <- which(is.na(survival$response))
survival <- survival[-index, ]
dataset <- dataset[, -index]

chr2q <-dataset["2q", ] %>% as.numeric

response <- survival$response
ggplot(survival, aes(x = as.factor(response), y = chr2q, group = response)) +
  geom_boxplot() + theme_bw() + scale_x_discrete(labels = c("0" = "No Response", "1" = "Response")) +
  labs(x = "Treatment Response", y = "Chromosome Arm-Level Aberration", 
       title = "Treatment Response by 2q") 
