# library 
library(readxl)
library(dplyr)
library(survival)
library(magrittr)
library(survminer)
library(ggsurvfit)
library(forestploter)
library(pheatmap)
library(tidyverse)

# functions 
source("~/MRes_project_1/codes/5_plots/script/function.R")

# dataframe 
lung.mtx <- read.table("~/MRes_project_1/codes/5_plots/file_final/lung_mtx.txt", sep = "\t", header = T, row.names = 1)
hh_ova.mtx <- read.table("~/MRes_project_1/codes/5_plots/file_final/hh_ova_mtx.txt", sep = "\t", header = T, row.names = 1)
ova.mtx <- read.table("~/MRes_project_1/codes/5_plots/file_final/ova_mtx.txt", sep = "\t", header = T, row.names = 1)

## plot 1 Kaplan Meier Plot ----------------------------------------------------------
conf_int = TRUE
name = "_"

# ovarian cancer 
ova.mtx.filt <- ova.mtx %>% dplyr::filter(os.time < 5000) 
#ova.mtx.filt <- ova.mtx.filt[-which(row.names(ova.mtx.filt) %in% c("TCGA.13.0886.01", "TCGA.24.2020.01")), ]
fit <- survfit(Surv(os.time, os.event) ~ as.factor(chr2qStatus), data = ova.mtx.filt)
ggsurv <- ggsurvplot_customised(dataframe = ova.mtx.filt, pval_coord = c(200, 0.25), break_time_by = 500, conf_int = conf_int, 
                                title = "Kaplan-Meier Curve for Ovarian Cancer", legend_title = "status", 
                                subtitle = "OS, hh_ova + tcga_ova (670), surv days < 5000") 
ggsurv$plot + geom_vline(xintercept = c(1826), size = 0.3, alpha = 0.5)


ggsurv

## plot 2 forest plot ----------------------------------------------------------------
# test: coxph test 
# method: continuous cna 
# direction: NA
## calculate correlation 
if(!file.exists("~/MRes_project_1/codes/5_plots/file_final/ova_survival.csv")){
  chr_cols <- colnames(hh_ova.mtx)[12:54]
  
  ova_survival <- matrix(NA, nrow = length(chr_cols), ncol = 5) %>% as.data.frame()
  rownames(ova_survival) <- chr_cols
  colnames(ova_survival) <- c("HR", "lower_CI", "upper_CI", "p_values", "FDR")
  
  for(arms in chr_cols){
    coxphmodel <- coxph(Surv(os.time, os.event)~hh_ova.mtx[, arms]+age+stage+TMB, data = hh_ova.mtx)
    
    ova_survival[arms, "HR"] <- summary(coxphmodel)$coef[1, 2]
    ova_survival[arms, "lower_CI"] <- summary(coxphmodel)$conf.int[1, 3]
    ova_survival[arms, "upper_CI"] <- summary(coxphmodel)$conf.int[1, 4]
    ova_survival[arms, "p_values"] <- summary(coxphmodel)$coef[1, 5]
  }
  ova_survival$FDR <- p.adjust(ova_survival$p_values, method = "fdr")
  write.csv(ova_survival, "~/MRes_project_1/codes/5_plots/file_final/ova_survival.csv")
} else {
  ova_survival <- read.csv("~/MRes_project_1/codes/5_plots/file_final/ova_survival.csv", row.names = 1)
}

## plot 
ova_survival$Index <- factor(rownames(ova_survival), levels = rownames(ova_survival))
#text_data <- data.frame(Index = ova_survival$Index, HR = ova_survival$HR, 
#                        label = paste0("p = ", format.pval(round(ova_survival$p_values, digits = 4))))
x_max <- max(log10(c(ova_survival$upper_CI, abs(ova_survival$lower_CI))))

chr_choice <- rev(c("chr2p", "chr2q", "chr4q", "chr8p", "chr9p", "chr12p", "chr12q", 
                    "chr13q", "chr17p", "chr20p", "chr20q", "chr22q"))
temp <- ova_survival[chr_choice, ]
temp$FDR <- p.adjust(temp$p_values, method = "fdr")
#text_data <- data.frame(Index = temp$Index, HR = temp$HR, 
#                        label = paste0("FDR = ", format.pval(round(temp$FDR, digits = 3))))

p <- ggplot(temp, aes(x = Index, y = HR)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.2) +
  coord_flip() + # Flip coordinates to make it a traditional forest plot layout
  theme_bw() +
  labs(y = "log10 Hazard Ratio (HR)", x = "Chromosome Arms") +
  geom_hline(yintercept = 1) + 
  ggtitle("Ovarian Cancer") + ylim(-5, 20) + 
  #geom_text(data = text_data, aes(x = Index, y = -3 ,label = label), check_overlap = TRUE, size = 3) + 
  theme(title = element_text(size = 9))
ggsave(filename = "~/MRes_project_1/codes/5_plots/plots/plot_2/forest_plot_ovarian.png", 
       plot = p, height = 5, width = 3, units = "in", dpi = 600)










# violin plot -------------------------------------------------------------------------------
long_data <- hh_ova.mtx %>% 
  pivot_longer(cols = "chr2q", names_to = "chromosome_arm", values_to = "cna") %>% na.omit()
cna_group0 <- long_data$cna[long_data$resp.rate == 0]
cna_group1 <- long_data$cna[long_data$resp.rate == 1]
wilcox.test(cna_group0, cna_group1)

# Create the boxplot for a specific chromosome arm, e.g., '15q'
p <- ggplot(long_data, aes(x = as.factor(resp.rate), y = chr2p)) +
  geom_violin() +
  labs(x = "Response Rate", y = "Chromosome 2q CNA", 
       caption = paste0("p-value = ", wilcox.test(cna_group0, cna_group1)$p.value)) +
  ggtitle("Ovarian Cancer") +
  theme_bw() +
  theme(title = element_text(size = 9),
        plot.caption = element_text(hjust = 0, size = 8)) 

ggsave(filename = "~/MRes_project_1/codes/5_plots/plots/plot_2/violin_plot_ova.png", plot = p, height = 2.5, width = 2.5, units = "in", dpi = 600)



# heatmap plot -------------------------------------------------------------------------------
# test: wilcox test 
# method: dichotimised chr 
# direction: amp
hh_ova <- read.table("~/MRes_project_1/docs/00_mitéra/clinical/hh_ova.txt", row.names = 1, sep = "\t", header = T)

#chr_cols <- colnames(hh_ova)[grep(pattern = "^chr.*Amp$", colnames(hh_ova))]
chr_cols <- chr_choice
immune_cols <- c("CD20", "CD4", "CD8", "PD1", "CD4_FOXP3", "CD8_PD1")

result <- matrix(NA, ncol = length(immune_cols), nrow = length(chr_cols)) %>% as.data.frame
colnames(result) <- immune_cols
rownames(result) <- chr_cols

significance <- matrix(NA, ncol = length(immune_cols), nrow = length(chr_cols)) %>% as.data.frame
colnames(significance) <- immune_cols
rownames(significance) <- chr_cols

for(i in immune_cols){
  for(j in chr_cols){
    ## direction 
    test <- lm(hh_ova[, i]~hh_ova[, j]+stage+age+sex, data = hh_ova)
    test <- summary(test)
    result[j, i] <- test$coefficients[2, 1]
    
    ## p value 
    j <- paste0(j, "Amp")
    test2 <- wilcox.test(hh_ova[, i]~hh_ova[, j], data = hh_ova, na.action = "na.omit")
    significance[j, i] <- test2$p.value 
  }
}


FDR <- matrix(NA, ncol = length(immune_cols), nrow = length(chr_cols)) %>% as.data.frame
colnames(FDR) <- immune_cols
rownames(FDR) <- chr_cols

FDR$CD20 <- p.adjust(significance$CD20, method = "fdr")
FDR$CD4 <- p.adjust(significance$CD4, method = "fdr")
FDR$CD8 <- p.adjust(significance$CD8, method = "fdr")
FDR$PD1 <- p.adjust(significance$PD1, method = "fdr")
FDR$CD4_FOXP3 <- p.adjust(significance$CD4_FOXP3, method = "fdr")
FDR$CD8_PD1 <- p.adjust(significance$CD8_PD1, method = "fdr")

color <- 1-FDR+0.4

pheatmap(t(result), number_format = "%.2f", #display_numbers = T, 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c("#1B3C73", "white", "#FF407D"))(100),
         breaks = seq(-1, 1, length.out = 101),
         annotation_legend = TRUE) 


pheatmap(t(FDR), number_format = "%.2f", #display_numbers = T, 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c("#FF407D", "white"))(100),
         breaks = seq(0, 1, length.out = 101),
         annotation_legend = TRUE) 









