## library packages 
library(dplyr)
library(ggplot2)
library(magrittr)
library(readxl)
library(MCPcounter)
library(gtools)
library(tidyr)


## Plot 1: TCGA and Immune Features Dot Plot =====================================================
## data ------------------------------------------------------------------------------------------
## function to create the file if it does not exist
create_df <- function(table_names, cancer_type, unwanted_cols=NULL, fdr=0.20, select_cols, estimate_thres=0){
  final_df <- data.frame(cancer_type=character(),
                         chr_arm=character(),
                         phenotypes=character(),
                         estimate=numeric(),
                         fdr=numeric(),
                         stringsAsFactors=FALSE)
  colnames(final_df) <- select_cols
  
  
  cancer_type <- cancer_type
  table_names <- table_names
  
  for(table_name in table_names){
    table <- get(table_name, envir=.GlobalEnv)
    
    for (phenotype in names(table)) {
      current_df <- table[[phenotype]]
      
      if (!'chr_arm' %in% colnames(current_df)) {
        current_df$chr_arm <- rownames(current_df)
      }
      
      current_df$phenotypes <- phenotype
      current_df$cancer_type <- tail(str_split(table_name, "_")[[1]], 1)
      current_df <- current_df %>% dplyr::select(all_of(select_cols))
      
      final_df <- rbind(final_df, current_df)
    }
  }
  
  ## filter 
  final_df <- final_df %>% dplyr::filter(!phenotypes %in% c("t.cell", "cd8_t.cell", "Fraction Genome Altered", 
                                                            "Endothelial cells", "Fibroblasts", "Neutrophils", 
                                                            "Myeloid dendritic cells", "CD8 T cells", "NK cells", 
                                                            "Cytotoxic lymphocytes"))
  
  ## filtering and add colors 
  final_df <- final_df[which(final_df[, 5] < fdr), ]
  final_df$color <- ifelse(final_df[, 4] > estimate_thres, "#FF407D", "#1B3C73")
  final_df[, 4] <- abs(final_df[, 4])
  final_df$alpha <- ifelse(final_df[, 5] < fdr, 1-final_df[, 5], 0)
  final_df$alpha_norm <- (final_df$alpha-min(final_df$alpha))/(max(final_df$alpha)-min(final_df$alpha))
  
  ## save file 
  return(final_df)
}

## data from previous analysis 
result_imm_lm_b_ova <- readRDS("~/MRes_project_1/codes/5_plots/files/result_imm_lm_b_ova.rds")
result_imm_lm_b_lung <- readRDS("~/MRes_project_1/codes/5_plots/files/result_imm_lm_b_lung.rds")
immuno_ova <- readRDS("~/MRes_project_1/codes/5_plots/files/immuno_ova.rds")
immuno_lung <- readRDS("~/MRes_project_1/codes/5_plots/files/immuno_lung.rds")

## apply function 
table_names <- c("result_imm_lm_b_ova", "result_imm_lm_b_lung", "immuno_ova", "immuno_lung")
cancer_type = c("lung", "ova")
unwanted_cols <- c("t.cell", "cd8_t.cell", "Fraction Genome Altered", 
                   "Endothelial cells", "Fibroblasts", "Neutrophils", 
                   "Myeloid dendritic cells", "CD8 T cells", "NK cells", 
                   "Cytotoxic lymphocytes")
select_cols <- c('cancer_type', 'chr_arm', 'phenotypes', 'estimate', 'FDR')
final_df <- create_df(table_names = table_names, cancer_type = cancer_type, 
                      unwanted_cols = unwanted_cols, select_cols = select_cols)

## small tweaks 
final_df$phenotypes[final_df$phenotypes=="tls_rna_seq$TLS_mean"] <- "TLS"
final_df$phenotypes[final_df$phenotypes=="TMB (nonsynonymous)"] <- "TMB"

## save file 
write.csv(final_df, "~/MRes_project_1/codes/5_plots/files/tcgaDotPlot.csv")

## plot  -------------------------------------------------------------------------------------
final_df <- read.csv("~/MRes_project_1/codes/5_plots/files/tcgaDotPlot.csv")
final_df <- final_df %>% 
  dplyr::filter(!chr_arm %in% c("Xp", "Xq")) %>% 
  dplyr::mutate(phenotypes = case_when(phenotypes == "B lineage" ~ "B cells", 
                                       phenotypes == "Monocytic lineage" ~ "Monocytes", 
                                       TRUE ~ phenotypes))
final_df$chr <- gsub("p|q", "", final_df$chr_arm) %>% as.numeric 
final_df$arm <- str_sub(final_df$chr_arm,-1,-1)
final_df <- final_df %>% group_by(cancer_type) %>% arrange(cancer_type, chr, arm) %>% 
  as.data.frame()
final_df$alpha_norm <- cut(final_df$alpha, 
                           breaks = c(0.80, 0.90, 0.95, 0.9999, 1), 
                           labels = c("0.20", "0.10", "0.05", "1e-4"), 
                           include.lowest = TRUE)
x_levels <- c("1p", "1q", "2p", "2q", "3p", "3q", "4p", "4q", "5p", "5q", 
  "6p", "6q", "7p", "7q", "8p", "8q", "9p", "9q", "10p", "10q", 
  "11p", "11q", "12p", "12q", "13p", "13q", "14p", "14q", "15p", "15q", 
  "16p", "16q", "17p", "17q", "18p", "18q", "19p", "19q", "20p", "20q", 
  "21p", "21q", "22p", "22q") %>% rev()
final_df$chr_arm <- factor(final_df$chr_arm, levels = x_levels)
y_levels <- c("PRF1", "MHCII", "MHCI", "GZMB", "GNLY", "CD274", "TMB", "TLS", "T cells", "Monocytes", "B cells")
final_df$phenotypes <- factor(final_df$phenotypes, levels = y_levels)

final_df$cancer_type <- factor(final_df$cancer_type,
                           levels = c("lung", "ova"),
                           labels = c("Lung", "Ovarian"))
final_df$match <- paste0(final_df$cancer_type, "_", final_df$chr_arm)

p <- ggplot(final_df, aes(x = chr_arm, y = phenotypes, size = estimate, color = color, alpha = alpha_norm)) + #, alpha = alpha_norm
  geom_point() + theme_bw() + coord_flip() +
  scale_size_continuous(range = c(0.1, 3.5)) + 
  scale_alpha_discrete(name = "FDR") + 
  scale_color_identity() + 
  labs(x = "Chromosome Arm", y = "Immune Phenotypes") +
  facet_grid(. ~ cancer_type, scales = "free_x", space = "free_x") +  # Adjust for your specific factor variable
  theme(axis.text.x = element_text(angle = 45, hjust = 1), #, hjust = 1, vjust = 0.5
        panel.spacing = unit(0, "lines"), 
        legend.position= "right", 
        legend.box = "vertical", 
        strip.background = element_blank(),
        strip.placement = "outside", 
        text = element_text(size = 6))
ggsave("~/MRes_project_1/codes/5_plots/plots/plot_1/Tcga_dotplots.png", p, width=3.6, height=5.43, units = "in", dpi = 600)



## plot 2 =======================================================================================
tcga_ova <- read.table("~/MRes_project_1/docs/00_mitéra/clinical/processed/tcga_ova.txt", sep = "\t", header = TRUE, row.names = 1)
os_ova <- read.table("~/MRes_project_1/codes/5_plots/files/tcga_ova_os.txt", sep = "\t", header = TRUE, row.names = 1)

# prepare 
os_ova <- os_ova %>% filter(p_values < 0.05)
os_ova$Index <- rownames(os_ova) #factor(rownames(os_ova), levels = rownames(os_ova))
text_data <- data.frame(Index = os_ova$Index, HR = os_ova$HR, 
                        label = paste0("p = ", format.pval(round(os_ova$p_values, digits = 4))))
x_max <- max(log10(c(os_ova$upper_CI, abs(os_ova$lower_CI))))

# plot
p <- ggplot(os_ova, aes(x = Index, y = HR)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.2) +
  coord_flip() + # Flip coordinates to make it a traditional forest plot layout 
  scale_y_reverse() + 
  theme_bw() +
  labs(y = "Hazard Ratio (HR)", x = "Chromosome Arms") +
  geom_hline(yintercept = 1) + 
  #ggtitle("Ovarian Cancer") + 
  ylim(0, 2) + 
  geom_text(data = text_data, aes(x = Index, y = 0.25, label = label), check_overlap = TRUE, size = 2) + 
  theme(text = element_text(size = 6)) 
ggsave(filename = "~/MRes_project_1/codes/5_plots/plots/plot_1/forest_plot_ovarian.png", 
       plot = p, height = 1.5, width = 3.10, units = "in", dpi = 600)





## plot 3 =======================================================================================

long_data <- tcga_ova %>% 
  pivot_longer(cols = "chr2q", names_to = "chromosome_arm", values_to = "cna") %>% na.omit()
cna_group0 <- long_data$cna[long_data$resp.rate == 0]
cna_group1 <- long_data$cna[long_data$resp.rate == 1]
wilcox.test(cna_group0, cna_group1)
## significant 


# Create the boxplot for a specific chromosome arm, e.g., '15q'
p <- ggplot(tcga_ova, aes(x = MSI, y = chr2p))+
  #long_data %>% filter(chromosome_arm == "chr2q"), aes(x = as.factor(resp.rate), y = cna)) +
  #geom_violin() +
  geom_point()+
  labs(x = "Response Rate", y = "Chromosome 2q CNA", caption = "p-value = 0.04519") +
  theme_bw() + #ggtitle("Lung Cancer") + 
  theme(text = element_text(size = 6), 
        plot.caption = element_text(hjust = 0.5, vjust = 120, size = 6)) 
p
ggsave(filename = "~/MRes_project_1/codes/5_plots/plots/plot_1/violin_plot_ova.png", plot = p, height = 2.5, width = 2.5, units = "in", dpi = 600)





## plot 3 =======================================================================================
tcga_lung <- read.table("~/MRes_project_1/docs/00_mitéra/clinical/processed/tcga_lung_mtx2.txt", sep = "\t", header = TRUE, row.names = 1)


coxph(Surv(os.time, os.event)~tcga_lung[, i]+TMB+stage.code+sex+sample.type+age, data = tcga_lung)


chromosome_cols <- colnames(tcga_lung)[grep(pattern = "^chr", colnames(tcga_lung))]

os_result <- matrix(data = NA, nrow = length(chromosome_cols), ncol = 5) %>% as.data.frame()
rownames(os_result) <- chromosome_cols
colnames(os_result) <- c("HR", "upper_CI", "lower_CI", "p_values", "FDR")

for(i in chromosome_cols){
  test <- coxph(Surv(os.time, os.event)~tcga_lung[, i]+TMB+stage.code+sex+sample.type+age, data = tcga_lung)
  test <- summary(test)
  
  os_result[i, "HR"] <- test$coefficients["tcga_lung[, i]", "exp(coef)"]
  os_result[i, "p_values"] <- test$coefficients["tcga_lung[, i]", "Pr(>|z|)"]
  os_result[i, "upper_CI"] <- test$conf.int["tcga_lung[, i]", "upper .95"]
  os_result[i, "lower_CI"] <- test$conf.int["tcga_lung[, i]", "lower .95"]
}

os_result$FDR <- p.adjust(os_result$p_values, method = "fdr")
rownames(os_result)[rownames(os_result) == "chrq"] <- "chrXq"
rownames(os_result)[rownames(os_result) == "chrp"] <- "chrXp"
#write.table(os_result, "~/MRes_project_1/codes/5_plots/files/tcga_lung_os.txt", sep = "\t", quote = F, row.names = T, col.names = T)



# prepare 
os_result <- read.table("~/MRes_project_1/codes/5_plots/files/tcga_lung_os.txt", sep = "\t", header = T, row.names = 1)
os_result <- os_result %>% filter(p_values < 0.10)
os_result$Index <- rownames(os_result) #factor(rownames(os_result), levels = rownames(os_result))
text_data <- data.frame(Index = os_result$Index, HR = os_result$HR, 
                        label = paste0("p = ", format.pval(round(os_result$p_values, digits = 4))))
x_max <- max(log10(c(os_result$upper_CI, abs(os_result$lower_CI))))

# plot
p <- ggplot(os_result, aes(x = Index, y = HR)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.2) +
  coord_flip() + # Flip coordinates to make it a traditional forest plot layout 
  scale_y_reverse() + 
  theme_bw() +
  labs(y = "Hazard Ratio (HR)", x = "Chromosome Arms") +
  geom_hline(yintercept = 1) + 
  #ggtitle("Ovarian Cancer") + 
  ylim(0, 2) + 
  geom_text(data = text_data, aes(x = Index, y = 0.25, label = label), check_overlap = TRUE, size = 2) + 
  theme(text = element_text(size = 6)) 
ggsave(filename = "~/MRes_project_1/codes/5_plots/plots/plot_1/forest_plot_lung.png", plot = p, height = 1.5, width = 3.10, units = "in", dpi = 600)










