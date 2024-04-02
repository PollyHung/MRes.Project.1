## Section 0: Purpose: 
## In this code we will create several list-like metadata each consisting of 3 layers and representing 1 feature. 
## For instance, survival list comparing the os and pfs features of lung, melanoma, and ovarian 

## Section 1: Environment ====================================================================================================
CODE <- "~/MRes_project_1/Codes/4_gistic_analysis/"
SAVE <- "~/MRes_project_1/0_RESULTS/"

## Library packages 
library(dplyr)
library(ggplot2)
library(gridExtra)
library(scatterplot3d)




## Section 2: Overall survival and Progression Free Survival =================================================================
## 2.1: Lung -----------------------------------------------------------------------------------------------------------------
tcga_lung <- read.table("~/MRes_project_1/GISTIC2/output/tcga_lung/lung_os.txt")
hh_lung <- read.table("~/MRes_project_1/GISTIC2/output/hh_lung/broad_l_os_surv_1.txt")

lung <- merge(tcga_lung, hh_lung, by="row.names")
rownames(lung) <- lung$Row.names
lung$Row.names <- NULL
lung <- merge(lung, cbioportal_lung, by="row.names")

p1 <- ggplot(lung, aes(x=p_values.x, y=p_values.y, label=Row.names)) + 
  geom_point(size=2, shape=20) + 
  geom_text(aes(label = Row.names), hjust=0, vjust=0, check_overlap = TRUE) + 
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "blue") +  # Vertical line at x = 0.05
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "blue") + 
  geom_smooth(method = "lm", se = FALSE, color = "red") + 
  ggtitle("Lung overall survival p_values") +
  xlab("TCGA lung") + 
  ylab("HH lung") 

## 2.2: Ovarian --------------------------------------------------------------------------------------------------------------
tcga_ova <- read.table("~/MRes_project_1/GISTIC2/output/tcga_ova/ova_os.txt")
hh_ova <- read.table("~/MRes_project_1/GISTIC2/output/hh_ova/broad_o_os_surv_1.txt")

ova <- merge(tcga_ova, hh_ova, by="row.names")

p3 <- ggplot(ova, aes(x=p_values.x, y=p_values.y, label=Row.names)) + 
  geom_point(size=2, shape=20) + 
  geom_text(aes(label = Row.names), hjust=0, vjust=0, check_overlap = TRUE) + 
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "blue") +  # Vertical line at x = 0.05
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "blue") + 
  geom_smooth(method = "lm", se = FALSE, color = "red") + 
  ggtitle("Ovarian overall survival p_values") +
  xlab("TCGA ovarian") + 
  ylab("HH ovarian") 

## 2.3: melanoma -------------------------------------------------------------------------------------------------------------
tcga_mel <- read.table("~/MRes_project_1/GISTIC2/output/tcga_mel/mel_os.txt")
cbioportal_mel <- read.table("~/MRes_project_1/GISTIC2/output/cbioportal_mel/broad_l_os_surv_1.txt")

mel <- merge(tcga_mel, cbioportal_mel, by="row.names")
p5 <- ggplot(mel, aes(x=p_values.x, y=p_values.y, label=Row.names)) + 
  geom_point(size=2, shape=20) + 
  geom_text(aes(label = Row.names), hjust=0, vjust=0, check_overlap = TRUE) + 
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "blue") +  # Vertical line at x = 0.05
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "blue") + 
  geom_smooth(method = "lm", se = FALSE, color = "red") + 
  ggtitle("Melanoma overall survival p_values") +
  xlab("TCGA melanoma") + 
  ylab("HH melanoma")

## 2.4: Lung -----------------------------------------------------------------------------------------------------------------
tcga_lung <- read.table("~/MRes_project_1/GISTIC2/output/tcga_lung/lung_pfi.txt")
hh_lung <- read.table("~/MRes_project_1/GISTIC2/output/hh_lung/broad_l_pfs_surv_1.txt")
cbioportal_lung <- read.table("~/MRes_project_1/GISTIC2/output/cbioportal_lung/broad_l_pfs_surv_1.txt")

lung <- merge(tcga_lung, hh_lung, by="row.names")
rownames(lung) <- lung$Row.names
lung$Row.names <- NULL
lung <- merge(lung, cbioportal_lung, by="row.names")

p2 <- ggplot(lung, aes(x=p_values.x, y=p_values.y, label=Row.names)) + 
  geom_point(size=2, shape=20) + 
  geom_text(aes(label = Row.names), hjust=0, vjust=0, check_overlap = TRUE) + 
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "blue") +  # Vertical line at x = 0.05
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "blue") + 
  geom_smooth(method = "lm", se = FALSE, color = "red") + 
  ggtitle("Lung progression free survival p_values") +
  xlab("TCGA ovarian") + 
  ylab("HH ovarian") 

p7 <- ggplot(lung, aes(x=p_values, y=p_values.y, label=Row.names)) + 
  geom_point(size=2, shape=20) + 
  geom_text(aes(label = Row.names), hjust=0, vjust=0, check_overlap = TRUE) + 
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "blue") +  # Vertical line at x = 0.05
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "blue") + 
  geom_smooth(method = "lm", se = FALSE, color = "red") + 
  ggtitle("Lung progression free survival p_values") +
  xlab("cBioPortal ovarian") + 
  ylab("HH ovarian") 

## 2.5: Ovarian --------------------------------------------------------------------------------------------------------------
tcga_ova <- read.table("~/MRes_project_1/GISTIC2/output/tcga_ova/ova_pfi.txt")
hh_ova <- read.table("~/MRes_project_1/GISTIC2/output/hh_ova/broad_o_pfs_surv_1.txt")

ova <- merge(tcga_ova, hh_ova, by="row.names")

p4 <- ggplot(ova, aes(x=p_values.x, y=p_values.y, label=Row.names)) + 
  geom_point(size=2, shape=20) + 
  geom_text(aes(label = Row.names), hjust=0, vjust=0, check_overlap = TRUE) + 
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "blue") +  # Vertical line at x = 0.05
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "blue") + 
  geom_smooth(method = "lm", se = FALSE, color = "red") + 
  ggtitle("Ovarian progression free survival p_values") +
  xlab("TCGA ovarian") + 
  ylab("HH ovarian") 

## 2.6: melanoma -------------------------------------------------------------------------------------------------------------
tcga_mel <- read.table("~/MRes_project_1/GISTIC2/output/tcga_mel/mel_pfi.txt")
cbioportal_mel <- read.table("~/MRes_project_1/GISTIC2/output/cbioportal_mel/broad_l_pfs_surv_1.txt")

mel <- merge(tcga_mel, cbioportal_mel, by="row.names")
p6 <- ggplot(mel, aes(x=p_values.x, y=p_values.y)) + 
  geom_point(size=2, shape=20) + 
  geom_text(aes(label = Row.names), hjust=0, vjust=0, check_overlap = TRUE) + 
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "blue") +  # Vertical line at x = 0.05
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "blue") + 
  geom_smooth(method = "lm", se = FALSE, color = "red") + 
  ggtitle("Melanoma progression free survival p_values") +
  xlab("TCGA melanoma") + 
  ylab("HH melanoma")


## 2.7: plotting -------------------------------------------------------------------------------------------------------------
pdf(paste0(SAVE, "plots/survival.pdf"), width = 8.3, height = 11.7)
grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 2, nrow = 3)
grid.arrange(p7, ncol = 2, nrow = 3)
dev.off()


## Section 3: treatment response =============================================================================================
## 3.1 Therapy Response ------------------------------------------------------------------------------------------------------
## 3.1.1 Lung 
tcga_lung <- read.table("~/MRes_project_1/GISTIC2/output/tcga_lung/lung_therapy_resp.txt")
hh_lung <- read.table("~/MRes_project_1/GISTIC2/output/hh_lung/broad_l_best_response_1.txt")
cbioportal_lung <- read.table("~/MRes_project_1/GISTIC2/output/cbioportal_lung/broad_l_clinical_benefit_1.txt")

lung <- merge(tcga_lung, hh_lung, by="row.names")
lung_2 <- merge(tcga_lung, cbioportal_lung, by="row.names")

p8 <- ggplot(lung, aes(x = p_values.x, y = p_values.y)) + 
  geom_point(size = 2, shape = 20) + 
  geom_text(aes(label = Row.names), hjust = 0, vjust = 0, check_overlap = TRUE) +  # Include label aesthetic here
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "blue") +  # Vertical line at x = 0.05
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "blue") +  # Horizontal line at y = 0.05
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add linear model fit line
  ggtitle("Lung therapy response") +
  xlab("TCGA Lung") + 
  ylab("HH Lung")

p9 <- ggplot(lung_2, aes(x = p_values.x, y = p_values.y)) + 
  geom_point(size = 2, shape = 20) + 
  geom_text(aes(label = Row.names), hjust = 0, vjust = 0, check_overlap = TRUE) +  # Include label aesthetic here
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "blue") +  # Vertical line at x = 0.05
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "blue") +  # Horizontal line at y = 0.05
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add linear model fit line
  ggtitle("Lung therapy response p_values") +
  xlab("TCGA Lung") + 
  ylab("cBioPortal Lung")

## 3.1.2 Ovarian 
tcga_ova <- read.table("~/MRes_project_1/GISTIC2/output/tcga_ova/ova_therapy_resp.txt")
hh_ova <- read.table("~/MRes_project_1/GISTIC2/output/hh_ova/broad_o_primary_chemo_outcome_1.txt")

ova <- merge(tcga_ova, hh_ova, by="row.names")

p10 <- ggplot(ova, aes(x=p_values.x, y=p_values.y, label=Row.names)) + 
  geom_point(size=2, shape=20) + 
  geom_text(aes(label = Row.names), hjust=0, vjust=0, check_overlap = TRUE) + 
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "blue") +  # Vertical line at x = 0.05
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "blue") + 
  geom_smooth(method = "lm", se = FALSE, color = "red") + 
  ggtitle("Ovarian therapy response p_values") +
  xlab("TCGA ovarian") + 
  ylab("HH ovarian") 

## 3.2 Stage -----------------------------------------------------------------------------------------------------------------
## 3.2.1 Lung 
tcga_lung <- read.table("~/MRes_project_1/GISTIC2/output/tcga_lung/lung_stage.txt")
hh_lung <- read.table("~/MRes_project_1/GISTIC2/output/hh_lung/broad_l_stage_1.txt")

lung <- merge(tcga_lung, hh_lung, by="row.names")
lung_2 <- merge(tcga_lung, cbioportal_lung, by="row.names")

p11 <- ggplot(lung, aes(x = p_values.x, y = p_values.y)) + 
  geom_point(size = 2, shape = 20) + 
  geom_text(aes(label = Row.names), hjust = 0, vjust = 0, check_overlap = TRUE) +  # Include label aesthetic here
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "blue") +  # Vertical line at x = 0.05
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "blue") +  # Horizontal line at y = 0.05
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add linear model fit line
  ggtitle("Lung stage p_values") +
  xlab("TCGA Lung") + 
  ylab("HH Lung")

## 3.2.2 Ova 
tcga_ova <- read.table("~/MRes_project_1/GISTIC2/output/tcga_ova/ova_stage.txt")
hh_ova <- read.table("~/MRes_project_1/GISTIC2/output/hh_ova/broad_o_stage_1.txt")

ova <- merge(tcga_ova, hh_ova, by="row.names")

p12 <- ggplot(ova, aes(x=p_values.x, y=p_values.y, label=Row.names)) + 
  geom_point(size=2, shape=20) + 
  geom_text(aes(label = Row.names), hjust=0, vjust=0, check_overlap = TRUE) + 
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "blue") +  # Vertical line at x = 0.05
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "blue") + 
  geom_smooth(method = "lm", se = FALSE, color = "red") + 
  ggtitle("Ovarian stage p_values") +
  xlab("TCGA ovarian") + 
  ylab("HH ovarian") 

## 3.2.3 Melanoma 
tcga_mel <- read.table("~/MRes_project_1/GISTIC2/output/tcga_mel/mel_stage.txt")
cbioportal_mel <- read.table("~/MRes_project_1/GISTIC2/output/cbioportal_mel/broad_l_stage_1.txt")

mel <- merge(tcga_mel, cbioportal_mel, by="row.names")
p13 <- ggplot(mel, aes(x=p_values.x, y=p_values.y)) + 
  geom_point(size=2, shape=20) + 
  geom_text(aes(label = Row.names), hjust=0, vjust=0, check_overlap = TRUE) + 
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "blue") +  # Vertical line at x = 0.05
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "blue") + 
  geom_smooth(method = "lm", se = FALSE, color = "red") + 
  ggtitle("Melanoma stage p_values") +
  xlab("TCGA melanoma") + 
  ylab("HH melanoma")

## 3.3 Age -------------------------------------------------------------------------------------------------------------------
## 3.2.1 Lung 
tcga_lung <- read.table("~/MRes_project_1/GISTIC2/output/tcga_lung/lung_age.txt")
hh_lung <- read.table("~/MRes_project_1/GISTIC2/output/hh_lung/broad_l_age_1.txt")

lung <- merge(tcga_lung, hh_lung, by="row.names")
lung_2 <- merge(tcga_lung, cbioportal_lung, by="row.names")

p14 <- ggplot(lung, aes(x = p_values.x, y = p_values.y)) + 
  geom_point(size = 2, shape = 20) + 
  geom_text(aes(label = Row.names), hjust = 0, vjust = 0, check_overlap = TRUE) +  # Include label aesthetic here
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "blue") +  # Vertical line at x = 0.05
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "blue") +  # Horizontal line at y = 0.05
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add linear model fit line
  ggtitle("Lung age p_values") +
  xlab("TCGA Lung") + 
  ylab("HH Lung")

## 3.2.2 Ova 
tcga_ova <- read.table("~/MRes_project_1/GISTIC2/output/tcga_ova/ova_age.txt")
hh_ova <- read.table("~/MRes_project_1/GISTIC2/output/hh_ova/broad_o_age_1.txt")

ova <- merge(tcga_ova, hh_ova, by="row.names")

p15 <- ggplot(ova, aes(x=p_values.x, y=p_values.y, label=Row.names)) + 
  geom_point(size=2, shape=20) + 
  geom_text(aes(label = Row.names), hjust=0, vjust=0, check_overlap = TRUE) + 
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "blue") +  # Vertical line at x = 0.05
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "blue") + 
  geom_smooth(method = "lm", se = FALSE, color = "red") + 
  ggtitle("Ovarian age p_values") +
  xlab("TCGA ovarian") + 
  ylab("HH ovarian") 

## 3.2.3 Melanoma 
tcga_mel <- read.table("~/MRes_project_1/GISTIC2/output/tcga_mel/mel_age.txt")
cbioportal_mel <- read.table("~/MRes_project_1/GISTIC2/output/cbioportal_mel/broad_l_age_1.txt")

mel <- merge(tcga_mel, cbioportal_mel, by="row.names")
p16 <- ggplot(mel, aes(x=p_values.x, y=p_values.y)) + 
  geom_point(size=2, shape=20) + 
  geom_text(aes(label = Row.names), hjust=0, vjust=0, check_overlap = TRUE) + 
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "blue") +  # Vertical line at x = 0.05
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "blue") + 
  geom_smooth(method = "lm", se = FALSE, color = "red") + 
  ggtitle("Melanoma age p_values") +
  xlab("TCGA melanoma") + 
  ylab("HH melanoma")

## 3.2 plotting --------------------------------------------------------------------------------------------------------------
pdf(paste0(SAVE, "plots/therapy_response.pdf"), width = 8.3, height = 11.7)
grid.arrange(p8, p9, p10, p11, p12, p13, ncol = 2, nrow = 3)
grid.arrange(p14, p15, p16, ncol = 2, nrow = 3)
dev.off()


## Section 4: immune features ================================================================================================
## 4.0 read in rds 
tcga_lung_rds <- readRDS("~/MRes_project_1/GISTIC2/output/tcga_lung/result_imm_lm_b_lung.rds")
tcga_ova_rds <- readRDS("~/MRes_project_1/GISTIC2/output/tcga_ova/result_imm_lm_b_ova.rds")
tcga_mel_rds <- readRDS("~/MRes_project_1/GISTIC2/output/tcga_mel/result_imm_lm_b_mel.rds")

## 4.1 TMB -------------------------------------------------------------------------------------------------------------------
## 4.1.2 TMB
## lung
tcga_lung <- tcga_lung_rds$`TMB (nonsynonymous)`
#hh_lung <- read.table("~/MRes_project_1/GISTIC2/output/hh_lung/")
cbioportal_lung <- read.table("~/MRes_project_1/GISTIC2/output/cbioportal_lung/broad_l_TMB_1.txt")

lung <- merge(tcga_lung, cbioportal_lung, by="row.names")

p17 <- ggplot(lung, aes(x = p_values.x, y = p_values.y)) + 
  geom_point(size = 2, shape = 20) + 
  geom_text(aes(label = Row.names), hjust = 0, vjust = 0, check_overlap = TRUE) +  # Include label aesthetic here
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "blue") +  # Vertical line at x = 0.05
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "blue") +  # Horizontal line at y = 0.05
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add linear model fit line
  ggtitle("Lung TMB p_values") +
  xlab("TCGA Lung") + 
  ylab("cBioPortal Lung") 

## mel
tcga_mel <- tcga_mel_rds$`TMB (nonsynonymous)`
cbioportal_mel <- read.table("~/MRes_project_1/GISTIC2/output/cbioportal_mel/broad_m_TMB_1.txt")

mel <- merge(tcga_mel, cbioportal_mel, by="row.names")
p18 <- ggplot(mel, aes(x=p_values.x, y=p_values.y)) + 
  geom_point(size=2, shape=20) + 
  geom_text(aes(label = Row.names), hjust=0, vjust=0, check_overlap = TRUE) + 
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "blue") +  # Vertical line at x = 0.05
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "blue") + 
  geom_smooth(method = "lm", se = FALSE, color = "red") + 
  ggtitle("Melanoma TMB p_value") +
  xlab("TCGA melanoma") + 
  ylab("cBioPortal melanoma")







## 3.2.2 Ova 
tcga_ova <- read.table("~/MRes_project_1/GISTIC2/output/tcga_ova/ova_age.txt")
hh_ova <- read.table("~/MRes_project_1/GISTIC2/output/hh_ova/broad_o_age_1.txt")

ova <- merge(tcga_ova, hh_ova, by="row.names")

p15 <- ggplot(ova, aes(x=p_values.x, y=p_values.y, label=Row.names)) + 
  geom_point(size=2, shape=20) + 
  geom_text(aes(label = Row.names), hjust=0, vjust=0, check_overlap = TRUE) + 
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "blue") +  # Vertical line at x = 0.05
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "blue") + 
  geom_smooth(method = "lm", se = FALSE, color = "red") + 
  ggtitle("Ovarian age p_values") +
  xlab("TCGA ovarian") + 
  ylab("HH ovarian") 

## 3.2.3 Melanoma 
tcga_mel <- read.table("~/MRes_project_1/GISTIC2/output/tcga_mel/mel_age.txt")
cbioportal_mel <- read.table("~/MRes_project_1/GISTIC2/output/cbioportal_mel/broad_l_age_1.txt")

mel <- merge(tcga_mel, cbioportal_mel, by="row.names")
p16 <- ggplot(mel, aes(x=p_values.x, y=p_values.y)) + 
  geom_point(size=2, shape=20) + 
  geom_text(aes(label = Row.names), hjust=0, vjust=0, check_overlap = TRUE) + 
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "blue") +  # Vertical line at x = 0.05
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "blue") + 
  geom_smooth(method = "lm", se = FALSE, color = "red") + 
  ggtitle("Melanoma age p_values") +
  xlab("TCGA melanoma") + 
  ylab("HH melanoma")






















