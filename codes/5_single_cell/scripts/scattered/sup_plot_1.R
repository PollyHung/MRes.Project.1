## datasets 
tcga_lung <- read.table("~/MRes_project_1/docs/00_mitéra/clinical/processed/tcga_lung_mtx2.txt", header = TRUE, row.names = 1)
tcga_ovarian <- read.csv("~/MRes_project_1/docs/00_mitéra/clinical/processed/tcga_ova_mtx.csv")
hh_ova <- read.table("~/MRes_project_1/codes/5_plots/file_final/hh_ova_mtx.txt", sep = "\t", header = T, row.names = 1)
hh_lung <- read.csv("~/MRes_project_1/docs/HH_lung/clinical_data/HH_lung_clinical.csv")
cbp_lung <- read.csv("~/MRes_project_1/docs/00_mitéra/clinical/processed/cbp_lung.csv")
sc_seq_ova <- readxl::read_xlsx("~/MRes_project_1/docs/00_mitéra/clinical/processed/msk_scseq_ova.xlsx", sheet = 2)

## plot setting 
color <- brewer.pal(6, "Set1")
names(color) <- c("tcga_lung", "tcga_ova", "hh_ova", "hh_lung", "mskcc_lung", "msk_scseq_ova")

## age 
age <- data.frame(age = c(tcga_ovarian$age, tcga_lung$age, hh_ova$age, hh_lung$age, cbp_lung$age, sc_seq_ova$patient_age_at_diagnosis),
                  dataset = rep(c("tcga_ova", "tcga_lung", "hh_ova", "hh_lung", "mskcc_lung", "msk_scseq_ova"), 
                                times = c(nrow(tcga_ovarian), nrow(tcga_lung), nrow(hh_ova), nrow(hh_lung), nrow(cbp_lung), nrow(sc_seq_ova))))
p <- ggplot(age, aes(x = dataset, y = age, fill = dataset)) + 
  geom_violin(trim = FALSE, alpha = 0.5, size = 0.2) +
  stat_summary(fun = median, geom = "crossbar", width = 0.4, color = "red") +
  scale_fill_manual(values = color) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_text(size = 6), 
        axis.title.y = element_text(size = 6), 
        text = element_text(size = 6), 
        title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 6), 
        legend.key.size = unit(0.3, units = "cm")) + 
  #guides(fill = FALSE) +  
  ggtitle("Patient Age Distribution") + 
  labs(x = "Dataset", y = "Age Distribution")

ggsave("~/MRes_project_1/codes/5_plots/plots/sup_plot_1/age_distribution.png", p, width = 2, height = 2, units = "in", dpi = 600)



## patient stage 
stage <- data.frame(table(tcga_ovarian$stage), 
                    table(tcga_lung$stage.code), 
                    table(hh_ova$stage), 
                    table(sc_seq_ova$Stage))
stage <- stage %>% dplyr::select(Var1, Freq, Freq.1, Freq.2, Freq.3) %>% 
  dplyr::rename(stage = Var1, tcga_ovarian = Freq, tcga_lung = Freq.1, hh_ova = Freq.2, sc_seq_ova = Freq.3)
stage$sc_seq_ova[1] <- 0
stage$sc_seq_ova[2] <- 0


# Convert data to long format
data_long <- stage %>% 
  gather(key = "dataset", value = "number_of_patients", -stage) %>% 
  group_by(dataset) %>% 
  mutate(percentage = number_of_patients / sum(number_of_patients) * 100)

# Plot
p1 <- ggplot(data_long, aes(fill = as.factor(stage), y = percentage, x = dataset)) + 
  geom_bar(position = "fill", stat = "identity") + 
  labs(x = "Dataset", y = "Percentage of Patients (%)", fill = "Stage") +
  scale_y_continuous(labels = scales::percent_format()) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_text(size = 6), 
        axis.title.y = element_text(size = 6), 
        text = element_text(size = 6), 
        title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 6), 
        legend.key.size = unit(0.3, units = "cm")) + 
  ggtitle("Patient Stage Distribution")
ggsave("~/MRes_project_1/codes/5_plots/plots/sup_plot_1/stage_distribution.png", p1,  width = 1.5, height = 2, units = "in", dpi = 600)




## TMB distribution 
sc_seq_ova <- read.table("~/MRes_project_1/docs/00_mitéra/clinical/processed/msk_spectrum_tme_2022_clinical_data.tsv", sep = "\t", header = TRUE)
tmb <- data.frame(tmb = c(tcga_ovarian$TMB, tcga_lung$TMB, hh_ova$TMB, cbp_lung$TMB_NONSYNONYMOUS, sc_seq_ova$TMB..nonsynonymous.),
                  dataset = rep(c("tcga_ova", "tcga_lung", "hh_ova", "mskcc_lung", "msk_scseq_ova"), 
                                times = c(nrow(tcga_ovarian), nrow(tcga_lung), nrow(hh_ova), nrow(cbp_lung), nrow(sc_seq_ova))))
tmb$log10tmb <- log10((tmb$tmb + 1))
p2 <- ggplot(tmb, aes(x = log10tmb, fill = dataset)) + 
  geom_density(alpha = 0.5, size = 0.2) + 
  labs(x = "log10 Tumor Mutational Burden (log10(TMB+1))", y = "Density") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_text(size = 6), 
        axis.title.y = element_text(size = 6), 
        text = element_text(size = 6), 
        title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 6), 
        legend.key.size = unit(0.3, units = "cm")) +
  ggtitle("Patient TMB Distribution") + 
  scale_fill_manual(values = color) # Optional: Use a color palette for clarity
ggsave("~/MRes_project_1/codes/5_plots/plots/sup_plot_1/log10TMB_distribution.png", p2,  width = 3, height = 2, units = "in", dpi = 600)




## purity and ploidy 
facets_qc <- readxl::read_xlsx("~/MRes_project_1/codes/5_plots/files/Key Resource Table.xlsx", sheet = 9) %>% as.data.frame()
facets_qc <- facets_qc[1:13]
facets_qc <- facets_qc %>% dplyr::filter(run_type == "purity")


p3 <- ggplot(facets_qc, aes(x = as.numeric(purity), fill = datset)) + 
  geom_density(alpha = 0.5, size = 0.2) + 
  labs(x = "purity by FACETs", y = "Density") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_text(size = 6), 
        axis.title.y = element_text(size = 6), 
        text = element_text(size = 6), 
        title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 6), 
        legend.key.size = unit(0.3, units = "cm")) +
  ggtitle("Purity Distribution") + 
  scale_fill_manual(values = color) # Optional: Use a color palette for clarity

p4 <- ggplot(facets_qc, aes(x = as.numeric(ploidy), fill = datset)) + 
  geom_density(alpha = 0.5, size = 0.2) + 
  labs(x = "ploidy by FACETs", y = "Density") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_text(size = 6), 
        axis.title.y = element_text(size = 6), 
        text = element_text(size = 6), 
        title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 6), 
        legend.key.size = unit(0.3, units = "cm")) +
  ggtitle("Ploidy Distribution") + 
  scale_fill_manual(values = color) # Optional: Use a color palette for clarity

ggsave("~/MRes_project_1/codes/5_plots/plots/sup_plot_1/purity_ploidy.png", p3+p4, width = 5, height = 2, units = "in", dpi = 600)



