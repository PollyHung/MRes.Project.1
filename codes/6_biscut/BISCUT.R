library(dplyr)
library(magrittr)

## read in biscut ------------------------------------------------------------------------
#biscut <- readxl::read_xlsx("~/MRes_project_1/codes/7_biscut/BISCUT.xlsx", sheet = "pancan_2q_only")
biscut <- read.table("~/MRes_project_1/codes/7_biscut/BISCUT_ovarian/results_2023_12_05_0.95/all_BISCUT_results.txt", 
                     header = TRUE)
biscut <- biscut %>% dplyr::filter(arm == "2q")
biscut_genes <- biscut$Gene


## read in tcga data ------------------------------------------------------------------------
tcga_rna <- read.table("~/MRes_project_1/docs/00_mitéra/rna.seq/tcgaOvaRNA", sep = "\t", header= TRUE, row.names = 1)
tcga_arm <- read.table("~/MRes_project_1/docs/GISTIC2/tcga/ova/broad_values_by_arm.txt",  sep = "\t", header = TRUE, row.names = 1)
tcga_cna <- read.table("~/MRes_project_1/docs/GISTIC2/tcga/ova/all_data_by_genes.txt", sep = "\t", header = TRUE, row.names = 1, 
                       colClasses = c(Gene.ID = "NULL")) #, Cytoband = "NULL"
clinical_mtx <- read.table("~/MRes_project_1/docs/00_mitéra/clinical/raw/tcgaOvaPheno.txt", sep = "\t", header = TRUE)
clinical_mtx <- clinical_mtx %>% 
  dplyr::select(sampleID, age_at_initial_pathologic_diagnosis, gender, clinical_stage) %>% 
  dplyr::mutate(
    sampleID = gsub("-", ".", sampleID), # change - to . in the sample ID to match
    age_at_initial_pathologic_diagnosis = as.numeric(age_at_initial_pathologic_diagnosis), # as numeric age
    gender = ifelse(gender == "FEMALE", 1, NA), # female = 1, male = 0
    clinical_stage = case_when(
      clinical_stage %in% c("Stage IA", "Stage IB", "Stage IC") ~ 1,
      clinical_stage %in% c("Stage IIA", "Stage IIB", "Stage IIC") ~ 2,
      clinical_stage %in% c("Stage IIIA", "Stage IIIB", "Stage IIIC") ~ 3,
      clinical_stage == "Stage IV" ~ 4,
      TRUE ~ NA)) %>% 
  dplyr::rename(age = age_at_initial_pathologic_diagnosis, stage = clinical_stage)
rownames(clinical_mtx) <- clinical_mtx$sampleID
clinical_mtx$sampleID <- NULL


## sample and gene names? ------------------------------------------------------------------------
# only get chromosome 2q 
tcga_cna <- tcga_cna[grep("2q", tcga_cna$Cytoband), ]
tcga_cna$Cytoband <- NULL

# intersecting gene names and sample names 
sample_names <- intersect(colnames(tcga_cna), colnames(tcga_rna))
genes <- intersect(rownames(tcga_cna), rownames(tcga_rna)) #intersect(biscut_genes)



## data alignment ------------------------------------------------------------------------------------
tcga_cna <- tcga_cna[genes, sample_names] #genes
tcga_rna <- tcga_rna[genes, sample_names]
tcga_arm <- tcga_arm[, sample_names]
clinical_mtx <- clinical_mtx[sample_names, ]


## correlational test --------------------------------------------------------------------------------
#result_lm <- matrix(NA, nrow = nrow(tcga_rna), ncol=5) %>% as.data.frame
#colnames(result_lm) <- c("estimate", "st.error", "p.val", "adj.r.sq", "adj.p.val")
#rownames(result_lm) <- rownames(tcga_rna)

#model <- lm(rna~arm+clinical_mtx$age+clinical_mtx$gender+clinical_mtx$stage) %>% summary()
#result_lm[i, "estimate"] <- model$coefficients[2, "Estimate"]
#result_lm[i, "st.error"] <- model$coefficients[2, "Std. Error"]
#result_lm[i, "p.val"] <- model$coefficients[2, "Pr(>|t|)"]
#result_lm[i, "adj.r.sq"] <- model$adj.r.squared
#result_lm$adj.p.val <- p.adjust(result_lm$p.val)
#result_lm <- result_lm %>% dplyr::arrange(adj.p.val)

result_cor <- matrix(NA, nrow = nrow(tcga_cna), ncol=3) %>% as.data.frame
colnames(result_cor) <- c("estimate", "p.val", "adj.p.val")
rownames(result_cor) <- rownames(tcga_cna)

for(i in rownames(tcga_rna)){
  rna <- unlist(tcga_rna[i, ])
  cna <- unlist(tcga_cna[i, ])
  arm <- unlist(tcga_arm["2q", ])
  
  model2 <- cor.test(x = rna, y = cna, method = "pearson")
  result_cor[i, "estimate"] <- model2$estimate
  result_cor[i, "p.val"] <- model2$p.value
}

result_cor$adj.p.val <- p.adjust(result_cor$p.val)
result_cor <- result_cor %>% dplyr::arrange(adj.p.val)
result_cor_filt <- result_cor %>% dplyr::arrange(adj.p.val) %>% dplyr::filter(adj.p.val<0.05)

cor_geneByArm <- result_cor
cor_geneByArm.filt <- result_cor_filt

cor_geneByCna <- result_cor
cor_geneByCna.filt <- result_cor_filt





## crispr screening --------------------------------------------------------------------------------
# path 1 
file_path <- "/rds/general/user/ph323/home/MRes_project_1/codes/8_crispr_screen/chr2q.txt"
data <- stream_in(file(file_path))
data$genes <- sapply(strsplit(data$genetargets, "::"), `[`, 1)
crispr_screen <- data %>% dplyr::select(cas, screentype, chr, start, end, score, hit, 
                                        cellline, pubmed, symbol, effect, rc_initial, rc_final, 
                                        log2fc, sequence, condition)
crispr_screen_biscut <- crispr_screen %>% dplyr::filter(symbol %in% biscut_genes & hit == "true")
crispr_screen_tcga <- crispr_screen %>% dplyr::filter(symbol %in% rownames(result_cor) & hit == "true")

# path 2 
biogrid <- readxl::read_xlsx("~/MRes_project_1/codes/8_crispr_screen/CRISPR/biogrid_total.xlsx", sheet = 1)
matching <- readxl::read_xlsx("~/MRes_project_1/codes/8_crispr_screen/CRISPR/biogrid_total.xlsx", sheet = 2)
matching <- matching %>% dplyr::select(`#SCREEN_ID`, SOURCE_ID)

biogrid_new <- biogrid %>% dplyr::select(OFFICIAL_SYMBOL, contains("-HIT")) %>% 
  dplyr::mutate(across(ends_with("-HIT"), ~case_when(
    .x == "YES" ~ 1,
    .x == "NO" ~ 0,
    .x == "-" ~ 0))) %>% 
  dplyr::mutate(hits = rowSums(across(ends_with("-HIT")))) %>% rowwise() %>% 
  dplyr::mutate(negative_hits = sum(c_across(ends_with("-HIT")) == -1, na.rm = TRUE), 
                positive_hits = sum(c_across(ends_with("-HIT")) == 1, na.rm = TRUE))

biogrid_filt <- biogrid_new %>% dplyr::filter(hits > 0) ## filter genes that are hits 

# biscut in biogrid? 
biscut.biogrid <- intersect(biscut$Gene, biogrid_filt$OFFICIAL_SYMBOL)

# tcga in biogrid? 
tcga.biogrid <- intersect(rownames(cor_geneByArm), biogrid_filt$OFFICIAL_SYMBOL)
tcga_filt.biogrid <- intersect(rownames(cor_geneByArm.filt), biogrid_filt$OFFICIAL_SYMBOL)

tcga.biogrid.2 <- intersect(rownames(cor_geneByCna), biogrid_filt$OFFICIAL_SYMBOL)
tcga_filt.biogrid.2 <- intersect(rownames(cor_geneByCna.filt), biogrid_filt$OFFICIAL_SYMBOL)


save.image("/rds/general/user/ph323/ephemeral/biscut/biscut.RData")










