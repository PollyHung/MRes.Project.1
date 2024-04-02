temp1 <- read.table("~/MRes_project_1/docs/SC_SEQ/ova/raw/lambrechtslab/InferCNV/result/scrSOL004_alt/HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat", header = TRUE)
temp1$cell_group_name %>% unique

temp1$cell_name <- gsub("\\..*","", temp1$cell_group_name)
summarise_temp1 <- temp1 %>% group_by(cell_name, chr) %>% summarise(temp1_mean = mean(state),
                                                                    temp1_IQR = IQR(state), 
                                                                    temp1_max = max(state), 
                                                                    temp1_min = min(state))
write.table(summarise_temp1, "~/MRes_project_1/docs/SC_SEQ/ova/processed/lambrechtslab/inferCNV/summarise_temp1_alt.txt", 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
