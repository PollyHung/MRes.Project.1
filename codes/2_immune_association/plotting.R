# Assuming 'result_imm_lm_b_ova' is your list of data frames

# Initialize an empty data frame to store the results
final_df <- data.frame(cancer_type=character(),
                       chr_arm=character(),
                       immune_phenotypes=character(),
                       estimate=numeric(),
                       fdr=numeric(),
                       stringsAsFactors=FALSE)

# Iterate through each item in the list
for (phenotype in names(result_imm_lm_b_lung)) {
  # Extract the current data frame
  current_df <- result_imm_lm_b_lung[[phenotype]]
  
  # Check if the 'chr_arm' column is not present and create it based on rownames
  if (!'chr_arm' %in% colnames(current_df)) {
    current_df$chr_arm <- rownames(current_df)
  }
  
  # Add the phenotype name and cancer type to the current data frame
  current_df$immune_phenotypes <- phenotype
  current_df$cancer_type <- 'lung'
  
  # Select only the necessary columns
  current_df <- current_df[, c('cancer_type', 'chr_arm', 'immune_phenotypes', 'estimate', 'FDR')]
  
  # Rename the 'FDR' column to 'fdr'
  colnames(current_df)[colnames(current_df) == 'FDR'] <- 'fdr'
  
  # Bind the current data frame to the final data frame
  final_df <- rbind(final_df, current_df)
}

final_df$immune_phenotypes[which(final_df$immune_phenotypes=="tls_rna_seq$TLS_mean")] <- "tls"

write.csv(final_df, "~/MRes_project_1/codes/5_plots/files/tcgaDotPlot.csv", row.names = FALSE, quote=FALSE)


total <- total %>% dplyr::filter(fdr<0.05) %>% dplyr::filter(!immune_phenotypes %in% c("Fraction Genome Altered", 
                                                                                "cd8_t.cell", 
                                                                                "TMB (nonsynonymous)")) 
total$color <- ifelse(total$estimate > 0, "red", "blue") # if positive correlation, color red, else blue
total$estimate <- abs(total$estimate) # we need to make all cor value positive so that the plot could create size 
total$alpha <- ifelse(total$fdr < 0.05, 1-total$fdr, 0) ## # if FDR is less than 0.05, we will create a new value (0.05-p-value) and store in the alpha column, else, store 0
total$alpha_norm <- (total$alpha-min(total$alpha))/(max(total$alpha)-min(total$alpha))
total$color <- ifelse(total$alpha == 0, "white", total$color) # if the alpha is set to 0 (FDR > 0.05), then the point will be colored white 
total$label <- paste0(gsub(" cancer", "", total$cancer_type), ".", total$chr_arm)

p <- ggplot(total, aes(x = total$label, y = immune_phenotypes, size = estimate, color = color, alpha = alpha_norm)) +
  geom_point() +
  scale_size_continuous(range = c(1, 6)) + # Adjust the size range as needed
  scale_alpha_continuous(limits = c(0.95, 1)) + # Transparency from 0 (fully transparent) to 1 (opaque)
  scale_color_identity() + # Use the colors defined in the dataframe
  labs(x = "Chromosome Arm", y = "Immune Phenotype") +
  ggtitle(paste0("Level copy number aberration at different immune phenotypes for various cancers")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
ggsave("~/plot_1.png", p, width=11, height=3.5, units = "in")

p
