library(dplyr)
library(fuzzyjoin)


cytoband <- read.table("~/MRes_project_1/codes/6_single_cell/reference/cytoBand.txt", sep = "\t")
colnames(cytoband) <- c("chr", "start", "end", "arm", "unknown")
directory <- "~/MRes_project_1/docs/SC_SEQ/lung/raw/lambrechtslab/inferCNV/72_hr_runtime/"
setwd(directory)
filename <- "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_genes.dat"
#samples <- list.files()

temp1 <- read.table(paste(directory, filename, sep = "/"), header = TRUE)
#temp1$cell_name <- gsub("\\..*","", temp1$cell_group_name)
temp1$sample <- gsub("_.*", "", temp1$cell_group_name)
#temp1$cell_name <- sub(".*_(.*?)\\..*", "\\1", temp1$cell_group_name)

# Make sure that 'start' and 'end' columns are numeric
cytoband <- cytoband %>% mutate(start = as.numeric(start), end = as.numeric(end))
temp1 <- temp1 %>% mutate(start = as.numeric(start), end = as.numeric(end))

# Use fuzzyjoin to join the dataframes based on the gene position overlapping with cytoband regions
mapped_genes <- fuzzyjoin::genome_join(temp1, cytoband, by = c("chr", "start", "end"))

# Filter to get only the rows where gene start and end are within the arm start and end
mapped_genes <- mapped_genes %>% dplyr::filter(start.x >= start.y & end.x <= end.y)

# Select the relevant columns 
mapped_genes <- mapped_genes %>%
  dplyr::select(sample, gene, state, chr=chr.x, start=start.x, end=end.x, arm)

# add a new column called chr_arm 
mapped_genes$chr_arm <- paste0(mapped_genes$chr, str_extract(mapped_genes$arm, pattern = "p|q"))

# Calculate the mean 
#summarise_temp1 <- mapped_genes %>% group_by(sample, chr_arm) %>% summarise(
    #median = quantile(state, probs = 0.5), 
    #mean = mean(state),
    #pct_1 = (count(state==1)/length(state))*100,  
    #pct_2 = (count(state==2)/length(state))*100, 
    #pct_3 = (count(state==3)/length(state))*100,  
    #pct_4 = (count(state==4)/length(state))*100,  
    #pct_5 = (count(state==5)/length(state))*100,  
    #pct_6 = (count(state==6)/length(state))*100,  
  #)

#summarise_long <- summarise_temp1 %>%
  #dplyr::select(sample, chr_arm, starts_with("pct_")) %>%
  #pivot_longer(cols = starts_with("pct_"), names_to = "pct_state", values_to = "value") %>%
  #mutate(pct_state = factor(pct_state, levels = paste0("pct_", 1:6)))

# Plot the stacked bar plot
#plot.OvaStackBarAneuStatus_bySample <- ggplot(summarise_long, aes(x = chr_arm, y = value, fill = pct_state)) +
  #geom_bar(stat = "identity") +
  #facet_wrap(~sample, scales = "free_x") +
  #theme_bw() +
  #labs(x = "Chromosome Arm", y = "Percentage", fill = "Aneuploidy State") +
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))
#ggsave("~/MRes_project_1/docs/SC_SEQ/plots/LungStackBarAneuStatus_bySample.png", plot.OvaStackBarAneuStatus_bySample, 
       #width=20, height=12, units = "in")

#summarise_long <- summarise_temp1 %>%
  #dplyr::select(sample, chr_arm, starts_with("pct_")) %>%
  #pivot_longer(cols = starts_with("pct_"), names_to = "pct_state", values_to = "value") %>%
  #mutate(pct_state = factor(pct_state, levels = paste0("pct_", 1:6)))

# Plot the stacked bar plot with samples on the x-axis and split by chr_arm
#plot.OvaStackBarAneuStatus_byArm <- ggplot(summarise_long, aes(x = sample, y = value, fill = pct_state)) +
  #geom_bar(stat = "identity", position = "stack") +
  #facet_wrap(~chr_arm, scales = "free_x", nrow = 8) + # Adjust nrow for better layout if needed
  #theme_bw() +
  #labs(x = "Sample", y = "Percentage", fill = "Aneuploidy State") +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1),
        #strip.text.x = element_text(size = 8)) # Adjust text size for readability
#ggsave("~/MRes_project_1/docs/SC_SEQ/plots/LungStackBarAneuStatus_byArm.png", plot.OvaStackBarAneuStatus_byArm, 
       #width=10, height=16, units = "in")


summarise_temp2 <- mapped_genes %>% 
  #dplyr::filter(chr_arm %in% c("chr2q", "chr10p")) %>% 
  group_by(sample, chr_arm) %>% summarise(aneu_score = mean(state)) #%>% 
  #pivot_wider(names_from = chr_arm, values_from = aneu_score) %>% as.data.frame
#rownames(summarise_temp2) <- summarise_temp2$sample
#summarise_temp2$sample <- NULL
summarise_temp2$aneu_score <- summarise_temp2$aneu_score - 3

write.table(summarise_temp2, "~/MRes_project_1/docs/SC_SEQ/lung/processed/lambrechtslab/InferCNV/armChrCNA_longer.txt", 
            sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)





















