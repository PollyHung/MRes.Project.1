## For GSE180661
cnv <- read.table("~/MRes_project_1/docs/SC_SEQ/ova/GSE180661/CNVbyGISTIC/broad_values_by_arm.txt", sep="\t", header=TRUE, row.names = 1)
colnames(cnv) <- gsub("\\.", "-", colnames(cnv))

map <- read_tsv("~/MRes_project_1/docs/SC_SEQ/ova/GSE180661/clinical_data.tsv") %>% as.data.frame
colnames(map) <- gsub(" ", "_", colnames(map)) %>% tolower 

wgs_update <- map %>% dplyr::filter(gene_panel == "WGS") %>%  dplyr::select(sample_id, patient_display_name)
name_mapping <- setNames(wgs_update$patient_display_name, wgs_update$sample_id)
wgs_cnv <- cnv[, 1:39]
colnames(wgs_cnv) <- name_mapping[colnames(wgs_cnv)]

impact_update <- map %>% dplyr::filter(gene_panel %in% c("IMPACT468", "IMPACT505")) %>% dplyr::select(sample_id, patient_display_name)
name_mapping <- setNames(impact_update$patient_display_name, impact_update$sample_id)
impact_cnv <- cnv[, 40:ncol(cnv)]
colnames(impact_cnv) <- name_mapping[colnames(impact_cnv)]

write.table(wgs_cnv, "~/MRes_project_1/docs/SC_SEQ/ova/GSE180661/CNVbyGISTIC/wgs_broad_values_by_arm.txt", 
            sep="\t", quote=FALSE, row.names = TRUE, col.names = TRUE)
write.table(impact_cnv, "~/MRes_project_1/docs/SC_SEQ/ova/GSE180661/CNVbyGISTIC/impact_broad_values_by_arm.txt", 
            sep="\t", quote=FALSE, row.names = TRUE, col.names = TRUE)

wgs_cnv$id <- rownames(wgs_cnv)
df_long <- melt(wgs_cnv, id.vars = "id")

impact_cnv$id <- rownames(impact_cnv)
df_long_2 <- melt(impact_cnv, id.vars = "id")

pdf("~/MRes_project_1/docs/SC_SEQ/ova/GSE180661/plots/inferCNVboxPlot.pdf", width=8, height = 5)
print(ggplot(df_long, aes(x = factor(id), y = value)) + geom_boxplot() + 
        xlab("chromosome arm") + ylab("copy number variation (0=diploid)") + ylim(c(-3, 3)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
        geom_hline(yintercept=0, color = "red", linewidth=0.5, alpha=0.5) +
        ggtitle("Boxplot by chromosome arm for WGS cohort")) 
print(ggplot(df_long_2, aes(x = factor(id), y = value)) + geom_boxplot() + 
        xlab("chromosome arm") + ylab("copy number variation (0=diploid)") + ylim(c(-3, 3)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        geom_hline(yintercept=0, color = "red", linewidth=0.5, alpha=0.5) +
        ggtitle("Boxplot by chromosome arm for IMPACT cohort"))
dev.off()



## For lambrechtslab
library(dplyr)
library(fuzzyjoin)

for(i in c("ova", "lung")){
  ## read in cytobands 
  cytoband <- read.table("~/MRes_project_1/codes/6_single_cell/reference/cytoBand.txt", sep = "\t")
  colnames(cytoband) <- c("chr", "start", "end", "arm", "unknown")
  cytoband <- cytoband %>% mutate(start = as.numeric(start), end = as.numeric(end))
  
  ## read in pred_cnv_genes.dat 
  pred_cnv_genes <- read.table(paste0("~/MRes_project_1/docs/SC_SEQ/", i,
                                      "/lambrechtslab/inferCNV/72_hr_runtime/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_genes.dat"), header=TRUE)
  pred_cnv_genes$sample_id <- gsub("_.*", "", pred_cnv_genes$cell_group_name)
  pred_cnv_genes <- pred_cnv_genes %>% mutate(start = as.numeric(start), end = as.numeric(end))
  
  ## fuzzyjoin to combine dataframe based on gene position overlapping with cytoband regions and filter out genes that does not lie in regions
  mapped_genes <- fuzzyjoin::genome_join(pred_cnv_genes, cytoband, by=c("chr", "start", "end"))
  mapped_genes <- mapped_genes %>% 
    dplyr::filter(start.x >= start.y & end.x <= end.y) %>% 
    dplyr::select(sample_id, gene, state, chr=chr.x, start=start.x, end=end.x, arm) 
  mapped_genes$chr_arm <- paste0(mapped_genes$chr, str_extract(mapped_genes$arm, pattern = "p|q")) ## add chr_arm column 
  
  ## save the final output dataframe 
  broad_value_by_arm <- mapped_genes %>% group_by(sample_id, chr_arm) %>% summarise(aneu_score = mean(state)) %>% 
    pivot_wider(names_from = chr_arm, values_from = aneu_score) %>% as.data.frame
  rownames(broad_value_by_arm) <- broad_value_by_arm$sample_id
  broad_value_by_arm$sample_id <- NULL
  broad_value_by_arm <- broad_value_by_arm - 3
  
  ## write table out 
  write.table(broad_value_by_arm, paste0("~/MRes_project_1/docs/SC_SEQ/", i,
                                         "/lambrechtslab/inferCNV/72_hr_runtime/broad_values_by_arm.txt"), 
              sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)
  
  ## visualisation via boxplot 
  broad_value_by_arm <- t(broad_value_by_arm) %>% as.data.frame
  broad_value_by_arm$id <- rownames(broad_value_by_arm)
  df_long_3 <- melt(broad_value_by_arm, id.vars = "id")
  
  pdf(paste0("~/MRes_project_1/docs/SC_SEQ/", i, "/lambrechtslab/plots/inferCNVboxPlot.pdf"), width=8, height = 5)
  print(ggplot(df_long_3, aes(x = factor(id), y = value)) + geom_boxplot() + 
          xlab("chromosome arm") + ylab("copy number variation (0=diploid)") + ylim(c(-3, 3)) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
          geom_hline(yintercept=0, color = "red", linewidth=0.5, alpha=0.5) +
          ggtitle(paste0("Boxplot by chromosome arm for ", i, " lambrechtslab by inferCNV"))) 
  dev.off()
}