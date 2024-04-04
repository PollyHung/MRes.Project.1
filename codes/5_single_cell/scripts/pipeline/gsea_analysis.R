## library 
library(dplyr)
library(magrittr)
library(stringr)

## remove key words 
keywords <- c("parasite", "viral", "virus", "bacteria", "symbiont", "bacterium")

## function 
gsea_process <- function(gsea){
  gsea <- gsea %>% dplyr::select(NAME, SIZE, ES, NES, `FDR q-val`, `FWER p-val`, `RANK AT MAX`) %>% 
    dplyr::mutate(description = tolower(gsub("GOBP_|REACTOME_|KEGG_|WP_", "", NAME))) %>% 
    dplyr::mutate(description = gsub("_", " ", description)) %>% 
    dplyr::rename(name = NAME, size = SIZE, enrichmentScore = ES, normalizedEnrichmentScore = NES, FDR = `FDR q-val`, 
                  FWER = `FWER p-val`, rank = `RANK AT MAX`)
  gsea$Pathway <- sapply(str_split(gsea$name, pattern = "_"), `[`, 1)
  gsea$FDR <- as.numeric(gsea$FDR)
  gsea$log10FDR <- -log10(gsea$FDR)
  gsea$log10FDR[which(gsea$log10FDR == Inf)] <- runif(length(which(gsea$log10FDR == Inf)), min = 3, max = 4)
  gsea$label <- ifelse(gsea$FDR < 0.05, as.character(gsea$description), "") 
  gsea$label <- gsub("regulation", "reg.", gsea$label)
  gsea$label <- gsub("response", "resp.", gsea$label)
  gsea$label <- gsub("within", "w/in", gsea$label)
  gsea$label <- gsub("second", "2nd", gsea$label)
  gsea$label <- ifelse(sapply(gsea$label, function(x) any(sapply(keywords, grepl, x))), "", gsea$label)
  #gsea$alpha <- ifelse(gsea$FDR < 0.05, 1, ifelse(gsea$FDR < 0.20, 0.5, 0))
  
  return(gsea)
}


#setwd("~/MRes_project_1/codes/6_single_cell/RNotebook/subcluster_analysis/b_cells/")
#gsea_go <- readxl::read_xlsx("~/MRes_project_1/docs/SC_SEQ/ova/GSE180661/GSEA_result/GSEA_GO.xlsx", sheet = "b")
#gsea_path <- readxl::read_xlsx("~/MRes_project_1/docs/SC_SEQ/ova/GSE180661/GSEA_result/gsea_pathway.xlsx", sheet = "b")

#gsea_go <- gsea_process(gsea_go)
#gsea_path <- gsea_process(gsea_path)

#colnames(gsea_go) == colnames(gsea_path)

#gsea <- rbind(gsea_go, gsea_path)

gsea <- readxl::read_xlsx("~/MRes_project_1/codes/6_single_cell/gsea_analysis/GSEA_results.xlsx", sheet = 1)
gsea <- gsea_process(gsea)

p1 <- ggplot(gsea, aes(x = as.numeric(normalizedEnrichmentScore), y = log10FDR, color = as.factor(Pathway))) + 
  geom_point(shape = 16, size = 1) + theme_bw() + scale_color_viridis_d(option = "inferno", end = 0.5) + 
  geom_text_repel(aes(label = label), #, alpha = alpha
                  size = 2,           # Adjust text size
                  #box.padding = unit(0.05, "lines"), # Adjust padding around text
                  point.padding = unit(0.5, "lines"),# Adjust space between text and point
                  segment.color = 'grey50', # Color of the connecting line
                  segment.size = 0.2, 
                  max.overlaps = 10) +
  labs(x = "Normalized Enrichment Score", y = "-log10(FDR)", color = "Functional Analysis") +  
  ggtitle("differentially regulated pathways in myeloid cells with amp2q status (GSEA volcano plot)") + 
  xlim(-4, 4) + ylim(0, 4) + 
  theme(legend.position = "top",
        axis.text.x = element_text(size = 6),  # Increase x-axis label size
        axis.text.y = element_text(size = 6),  # Increase y-axis label size
        axis.title.x = element_text(size = 6),  # Increase x-axis title size
        axis.title.y = element_text(size = 6),  # Increase y-axis title size
        axis.ticks.length = unit(0.2, "cm"), 
        title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 6)) + 
  guides(alpha = FALSE) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "pink", linewidth = 0.5) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "pink", linewidth = 0.5) 

ggsave("cancer_DE_pathway.png", p1, width = 4, height = 4, units = "in")




#gsea <- readxl::read_xlsx("~/MRes_project_1/codes/6_single_cell/RNotebook/subcluster_analysis/myeloid/cancer_myeloid.xlsx")
#gsea <- gsea_process(gsea)
gsea <- gsea %>% dplyr::filter(label != "")
gsea <- gsea %>% mutate(Sign = ifelse(normalizedEnrichmentScore >= 0, "Positive", "Negative")) #%>% 
  #dplyr::filter(Pathway == "REACTOME")

gsea <- dplyr::filter(gsea, FDR < 0.20)

top_positive <- gsea %>% filter(Sign == "Positive") %>% top_n(10, normalizedEnrichmentScore)
top_negative <- gsea %>% filter(Sign == "Negative") %>% top_n(-10, normalizedEnrichmentScore)
top_gsea <- bind_rows(top_positive, top_negative)
top_gsea$label <- paste0("[", top_gsea$Pathway, "] ", top_gsea$description)
top_gsea$alpha <- ifelse(top_gsea$FDR < 0.05, 1, 0.65)

ggplot(top_gsea, aes(x = reorder(description, normalizedEnrichmentScore), y = normalizedEnrichmentScore, fill = Sign, alpha = alpha)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c("Positive" = "#FF407D", "Negative" = "#1B3C73"), 
                    labels = c("Downregulation", "Upregulation")) +
  scale_alpha_identity() +  
  labs(x = "Pathway", 
       y = "Normalized Enrichment Score", 
       title = "Top Enriched Reactome Pathways in chr2q amp T cells",
       fill = "Direction") + 
  theme_bw() +
  theme(legend.position = "top", 
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        axis.ticks.length = unit(0.1, "cm"), 
        #title = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 13)) + 
  guides(alpha = FALSE)





