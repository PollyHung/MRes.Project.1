targets <- read.table("~/MRes_project_1/codes/3_processBAM/4_docs/sureSelect/S07604514_Targets.txt")
colnames(targets) <- "Marker_Name"

library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", 
                      host = "useast.ensembl.org")

ENSG_start <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'chromosome_name', 
                                  'start_position', 'end_position', 'strand', 'ensembl_transcript_id'),
                   filters = 'ensembl_transcript_id', 
                   values = targets,
                   mart = ensembl)

left <- targets %>% filter(!targets$Marker_Name %in% ENSG_start$ensembl_transcript_id)

