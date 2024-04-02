library(dplyr)
library(magrittr)
library(Rsamtools)

DATA <- "/rds/general/user/ph323/ephemeral/HH_ova/Alignment_ov"

file_list <- list.files(DATA)  
file_list <- file_list[grep(".bam$", file_list)]
ids <- sub("(X[0-9]+)[NT]_sorted_nodup.bam", "\\1", file_list)
table_ids <- table(ids)

paired_ova_bam <- file_list[ids %in% names(table_ids[table(ids) == 2])]

bam_ids <- paired_ova_bam[grepl("N_", paired_ova_bam)]

set.seed(1234)
selected_ids <- paste(sample(bam_ids, 5), collapse = " ")

setwd(DATA) ## change to directory 
cmd <- paste("samtools merge pooled.bam", selected_ids) ## create command 
system(cmd) ## execute command in environment 