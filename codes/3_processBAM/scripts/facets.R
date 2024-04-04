library(dplyr)
library(magrittr)
library(httr)
library(facets)

sample_ids <- read.table("~/MRes_project_1/codes/3_processBAM/4_docs/sample_ids/lung_ids.txt") %>% unlist
sample_ids <- gsub("_sorted_nodup.bam", "", sample_ids) %>% unique

#"N_sorted_nodup.bam|T_sorted_nodup.bam"

## execute facets 
for(cancer in sample_ids){
  ## output file: 
  output <- paste0("~/MRes_project_1/docs/HH_lung/facets_original/dataframe/", cancer, ".csv")
  rds <- paste0("~/MRes_project_1/docs/HH_lung/facets_original/rds/", cancer, ".rds")
  
  ## if the output file does not exist 
  #if (!file.exists(output)) {
  
  datafile <- paste0("~/MRes_project_1/docs/HH_lung/facets/facet_input/", cancer, "_ordered.csv")
  
  ## preprocessing 
  set.seed(1234)
  rcmat <-  readSnpMatrix(datafile)
  xx <- preProcSample(rcmat)
  oo <- procSample(xx, dipLogR=xx$dipLogR) ## cval 100
  
  ## fitting 
  fit <- emcncf(oo)
  
  ## plotting the results and diagnostic plot 
  pdf(paste0("~/MRes_project_1/docs/HH_lung/facets_original/plots/", cancer, ".pdf"), 
      width = 10, height = 7)
  plotSample(x = oo, emfit = fit) %>% print
  logRlogORspider(oo$out, oo$dipLogR) %>% print
  dev.off()
  
  ## write the table and save rds object
  cncf <- fit$cncf
  write.csv(cncf, output, quote = FALSE, col.names = TRUE, row.names = FALSE)
  saveRDS(fit, rds)
  #}
  
  print(paste0("finish running facets on ", cancer))
}


