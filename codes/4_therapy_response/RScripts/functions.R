## Read files ============================================================================
read_files <- function(directory, folder, cancer, file_condition, return = FALSE){ 
  ## This function reads in gistic processed files from RDS directory, parameters are as following: 
  ## directory = GISTIC2
  ## folder = tcga, cbioportal, hh_lung etc...
  ## cancer = TCGA <- c("lung", "ova", "mel")
  ## file_condition = either "remote_|_cutoffs$" or "all_lesions|broad_values" or missing
  
  setwd(directory) ## set working directory 
  list <- list.files(paste0(directory, "/", folder,"/", cancer)) ## list files in each cancer type gistic2 result
  
  ## if file_condition parameter is not missing, we perform this 
  if(!missing(file_condition)){ list <- list[grepl(file_condition, list)] } 
  
  ## read in files 
  for(file in list){
    file_full <- paste0(directory, "/", folder,"/", cancer, "/", file) ## formulate the full file path 
    file_name <- make.names(paste0(cancer, "_", gsub("\\.txt$|\\.tsv$|\\.tab$", "", file))) # Create a file name 
    
    print(paste0("read in ", file_name))
    df <- read.delim(file_full, stringsAsFactors = FALSE)
    assign(file_name, df, envir = .GlobalEnv)
  }
  
  ## should the file be returned?
  if(return){
    return(df)
  }
}

## focal cleaning ============================================================================ 
focal_cleaning <- function(cancer){
  ## needs to have the patient list 
  
  ## get files 
  focal <- get(paste0(cancer, "_all_lesions.conf_95"), envir = .GlobalEnv)
  sample_name <- get(cancer, envir = .GlobalEnv)
  
  ## cleaning 
  focal <- focal[grepl("CN values$", focal$Unique.Name), ]  ## select the unthresholded values 
  focal$Unique.Name <- gsub("lification Peak|- CN values|etion Peak|\\s+", "", focal$Unique.Name) ## clean name 
  focal$Unique.Name <- paste0(focal$Unique.Name, "_", focal$Descriptor) %>% tolower ## create unique identifier 
  rownames(focal) <- focal$Unique.Name ## rename 
  colnames(focal) <- gsub("\\.", "-", colnames(focal)) ## change the sample_name to uniform format with -
  focal <- focal[, sample_name] ## select AND align the dataframe in the sequence of sample_name 
  focal <- as.matrix(focal) ## make it to matrix for subsequent coxph 
  
  return(focal)
}


## broad cleaning ============================================================================  
broad_cleaning <- function(cancer){
  ## get files 
  broad <- get(paste0(cancer, "_broad_values_by_arm"), envir = .GlobalEnv)
  sample_name <- get(cancer, envir = .GlobalEnv)
  
  ## cleaning 
  rownames(broad) <- broad$Chromosome.Arm ## add row names 
  colnames(broad) <- gsub("\\.", "-", colnames(broad)) ## change the sample_name to uniform format with -
  broad <- broad[, sample_name] ## select AND align the dataframe in the sequence of sample_name 
  broad <- as.matrix(broad) ## make it to matrix for subsequent coxph 
  
  return(broad)
}


## phenotype_cleaning ============================================================================  
phenotype_cleaning <- function(cancer, confounders, col_names){
  ## confounders = c("^age_", "gender", "_stage$", "^primary_therapy")
  ## col_names = c("age", "sex", "stage", "primary_outcome")
  
  ## get files 
  phenotype <- get(paste0(cancer, "_phenotype"), envir = .GlobalEnv)
  sample_name <- get(cancer, envir = .GlobalEnv)
  
  ## data cleaning general 
  rownames(phenotype) <- phenotype$sampleID
  phenotype <- phenotype[sample_name, unlist(lapply(confounders, function(patt) grep(patt, names(phenotype), value = TRUE)))]
  colnames(phenotype) <- col_names
  
  ## data cleaning for stage 
  phenotype <- phenotype %>%
    mutate(stage = case_when(
      grepl("^Stage I[AB]?$|^I/II NOS$", stage) ~ 1,
      grepl("^Stage II[ABC]?$", stage) ~ 2,
      grepl("^Stage III[ABC]?$", stage) ~ 3,
      stage == "Stage IV" ~ 4,
      TRUE ~ NA  # For anything else, or you can adjust this as needed
    ))
  phenotype$stage <- as.numeric(phenotype$stage)
  
  return(phenotype)
}


## survival ============================================================================  
survival <- function(cancer, surv_type){
  ## surv_type either OS, DSS, DFI, or PFI
  
  ## get files 
  sample_name <- get(cancer, envir = .GlobalEnv)
  survival <- get(paste0(cancer, "_survival"), envir = .GlobalEnv)
  rownames(survival) <- survival$sample
  survival <- survival[sample_name, ]
  
  ## calculate survival objects 
  surv_time <- paste0(surv_type, ".time")
  time <- as.numeric(survival[sample_name, surv_time]) ## this ensures that the survival object produced are aligned
  event <- as.numeric(survival[sample_name, surv_type])
  obj <- Surv(time, event)
  
  ## store it to environment 
  name <- make.names(paste0(cancer, ".", tolower(surv_type)))
  assign(name, obj, envir = .GlobalEnv) ## get the survival object 
}


## surv_analysis ============================================================================  
surv_analysis <- function(cancer, surv_type, confounders, col_names, focal = FALSE, broad = FALSE){
  ## confounders = c("^age_", "gender", "_stage$", "^primary_therapy")
  ## surv_type = either "os" or "pfi"
  ## col_names = c("age", "sex", "stage", "primary_outcome")
  ## focal = focal analysis: FALSE or TRUE
  ## broad = broad analysis: FALSE or TRUE
  
  ## get basic objects 
  survival <- get(paste0(cancer, ".", surv_type), envir = .GlobalEnv)
  phenotype <- get(paste0(cancer, "_phenotype"), envir = .GlobalEnv)
  sample_name <- get(cancer, envir = .GlobalEnv)
  
  ## do broad or focal analysis 
  if(focal){gistic <- focal_cleaning(cancer)} 
  if(broad){gistic <- broad_cleaning(cancer)}
  
  ## phenotype data cleaning, extract multivariate confoudners  
  phenotype <- phenotype_cleaning(cancer, confounders = confounders, col_names = col_names)
  age <- phenotype$age
  stage <- phenotype$stage 
  sex <- phenotype$sex
  
  ## This builds an empty dataframe to hold the results 
  results <- array(NA, c(nrow(gistic), 5))  ## dataframe size determined by how many genes there are 
  rownames(results) <- rownames(gistic)  ## rownames = genes 
  colnames(results) <- c("HR", "lower_CI", "upper_CI", "z_scores", "p_values")  
  results <- as.data.frame(results) ## set it as a dataframe 
  
  ## multivariate analysis accounting for age and stage 
  for(j in 1:nrow(gistic)){
    
    if(cancer %in% c("lung", "mel")){ ## if the cancer occurs in both sex, add sex as confounding variable 
      coxphmodel <- coxph(survival ~ gistic[j, ]+age+stage+sex) 
    } else { ## if the cancer does not occur in both sex (ovarian), do NOT add sex as confounding variable 
      coxphmodel <- coxph(survival ~ gistic[j, ]+age+stage) # survival ~ gistic[j, ]
    }
    
    results$HR[j] <- summary(coxphmodel)$coef[1, 2]
    results$lower_CI[j] <- summary(coxphmodel)$conf.int[1, 3]
    results$upper_CI[j] <- summary(coxphmodel)$conf.int[1, 4]
    results$z_scores[j] <- summary(coxphmodel)$coef[1, 4]
    results$p_values[j] <- summary(coxphmodel)$coef[1, 5]
  }
  
  ## Multiple-testing correction by False Discovery Rate and rank the dataframe 
  results <- as.data.frame(results)
  results$FDR <- p.adjust(results$p_values, method = "fdr")
  results <- results[order(results$FDR, decreasing = FALSE), ] 
  
  return(results)
}


## pheno_analysis  ============================================================================  
pheno_analysis <- function(cancer, event, confounders, col_names, focal = FALSE, broad = FALSE){
  ## cancer = what type of cancer you have 
  ## confounders = c("^age_", "gender", "_stage$", "^primary_therapy")
  ## surv_type = either "os" or "pfi"
  ## col_names = c("age", "sex", "stage", "primary_outcome")
  ## focal = focal analysis: FALSE or TRUE
  ## broad = broad analysis: FALSE or TRUE
  ## event = either age, or sex, stage, or primary outcome 
  
  ## get basic objects 
  phenotype <- get(paste0(cancer, "_phenotype"), envir = .GlobalEnv)
  sample_name <- get(cancer, envir = .GlobalEnv)
  
  ## do broad or focal analysis 
  if(focal){gistic <- focal_cleaning(cancer)} 
  if(broad){gistic <- broad_cleaning(cancer)}
  
  ## phenotype data cleaning, extract event for analyais 
  phenotype <- phenotype_cleaning(cancer, confounders = confounders, col_names = col_names)
  event <- phenotype[, which(colnames(phenotype) == event)]
  print(paste0(round(length(which(event == ""))/length(event)*100, digits = 2), "% of data is missing from ", cancer, " dataset."))
  
  ## This builds an empty dataframe to hold the results 
  results <- array(NA, c(nrow(gistic), 3))  
  rownames(results) <- rownames(gistic)  
  colnames(results) <- c("test", "cor.val", "p_values")  
  results <- as.data.frame(results) ## set it as a dataframe 
  
  ## if the event is continuous like age 
  for(j in 1:nrow(gistic)){
    if(is.numeric(event)){
      rough <- cor.test(event, gistic[j, ], method = "spearman")
      #neat <- lm(event ~ gistic[j, ], data = gistic)
      results$test <- "cor.test"
      results$cor.val[j] <- rough$estimate
    } else if(nlevels(as.factor(event)) == 2){
      event <- as.factor(event)
      rough <- wilcox.test(gistic[j, ] ~ event, data = gistic)
      results$test <- "wilcox.test"
    } else if(nlevels(as.factor(event)) > 2){
      event <- as.factor(event)
      rough <- kruskal.test(gistic[j, ] ~ event, data = gistic)
      results$test <- "kruskal.test"
      results$cor.val[j] <- rough$statistic
    }
    results$p_values[j] <- rough$p.value
  }
  
  ## Multiple-testing correction -------------
  results <- as.data.frame(results)
  results$FDR <- p.adjust(results$p_values, method = "fdr")
  results <- results[order(results$FDR, decreasing = FALSE), ] 
  
  return(results)
}


## immune_pheno  ============================================================================  
immune_pheno <- function(cancer, rna_seq, mcp_pat, tls_sig, tmb_col, gene_pat, write_each = FALSE, write_uni = TRUE){
  ## cancer = cancer type, usually given as "i" in a for (i in TCGA){} loop 
  ## rna_seq = the rna_seq dataframe 
  ## write_result = TRUE or FALSE (default to FALSE), determines if the result is write to environment. 
  ## tls_signature = c("CD79B", "CD1D", "CCR6", "LAT", "SKAP1", "CETP", "EIF1AY", "RBP5", "PTGDS", "CCL19", "CCL21", "CXCL13", "CCR7", "CXCR5", "SELL", "LAMP3")
  ## tmb_col = c("Diagnosis Age", "Sex", "Fraction Genome Altered", "Mutation Count", "TMB (nonsynonymous)") 
  ## gene_pat = "^HLA.*[ABC]$|^CD274$|^GZMB$|^GNLY$|^PRF1$"
  ## mcp_pat = c("t.cell", "cd8_t.cell", "cytotoxic_lymphocytes", "b_lineage", "nk_cells", "monocytic_lineage", "myeloid_dendritic_cells", "neutrophils", "endothelial_cells", "fibroblasts")
  
  ## Read in samples and data cleaning -----------------------
  sample_names <- get(cancer, envir = .GlobalEnv)
  rna_seq <- get(rna_seq, envir = .GlobalEnv)
  sample_names <- intersect(sample_names, colnames(rna_seq)) ## this will be our new sample_names 
  rna_seq <- rna_seq[, sample_names] ## we will only preserve cancer patients we care about, this is a DATAFRAME!!
  
  ## MCP counter ---------------------------- 
  mcp_results <- MCPcounter.estimate(
    expression = rna_seq, ## dataframe rna_seq 
    featuresType = "HUGO_symbols")
  rownames(mcp_results) <- c("t.cell", "cd8_t.cell", "cytotoxic_lymphocytes", "b_lineage", "nk_cells", 
                             "monocytic_lineage", "myeloid_dendritic_cells", "neutrophils", "endothelial_cells", 
                             "fibroblasts")
  mcp_results <- as.data.frame(t(mcp_results))
  mcp_results <- select(mcp_results, all_of(mcp_pat))
  
  ## TLS status -----------------------
  tls_signature <- tls_sig  ## from literature 
  tls_rna_seq <- rna_seq[tls_signature, ] ## isolate their rna_expression
  tls_rna_seq <- as.matrix(t(tls_rna_seq)) ## transform it into a matrix 
  tls_rna_seq <- scale(tls_rna_seq, center = F) ## scale it 
  tls_rna_seq <- as.data.frame(tls_rna_seq) ## transform it back into a dataframe 
  tls_rna_seq$TLS_mean <- apply(tls_rna_seq, 1, mean) ## calculate the TLS by mean 
  
  ## TMB calculation ----------------------
  tmb <- read_tsv(paste0(MUTATION, "/", cancer,"_tmb.tsv"), show_col_types = FALSE) ## the read in is a tibble 
  tmb <- as.data.frame(tmb)
  rownames(tmb) <- tmb$`Sample ID`
  sample_names <- intersect(sample_names, rownames(tmb)) ## there are fewer patients than the original patient list. 
  tmb <- tmb[sample_names, tmb_col]
  
  ## PDL1 exp, GZMB, GNLY, PRF1, MHC1 expression ----------------------
  genes <- which(grepl(gene_pat, rownames(RNA_seq)))
  spec_gene <- rna_seq[genes, ] 
  MHCI <- colMeans(spec_gene[rownames(spec_gene) %in% c("HLA-A", "HLA-B", "HLA-C"), ], na.rm = TRUE)
  MHCII <- colMeans(spec_gene[which(grepl("^HLA-D", rownames(spec_gene))), ], na.rm = TRUE)
  spec_gene <- rbind(spec_gene, MHCI = MHCI, MHCII = MHCII)
  spec_gene <- as.data.frame(t(spec_gene), )
  spec_gene <- spec_gene[, which(!grepl("^HLA", colnames(spec_gene)))]
  
  ## immuno_phenotype -------------------------
  ## this checks if the order of the rownames are identical because cbind does not bind by rownames 
  CHECK <- identical(rownames(mcp_results), rownames(tls_rna_seq)) && identical(rownames(mcp_results), rownames(spec_gene))
  if(CHECK){
    temp <- cbind(mcp_results, tls_rna_seq$TLS_mean, spec_gene)
    immuno_phenotype <- merge(temp, tmb, by = "row.names", all.x = TRUE)
    rownames(immuno_phenotype) <- immuno_phenotype$Row.names
    immuno_phenotype$Row.names <- NULL
    immuno_phenotype$y <- NULL
  }
  
  ## Write the unified table ------------------
  if(write_uni){
    print(paste0("This is the organised immuno_phenotype table for ", cancer, ". ", 
                 "Number of samples: ", nrow(immuno_phenotype)))
    print(immuno_phenotype[1:3, ])
    
    immuno_phenotype_name <- make.names(paste0(cancer, "_immuno_pheno"))
    assign(immuno_phenotype_name, immuno_phenotype, envir = .GlobalEnv)
    
    rna_seq_name <- make.names(paste0(cancer, "_rna_seq"))
    assign(rna_seq_name, rna_seq, envir = .GlobalEnv)
  }
  
  ## Write out each separated table ----------------
  if(write_each){
    print(paste0("This is immune phenotype for ", cancer, " cancer"))
    print(paste0("This is MCP, number of samples: ", nrow(mcp_results))); print(mcp_results[1:3, ])
    print(paste0("This is TLS, number of samples: ", nrow(tls_rna_seq))); print(tls_rna_seq[1:3, ])
    print(paste0("This is TMB, number of samples: ", nrow(tmb))); print(tmb[1:3, ])
    print(paste0("This is selected genes, number of samples: ", nrow(spec_gene))); print(spec_gene[1:3, ])
    
    rna_seq_name <- make.names(paste0(cancer, "_rna_seq"))
    assign(rna_seq_name, rna_seq, envir = .GlobalEnv)
    write.table(rna_seq, file = paste0(RNA_SEQ, "/", rna_seq_name, ".txt"), quote = FALSE, 
                row.names = TRUE, col.names = TRUE)
    
    mcp_name <- make.names(paste0(cancer, "_mcp_results"))
    assign(mcp_name, mcp_results, envir = .GlobalEnv)
    
    tls_rna_seq_name <- make.names(paste0(cancer, "_tls"))
    assign(tls_rna_seq_name, tls_rna_seq, envir = .GlobalEnv)
    
    tmb_name <- make.names(paste0(cancer, "_tmb"))
    assign(tmb_name, tmb, envir = .GlobalEnv)
    
    spec_gene_name <- make.names(paste0(cancer, "_spec_gene"))
    assign(spec_gene_name, spec_gene, envir = .GlobalEnv)
  }
}


## statistics_test ============================================================================  
statistics_test <- function(cancer, focal = FALSE, broad = FALSE, not_gistic = FALSE, test){
  ## cancer = lung, ova, or mel 
  ## test = either spearman or linear_regression
  
  ## build a list
  immuno <- list()
  
  ## get needed data 
  if(focal){cnv <- focal_cleaning(cancer)}
  if(broad){cnv <- broad_cleaning(cancer)}
  if(not_gistic){
    cnv <- get(paste0(cancer, "_scna_as"), envir = .GlobalEnv)
    cnv <- t(cnv) %>% as.data.frame()
  }
  immuno_data <- get(paste0(cancer, "_immuno_pheno"), envir = .GlobalEnv)
  pheno <- phenotype_cleaning(cancer, confounders = c("^age_", "_stage$"), col_names = c("age", "stage"))
  index <- intersect(rownames(immuno_data), colnames(cnv)) ## select out common samples 
  
  ## further process 
  pheno <- pheno[index, ]
  immuno_data <- as.matrix(immuno_data[index, ])
  cnv <- as.matrix(cnv[, index]) 
  col_names <- colnames(immuno_data)
  
  if(test == "spearman"){ ## perform spearman statistical test 
    ## run test 
    for(y in 1:ncol(immuno_data)){
      ## set up empty dataframe
      results <- array(NA, c(nrow(cnv), 2))
      rownames(results) <- rownames(cnv)
      colnames(results) <- c("cor.val", "p_values")
      results <- as.data.frame(results) ## set it as a dataframe
      
      for(x in 1:nrow(cnv)){ 
        temp <- cor.test(immuno_data[, y], as.numeric(cnv[x, ]), method = "spearman")
        results$cor.val[x] <- temp$estimate
        results$p_values[x] <- temp$p.value
      }
      
      ## multiple testing correlation by fdr? 
      results <- as.data.frame(results)
      results$FDR <- p.adjust(results$p_values, method = "fdr")
      results <- results[order(results$FDR, decreasing = FALSE), ]
      
      #print(paste0("finished analysing column ", y))
      immuno[[col_names[y]]] <- results
    }
    
    ## random test, let's use 1p 
    check_test <- cor.test(immuno_data[, 3], as.numeric(cnv[1, ]), method = "spearman")
  }
  
  if(test == "linear_regression"){ ## perform linear regression statistical test 
    ## run test 
    for(y in 1:ncol(immuno_data)){ #1:ncol(immuno_data) to replace 1:2
      ## set up empty dataframe
      results <- array(NA, c(nrow(cnv), 5))
      rownames(results) <- rownames(cnv)
      colnames(results) <- c("estimate", "std_error", "t_value", "p_values", "r_sq")
      results <- as.data.frame(results) ## set it as a dataframe
      
      for(x in 1:nrow(cnv)){ #1:nrow(cnv)
        temp <- lm(immuno_data[, y]~cnv[x, ]+pheno$age+pheno$stage, data = as.data.frame(cnv))
        temp <- summary(temp)
        results$estimate[x] <- temp$coefficients[2, 1]
        results$std_error[x] <- temp$coefficients[2, 2]
        results$t_value[x] <- temp$coefficients[2, 3]
        results$p_values[x] <- temp$coefficients[2, 4]
        results$r_sq[x] <- temp$adj.r.squared
      }
      
      ## multiple testing correlation by fdr? 
      results <- as.data.frame(results)
      results$FDR <- p.adjust(results$p_values, method = "fdr")
      results <- results[order(results$FDR, decreasing = FALSE), ]
      
      #print(paste0("finished analysing column ", y))
      immuno[[col_names[y]]] <- results
    }
    
    ## random test, let's use the first row of the dataframe against tls 
    check_test <- lm(immuno_data[, 3]~cnv[1, ]+pheno$age+pheno$stage, data = as.data.frame(cnv))
    check_test <- summary(check_test)
  }
  
  ## check test 
  print(paste0("this is a check test using first row of the cnv data and tls for ", cancer))
  print(check_test)
  
  return(immuno)
}


save.image(file = "/rds/general/user/ph323/home/MRes_project_1/Codes/function.RData")

