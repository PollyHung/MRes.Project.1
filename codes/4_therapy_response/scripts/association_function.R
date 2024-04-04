## the clinical data and immune data file must begin with a column named "sample"
## follwed by columns of interested features 
## example: 
## sample.ids.dir <- "~/MRes_project_1/Codes/3_therapy_resp/sample_ids.txt"
## clinical.data.dir <- "~/MRes_project_1/docs/HH_ova/clinical_data/ho_clinical.txt"
## survival.data.dir <- "~/MRes_project_1/docs/HH_ova/clinical_data/ho_survival.txt"
## focal.dir <- "~/MRes_project_1/docs/HH_ova/gistic/facets_seg/cval_50/hisens_adj/all_lesions.conf_95.txt"
## broad.dir <- "~/MRes_project_1/docs/HH_ova/gistic/facets_seg/cval_50/hisens_adj/broad_values_by_arm.txt" 

## small function ========================================================

mk_df <- function(ROW_NUM, ## c(nrow(gistic)
                  COL_NUM, ## 5
                  ROW_NAMES, ## rownames(gistic)
                  COL_NAMES){ ## c("estimate", "std_error", "t_value", "p_values", "r_sq")
  ## this small function makes the empty dataframe                 
  output <- matrix(NA, ROW_NUM, COL_NUM) %>% as.data.frame
  rownames(output) <- ROW_NAMES
  colnames(output) <- COL_NAMES
  
  return(output)
}


after_proc <- function(INPUT, ## input dataframe, usually output 
                       OUTPUT_NAME, ## paste0(i, "_", names(immune_data[j]))
                       #p.val, 
                       result.dir){ 
  ## this small function post-process the output dataframe 
  output <- as.data.frame(INPUT)
  #output <- filter(output, output$p_values < p.val)
  output$FDR <- p.adjust(output$p_values, method = "fdr")
  output <- output[order(output$FDR, output$p_values, decreasing = FALSE), ] 
  
  output_name <- make.names(OUTPUT_NAME)
  assign(output_name, output, envir = .GlobalEnv)
  write.table(output, paste0(result.dir, output_name, ".txt"), sep = "\t", 
              col.names = TRUE, row.names = TRUE, quote = FALSE)
}



## Association function ========================================================

association <- function(result.dir, ## path to your saving directory 
                        sample.ids.dir, ## path to your sample ids txt file 
                        clinical.data.dir, ## path to your clinical features (age, stage, sex etc.) txt file 
                        survival.data.dir, ## path to your survival features (survival) txt file 
                        focal.dir, ## path to your all_lesions.conf_95.txt file 
                        broad.dir){ ## path to your broad_values_by_arm.txt file
  
  
  ## Environment ===============================================================
  ## Library all the packages needed 
  library(dplyr)
  library(magrittr)
  library(tibble)
  library(readxl)
  library(stringr)
  library(survival)
  
  
  ## Prepare files =============================================================
  ## read in raw files 
  sample_ids <- read.table(sample.ids.dir) %>% unlist 
  clinical_data <- read.table(clinical.data.dir, sep = "\t", header = TRUE, row.names = 1) 
  survival_data <- read.table(survival.data.dir, sep = "\t", header = TRUE, row.names = 1)
  focal <- read.table(focal.dir, sep = "\t", header = TRUE)
  broad <- read.table(broad.dir, sep = "\t", header = TRUE)
  
  ## edit and reorganize the focal and broad data    
  focal <- focal[grepl("CN values$", focal$Unique.Name), ]
  focal$Unique.Name <- gsub("lification Peak|- CN values|etion Peak|\\s+", "", focal$Unique.Name) ## clean name 
  focal$Unique.Name <- paste0(focal$Unique.Name, "_", focal$Descriptor) %>% tolower ## create unique identifier 
  rownames(focal) <- focal$Unique.Name ## rename 
  rownames(broad) <- broad$Chromosome.Arm
  
  ## reorder dataframe and remove unwanted columns/rows 
  common_samples <- rownames(clinical_data) %>% intersect(rownames(survival_data)) %>% intersect(colnames(broad))
  clinical_data <- clinical_data[common_samples, ]
  survival_data <- survival_data[common_samples, ]
  focal <- focal[common_samples] %>% as.matrix()
  broad <- broad[common_samples] %>% as.matrix()
  
  
  ## Perform statistics test according to data type ============================
  for(i in c("focal", "broad")){
     for(j in 1:ncol(clinical_data)) {
      ## test if the column is qualitative or quantitative 
      
      ## if the column only has 2 layers 
      if(length(unique(na.omit(clinical_data[, j]))) == 2) {
        gistic <- get(i)
        output <- mk_df(nrow(gistic), 2, rownames(gistic), 
                        c("estimate", "p_values")) ## empty df 
        ## wilcox test 
        for (z in 1:nrow(gistic)) {
          group1 <- gistic[z, which(clinical_data[, j] == unique(na.omit(clinical_data[, j]))[1])]
          group2 <- gistic[z, which(clinical_data[, j] == unique(na.omit(clinical_data[, j]))[2])]
          
          wilcox <- wilcox.test(group1, group2) 
          output$estimate[z] <- wilcox$statistic %>% as.numeric
          output$p_values[z] <- wilcox$p.value %>% as.numeric  
        } 
        after_proc(output, paste(i, names(clinical_data[j]), sep = "."), result.dir) ## post-process 
        
        ## if the column has > 2 columns, then use linear model 
      } else {
        gistic <- get(i)
        output <- mk_df(nrow(gistic), 5, rownames(gistic), 
                        c("estimate", "std_error", "t_value", "p_values", "r_sq")) ## empty df 
        for (z in 1:nrow(gistic)) { 
          coxphmodel <- lm(clinical_data[, j]~gistic[z, ]+clinical_data$age+clinical_data$stage+clinical_data$sex, 
                           data = as.data.frame(gistic)) %>% summary
          output[z, ] <- c(coxphmodel$coefficients[2, 1:4], coxphmodel$adj.r.squared)
        } 
        after_proc(output, paste0(i, names(clinical_data[j]), sep = "."), result.dir) ## post-process 
      }
    }
  }
  
  ## Perform survival associated analysis ======================================
  ## isolate survival time and events for os and pfs 
  os_time <- as.numeric(survival_data$os_time)
  os_event <- as.numeric(survival_data$os_event)
  pfs_time <- as.numeric(survival_data$pfs_time)
  pfs_event <- as.numeric(survival_data$pfs_event)
  
  ## Make survival object 
  os_surv <- Surv(os_time, os_event)
  pfs_surv <- Surv(pfs_time, pfs_event)
  
  for(obj in c("os_surv", "pfs_surv")){ ## loop though the os and pfs 
    ## get the survival object from environment 
    surv_obj <- get(obj)
    
    for(i in c("focal", "broad")){  ## loop though the focal and broad 
      ## get the focal or broad object from environment 
      gistic <- get(i)
      
      ## build an empty dataframe 
      output <- mk_df(nrow(gistic), 5, rownames(gistic), 
                      c("HR", "lower_CI", "upper_CI", "z_scores", "p_values")) ## empty df 
      ## perform association 
      for(j in 1:nrow(gistic)){
        coxphmodel <- coxph(surv_obj~gistic[j, ]+clinical_data$age+clinical_data$stage+clinical_data$sex) %>% summary
        output[j,] <- c(coxphmodel$coef[1, 1:5])
      }
      after_proc(output, paste0(i, obj, sep = "."), result.dir) ## post-process 
    }
  }

  ## Get the list of all objects in the environment
  all_objects <- ls()
  #objects_to_save <- grep("^(focal|broad)", all_objects, value = TRUE)
  save(list = all_objects, file = paste0(result.dir, "association.RData"))
}














