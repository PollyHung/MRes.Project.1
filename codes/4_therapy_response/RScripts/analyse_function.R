## the clinical data and immune data file must begin with a column named "sample"
## follwed by columns of interested features 
## example: 
## sample.ids.dir <- "~/MRes_project_1/Codes/3_therapy_resp/sample_ids.txt"
## clinical.data.dir <- "~/MRes_project_1/docs/HH_ova/clinical_data/hh_ova.txt"
## survival.data.dir <- "~/MRes_project_1/docs/HH_ova/clinical_data/hh_ova_surv.txt"
## focal.dir <- "~/MRes_project_1/docs/HH_ova/gistic/facets_seg/cval_50/hisens_adj/all_lesions.conf_95.txt"
## broad.dir <- "~/MRes_project_1/docs/HH_ova/gistic/facets_seg/cval_50/hisens_adj/broad_values_by_arm.txt" 

immune_association <- function(prefix, ## the prefix you wish to add 
                               result.dir, ## path to your saving directory 
                               sample.ids.dir, ## path to your sample ids txt file 
                               clinical.data.dir, ## path to your clinical features (age, stage, sex etc.) txt file 
                               survival.data.dir, ## path to your survival features (survival) txt file 
                               focal.dir, ## path to your all_lesions.conf_95.txt file 
                               broad.dir, ## path to your broad_values_by_arm.txt file 
                               p.val
                               ){
  
  
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
  immune_data <- read.table(immune.data.dir, sep = "\t", header = TRUE, row.names = 1) 
  focal <- read.table(focal.dir, sep = "\t", header = TRUE)
  broad <- read.table(broad.dir, sep = "\t", header = TRUE)
  
  ## edit and reorganize the focal and broad data    
  focal <- focal[grepl("CN values$", focal$Unique.Name), ]
  focal$Unique.Name <- gsub("lification Peak|- CN values|etion Peak|\\s+", "", focal$Unique.Name) ## clean name 
  focal$Unique.Name <- paste0(focal$Unique.Name, "_", focal$Descriptor) %>% tolower ## create unique identifier 
  rownames(focal) <- focal$Unique.Name ## rename 
  rownames(broad) <- broad$Chromosome.Arm
  
  ## reorder dataframe and remove unwanted columns/rows 
  common_samples <- rownames(clinical_data) %>% intersect(rownames(survival_data)) %>% 
    intersect(rownames(immune_data)) %>% intersect(colnames(broad))
  clinical_data <- clinical_data[common_samples, ]
  survival_data <- survival_data[common_samples, ]
  immune_data <- immune_data[common_samples, ] 
  focal <- focal[common_samples] %>% as.matrix()
  broad <- broad[common_samples] %>% as.matrix()

  
  
  ## Perform immune association analysis =======================================
  for(i in c("focal", "broad")){
    ## get the focal or broad object from environment 
    gistic <- get(i)
    
    ## perform immune association analysis 
    for(j in 1:ncol(immune_data)){
      ## build an empty dataframe 
      output <- array(NA, c(nrow(gistic), 5))
      rownames(output) <- rownames(gistic)
      colnames(output) <- c("estimate", "std_error", "t_value", "p_values", "r_sq")
      output <- output %>% as.data.frame
      
      ## multivariate analysis accounting for age and stage 
      for (z in 1:nrow(gistic)) {
        coxphmodel <- lm(immune_data[, j]~gistic[z, ]+clinical_data$age+clinical_data$stage, data = as.data.frame(gistic))
        temp <- summary(coxphmodel)
        output$estimate[z] <- temp$coefficients[2, 1]
        output$std_error[z] <- temp$coefficients[2, 2]
        output$t_value[z] <- temp$coefficients[2, 3]
        output$p_values[z] <- temp$coefficients[2, 4]
        output$r_sq[z] <- temp$adj.r.squared
      }
      
      ## first save the file as dataframe 
      output <- as.data.frame(output)
      
      ## remove unsignificant p-values 
      output <- filter(output, output$p_values < p.val)
      
      ## perform multiple testing correction 
      output$FDR <- p.adjust(output$p_values, method = "fdr")
      output <- output[order(output$FDR, output$p_values, decreasing = FALSE), ] 
      
      output_name <- make.names(paste0(i, "_", prefix, "_", names(immune_data[j])))
      assign(output_name, output, envir = .GlobalEnv)
    }
  }
  
  

  ## perform survival associated analysis ======================================
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
      output <- array(NA, c(nrow(gistic), 5))
      rownames(output) <- rownames(gistic)
      colnames(output) <- c("HR", "lower_CI", "upper_CI", "z_scores", "p_values")  
      output <- output %>% as.data.frame
      
      ## perform association 
      for(j in 1:nrow(gistic)){
        coxphmodel <- coxph(surv_obj~gistic[j, ]+clinical_data$age+clinical_data$stage)
        output$HR[j] <- summary(coxphmodel)$coef[1, 2]
        output$lower_CI[j] <- summary(coxphmodel)$conf.int[1, 3]
        output$upper_CI[j] <- summary(coxphmodel)$conf.int[1, 4]
        output$z_scores[j] <- summary(coxphmodel)$coef[1, 4]
        output$p_values[j] <- summary(coxphmodel)$coef[1, 5]
      }
      
      ## multiple testing analysis 
      output <- as.data.frame(output)
      
      ## remove unsignificant p-values 
      output <- filter(output, output$p_values < p.val)
      
      ## perform multiple testing correction 
      output$FDR <- p.adjust(output$p_values, method = "fdr")
      output <- output[order(output$FDR, output$p_values, decreasing = FALSE), ] 
      
      output_name <- make.names(paste0(i, "_", prefix, "_", obj))
      assign(output_name, output, envir = .GlobalEnv)
    }
  }
  
  
  
  ## Perform the clinical features assocation analysis =========================
  for(i in c("focal", "broad")){
    ## get the focal or broad object from environment 
    gistic <- get(i)
    
    ## perform immune association analysis 
    for(j in 1:ncol(clinical_data)){
      ## build an empty dataframe 
      output <- array(NA, c(nrow(gistic), 5))
      rownames(output) <- rownames(gistic)
      colnames(output) <- c("estimate", "std_error", "t_value", "p_values", "r_sq")
      output <- output %>% as.data.frame
      
      ## multivariate analysis accounting for age and stage 
      for (z in 1:nrow(gistic)) {
        
        ## if the data is two level --> wilcox test, else perform linear regression 
        if(length(unique(na.omit(clinical_data[, j]))) == 2){
          group1 <- gistic[z, which(clinical_data[, j] == unique(na.omit(clinical_data[, 4]))[1])]
          group2 <- gistic[z, which(clinical_data[, j] == unique(na.omit(clinical_data[, 4]))[2])]
          
          wilcox <- wilcox.test(group1, group2) 
          output$estimate[z] <- wilcox$statistic %>% as.numeric
          output$p_values[z] <- wilcox$p.value %>% as.numeric
          
          ## if it is age and stage column, then linear model without addressing age and stage 
        } else if(names(clinical_data[j]) %in% c("age", "stage")) {
          coxphmodel <- lm(clinical_data[, j]~gistic[z, ], data = as.data.frame(gistic))
          temp <- summary(coxphmodel)
          output$estimate[z] <- temp$coefficients[2, 1]
          output$std_error[z] <- temp$coefficients[2, 2]
          output$t_value[z] <- temp$coefficients[2, 3]
          output$p_values[z] <- temp$coefficients[2, 4]
          output$r_sq[z] <- temp$adj.r.squared
          
          ## for everything else, use the formula linear regression + age + stage. 
        } else {
          coxphmodel <- lm(clinical_data[, j]~gistic[z, ]+clinical_data$age+clinical_data$stage, data = as.data.frame(gistic))
          temp <- summary(coxphmodel)
          output$estimate[z] <- temp$coefficients[2, 1]
          output$std_error[z] <- temp$coefficients[2, 2]
          output$t_value[z] <- temp$coefficients[2, 3]
          output$p_values[z] <- temp$coefficients[2, 4]
          output$r_sq[z] <- temp$adj.r.squared
          
        }
      }
      
      ## first save the file as dataframe 
      output <- as.data.frame(output)
      
      ## remove unsignificant p-values 
      output <- filter(output, output$p_values < p.val)
      
      ## perform multiple testing correction 
      output$FDR <- p.adjust(output$p_values, method = "fdr")
      output <- output[order(output$FDR, output$p_values, decreasing = FALSE), ] 
      
      output_name <- make.names(paste0(i, "_", prefix, "_", names(clinical_data[j])))
      assign(output_name, output, envir = .GlobalEnv)
    }
  }  
  
  
  
  ## Get the list of all objects in the environment
  all_objects <- ls()
  objects_to_save <- grep("^(focal|broad)", all_objects, value = TRUE)
  save(list = objects_to_save, file = result.dir)
}














