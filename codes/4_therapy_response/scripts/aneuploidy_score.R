## This is a script calculating the aneuploidy score provided by methods detailed in https://doi.org/10.1016/j.ccell.2018.03.007


## To calculate arm level events and aneuploidy score using FACETS algorithm. For each sample, we determined the likeliest ploidy and 
##    absolute total copy number of each segment in the genome. 
## Each segment was designated as amplified, deleted, or neutral based on whether its copy number was greater, smaller, or equal to the sample's 
##    ploidy (round to nearest integer). 
## For amplifications and deletions separately, segments were joined until either the entire chromosome was considered altered or more than 20% 
##    of the genomic length between the start and ends were not altered in the same direction 
## For each combination of arm/chromosome and direction of alteration within each dataset, the start coordinates, end coordinates, and percentage 
##    length of the longest joined segment were clustered across samples using a Gaussian Mixture Model. 
## The optimal clustering solution between 2-9 clusters inclusive was chosen based on the lowest BIC (Bayesian information criterion).
## Tumors in clusters whose mean fraction altered in either specific direction was >=80% were considered ‘‘aneuploid.’’ Tumors altered 


## Section 0: Environment ==========================================================
setwd("~/MRes_project_1/docs/")

library(mclust)
library(dplyr)


DATASET = "HH_lung" ## or HH_ova



## Section 1: Load in the FACETS ploidy and purity scores.===========================  
sample_ids = list.files("~/MRes_project_1/docs/HH_lung/facets/facet_cval_50/")
sample = read.table("~/MRes_project_1/docs/HH_lung/facets/facet_cval_50/MDL15-1352/MDL15-1352.txt", 
                    sep = "\t", header = TRUE)

for(i in 2:length(sample_ids)){
    ## set working directory 
    setwd("~/MRes_project_1/docs/HH_lung/facets/facet_cval_50/")
    
    ## read in the txt 
    temp = read.table(paste0(sample_ids[i], "/", sample_ids[i], ".txt"), 
                      sep = "\t", header = TRUE)
    sample = rbind(sample, temp)
}

hisens <- sample %>% filter(run_type %in% "hisens")
hisens <- hisens[c("sample", "purity", "ploidy", "dipLogR", "fraction_cna")]


segmentation <- read.table("~/MRes_project_1/docs/HH_lung/facets/seg_cval_50/_hisens_diplogR.adjusted.seg", 
                           sep = "\t")
colnames(segmentation) <- c("sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")
segmentation <- segmentation %>%
  left_join(hisens, by = "sample") %>%
  mutate(fraction = fraction_cna)


## Section 2: 
X <- as.matrix(data[c("fraction", "start_coordinate", "end_coordinate")])

# Find the optimal number of clusters using BIC
BIC_values <- mclustBIC(X)
optimal_model <- Mclust(X, BIC_values)
n <- optimal_model$G # Optimal number of components

# Fit the GMM model
g <- Mclust(X, G=n)

# Assign (cluster, probability) to each tumor based on max likelihood
data$GMM_cluster <- paste("k", g$classification, sep="")
data$clusterprob <- apply(g$posterior, 1, max)

# Calculate means for each cluster
means <- tapply(data$fraction, g$classification, mean)
data$fractionmean <- data$GMM_cluster %>% sapply(function(k) means[as.numeric(sub("k", "", k))])

# Similarly calculate means for start and end
means_start <- tapply(data$start_coordinate, g$classification, mean)
data$startmean <- data$GMM_cluster %>% sapply(function(k) means_start[as.numeric(sub("k", "", k))])

means_end <- tapply(data$end_coordinate, g$classification, mean)
data$endmean <- data$GMM_cluster %>% sapply(function(k) means_end[as.numeric(sub("k", "", k))])

# Define aneuploidy
data$aneu <- ifelse(data$fractionmean >= 0.8, 1, ifelse(data$fraction <= 0.2, 0, NA))

# 'data' now contains the additional columns with cluster assignments and means









