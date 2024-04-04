## HH_ova
setwd("~/MRes_project_1/docs/HH_ova/facets/facet_cval_50/")
samples <- list.files()

samples <- samples[ !samples %in% c('X1159', 'X1876')]

temp1 <- readRDS("X991/X991_hisens.rds")
temp1 <- temp1$segs
temp1$sample <- "X991"

dataframe <- matrix(NA, ncol = 18, nrow = 1) %>% as.data.frame()
colnames(dataframe) <- colnames(temp1)

for(i in samples){ 
  object <- readRDS(paste0(i, "/", i, "_hisens.rds"))
  object_1 <- object$segs
  object_1$sample <- i
  
  dataframe <- rbind(dataframe, object_1)
}

dataframe <- dataframe[2:nrow(dataframe), ]
write.table(dataframe, "~/MRes_project_1/docs/HH_ova/facets/complete_seg.txt", 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)



## HH_lung
setwd("~/MRes_project_1/docs/HH_lung/facets/facet_cval_50/")
samples <- list.files()

dataframe2 <- matrix(NA, ncol = 18, nrow = 1) %>% as.data.frame()
colnames(dataframe2) <- colnames(temp1)

for(i in samples){ 
  object <- readRDS(paste0(i, "/", i, "_hisens.rds"))
  object_1 <- object$segs
  object_1$sample <- i
  
  dataframe2 <- rbind(dataframe2, object_1)
}

dataframe2 <- dataframe2[2:nrow(dataframe2), ]
write.table(dataframe2, "~/MRes_project_1/docs/HH_lung/facets/complete_seg.txt", 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)



