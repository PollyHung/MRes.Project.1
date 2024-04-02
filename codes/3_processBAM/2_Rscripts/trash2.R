library(dplyr)
library(magrittr)
library(httr)
library(Rsamtools)
library(facets)
library(ASCAT)
library(data.tree)
library(ABSOLUTE)
library(caroline)


PRINT <- TRUE
EXECUTE <- FALSE
DATA <- "/rds/general/user/ph323/ephemeral/HH_ova/Alignment_ov"
RESULT <- "~/MRes_project_1/docs/HH_ova/facet_input"
EXAMPLE <- "~/MRes_project_1/docs/HH_ova/example"
GATK <- "~/MRes_project_1/otherCodes/gatk-4.2.5.0/gatk" ## this is also a file exec
SNP_PILEUP <- "~/MRes_project_1/Codes/3_therapy_resp/pileup_exe/snp-pileup" ## this is a file exec
VCF_FILE <- "/MRes_project_1/Codes/3_therapy_resp/vcf_file"

## extract the files that is not paired 
# Assuming you're working in the directory containing the files:
file_list <- list.files(DATA)  
file_list <- file_list[grep(".bam$", file_list)]
ids <- sub("(X[0-9]+)[NT]_sorted_nodup.bam", "\\1", file_list)
table_ids <- table(ids)

## unpaired files 
unpaired_ova_bam <- file_list[ids %in% names(table_ids[table(ids) == 1])]
## paired files 
paired_ova_bam <- file_list[ids %in% names(table_ids[table(ids) == 2])]

## Process the bam files 
paired_tumour <- paired_ova_bam[grep("T_sorted_nodup.bam$", paired_ova_bam)] %>% sort ## get all the tumour samples 
paired_normal <- paired_ova_bam[grep("N_sorted_nodup.bam$", paired_ova_bam)] %>% sort ## get all the normal samples 
unpaired_tumour <- unpaired_ova_bam[grep("sorted_nodup.bam", unpaired_ova_bam)]
unpaired_normal <- rep("X991N_sorted_nodup.bam", length(unpaired_tumour))

## combine the files 
normal <- c(unpaired_normal, paired_normal) ## the first 21 samples will not be matched 
tumour <- c(unpaired_tumour, paired_tumour) 

## check if the paired bam files are matched up by randomly sample a few  
set.seed(1234)
random <- sample(1:119, 12, replace = FALSE)
for(i in random){
  print(c(i, normal[i], tumour[i]))
}


## reorder facet input -------------------------
## sample ids 
sample_ids <- list.files("~/MRes_project_1/docs/HH_ova/facets/facet_input")
sample_ids <- grep("_order", sample_ids, value = TRUE, invert = TRUE)
sample_ids <- gsub(".csv", "", sample_ids)

list.files(path = "~/MRes_project_1/docs/HH_ova/facets/facet_input", pattern = "_ordered.csv") %>% length %>% print

## run facets ------------------
unpaired <- gsub("T_sorted_nodup.bam", "", unpaired_tumour)

## execute facets 
setwd("~/MRes_project_1/docs/HH_ova/example/")
datafile <- paste0("~/MRes_project_1/docs/HH_ova/facets/facet_input/X2316_ordered.csv")
## preprocessing 
set.seed(1234)
rcmat <-  readSnpMatrix(datafile)
xx <- preProcSample(rcmat, unmatched=TRUE, gbuild = "hg38")
oo <- procSample(xx, cval = 400)

## fitting 
fit <- emcncf(oo)

## plotting the results and diagnostic plot 
plotSample(x = oo, emfit = fit) ## result 
logRlogORspider(oo$out, oo$dipLogR) ## diagnostic plot 

write.csv(fit$cncf, "X2316_400.csv", quote = FALSE, row.names = FALSE)
save.image("X2316.FACETs.RData")

## ABSOLUTE ---------------------------------
absolute <- read.table("~/MRes_project_1/docs/HH_ova/facets/facet_cval_50/X1032/X1032_hisens_diplogR.adjusted.seg", 
                       header = TRUE) ## read in table 

## edit and write the table 
absolute <- absolute[2:6]
names(absolute) <- c("Chromosome","Start","End","Num_Probes","Segment_Mean")
absolute <- absolute[absolute[ , 1] %in% c(1:22),]
write.table(absolute, "~/MRes_project_1/docs/HH_ova/example/absolute", sep = "\t", row.names = FALSE)

setwd("~/MRes_project_1/docs/HH_ova/example/")
## run absolute
RunAbsolute("~/MRes_project_1/docs/HH_ova/example/absolute", 
            sigma.p=0, 
            max.sigma.h=0.015, min.ploidy=0.95, max.ploidy=10, 
            primary.disease="ov", platform="Illumina_WES", 
            sample.name= "X1032", 
            results.dir="~/MRes_project_1/docs/HH_ova/example/X1032", 
            max.as.seg.count=5000, max.non.clonal=0.05, max.neg.genome=0.005, 
            copy_num_type="total", 
            maf.fn=NULL, min.mut.af=NULL, 
            output.fn.base=NULL, verbose=FALSE)

CreateReviewObject( 
  obj.name = "X1032", 
  absolute.files = "X991_50.ABSOLUTE.RData", 
  indv.results.dir = "~/MRes_project_1/docs/HH_ova/example/X1032/ploidy", 
  copy_num_type = "total", 
  plot.modes = TRUE
)


## set working directory 
load("~/MRes_project_1/docs/HH_ova/example/X1032/X1032.ABSOLUTE.RData")
write.csv(seg.dat$segtab, file = "X991_50_abs.csv")

Seg.CN <- seg.dat$segtab 
Seg.CN$Sample <- "X991"
Seg.CN$Seg.CN <- log2(Seg.CN$copy_num)-1
Seg.CN <- Seg.CN[c("Sample", "Chromosome", "Start.bp", "End.bp", "n_probes", "Seg.CN")]
colnames(Seg.CN) <- c("Sample", "Chromosome", "Start_position", "End_position", "Num_markers", "Seg.CN")


## CreateReviewObject ---------------------
list.files("~/MRes_project_1/docs/HH_ova/example/")
CreateReviewObject( 
  obj.name = "DoAbsolute", 
  absolute.files = "98.ABSOLUTE.RData", 
  indv.results.dir = "~/MRes_project_1/docs/HH_ova/example/X991_50", 
  copy_num_type = "total", 
  plot.modes = TRUE
)

ExtractReviewedResults( 
  reviewed.pp.calls.fn = "DoAbsolute.PP-calls_tab.txt", 
  analyst.id = "wsx", 
  modes.fn = "DoAbsolute.PP-modes.data.RData", 
  out.dir.base = "~/MRes_project_1/docs/HH_ova/example/X991_50", 
  obj.name = "DoAbsolute", 
  copy_num_type = "total", 
)

path_to_file = system.file("extdata", "ABSOLUTE_1.0.6.tar.gz", package = "DoAbsolute", mustWork = T, 
                           lib.loc = "/rds/general/user/ph323/home/anaconda3/lib/R/library")
install.packages(path_to_file, repos = NULL, type="source")




