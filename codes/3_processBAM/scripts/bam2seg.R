## this is a script that compiles the pipeline of executing the process of processing BAM files for HH lung and HH ova BAM files.  
## Section 0 environment =====================================================================================
## Library packages ------------------------------------------------------------------------------------------
library(dplyr)
library(magrittr)

## Setting paths ---------------------------------------------------------------------------------------------
BASH.DIR = "~/MRes_project_1/codes/3_processBAM/1_bash/"
RCODE.DIR = "~/MRes_project_1/codes/3_processBAM/2_Rscripts/"
REF.DIR = "~/MRes_project_1/codes/3_processBAM/4_docs/"
RESULT.DIR = "~/MRes_project_1/codes/3_processBAM/6_results/"


## Section 1 quality control =================================================================================
## compile customised fasta files ----------------------------------------------------------------------------
cmd = paste0("qsub ", BASH.DIR, "compile_fasta.sh")
system(cmd) 

## create fasta extension files ------------------------------------------------------------------------------
cmd = paste0("qsub ", BASH.DIR, "fasta_ext_files.sh")
system(cmd)

## create interval lists 
cmd = paste0("qsub ", BASH.DIR, "interval_lists.sh")
system(cmd)

## determining the mapping quality of bam files --------------------------------------------------------------
cmd = paste0("qsub ", BASH.DIR, "map_quality.sh")
system(cmd)

## determining the WES quality by hs metrics ----------------------------------------------------------------- 
cmd = paste0("qsub", BASH.DIR, "hs_metrics.sh")
system(cmd)

## collect GCbias --------------------------------------------------------------------------------------------
cmd = paste0("qsub", BASH.DIR, "GcBias.sh")
system(cmd)

## Section 2 snp pile up =====================================================================================
## check vcf file order --------------------------------------------------------------------------------------
cmd = paste0("qsub", BASH.DIR, "get_order_files.sh")
system(cmd)

## rewrite a new vcf order -----------------------------------------------------------------------------------
tumour_order <- read.table("~/MRes_project_1/codes/3_processBAM/4_docs/tumourChrIDs_order.txt")
vcf_new_order <- order$X970N[c(1:22, 24:25)]
write.table(vcf_new_order, file = "~/MRes_project_1/codes/3_processBAM/4_docs/vcf_file/vcf_new_order.txt", 
            row.names = FALSE, quote = FALSE, col.names = FALSE)

## order vcf file --------------------------------------------------------------------------------------------
cmd = paste0("qsub", BASH.DIR, "resort_vcf_file.sh")
system(cmd)

## perform snp pileup for hh_lung ----------------------------------------------------------------------------
## The snp-pileup programme is already done on local computer and then uploaded to HPC following the command 
## provided by FACETS github page. 
## snp-pileup <vcf file> <output file> <sequence files...>
## Usage of snp-pileup requires a VCF file and one (or multiple) sequence files containing DNA. 
## The sequence files should be in the BAM format, and both the VCF and all sequence files must be sorted. 
## A suitable option for VCF is one from NCBI. ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/VCF/00-common_all.vcf.gz     
## Use a VCF file that is consistent with the genome build used for aligning the sequencing data. 
## Available snp/genome build versions are in directories human_9606_b149_GRCh38p7, human_9606_b150_GRCh38p7.
cmd = paste0("qsub", BASH.DIR, "snp_pileup.sh")
system(cmd)

## reorder the snp-pileup outputs before executing facets -----------------------------------------------------
setwd("~/MRes_project_1/codes/3_processBAM/4_docs/")
sample_ids = read.table("sample_ids.txt") %>% unlist
cancer = c("ova", "lung")

for(i in cancer){
  for(sample in sample_ids){
    ## output file: 
    setwd("~/MRes_project_1/docs/HH_", i,"/facets/facet_input") ## in this case is for ovarian i = ova or lung
    output <- paste0(sample, "_ordered.csv")
    gzip_output <- paste0(sample, "_ordered.csv.gz")
    
    if (!file.exists(output)) {
      ## reorder the snp_pileup files to 1, 2, 3... 22, X, Y format 
      temp <- read.csv(paste0("~/MRes_project_1/docs/HH_ova/facets/facet_input/", sample, ".csv")) ## read in files 
      chrom_order <- c(as.character(1:22), "X", "Y")  ## set up the order 
      temp <- temp[order(match(temp$Chromosome, chrom_order), temp$Position), ] ## reorder 
      
      write.csv(temp, paste0("~/MRes_project_1/docs/HH_ova/facets/facet_input/", sample, "_ordered.csv"), 
                quote = FALSE, col.names = TRUE, row.names = FALSE)  ## save the files 
    } else if (!file.exists(gzip_output)){
      temp <- read.csv(paste0("~/MRes_project_1/docs/HH_ova/facets/facet_input/", sample, "_ordered.csv"))
      write.csv(temp, gzfile(paste0("~/MRes_project_1/docs/HH_ova/facets/facet_input/", sample, "_ordered.csv.gz")), 
                row.names = FALSE)
    }
  }
}

## perform facets --------------------------------------------------------------------------------------------
## we have the choice to either run the original FACETS package or to run the FACETS wrapper script package 
## in this analysis, we will run the facet_wrapper script 
cmd = paste0("qsub", BASH.DIR, "facet_wrapper.sh")
system(cmd)

## this will result in a list of segmentation files 
print("Processing from sorted duplication removed BAM files to individual segmentation file done!")