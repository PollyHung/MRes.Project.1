## extract the files that is not paired 
# Assuming you're working in the directory containing the files:
file_list <- list.files("/rds/general/user/ph323/ephemeral/HH_ova/Alignment_ov")  
file_list <- file_list[grep(".bam$", file_list)]
ids <- sub("(X[0-9]+)[NT]_sorted_nodup.bam", "\\1", file_list)
table_ids <- table(ids)

## Process the bam files 
# Assuming you have two vectors with matched normal and tumor bam file names
normal_bams <- list.files(pattern = "N_sorted_nodup.bam$", path = "/rds/general/user/ph323/ephemeral/HH_ova/Alignment_ov")
tumor_bams <- list.files(pattern = "T_sorted_nodup.bam$", path = "/rds/general/user/ph323/ephemeral/HH_ova/Alignment_ov")
# remove the single bam files 
single <- file_list[ids %in% names(table_ids[table(ids) == 1])]
tumor_bams <- tumor_bams[!tumor_bams %in% single]
# sort the bam files to match them up 
normal_bams <- sort(normal_bams)
tumor_bams <- sort(tumor_bams)

## other parameters: 
# path to snp-pileup 
snp_pileup_path <- "/rds/general/user/ph323/home/MRes_project_1/docs/snp_pileup/codes/snp-pileup"
# snp-pileup parameters 
params <- "-q15 -Q20 -P100 -r25,0"
# vcf file 
vcf <- "/rds/general/user/ph323/home/MRes_project_1/docs/snp_pileup/sorted_vcf_file.vcf.gz"

# Loop over the BAM files
for (i in 59:98) { ## seq_along(normal_bams)
  ## set working directory 
  setwd("/rds/general/user/ph323/ephemeral/HH_ova/Alignment_ov")
  
  ## index the correct bam file pairs 
  normal_bam <- normal_bams[i]
  tumor_bam <- tumor_bams[i]
  
  ## define the output file name based on the BAM file names
  output_file <- paste0("~/MRes_project_1/docs/HH_ova/facet_input/", 
                        sub("N_sorted_nodup.bam$", "", normal_bams[i]), ".csv")
  
  ## construct the command
  cmd <- paste(snp_pileup_path, params, vcf, output_file, normal_bam, tumor_bam)
  
  #cmd <- sprintf('%s %s %s %s %s %s', snp_pileup_path, params, vcf, output_file, normal_bams[i], tumor_bams[i])
  
  # Run the command
  system(cmd, intern = TRUE)
  
  print(paste0("finish piling ", sub("N_sorted_nodup.bam$", "", normal_bams[i])))
}

