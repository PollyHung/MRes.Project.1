## This is a script from a series of scripts to test all kinds of possibilities for facets package 
## This script focuses on cval 















## calculation 
for(cancer in sample_ids){
  
  datafile <- paste0(cancer, "_ordered.csv")
  
  ## preprocessing 
  set.seed(1234)
  rcmat <-  readSnpMatrix(datafile)
  xx <- preProcSample(rcmat, gbuild = "hg38")
  oo <- procSample(xx, cval = CVAL, dipLogR = xx$dipLogR)
  
  ## fitting 
  fit <- emcncf(oo)
  
  ## plotting the results and diagnostic plot 
  plotSample(x = oo, emfit = fit) ## result 
  logRlogORspider(oo$out, oo$dipLogR) ## diagnostic plot 
  
  ## extract segmentation file and assign them to environment 
  cncf <- fit$cncf
  
  ## assign the data 
  xx_name <- paste0("preproc_", cancer, "_", CVAL)
  oo_name <- paste0("proc_", cancer, "_", CVAL)
  fit_name <- paste0("fit_", cancer, "_", CVAL)
  cncf_name <- paste0("segFile_", cancer, "_", CVAL)
  
  assign(xx_name, xx, envir = .GlobalEnv)
  assign(oo_name, oo, envir = .GlobalEnv)
  assign(fit_name, fit, envir = .GlobalEnv)
  assign(cncf_name, cncf, envir = .GlobalEnv)
}






library(ABSOLUTE)
filelist=list.files(pattern="*_sorted_nodup.cns")

for(i in filelist) {
  setwd("I:/BRC_backup/cnv_ov")
  x=read.delim(i,sep="\t")
  x=x[,c(1,2,3,7,5)]
  names(x)=c("Chromosome","Start","End","Num_Probes","Segment_Mean")
  x=x[x[,1] %in% c(1:22),]
  write.table(x, "x", sep = "\t", row.names = FALSE)
  RunAbsolute("x", sigma.p=0, max.sigma.h=0.015, min.ploidy=0.95, max.ploidy=10, primary.disease="ov", platform="Illumina_WES", sample.name= i, results.dir="C:/Users/OPCML/Downloads/CNV_ABSOLUTE", max.as.seg.count=1500, max.non.clonal=0.05, max.neg.genome=0.005, copy_num_type="total", maf.fn=NULL, min.mut.af=NULL, output.fn.base=NULL, verbose=FALSE)
  setwd("C:/Users/OPCML/Downloads/CNV_ABSOLUTE")
  load(paste(i,".ABSOLUTE.RData",sep=""))
  write.csv(seg.dat$segtab,file=paste(i,".csv",sep=""))
}












