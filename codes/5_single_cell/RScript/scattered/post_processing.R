## this is to do post-processing on inferCNV object using infercnv_postprocessing_v2.py 
## from https://github.com/ding-lab/infer_cnv_postprocesssing/blob/main/infercnv_postprocessing_v2.py

## Base Controls 
CANCER <- "ova" ## could also be lung
STUDY <- "lambrechtslab" ## other studies such as bischoff
SAMPLE <- "scrSOL001"
BASE <- "~/MRes_project_1/docs/SC_SEQ"

raw.dir <- paste(BASE, CANCER, "raw", STUDY, "InferCNV/",sep = "/")
result.dir <- paste(BASE, CANCER, "processed", STUDY, "InferCNV/",sep = "/")

if(!dir.exists(result.dir)){ ## if directory does not exist, make directory 
  dir.create(result.dir)
}

## file paths: 
setwd(result.dir)
py <- "~/MRes_project_1/codes/6_single_cell/Python/infercnv_postprocessing_v2.py"
observation <- paste0(INFERCNV, "result/", SAMPLE, "/infercnv.observations.txt")
reference <- paste0(INFERCNV, "result/", SAMPLE, "/infercnv.references.txt")

cmd <- paste("python", py, observation, reference, SAMPLE)
system(cmd)













