## This code analyses the gistic output   
rm(list = ls())
source("~/MRes_project_1/codes/4_therapy_response/association_function.R")
association(result.dir = "~/MRes_project_1/docs/HH_ova/gistic/output_2/", 
            sample.ids.dir = "~/MRes_project_1/codes/3_processBAM/4_docs/sample_ids/ovarian_ids.txt", 
            clinical.data.dir = "~/MRes_project_1/docs/HH_ova/clinical_data/ho_clinical.txt", 
            survival.data.dir <- "~/MRes_project_1/docs/HH_ova/clinical_data/ho_survival.txt", 
            focal.dir = "~/MRes_project_1/docs/HH_ova/gistic/input_2/all_lesions.conf_99.txt", 
            broad.dir = "~/MRes_project_1/docs/HH_ova/gistic/input_2/broad_values_by_arm.txt")

rm(list = ls())
source("~/MRes_project_1/Codes/3_therapy_response/HH_ova/immune_association/association_function.R")
association(prefix = "l", 
            result.dir = "~/MRes_project_1/docs/HH_lung/gistic/output/", 
            sample.ids.dir = "~/MRes_project_1/Codes/3_therapy_response/HH_lung/facets/sample_ids_2.txt", 
            clinical.data.dir = "~/MRes_project_1/docs/HH_ova/clinical_data/hl_clinical.txt", 
            survival.data.dir <- "~/MRes_project_1/docs/HH_ova/clinical_data/hl_survival.txt", 
            focal.dir = "~/MRes_project_1/docs/HH_lung/gistic/gistic_output/hisens_adj/all_lesions.conf_95.txt", 
            broad.dir = "~/MRes_project_1/docs/HH_lung/gistic/gistic_output/hisens_adj/broad_values_by_arm.txt", 
            p.val = 1)

rm(list = ls())
source("~/MRes_project_1/codes/4_therapy_response/association_function.R")
association(result.dir = "~/MRes_project_1/docs/HH_lung/gistic/output_2/", 
            sample.ids.dir = "~/MRes_project_1/codes/3_processBAM/4_docs/sample_ids/lung_ids.txt", 
            clinical.data.dir = "~/MRes_project_1/docs/HH_ova/clinical_data/hl_clinical.txt", 
            survival.data.dir = "~/MRes_project_1/docs/HH_ova/clinical_data/hl_survival.txt", 
            focal.dir = "~/MRes_project_1/docs/HH_lung/gistic/input_2/all_lesions.conf_99.txt", 
            broad.dir = "~/MRes_project_1/docs/HH_lung/gistic/input_2/broad_values_by_arm.txt")

temp <- read.table("~/MRes_project_1/docs/HH_lung/gistic/input_2/broad_values_by_arm.txt", sep = "\t",header = TRUE)

