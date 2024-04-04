## This code analyses the gistic output   
rm(list = ls())
source("~/MRes_project_1/Codes/3_therapy_response/HH_ova/immune_association/association_function.R")
association(prefix = "o", 
            result.dir = "~/MRes_project_1/GISTIC2/output/hh_ova/", 
            sample.ids.dir = "~/MRes_project_1/Codes/3_therapy_response/HH_ova/immune_association/sample_ids.txt", 
            clinical.data.dir = "~/MRes_project_1/docs/HH_ova/clinical_data/ho_clinical.txt", 
            survival.data.dir <- "~/MRes_project_1/docs/HH_ova/clinical_data/ho_survival.txt", 
            focal.dir = "~/MRes_project_1/docs/HH_ova/gistic/facets_seg/cval_50/hisens_adj/all_lesions.conf_95.txt", 
            broad.dir = "~/MRes_project_1/docs/HH_ova/gistic/facets_seg/cval_50/hisens_adj/broad_values_by_arm.txt", 
            p.val = 1)

rm(list = ls())
source("~/MRes_project_1/Codes/3_therapy_response/HH_ova/immune_association/association_function.R")
association(prefix = "l", 
            result.dir = "~/MRes_project_1/GISTIC2/output/hh_lung/", 
            sample.ids.dir = "~/MRes_project_1/Codes/3_therapy_response/HH_lung/facets/sample_ids_2.txt", 
            clinical.data.dir = "~/MRes_project_1/docs/HH_ova/clinical_data/hl_clinical.txt", 
            survival.data.dir <- "~/MRes_project_1/docs/HH_ova/clinical_data/hl_survival.txt", 
            focal.dir = "~/MRes_project_1/docs/HH_lung/gistic/gistic_output/hisens_adj/all_lesions.conf_95.txt", 
            broad.dir = "~/MRes_project_1/docs/HH_lung/gistic/gistic_output/hisens_adj/broad_values_by_arm.txt", 
            p.val = 1)

rm(list = ls())
source("~/MRes_project_1/Codes/3_therapy_response/HH_ova/immune_association/association_function.R")
association(prefix = "l", 
            result.dir = "~/MRes_project_1/GISTIC2/output/cbioportal_mel/", 
            sample.ids.dir = "~/MRes_project_1/GISTIC2/sample_ids/cbioportal_mel_ids.txt", 
            clinical.data.dir = "~/MRes_project_1/GISTIC2/clinical_docs/cbioportal_mel_pheno.txt", 
            survival.data.dir <- "~/MRes_project_1/GISTIC2/clinical_docs/cbioportal_mel_surv.txt", 
            focal.dir = "~/MRes_project_1/GISTIC2/cbioportal/mel/all_lesions.conf_95.txt", 
            broad.dir = "~/MRes_project_1/GISTIC2/cbioportal/mel/broad_values_by_arm.txt", 
            p.val = 1)


rm(list = ls())
source("~/MRes_project_1/Codes/3_therapy_response/HH_ova/immune_association/association_function.R")
association(prefix = "l", 
            result.dir = "~/MRes_project_1/GISTIC2/output/cbioportal_lung/", 
            sample.ids.dir = "~/MRes_project_1/GISTIC2/sample_ids/cbioportal_lung_ids.txt", 
            clinical.data.dir = "~/MRes_project_1/GISTIC2/clinical_docs/cbioportal_lung_pheno.txt", 
            survival.data.dir <- "~/MRes_project_1/GISTIC2/clinical_docs/cbioportal_lung_surv.txt", 
            focal.dir = "~/MRes_project_1/GISTIC2/cbioportal/lung/all_lesions.conf_95.txt", 
            broad.dir = "~/MRes_project_1/GISTIC2/cbioportal/lung/broad_values_by_arm.txt", 
            p.val = 1)


