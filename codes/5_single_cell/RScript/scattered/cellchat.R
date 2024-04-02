cellchat_amp <- readRDS("~/MRes_project_1/docs/SC_SEQ/ova/GSE180661/cellChat/cancer_myeloid/cellchat_amp.rds")
cellchat_wt <- readRDS("~/MRes_project_1/docs/SC_SEQ/ova/GSE180661/cellChat/cancer_myeloid/cellchat_wt.rds")
## cancer communication with myeloid 
cellchat_amp@netP$pathways
cellchat_wt@netP$pathways


cellchat_amp_2 <- readRDS("~/MRes_project_1/docs/SC_SEQ/ova/GSE180661/cellChat/myeloid_tcells/cellchat_amp.rds")
cellchat_wt_2 <- readRDS("~/MRes_project_1/docs/SC_SEQ/ova/GSE180661/cellChat/myeloid_tcells/cellchat_wt.rds")
## cancer communication with t cells 
cellchat_amp_2@netP$pathways
cellchat_wt_2@netP$pathways
