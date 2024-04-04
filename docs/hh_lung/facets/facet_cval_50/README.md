cell_types <- list(
    ## skin cells 
    keratinocytes = c("KRT14", "KRT1", "DMKN", "KRT10"),
    melanocytes = c("DCT", "TYRP1", "PMEL", "MLANA"),
    eccrine_gland_cells = c("PIP", "DCD", "MUCL1"),
    fibroblasts = c("DCN", "COL1A1", "COL1A2"), 
    
    ## nerve cells 
    nerve_cells = c("MPZ", "PLP1", "S100B"),
    schwann = c("SOX10", "S100B", "MPZ"), 
    
    ## immune cells 
    #t_cells = c("PTPRC","CD3G", "CD3E", "CD4"), 
    #cd4_t = c("PTPRC","CD3G", "CD3E", "CD4"), 
    #cd8_t = c("CD8A", "CD8B", "CTLA4"), 
    cd4_t_rest = c("CCR7", "SELL", "TCF7"), 
    cd4_t_act = c("IL2", "TNF", "IL4R"), 
    t_resident_mem = c("CXCR6", "ITGA113"),
    t_reg = c("FOXP3", "IL2RA", "CTLA4"), 
    cd8_t_TRM.TEM = c("CCL5", "GZMB", "GZMK", "ITGA1"), 
    cd8_t_TRM.TEM.act = c("IFNG", "CCL4", "CCL3"), 
    cd8_t_effector_cells = c("PRF1", "NKG7"), 
    dendritic = c("CD86", "HLA-DR", "CD45", "CD11c"), 
    mast_cells = c("CPA3", "TPSAB1", "CTSG"), 
    macrophage = c("CD68", "IL1B", "LGMN", "MS4A4A"), ## "IL1B", "CCL4", "THBD", "CD11C"
    #macrophage_2 = c("CD68", "CD163", "LGMN", "MS4A4A"), ## "F13A1", "FOLR2", "CXCL8", "CD11B"
    myeloid_cells = c("CD74", "HLA-DRA", "HLA-DPB1"),
    #lactoglobulin = c("CD207", "CLDN1", "PLEK2", "CD1B", "CD1A", "CD1C", "CXCL8", "CD11C"), 

    ## others 
    endothelial_cells = c("PECAM1", "CDH5", "CLDN5", "CD34"),
    mesenchymal = c("PDGFRB", "LUM", "COL1A1"))