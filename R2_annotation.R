rm(list=ls())
pkgs <- c('Seurat','tidyr','plyr','dplyr','stringr','SingleR','ggsci','dior','tibble','qs','BiocParallel','scGate','Matrix','SingleCellExperiment','scran','parallel','scGate','burgertools')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)
options(max.print = 10000)

# Loading dataset
SKCM_Becker <- qread('./PanCancer_ICI/data/SKCM_Becker/seu_r1.qs')
BRCA_Bassez1 <- qread('./PanCancer_ICI/data/BRCA_Bassez1/seu_r1.qs')
BRCA_Bassez2 <- qread('./PanCancer_ICI/data/BRCA_Bassez2/seu_r1.qs')
TNBC_Zhang <- qread('./PanCancer_ICI/data/TNBC_Zhang/seu_r1.qs')
BCC_Yost <- qread('./PanCancer_ICI/data/BCC_Yost/seu_r1.qs')
BCC_Yost1 <- qread('./PanCancer_ICI/data/BCC_Yost/seu_r1.qs') %>% 
  subset(subset = patient %in% c('BCC/SCC_Yost_su009', 'BCC/SCC_Yost_su012'), invert = T)
BCC_Yost2 <- qread('./PanCancer_ICI/data/BCC_Yost/seu_r1.qs') %>% 
  subset(subset = patient %in% c('BCC/SCC_Yost_su009', 'BCC/SCC_Yost_su012')) 
SCC_Yost <- qread('./PanCancer_ICI/data/SCC_Yost/seu_r1.qs')
SCC_Yost$dataset <- 'SCC_Yost'
HNSC_IMCISION <- qread('./PanCancer_ICI/data/HNSC_IMCISION/seu_r1.qs')
HNSC_Luoma <- qread('./PanCancer_ICI/data/HNSC_Luoma/seu_r1.qs')
NSCLC_Liu <- qread('./PanCancer_ICI/data/NSCLC_Liu/seu_r1.qs')
CRC_Li <- qread('./PanCancer_ICI/data/CRC_Li/seu_r1.qs')
PCa_Hawley <- qread('./PanCancer_ICI/data/PCa_Hawley/seu_r1.qs')
TNBC_Shiao <- qread('./PanCancer_ICI/data/TNBC_Shiao/seu_r1.qs')
HNSC_Franken <- qread('./PanCancer_ICI/data/HNSC_Franken/seu_r1.qs')
CRC_Chen <- qread('./PanCancer_ICI/data/CRC_Chen/seu_r1.qs')

# Run SingleR with mutiple reference datasets
# arguments:
# seu: seurat object 
# ref_list: reference (list)
# major: major cell types (T/NK, Myeloids)
# type_name: for file name of the predicjtion results
# datset: for file name of the prediction results
SingleRMultiRef <- function(seu, ref_list, if.subset = F, major, type_name, n_cores = 50){
  label_list <- lapply(ref_list, function(x){
    label <- x$label
    return(label)
  })
  # Create the formatted input string
  # Set the number of combinations
  num_combinations <- length(ref_list)
  # Generate random combinations
  random_combinations <- replicate(num_combinations, paste(sample(letters, 4), collapse = ""))
  input_string <- paste0(random_combinations, "=", "ref_list[[", seq_along(ref_list), "]]", collapse = ", ")
  final_input <- paste("ref = list(", input_string, ")", sep = "")
  if (if.subset == T){
    seu <- subset(seu, subset = celltype_major %in% major)
  }
  pred <- seu %>%
    as.SingleCellExperiment() %>%
    logNormCounts() %>%
    SingleR(ref = eval(parse(text = final_input)), labels = label_list, assay.type.test=1, BPPARAM=MulticoreParam(n_cores))
  print('Predictiion done!')
  qsave(pred, paste0("./PanCancer_ICI/data/", unique(seu$dataset), "/", type_name, ".qs"))
  return(pred)
}
# run SingleR with one reference dataset
SingleROneRef <- function(seu, ref, label, major, type_name, n_cores = 50){
  pred <- seu %>% 
    subset(subset = celltype_major %in% major) %>% 
    as.SingleCellExperiment() %>% 
    logNormCounts() %>% 
    SingleR(ref = ref, labels = label, assay.type.test=1, BPPARAM=MulticoreParam(50))
  qsave(pred, paste0("./PanCancer_ICI/data/", unique(seu$dataset), "/", type_name, ".qs"))
}

# datasets <- list(SKCM_Becker, BRCA_Bassez1, BRCA_Bassez2, TNBC_Zhang, BCC_Yost1, HNSC_IMCISION, HNSC_Luoma, NSCLC_Liu, CRC_Li, PCa_Hawley)

# T_NK
ref_T <- qread("./PanCancer_ICI/data/Ref_SingleR/T_ref.qs")
ref_CD4 <- qread("./PanCancer_ICI/data/Ref_SingleR/T_CD4_ref.qs")
ref_CD8 <- qread("./PanCancer_ICI/data/Ref_SingleR/T_CD8_ref.qs")
ref_main_NK <- qread("./PanCancer_ICI/data/Ref_SingleR/NK_main_ref.qs")
ref_fine_NK <- qread("./PanCancer_ICI/data/Ref_SingleR/NK_fine_ref.qs")

hierSingleRMultiRef_T_NK <- function(seu, main_ref_list, ref_list_cd4, ref_list_cd8, ref_list_nk_main, ref_list_nk_fine, major = c("CD8+ T-cells","CD4+ T-cells","NK cells")){
  print(unique(seu$dataset))
  seu <- subset(seu, subset = celltype_major %in% major) 
  seu <- NormalizeData(seu) %>%
    FindVariableFeatures()%>%
    ScaleData() %>%
    RunPCA() %>% 
    RunUMAP(dims = 1:20) %>%
    FindNeighbors(dims = 1:20) 
  print('Gating NK cells')
  NK_KLRD1 <- gating_model(name = "NK_KLRD1", signature = c("KLRD1","CD3D-"))
  NK_NCAM1 <- gating_model(name = "NK_NCAM1", signature = c("NCAM1","CD3D-"))
  seu_1<- scGate(data = seu, model = NK_KLRD1)
  seu_2<- scGate(data = seu, model = NK_NCAM1)
  barcode_NK <- unique(c(colnames(seu_1)[seu_1$is.pure == 'Pure'],
                         colnames(seu_2)[seu_2$is.pure == 'Pure']))
  print(paste0(length(barcode_NK), ' NK cells were detected.'))
  write.csv(barcode_NK, paste0('./PanCancer_ICI/data/', unique(seu$dataset), '/barcode_NK.csv'))
  print('Gating gdT cells')
  gdT <- gating_model(name = "gdT", signature = c("CD3D","TRGC2","TRDV2","TRGV9","TRGV10","TRDC","CD8A-","CD4-"))
  seu_gdT<- scGate(data = seu, model = gdT)
  barcode_gdT <- colnames(seu_gdT)[seu_gdT$is.pure == 'Pure']
  # barcode_gdT <- setdiff(barcode_gdT, barcode_NK)
  print(paste0(length(barcode_gdT), ' gdT cells were detected'))
  write.csv(barcode_gdT, paste0('./PanCancer_ICI/data/', unique(seu$dataset), '/barcode_gdT.csv'))
  scGate_models_DB <- get_scGateDB()
  seu_cd4cd8 <- seu[, !colnames(seu) %in% c(barcode_NK, barcode_gdT)]
  print('Gating CD4+T cells')
  seu_cd4 <- scGate(seu_cd4cd8, model = scGate_models_DB$human$generic$CD4T)
  barcode_cd4 <- colnames(seu_cd4)[seu_cd4$is.pure == 'Pure']
  print('Gating CD8+T cells')
  seu_cd8 <- scGate(seu_cd4cd8, model = scGate_models_DB$human$generic$CD8T)
  barcode_cd8 <- colnames(seu_cd8)[seu_cd8$is.pure == 'Pure']
  dual <- intersect(barcode_cd4, barcode_cd8)
  barcode_cd4 <- setdiff(barcode_cd4, dual)
  barcode_cd8 <- setdiff(barcode_cd8, dual)
  print('CD4/CD8')
  # Unresolved CD4/CD8 cells
  pred_T <- SingleRMultiRef(seu[, !colnames(seu) %in% c(barcode_NK, barcode_gdT, barcode_cd4, barcode_cd8)], ref_list = main_ref_list, type_name = 'T')
  barcode_CD4 <- c(rownames(pred_T)[pred_T$pruned.labels == 'CD4'], barcode_cd4)
  barcode_CD8 <- c(rownames(pred_T)[pred_T$pruned.labels == 'CD8'], barcode_cd8)
  print('CD4')
  pred_CD4 <- SingleRMultiRef(seu = seu[, barcode_CD4], ref_list = ref_list_cd4, type_name = 'CD4')
  print('CD8')
  pred_CD8 <- SingleRMultiRef(seu = seu[, barcode_CD8], ref_list = ref_list_cd8, type_name = 'CD8')
  print('NK main')
  pred_NK_main <- SingleRMultiRef(seu = seu[, barcode_NK], ref_list = ref_list_nk_main, type_name = 'NK_main')
  print('NK fine')
  pred_NK_fine <- SingleRMultiRef(seu = seu[, barcode_NK], ref_list = ref_list_nk_fine, type_name = 'NK_fine')
}
datasets <- list(SKCM_Becker, BRCA_Bassez1, BRCA_Bassez2, TNBC_Zhang, BCC_Yost1, HNSC_IMCISION, HNSC_Luoma, NSCLC_Liu, CRC_Li, PCa_Hawley)
mclapply(datasets, function(seu){
  hierSingleRMultiRef_T_NK(seu, main_ref_list = ref_T, 
                          ref_list_cd4 = ref_CD4, ref_list_cd8 = ref_CD8, 
                          ref_list_nk_main = ref_main_NK, ref_list_nk_fine = ref_fine_NK,
                          major = c("CD8+ T-cells","CD4+ T-cells","NK cells"))}, mc.cores = 100)
hierSingleRMultiRef_T_NK(CRC_Chen, main_ref_list = ref_T, 
                         ref_list_cd4 = ref_CD4, ref_list_cd8 = ref_CD8, 
                         ref_list_nk_main = ref_main_NK, ref_list_nk_fine = ref_fine_NK,
                         major = c("CD8+ T-cells","CD4+ T-cells","NK cells"))

# T cells only (CD45+CD3+ sorted)
hierSingleRMultiRef_T <- function(seu, main_ref_list, ref_list_cd4, ref_list_cd8){
  print(unique(seu$dataset))
  seu <- NormalizeData(seu) %>%
    FindVariableFeatures()%>%
    ScaleData() %>%
    RunPCA() %>% 
    RunUMAP(dims = 1:20) %>%
    FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 1)
  print('Gating gdT cells')
  gdT <- gating_model(name = "gdT", signature = c("CD3D","TRGC2","TRDV2","TRGV9","TRGV10","TRDC","CD8A-","CD4-"))
  seu_gdT<- scGate(data = seu, model = gdT)
  barcode_gdT <- colnames(seu_gdT)[seu_gdT$is.pure == 'Pure']
  # barcode_gdT <- setdiff(barcode_gdT, barcode_NK)
  print(paste0(length(barcode_gdT), ' gdT cells were detected'))
  write.csv(barcode_gdT, paste0('./PanCancer_ICI/data/', unique(seu$dataset), '/barcode_gdT_sorted.csv'))
  scGate_models_DB <- get_scGateDB()
  seu_cd4cd8 <- seu[,!colnames(seu) %in% barcode_gdT]
  print('Gating CD4+T cells')
  seu_cd4 <- scGate(seu_cd4cd8, model = scGate_models_DB$human$generic$CD4T)
  barcode_cd4 <- colnames(seu_cd4)[seu_cd4$is.pure == 'Pure']
  print('Gating CD8+T cells')
  seu_cd8 <- scGate(seu_cd4cd8, model = scGate_models_DB$human$generic$CD8T)
  barcode_cd8 <- colnames(seu_cd8)[seu_cd8$is.pure == 'Pure']
  dual <- intersect(barcode_cd4, barcode_cd8)
  barcode_cd4 <- setdiff(barcode_cd4, dual)
  barcode_cd8 <- setdiff(barcode_cd8, dual)
  print('CD4/CD8')
  # Unresolved CD4/CD8 cells
  pred_T <- SingleRMultiRef(seu[, !colnames(seu) %in% c(barcode_gdT, barcode_cd4, barcode_cd8)], ref_list = main_ref_list, type_name = 'T') # For CD45CD3 sorted (SCC_Yost)
  barcode_CD4 <- c(rownames(pred_T)[pred_T$pruned.labels == 'CD4'], barcode_cd4)
  barcode_CD8 <- c(rownames(pred_T)[pred_T$pruned.labels == 'CD8'], barcode_cd8)
  print('CD4')
  pred_CD4 <- SingleRMultiRef(seu = seu[, barcode_CD4], ref_list = ref_list_cd4, type_name = 'sorted_CD4')
  print('CD8')
  pred_CD8 <- SingleRMultiRef(seu = seu[, barcode_CD8], ref_list = ref_list_cd8, type_name = 'sorted_CD8')
}
hierSingleRMultiRef_T(BCC_Yost2, main_ref_list = ref_T, ref_list_cd4 = ref_CD4, ref_list_cd8 = ref_CD8)
hierSingleRMultiRef_T(SCC_Yost, main_ref_list = ref_T, ref_list_cd4 = ref_CD4, ref_list_cd8 = ref_CD8)

# Myeloids (doi.org/10.1016/j.cell.2021.01.010)
ref_Myeloids <- qread("./PanCancer_ICI/data/Ref_SingleR/Myeloids_major_ref.qs")
ref_cdc2 <- qread("./PanCancer_ICI/data/Ref_SingleR/Myeloids_cdc2_ref.qs")
ref_mono <- qread("./PanCancer_ICI/data/Ref_SingleR/Myeloids_mono_ref.qs")
ref_macro <- qread("./PanCancer_ICI/data/Ref_SingleR/Myeloids_macro_ref.qs")
datasets <- list(SKCM_Becker, BRCA_Bassez1, BRCA_Bassez2, TNBC_Zhang, TNBC_Shiao, BCC_Yost1, HNSC_IMCISION, HNSC_Luoma, CRC_Li, PCa_Hawley)

hierSingleRMultiRef_Myeloids <- function(seu, main_ref_list, ref_list_cdc2, ref_list_mono, ref_list_macro){
  pred <- SingleRMultiRef(seu = seu, ref_list = main_ref_list, if.subset = T, major = c('DC', 'Monocytes', 'Macrophages'), type_name = 'Myeloids_major')
  print('Major done!')
  barcode_cdc2 <- rownames(pred)[pred$pruned.labels == 'cDC2']
  barcode_mono <- rownames(pred)[pred$pruned.labels == 'Mono']
  barcode_macro <- rownames(pred)[pred$pruned.labels == 'Macro']
  print('cDC2s')
  pred_cdc2 <- SingleRMultiRef(seu = seu[, barcode_cdc2], ref_list = ref_list_cdc2, type_name = 'cDC2')
  print('Monocytes')
  pred_mono <- SingleRMultiRef(seu = seu[, barcode_mono], ref_list = ref_list_mono, type_name = 'Mono')
  print('Macrophages')
  pred_macro <- SingleRMultiRef(seu = seu[, barcode_macro], ref_list = ref_list_macro, type_name = 'Macro')
}
mclapply(datasets, function(seu){
  hierSingleRMultiRef_Myeloids(seu, main_ref_list = ref_Myeloids, ref_list_cdc2 = ref_cdc2, 
                               ref_list_mono = ref_mono, ref_list_macro = ref_macro)}, mc.cores = 100)
hierSingleRMultiRef_Myeloids(CRC_Chen, main_ref_list = ref_Myeloids, 
                             ref_list_mono = ref_mono, ref_list_macro = ref_macro, ref_list_cdc2 = ref_cdc2)

# pan-cancer B (doi.org/10.1126/science.adj4857)
ref_Bplasma <- qread("./PanCancer_ICI/data/Ref_SingleR/Pan_B_major.qs")
ref_b <- qread("./PanCancer_ICI/data/Ref_SingleR/B.qs")
ref_gcb <- qread("./PanCancer_ICI/data/Ref_SingleR/GCB.qs")
ref_plasma <- qread("./PanCancer_ICI/data/Ref_SingleR/Plasma.qs")
datasets <- list(SKCM_Becker, BRCA_Bassez1, BRCA_Bassez2, TNBC_Zhang, TNBC_Shiao, BCC_Yost1, HNSC_IMCISION, HNSC_Luoma, HNSC_Franken, CRC_Li, PCa_Hawley)
hierSingleRMultiRef_Bplasma <- function(seu, main_ref_list, ref_list_b, ref_list_gcb, ref_list_plasma){
  pred <- SingleRMultiRef(seu = seu, ref_list = main_ref_list, if.subset = T, major = c('B-cells', 'Plasma'), type_name = 'Bplasma_major')
  print('Major done!')
  pred_b <- rownames(pred)[pred$pruned.labels == 'B']
  pred_gcb <- rownames(pred)[pred$pruned.labels == 'GCB']
  pred_plasma <- rownames(pred)[pred$pruned.labels == 'Plasma']
  print('B')
  pred_b <- SingleRMultiRef(seu = seu[, pred_b], ref_list = ref_list_b, type_name = 'pan-B')
  if (length(pred_gcb)>0){
    print('GCB')
    pred_gcb <- SingleRMultiRef(seu = seu[, pred_gcb], ref_list = ref_list_gcb, type_name = 'GCB')
  }
  if (length(pred_plasma)>0){
  print('Plasma')
  pred_plasma <- SingleRMultiRef(seu = seu[, pred_plasma], ref_list = ref_list_plasma, type_name = 'Plasma')
  }
}
mclapply(datasets, function(seu){
  hierSingleRMultiRef_Bplasma(seu, main_ref_list = ref_Bplasma, ref_list_b = ref_b, 
                               ref_list_gcb = ref_gcb, ref_list_plasma = ref_plasma)}, mc.cores = 100)
hierSingleRMultiRef_Bplasma(CRC_Chen, main_ref_list = ref_Bplasma, ref_list_b = ref_b, 
                             ref_list_gcb = ref_gcb, ref_list_plasma = ref_plasma)

# Pan-cancer Endothelial cells (doi.org/10.1093/nsr/nwae231)
ref_endo_major <- qread('./PanCancer_ICI/data/Ref_SingleR/Endo_major.qs')
ref_endo_vascular <- qread('./PanCancer_ICI/data/Ref_SingleR/Endo_vascular.qs')
hierSingleRMultiRef_Endo <- function(seu, main_ref_list, ref_list_vas){
  pred <- SingleRMultiRef(seu = seu, ref_list = main_ref_list, if.subset = T, major = c('Endothelial cells'), type_name = 'pan-Endo')
  print('Major done!')
  pred_vas <- rownames(pred)[pred$pruned.labels == 'vascular']
  print('Vascular')
  pred_vas <- SingleRMultiRef(seu = seu[, pred_vas], ref_list = ref_list_vas, type_name = 'pan-Endo_vas')
}
datasets <- list(SKCM_Becker, BRCA_Bassez1, BRCA_Bassez2, BCC_Yost1, CRC_Li, CRC_Chen, PCa_Hawley)
mclapply(datasets, function(seu){
  hierSingleRMultiRef_Endo(seu, main_ref_list = ref_endo_major, ref_list_vas = ref_endo_vascular)}, 
  mc.cores = 100)
hierSingleRMultiRef_Endo(CRC_Chen, main_ref_list = ref_endo_major, ref_list_vas = ref_endo_vascular)
# Markers to validate
genes_to_check = c('GJA4','GJA5','FBLN5',
                   'CA4','CD36','RGCC',
                   'PROX1','LYVE1','CCL2',
                   'COL4A1','KDR','ESM1',
                   'ACKR1','SELP','CLU')

# Endothelial cells
ref_Endo <- qread("./PanCancer_ICI/data/Ref_SingleR/Endo_ref.qs")
datasets <- list(SKCM_Becker, BRCA_Bassez1, BRCA_Bassez2, BCC_Yost1, CRC_Li, CRC_Chen, PCa_Hawley)
lapply(datasets, function(seu){SingleROneRef(seu, ref = ref_Endo, label = ref_Endo$label, major = 'Endothelial cells', type_name = 'Endo')})
SingleROneRef(HNSC_Franken, ref = ref_Endo, label = ref_Endo$label, major = 'Endothelial cells', type_name = 'Endo')

# CAF (doi.org/10.1038/s41467-024-48310-4)
ref_CAF <- qread("./PanCancer_ICI/data/Ref_SingleR/CAF_ref.qs")
datasets <- list(SKCM_Becker, BRCA_Bassez1, BRCA_Bassez2, BCC_Yost1, CRC_Li, PCa_Hawley)
lapply(datasets, function(seu){
  print(unique(seu$dataset))
  SingleRMultiRef(seu, ref_list = ref_CAF, major = c('Fibroblasts','Myocytes'), if.subset = T, type_name = 'CAF')})
SingleRMultiRef(CRC_Chen, ref_list = ref_CAF, major = c('Fibroblasts','Myocytes'), if.subset = T, type_name = 'CAF')

add_to_seu_r2 <- function(cd45_sorted = F, dataset, majortype, celltype_major_remove, label_overlapped, celltype_r2_remove){
  seu <- qread(paste0('./PanCancer_ICI/data/', dataset,'/seu_r1.qs'))
  seu$celltype_major <- as.character(seu$celltype_major)
  seu <- subset(seu, subset = celltype_major %in% celltype_major_remove, invert = T)
  pred_list <- lapply(paste0('./PanCancer_ICI/data/', dataset, '/', majortype, ".qs"), qread)
  names(pred_list) <- majortype
  pred_concat <- do.call(rbind, lapply(pred_list, function(df) {
    data.frame(cell = rownames(df), celltype_r2 = df$pruned.labels)
  }))
  barcode_gdT <- read.csv(paste0('./PanCancer_ICI/data/', dataset, '/barcode_gdT.csv')) 
  pred_gdT <- data.frame(cell = barcode_gdT[,'x'], celltype_r2 = 'gdT')
  pred_concat <- rbind(pred_concat, pred_gdT)
  pred_concat <- filter(pred_concat, !celltype_r2 %in% label_overlapped)
  seu$celltype_r2 <- pred_concat$celltype_r2[match(colnames(seu), pred_concat$cell)]
  seu$celltype_r2[seu$celltype_major == 'pDC'] <- 'pDC'
  seu$celltype_r2[seu$celltype_major == 'Mast'] <- 'Mast'
  seu$celltype_major[seu$celltype_r2 %in% c(unique(pred_list$CD4$labels), unique(pred_list$CD8$labels), 'gdT')] <- 'T_cells'
  seu$celltype_major[seu$celltype_r2 %in% c(unique(pred_list$NK_main$labels))] <- 'NK_cells'
  print(table(seu$celltype_r2, seu$celltype_major, useNA = 'ifany'))
  seu <- seu[,!is.na(seu$celltype_r2)]
  seu$celltype_main <- 'celltype'
  seu$celltype_main[str_detect(seu$celltype_r2, 'CD4')] <- 'CD4+T'
  seu$celltype_main[str_detect(seu$celltype_r2, 'CD8')] <- 'CD8+T'
  seu$celltype_main[str_detect(seu$celltype_r2, 'gdT')] <- 'CD8+T'
  seu$celltype_main[seu$celltype_r2 %in% c(unique(pred_list$NK_main$pruned.labels))] <- 'NK'
  seu$celltype_main[seu$celltype_r2 %in% c(unique(pred_list$`pan-B`$pruned.labels), unique(pred_list$GCB$pruned.labels), unique(pred_list$Plasma$pruned.labels))] <- 'B'
  seu$celltype_main[str_detect(seu$celltype_r2, 'Plasma')] <- 'Plasma'
  seu$celltype_main[str_detect(seu$celltype_major, 'pDC')] <- 'pDC'
  seu$celltype_main[str_detect(seu$celltype_r2, 'cDC')] <- 'cDC'
  seu$celltype_main[str_detect(seu$celltype_r2, 'Mono')] <- 'Mono'
  seu$celltype_main[str_detect(seu$celltype_r2, 'Macro')] <- 'Macro'
  seu$celltype_main[seu$celltype_r2 == 'Mast'] <- 'Mast'
  if (cd45_sorted == F){
    seu$celltype_main[seu$celltype_r2 %in% c(unique(pred_list$`pan-Endo`$pruned.labels), unique(pred_list$`pan-Endo_vas`$pruned.labels))] <- 'Endo'
    seu$celltype_main[str_detect(seu$celltype_r2, 'CAF') | str_detect(seu$celltype_r2, 'Myofibroblast')] <- 'CAF'
  }
  print(table(seu$celltype_main, seu$celltype_major, useNA = 'ifany'))
  seu <- subset(seu, subset = celltype_r2 == celltype_r2_remove, invert = T)
  return(seu)
}

# TME
## SKCM_Becker BRCA_Bassez1 BRCA_Bassez2 HNSC_Franken CRC_Li CRC_Chen PCa_Hawley
dataset <- 'SKCM_Becker'
seu <- add_to_seu_r2(dataset = dataset, 
                     majortype = c('CD4','CD8','NK_main','pan-B','GCB','Plasma','Myeloids_major','cDC2','Mono','Macro','pan-Endo','pan-Endo_vas','CAF'),
                     celltype_major_remove = c('Melanoma', 'Epithelial cells'),
                     label_overlapped = c('cDC2','Macro', 'Mono','vascular'),
                     celltype_r2_remove = 'CD56highCD16high')
# keep only site-matched samples
seu <- subset(seu, subset = patient %in% c('SKCM_Becker_P2', 'SKCM_Becker_P5', 'SKCM_Becker_P6','SKCM_Becker_P7','SKCM_Becker_P9','SKCM_Becker_P10','SKCM_Becker_P11','SKCM_Becker_P12'))
qsave(seu, paste0('./PanCancer_ICI/data/', dataset, '/seu_r2.qs'))

dataset <- 'BRCA_Bassez1'
seu <- add_to_seu_r2(dataset = dataset, 
                     majortype = c('CD4','CD8','NK_main','pan-B','GCB','Plasma','Myeloids_major','cDC2','Mono','Macro','pan-Endo','pan-Endo_vas','CAF'),
                     celltype_major_remove = c('Melanoma', 'Epithelial cells'),
                     label_overlapped = c('cDC2','Macro', 'Mono','vascular'),
                     celltype_r2_remove = 'CD56highCD16high')
qsave(seu, paste0('./PanCancer_ICI/data/', dataset, '/seu_r2.qs'))

dataset <- 'BRCA_Bassez2'
seu <- add_to_seu_r2(dataset = dataset, 
                     majortype = c('CD4','CD8','NK_main','pan-B','GCB','Plasma','Myeloids_major','cDC2','Mono','Macro','pan-Endo','pan-Endo_vas','CAF'),
                     celltype_major_remove = c('Melanoma', 'Epithelial cells'),
                     label_overlapped = c('cDC2','Macro', 'Mono','vascular'),
                     celltype_r2_remove = 'CD56highCD16high')
qsave(seu, paste0('./PanCancer_ICI/data/', dataset, '/seu_r2.qs'))

dataset <- 'HNSC_Franken'
seu <- add_to_seu_r2(dataset = dataset, 
                     majortype = c('CD4','CD8','NK_main','pan-B','GCB','Plasma','Myeloids_major','cDC2','Mono','Macro','pan-Endo','pan-Endo_vas','CAF'),
                     celltype_major_remove = c('Melanoma', 'Epithelial cells'),
                     label_overlapped = c('cDC2','Macro', 'Mono','vascular'),
                     celltype_r2_remove = 'CD56highCD16high')
qsave(seu, paste0('./PanCancer_ICI/data/', dataset, '/seu_r2.qs'))

dataset <- 'CRC_Li'
seu <- add_to_seu_r2(dataset = dataset, 
                     majortype = c('CD4','CD8','NK_main','pan-B','GCB','Plasma','Myeloids_major','cDC2','Mono','Macro','pan-Endo','pan-Endo_vas','CAF'),
                     celltype_major_remove = c('Melanoma', 'Epithelial cells'),
                     label_overlapped = c('cDC2','Macro', 'Mono','vascular'),
                     celltype_r2_remove = 'CD56highCD16high')
qsave(seu, paste0('./PanCancer_ICI/data/', dataset, '/seu_r2.qs'))

dataset <- 'CRC_Chen'
seu <- add_to_seu_r2(dataset = dataset, 
                     majortype = c('CD4','CD8','NK_main','pan-B','GCB','Plasma','Myeloids_major','cDC2','Mono','Macro','pan-Endo','pan-Endo_vas','CAF'),
                     celltype_major_remove = c('Melanoma', 'Epithelial cells'),
                     label_overlapped = c('cDC2','Macro', 'Mono','vascular'),
                     celltype_r2_remove = 'CD56highCD16high')
qsave(seu, paste0('./PanCancer_ICI/data/', dataset, '/seu_r2.qs'))

dataset <- 'PCa_Hawley'
cell_subtype <- c('CD4','CD8','NK_main','Myeloids_major','cDC2','Mono','Macro','pan-Endo','pan-Endo_vas','CAF','pan-B','Plasma')
PCa_Hawley$celltype_major <- as.character(PCa_Hawley$celltype_major)
PCa_Hawley <- subset(PCa_Hawley, subset = celltype_major %in% c('Epithelial cells'), invert = T)
pred_list <- lapply(paste0('./PanCancer_ICI/data/', dataset, '/', cell_subtype, ".qs"), qread)
names(pred_list) <- cell_subtype
PCa_Hawley$celltype_r2 <- 'Unresolved'
pred_concat <- do.call(rbind, lapply(pred_list, function(df) {
  data.frame(cell = rownames(df), celltype_r2 = df$pruned.labels)
}))
pred_concat <- filter(pred_concat, !celltype_r2 %in% c('cDC2', 'Macro', 'Mono','vascular'))
PCa_Hawley$celltype_r2 <- pred_concat$celltype_r2[match(colnames(PCa_Hawley), pred_concat$cell)]
PCa_Hawley$celltype_r2[PCa_Hawley$celltype_major == 'pDC'] <- 'pDC'
PCa_Hawley$celltype_r2[PCa_Hawley$celltype_major == 'Mast'] <- 'Mast'
PCa_Hawley$celltype_major[PCa_Hawley$celltype_r2 %in% c(unique(pred_list$CD4$labels), unique(pred_list$CD8$labels), 'gdT')] <- 'T_cells'
PCa_Hawley$celltype_major[PCa_Hawley$celltype_r2 %in% c(unique(pred_list$NK_main$labels))] <- 'NK_cells'
table(PCa_Hawley$celltype_r2, PCa_Hawley$celltype_major, useNA = 'ifany')
PCa_Hawley <- PCa_Hawley[,!is.na(PCa_Hawley$celltype_r2)]
PCa_Hawley$celltype_main <- 'celltype'
PCa_Hawley$celltype_main[str_detect(PCa_Hawley$celltype_r2, 'CD4')] <- 'CD4+T'
PCa_Hawley$celltype_main[str_detect(PCa_Hawley$celltype_r2, 'CD8')] <- 'CD8+T'
PCa_Hawley$celltype_main[str_detect(PCa_Hawley$celltype_r2, 'gdT')] <- 'CD8+T'
PCa_Hawley$celltype_main[str_detect(PCa_Hawley$celltype_major, 'NK')] <- 'NK'
PCa_Hawley$celltype_main[PCa_Hawley$celltype_r2 %in% c(unique(pred_list$`pan-B`$pruned.labels), unique(pred_list$GCB$pruned.labels), unique(pred_list$Plasma$pruned.labels))] <- 'B'
PCa_Hawley$celltype_main[str_detect(PCa_Hawley$celltype_r2, 'Plasma')] <- 'Plasma'
PCa_Hawley$celltype_main[str_detect(PCa_Hawley$celltype_major, 'pDC')] <- 'pDC'
PCa_Hawley$celltype_main[str_detect(PCa_Hawley$celltype_r2, 'cDC')] <- 'cDC'
PCa_Hawley$celltype_main[str_detect(PCa_Hawley$celltype_r2, 'Mono')] <- 'Mono'
PCa_Hawley$celltype_main[str_detect(PCa_Hawley$celltype_r2, 'Macro')] <- 'Macro'
PCa_Hawley$celltype_main[PCa_Hawley$celltype_r2 %in% c(unique(pred_list$`pan-Endo`$pruned.labels), unique(pred_list$`pan-Endo_vas`$pruned.labels))] <- 'Endo'
PCa_Hawley$celltype_main[str_detect(PCa_Hawley$celltype_r2, 'CAF') | str_detect(PCa_Hawley$celltype_r2, 'Myofibroblast')] <- 'CAF'
PCa_Hawley$celltype_main[str_detect(PCa_Hawley$celltype_major, 'Mast')] <- 'Mast'
table(PCa_Hawley$celltype_main, PCa_Hawley$celltype_major, useNA = 'ifany')
PCa_Hawley <- subset(PCa_Hawley, subset = celltype_r2 =='CD56highCD16high', invert = T)
qsave(PCa_Hawley, './PanCancer_ICI/data/PCa_Hawley/seu_r2.qs')

# CD45+Sorted 
# TNBC_Zhang TNBC_Shiao HNSC_IMCISION HNSC_Luoma 
dataset <- 'TNBC_Zhang'
seu <- add_to_seu_r2(dataset = dataset, 
                     majortype = c('CD4','CD8','NK_main','Myeloids_major','cDC2','Mono','Macro','pan-B','Plasma','GCB'),
                     celltype_major_remove = NA,
                     label_overlapped = c('cDC2','Macro', 'Mono'),
                     celltype_r2_remove = 'CD56highCD16high')
qsave(seu, paste0('./PanCancer_ICI/data/', dataset, '/seu_r2.qs'))

dataset <- 'TNBC_Shiao'
seu <- add_to_seu_r2(dataset = dataset, 
                     majortype = c('CD4','CD8','NK_main','Myeloids_major','cDC2','Mono','Macro','pan-B','Plasma','GCB'),
                     celltype_major_remove = NA,
                     label_overlapped = c('cDC2','Macro', 'Mono'),
                     celltype_r2_remove = 'CD56highCD16high')
qsave(seu, paste0('./PanCancer_ICI/data/', dataset, '/seu_r2.qs'))

dataset <- 'HNSC_IMCISION'
seu <- add_to_seu_r2(dataset = dataset, 
                     majortype = c('CD4','CD8','NK_main','Myeloids_major','cDC2','Mono','Macro','pan-B','Plasma','GCB'),
                     celltype_major_remove = NA,
                     label_overlapped = c('cDC2','Macro', 'Mono'),
                     celltype_r2_remove = 'CD56highCD16high')
qsave(seu, paste0('./PanCancer_ICI/data/', dataset, '/seu_r2.qs'))

dataset <- 'HNSC_Luoma'
seu <- add_to_seu_r2(dataset = dataset, 
                     majortype = c('CD4','CD8','NK_main','Myeloids_major','cDC2','Mono','Macro','pan-B','Plasma','GCB'),
                     celltype_major_remove = NA,
                     label_overlapped = c('cDC2','Macro', 'Mono'),
                     celltype_r2_remove = 'CD56highCD16high')
qsave(seu, paste0('./PanCancer_ICI/data/', dataset, '/seu_r2.qs'))

# Mixed BCC_Yost
dataset <- 'BCC_Yost'
cell_subtype <- c('CD4','CD8','NK_main','Myeloids_major','cDC2','Mono','Macro','pan-Endo','pan-Endo_vas','CAF','pan-B','Plasma','GCB')
BCC_Yost$celltype_major <- as.character(BCC_Yost$celltype_major)
BCC_Yost <- subset(BCC_Yost, subset = celltype_major %in% c('Epithelial cells', 'Melanocytes'), invert = T)
pred_list <- lapply(paste0('./PanCancer_ICI/data/', dataset, '/', cell_subtype, ".qs"), qread)
names(pred_list) <- cell_subtype
pred_concat <- do.call(rbind, lapply(pred_list, function(df) {
  data.frame(cell = rownames(df), celltype_r2 = df$pruned.labels)
}))
barcode_gdT <- read.csv(paste0('./PanCancer_ICI/data/', dataset, '/barcode_gdT.csv')) 
pred_gdT <- data.frame(cell = barcode_gdT[,'x'], celltype_r2 = 'gdT')
pred_concat <- rbind(pred_concat, pred_gdT)
pred_concat <- filter(pred_concat, !celltype_r2 %in% c('cDC2', 'Macro', 'Mono','vascular'))
# su009 su012
cell_subtype <- c('sorted_CD4','sorted_CD8')
pred_list2 <- lapply(paste0('./PanCancer_ICI/data/', dataset, '/', cell_subtype, ".qs"), qread)
names(pred_list2) <- cell_subtype
pred_concat2 <- do.call(rbind, lapply(pred_list2, function(df) {
  data.frame(cell = rownames(df), celltype_r2 = df$pruned.labels)
}))
barcode_gdT2 <- read.csv(paste0('./PanCancer_ICI/data/', dataset, '/barcode_gdT_sorted.csv')) 
pred_gdT2 <- data.frame(cell = barcode_gdT2[,'x'], celltype_r2 = 'gdT')
pred_concat <- purrr::reduce(list(pred_concat, pred_gdT2, pred_concat2), rbind)
BCC_Yost$celltype_r2 <- pred_concat$celltype_r2[match(colnames(BCC_Yost), pred_concat$cell)]
BCC_Yost$celltype_r2[BCC_Yost$celltype_major == 'pDC'] <- 'pDC'
BCC_Yost$celltype_r2[BCC_Yost$celltype_major == 'Mast'] <- 'Mast'
BCC_Yost$celltype_major[BCC_Yost$celltype_r2 %in% c(unique(pred_list$CD4$labels), unique(pred_list$CD8$labels), unique(pred_list2[[1]]$pruned.labels), unique(pred_list2[[2]]$pruned.labels), 'gdT')] <- 'T_cells'
BCC_Yost$celltype_major[BCC_Yost$celltype_r2 %in% c(unique(pred_list$NK_main$labels))] <- 'NK_cells'
table(BCC_Yost$celltype_r2, BCC_Yost$celltype_major, useNA = 'ifany')
BCC_Yost <- BCC_Yost[,!is.na(BCC_Yost$celltype_r2)]
BCC_Yost$celltype_main <- 'celltype'
BCC_Yost$celltype_main[str_detect(BCC_Yost$celltype_r2, 'CD4')] <- 'CD4+T'
BCC_Yost$celltype_main[str_detect(BCC_Yost$celltype_r2, 'CD8')] <- 'CD8+T'
BCC_Yost$celltype_main[str_detect(BCC_Yost$celltype_r2, 'gdT')] <- 'CD8+T'
BCC_Yost$celltype_main[str_detect(BCC_Yost$celltype_major, 'NK')] <- 'NK'
BCC_Yost$celltype_main[BCC_Yost$celltype_r2 %in% c(unique(pred_list$`pan-B`$pruned.labels), unique(pred_list$GCB$pruned.labels), unique(pred_list$Plasma$pruned.labels))] <- 'B'
BCC_Yost$celltype_main[str_detect(BCC_Yost$celltype_r2, 'Plasma') & !BCC_Yost$patient %in% c('su009', 'su012')] <- 'Plasma'
BCC_Yost$celltype_main[str_detect(BCC_Yost$celltype_major, 'pDC')] <- 'pDC'
BCC_Yost$celltype_main[str_detect(BCC_Yost$celltype_r2, 'cDC')] <- 'cDC'
BCC_Yost$celltype_main[str_detect(BCC_Yost$celltype_r2, 'Mono')] <- 'Mono'
BCC_Yost$celltype_main[str_detect(BCC_Yost$celltype_r2, 'Macro')] <- 'Macro'
BCC_Yost$celltype_main[BCC_Yost$celltype_r2 %in% c(unique(pred_list$`pan-Endo`$pruned.labels), unique(pred_list$`pan-Endo_vas`$pruned.labels))] <- 'Endo'
BCC_Yost$celltype_main[str_detect(BCC_Yost$celltype_major, 'Mast')] <- 'Mast'
BCC_Yost$celltype_main[str_detect(BCC_Yost$celltype_r2, 'CAF') | str_detect(BCC_Yost$celltype_r2, 'Myofibroblast')] <- 'CAF'
table(BCC_Yost$patient, BCC_Yost$celltype_main, useNA = 'ifany')
BCC_Yost <- BCC_Yost[,!BCC_Yost$celltype_main == 'celltype']
BCC_Yost <- subset(BCC_Yost, subset = celltype_r2 == 'CD56highCD16high', invert = T)
BCC_Yost$dataset <- 'BCC_Yost'
qsave(BCC_Yost, './PanCancer_ICI/data/BCC_Yost/seu_r2.qs')

# SCC_Yost
dataset <- 'SCC_Yost'
cell_subtype <- c('CD4','CD8')
SCC_Yost$celltype_major <- as.character(SCC_Yost$celltype_major)
pred_list <- lapply(paste0('./PanCancer_ICI/data/', dataset, '/sorted_', cell_subtype, ".qs"), qread)
names(pred_list) <- cell_subtype
SCC_Yost$celltype_r2 <- 'Unresolved'
pred_concat <- do.call(rbind, lapply(pred_list, function(df) {
  data.frame(cell = rownames(df), celltype_r2 = df$pruned.labels)
}))
barcode_gdT <- read.csv(paste0('./PanCancer_ICI/data/', dataset, '/barcode_gdT_sorted.csv')) 
pred_gdT <- data.frame(cell = barcode_gdT[,'x'], celltype_r2 = 'gdT')
pred_concat <- rbind(pred_concat, pred_gdT)
mapping <- setNames(pred_concat$celltype_r2, pred_concat$cell)
SCC_Yost$celltype_r2 <- mapping[colnames(SCC_Yost)]
table(SCC_Yost$celltype_r2, SCC_Yost$celltype_major, useNA = 'ifany')
SCC_Yost <- SCC_Yost[,!is.na(SCC_Yost$celltype_r2)]
SCC_Yost$celltype_main <- 'celltype'
SCC_Yost$celltype_main[str_detect(SCC_Yost$celltype_r2, 'CD4')] <- 'CD4+T'
SCC_Yost$celltype_main[str_detect(SCC_Yost$celltype_r2, 'CD8')] <- 'CD8+T'
SCC_Yost$celltype_main[str_detect(SCC_Yost$celltype_r2, 'gdT')] <- 'CD8+T'
table(SCC_Yost$celltype_r2, SCC_Yost$celltype_main, useNA = 'ifany')
qsave(SCC_Yost, './PanCancer_ICI/data/SCC_Yost/seu_r2.qs')

dataset <- 'NSCLC_Liu'
cell_subtype <- c('CD4','CD8','NK_main')
NSCLC_Liu$celltype_major <- as.character(NSCLC_Liu$celltype_major)
pred_list <- lapply(paste0('./PanCancer_ICI/data/', dataset, '/', cell_subtype, ".qs"), qread)
names(pred_list) <- cell_subtype
NSCLC_Liu$celltype_r2 <- 'Unresolved'
pred_concat <- do.call(rbind, lapply(pred_list, function(df) {
  data.frame(cell = rownames(df), celltype_r2 = df$pruned.labels)
}))
mapping <- setNames(pred_concat$celltype_r2, pred_concat$cell)
NSCLC_Liu$celltype_r2 <- mapping[colnames(NSCLC_Liu)]
NSCLC_Liu$celltype_major[NSCLC_Liu$celltype_r2 %in% c(unique(pred_list[[1]]$pruned.labels), unique(pred_list[[2]]$pruned.labels))] <- 'T_cells'
NSCLC_Liu$celltype_major[NSCLC_Liu$celltype_r2 %in% c(unique(pred_list[[3]]$pruned.labels))] <- 'NK_cells'
table(NSCLC_Liu$celltype_r2, NSCLC_Liu$celltype_major, useNA = 'ifany')
NSCLC_Liu <- NSCLC_Liu[,!is.na(NSCLC_Liu$celltype_r2)]
NSCLC_Liu <- subset(NSCLC_Liu, subset = celltype_major == 'NK_cells', invert = T)
NSCLC_Liu$celltype_main <- 'celltype'
NSCLC_Liu$celltype_main[str_detect(NSCLC_Liu$celltype_r2, 'CD4')] <- 'CD4+T'
NSCLC_Liu$celltype_main[str_detect(NSCLC_Liu$celltype_r2, 'CD8')] <- 'CD8+T'
table(NSCLC_Liu$celltype_main, NSCLC_Liu$celltype_major, useNA = 'ifany')
qsave(NSCLC_Liu, './PanCancer_ICI/data/NSCLC_Liu/seu_r2.qs')

# Adjust labels
datasets <- c('SKCM_Becker', 'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Zhang', 'BCC_Yost', 'SCC_Yost', 'HNSC_IMCISION', 'HNSC_Luoma', 'NSCLC_Liu', 'CRC_Li', 'CRC_Chen', 'PCa_Hawley', 'TNBC_Shiao', 'HNSC_Franken')
lapply(datasets, function(dataset){
  print(dataset)
  seu <- qread(paste0('./PanCancer_ICI/data/', dataset, '/seu_r2.qs')) 
  # Modification
  # CD4+T
  seu$celltype_r2[seu$celltype_r2 == 'CD4_Tn_IL7R-'] <- 'CD4_Prolif'
  seu$celltype_r2[seu$celltype_r2 == 'CD4_Tn_ADSL'] <- 'CD4_Naive'
  seu$celltype_r2[seu$celltype_r2 == 'CD4_pre-Tfh_CXCR5+'] <- 'CD4_pre-Tfh_CXCR5'
  seu$celltype_r2[seu$celltype_r2 == 'CD4_Tm_CAPG+CREM-'] <- 'CD4_Tm_CAPG'
  seu$celltype_r2[seu$celltype_r2 == 'CD4_Tm_TNF'] <- 'CD4_Tm_CREM-'
  seu$celltype_r2[seu$celltype_r2 %in% c('CD4_Treg_TNFRSF9-', 'CD4_Treg_S1PR1')] <- 'CD4_Treg_Early'
  seu$celltype_r2[seu$celltype_r2 == 'CD4_Th_ISG'] <- 'CD4_Th_ISG15'
  seu$celltype_r2[seu$celltype_r2 == 'CD4_Treg_ISG'] <- 'CD4_Treg_ISG15'
  # CD8
  seu$celltype_r2[seu$celltype_r2 == 'CD8_Tm_IL7R'] <- 'CD8_Tcm_IL7R'
  seu$celltype_r2[seu$celltype_r2 == 'CD8_Tex_TCF7'] <- 'CD8_Tpex_TCF7'
  seu$celltype_r2[seu$celltype_r2 == 'CD8_Tm_NME1'] <- 'CD8_Prolif'
  seu$celltype_r2[seu$celltype_r2 == 'CD8_Tm_ZNF683'] <- 'CD8_Trm_ZNF683'
  seu$celltype_r2[seu$celltype_r2 %in% c('CD8_NK-like_EOMES', 'CD8_NK-like_TXK')] <- 'CD8_NK-like'
  seu$celltype_r2[seu$celltype_r2 == 'CD8_Tm_ZNF683'] <- 'CD8_Trm_ZNF683'
  seu$celltype_r2[seu$celltype_r2 == 'CD8_MAIT_SLC4A10'] <- 'MAIT'
  seu$celltype_r2[seu$celltype_r2 == 'CD8_Tex_OXPHOS-'] <- 'CD8_Tex_CXCL13'
  seu$celltype_r2[seu$celltype_r2 == 'CD8_ISG'] <- 'CD8_Tex_ISG15'
  # NK
  seu$celltype_r2[seu$celltype_r2 == 'CD56highCD16low'] <- 'NK_CD56hiCD16lo'
  seu$celltype_r2[seu$celltype_r2 == 'CD56lowCD16high'] <- 'NK_CD56loCD16hi'
  # DC
  seu$celltype_r2[seu$celltype_r2 == 'cDC2_CD1A'] <- 'DC_LC-like'
  seu$celltype_r2[seu$celltype_r2 == 'cDC2_FCN1'] <- 'MoDC'
  seu$celltype_r2[seu$celltype_r2 == 'cDC2_CXCR4hi'] <- 'cDC2_CD1C'
  seu$celltype_r2[seu$celltype_r2 == 'cDC3'] <- 'MigrDC'
  # Macro
  seu$celltype_r2[seu$celltype_r2 == 'Macro_GPNMB'] <- 'Macro_TREM2'
  seu$celltype_r2[seu$celltype_r2 == 'Macro_IL1B'] <- 'Macro_TNF'
  seu$celltype_r2[seu$celltype_r2 == 'Macro_NLRP3'] <- 'Macro_IL1B'
  # B/Plasma
  seu$celltype_r2[seu$celltype_r2 == 'B.01.TCL1A+naiveB'] <- 'B_Naive'
  seu$celltype_r2[seu$celltype_r2 == 'B.02.IFIT3+B'] <- 'B_ISG15'
  seu$celltype_r2[seu$celltype_r2 == 'B.03.HSP+B'] <- 'B_HSP'
  seu$celltype_r2[seu$celltype_r2 == 'B.04.MT1X+B'] <- 'B_MT2A'
  seu$celltype_r2[seu$celltype_r2 == 'B.05.EGR1+ACB'] <- 'ACB_EGR1'
  seu$celltype_r2[seu$celltype_r2 == 'B.06.NR4A2+ACB2'] <- 'ACB_NR4A2'
  seu$celltype_r2[seu$celltype_r2 == 'B.07.CCR7+ACB3'] <- 'ACB_CCR7'
  seu$celltype_r2[seu$celltype_r2 == 'B.08.ITGB1+SwBm'] <- 'B_Memory'
  seu$celltype_r2[seu$celltype_r2 == 'B.09.DUSP4+AtM'] <- 'B_AtM'
  seu$celltype_r2[seu$celltype_r2 == 'B.10.ENO1+Pre_GCB'] <- 'GCB_Pre'
  seu$celltype_r2[seu$celltype_r2 == 'B.11.SUGCT+DZ_GCB'] <- 'GCB_SUGCT'
  seu$celltype_r2[seu$celltype_r2 == 'B.12.LMO2+LZ_GCB'] <- 'GCB_LMO2'
  seu$celltype_r2[seu$celltype_r2 == 'B.13.Cycling_GCB'] <- 'GCB_Prolif'
  seu$celltype_r2[seu$celltype_r2 == 'B.14.Plasmablast'] <- 'Plasmablast'
  seu$celltype_r2[seu$celltype_r2 == 'B.15.Plasma cell'] <- 'Plasma_cell'
  # CAF
  seu$celltype_r2[str_detect(seu$celltype_r2, 'EndMT')] <- 'EndMT'
  # Endo
  seu$celltype_r2[seu$celltype_r2 == 'lymphatics'] <- 'Endo_lymphatic'
  seu$celltype_r2[seu$celltype_r2 == 'arteries'] <- 'Endo_artery'
  seu$celltype_r2[seu$celltype_r2 == 'capillaries'] <- 'Endo_capillary'
  seu$celltype_r2[seu$celltype_r2 == 'tip cell'] <- 'Endo_tip'
  seu$celltype_r2[seu$celltype_r2 == 'veins'] <- 'Endo_vein'
  
  seu$time_point <- ifelse(seu$time_point == 'Post', 'On', seu$time_point)
  
  qsave(seu, paste0('./PanCancer_ICI/data/', dataset, '/seu_r2.qs')) 
  if (!dir.exists(paste0('./PanCancer_ICI/data/', dataset, '/seu_r2/'))){
    dir.create(paste0('./PanCancer_ICI/data/', dataset, '/seu_r2/'))
  }
  Export10X(seu, dir =paste0('./PanCancer_ICI/data/', dataset, '/seu_r2/'), 
            append_reductions = NULL, gzip = F)
})





