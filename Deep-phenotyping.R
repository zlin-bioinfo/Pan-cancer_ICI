rm(list=ls())
pkgs <- c('Seurat','tidyr','plyr','dplyr','stringr','SingleR','ggsci','dior','tibble','qs','BiocParallel','scGate','Matrix','SingleCellExperiment','scran','parallel','scGate','burgertools')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)

# Loading dataset
SKCM_Becker <- qread('./SKCM_Becker/seu_r1.qs')
BRCA_Bassez1 <- qread('./BRCA_Bassez1/seu_r1.qs')
BRCA_Bassez2 <- qread('./BRCA_Bassez2/seu_r1.qs')
TNBC_Zhang <- qread('./TNBC_Zhang/seu_r1.qs')
BCC_Yost <- qread('./BCC_Yost/seu_r1.qs')
BCC_Yost1 <- qread('./BCC_Yost/seu_r1.qs') %>% 
  subset(subset = patient %in% c('BCC/SCC_Yost_su009', 'BCC/SCC_Yost_su012'), invert = T)
BCC_Yost1$dataset <- 'BCC_Yost'
BCC_Yost2 <- qread('./BCC_Yost/seu_r1.qs') %>% 
  subset(subset = patient %in% c('BCC/SCC_Yost_su009', 'BCC/SCC_Yost_su012')) 
BCC_Yost2$dataset <- 'BCC_Yost'
SCC_Yost <- qread('./SCC_Yost/seu_r1.qs')
SCC_Yost$dataset <- 'SCC_Yost'
HNSC_IMCISION <- qread('./HNSC_IMCISION/seu_r1.qs')
HNSC_Luoma <- qread('./HNSC_Luoma/seu_r1.qs')
NSCLC_Liu <- qread('./NSCLC_Liu/seu_r1.qs')
CRC_Li <- qread('./CRC_Li/seu_r1.qs')
PCa_Hawley <- qread('./PCa_Hawley/seu_r1.qs')
TNBC_Shiao <- qread('./TNBC_Shiao/seu_r1.qs')
HNSC_Franken <- qread('./HNSC_Franken/seu_r1.qs')

# datasets <- list(SKCM_Becker, BRCA_Bassez1, BRCA_Bassez2, TNBC_Zhang, BCC_Yost1, HNSC_IMCISION, HNSC_Luoma, NSCLC_Liu, CRC_Li, PCa_Hawley)

# T_NK 
ref_T <- qread("./Ref_SingleR/T_ref.qs")
ref_CD4 <- qread("./Ref_SingleR/T_CD4_ref.qs")
ref_CD8 <- qread("./Ref_SingleR/T_CD8_ref.qs")
ref_main_NK <- qread("./Ref_SingleR/NK_main_ref.qs")
ref_fine_NK <- qread("./Ref_SingleR/NK_fine_ref.qs")

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
  qsave(pred, paste0("./", unique(seu$dataset), "/", type_name, ".qs"))
  return(pred)
}

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
  write.csv(barcode_NK, paste0('./', unique(seu$dataset), '/barcode_NK.csv'))
  print('Gating gdT cells')
  gdT <- gating_model(name = "gdT", signature = c("CD3D","TRGC2","TRDV2","TRGV9","TRGV10","TRDC","CD8A-","CD4-"))
  seu_gdT<- scGate(data = seu, model = gdT)
  barcode_gdT <- colnames(seu_gdT)[seu_gdT$is.pure == 'Pure']
  # barcode_gdT <- setdiff(barcode_gdT, barcode_NK)
  print(paste0(length(barcode_gdT), ' gdT cells were detected'))
  write.csv(barcode_gdT, paste0('./', unique(seu$dataset), '/barcode_gdT.csv'))
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
hierSingleRMultiRef_T_NK(HNSC_Franken, main_ref_list = ref_T, 
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
  write.csv(barcode_gdT, paste0('./', unique(seu$dataset), '/barcode_gdT_sorted.csv'))
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
# Myeloids
ref_Myeloids <- qread("./Ref_SingleR/Myeloids_major_ref.qs")
ref_cdc2 <- qread("./Ref_SingleR/Myeloids_cdc2_ref.qs")
ref_mono <- qread("./Ref_SingleR/Myeloids_mono_ref.qs")
ref_macro <- qread("./Ref_SingleR/Myeloids_macro_ref.qs")
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
hierSingleRMultiRef_Myeloids(HNSC_Franken, main_ref_list = ref_Myeloids, 
                             ref_list_mono = ref_mono, ref_list_macro = ref_macro, ref_list_cdc2 = ref_cdc2)

# run SingleR with one reference dataset
SingleROneRef <- function(seu, ref, label, major, type_name, n_cores = 50){
  pred <- seu %>% 
    subset(subset = celltype_major %in% major) %>% 
    as.SingleCellExperiment() %>% 
    logNormCounts() %>% 
    SingleR(ref = ref, labels = label, assay.type.test=1, BPPARAM=MulticoreParam(50))
  qsave(pred, paste0("./", unique(seu$dataset), "/", type_name, ".qs"))
}

# pan-cancer B
ref_Bplasma <- qread("./Ref_SingleR/Pan_B_major.qs")
ref_b <- qread("./Ref_SingleR/B.qs")
ref_gcb <- qread("./Ref_SingleR/GCB.qs")
ref_plasma <- qread("./Ref_SingleR/Plasma.qs")
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
hierSingleRMultiRef_Bplasma(PCa_Hawley, main_ref_list = ref_Bplasma, ref_list_b = ref_b, 
                             ref_list_gcb = ref_gcb, ref_list_plasma = ref_plasma)
# Endothelial cells
ref_Endo <- qread("./Ref_SingleR/Endo_ref.qs")
datasets <- list(SKCM_Becker, BRCA_Bassez1, BRCA_Bassez2, BCC_Yost1, CRC_Li, PCa_Hawley)
lapply(datasets, function(seu){SingleROneRef(seu, ref = ref_Endo, label = ref_Endo$label, major = 'Endothelial cells', type_name = 'Endo')})
SingleROneRef(HNSC_Franken, ref = ref_Endo, label = ref_Endo$label, major = 'Endothelial cells', type_name = 'Endo')

# CAF
ref_CAF <- qread("./Ref_SingleR/CAF_ref.qs")
datasets <- list(SKCM_Becker, BRCA_Bassez1, BRCA_Bassez2, BCC_Yost1, CRC_Li, PCa_Hawley)
lapply(datasets, function(seu){
  print(unique(seu$dataset))
  SingleRMultiRef(seu, ref_list = ref_CAF, major = c('Fibroblasts','Myocytes'), if.subset = T, type_name = 'CAF')})
SingleRMultiRef(HNSC_Franken, ref_list = ref_CAF, major = c('Fibroblasts','Myocytes'), if.subset = T, type_name = 'CAF')



