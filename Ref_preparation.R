rm(list=ls())
pkgs <- c('Seurat','tidyr','plyr','dplyr','stringr','SingleR','ggsci','dior','tibble','qs','BiocParallel','Matrix','SingleCellExperiment','scran','parallel','scGate')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)

# Reference Preparation
# Myeloids (normalized)
directory <- './PanCancer_ICI/data/GSE154763/'
file_list <- list.files(directory, pattern="expression.csv$")[-1]
expression_list <- lapply(file_list, function(x){t(read.csv(paste0(directory, x), row.names = 1))})
meta_list <- list.files(directory, pattern="metadata.csv$")[-1]
metadata_list <- lapply(meta_list, function(x){read.csv(paste0(directory, x),row.names = 1)})

ref_list <- c()
for (i in 1:length(file_list)){
  print(i)
  sce <- SingleCellExperiment(assays = list(logcounts = expression_list[[i]]), colData = metadata_list[[i]])
  sce$subtype <- sapply(strsplit(sce$MajorCluster,'_'), function(x){paste0(x[2],'_',x[3])})
  sce$majortype <- sapply(strsplit(sce$MajorCluster,'_'), function(x){x[2]})
  sce <- sce[, !sce$subtype %in% c('Mast_KIT','pDC_LILRA4')]
  ref_list <- c(ref_list, sce)
}
ref_major <- lapply(ref_list, function(x){aggregateReference(x, x$majortype)})
qsave(ref_major, "./PanCancer_ICI/data/Ref_SingleR/Myeloids_major_ref.qs")
ref_mono <- list()
for (i in 1:length(ref_list)){
  sce <- ref_list[[i]]
  sce <- sce[, str_split(sce$subtype, '_', simplify = T)[, 1] == 'Mono']
  sce <- aggregateReference(sce, sce$subtype)
  ref_mono <- c(ref_mono, sce)
}
qsave(ref_mono, "./PanCancer_ICI/data/Ref_SingleR/Myeloids_mono_ref.qs")
ref_macro <- list()
for (i in 1:length(ref_list)){
  sce <- ref_list[[i]]
  sce <- sce[, str_split(sce$subtype, '_', simplify = T)[, 1] == 'Macro']
  sce <- aggregateReference(sce, sce$subtype)
  ref_macro <- c(ref_macro, sce)
}
qsave(ref_macro, "./PanCancer_ICI/data/Ref_SingleR/Myeloids_macro_ref.qs")
# cDC2
expr <- list.files(directory, pattern="expression.csv$")[1]
expr <- t(read.csv(paste0(directory, expr), row.names = 1))
meta <- list.files(directory, pattern="metadata.csv$")[1]
meta <- read.csv(paste0(directory, meta),row.names = 1)
sce <- SingleCellExperiment(assays = list(logcounts = expr), colData = meta) 
sce$subtype <- sapply(strsplit(sce$MajorCluster,'_'), function(x){paste0(x[2],'_',x[3])})
ref_cdc2 <- list()
for (i in 1:length(unique(sce$cancer))){
  sce_sub <- sce[, sce$cancer == unique(sce$cancer)[i]]
  sce_sub <- aggregateReference(sce_sub, sce_sub$subtype)
  ref_cdc2 <- c(ref_cdc2, sce_sub)
}
qsave(ref_cdc2, "./PanCancer_ICI/data/Ref_SingleR/Myeloids_cdc2_ref.qs")

#T/NK (counts)
file_list <- list.files("./PanCancer_ICI/data/GSE156728(ref_Tcell)/", pattern = 'counts', full.names = T)
cancertype <- unique(str_split(file_list,'_',simplify = T)[,4])
df_metadata <- data.table::fread('./PanCancer_ICI/data/GSE156728(ref_Tcell)/GSE156728_metadata.txt.gz') %>%
  filter(cancerType %in% cancertype, loc=='T') %>% 
  column_to_rownames(var = 'cellID')
df_metadata$meta.cluster <- str_replace_all(df_metadata$meta.cluster, '\\.', '_')
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD4') & df_metadata$meta.cluster == 'CD4_c01_Tn_TCF7'] <- 'CD4_c01_Naive'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD4') & df_metadata$meta.cluster == 'CD4_c02_Tn_PASK'] <- 'CD4_c02_pre-Tfh_CXCR5+'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD4') & df_metadata$meta.cluster == 'CD4_c03_Tn_ADSL'] <- 'CD4_c03_Tn_ADSL'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD4') & df_metadata$meta.cluster == 'CD4_c04_Tn_il7r'] <- 'CD4_c04_Tn_IL7R-'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD4') & df_metadata$meta.cluster == 'CD4_c06_Tm_ANXA1'] <- 'CD4_c06_Tm_AREG'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD4') & df_metadata$meta.cluster == 'CD4_c07_Tm_ANXA2'] <- 'CD4_c07_Tm_TIMP1'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD4') & df_metadata$meta.cluster == 'CD4_c11_Tm_GZMA'] <- 'CD4_c11_Tm_CAPG+CREM-'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD4') & df_metadata$meta.cluster == 'CD4_c14_Th17_SLC4A10'] <- 'CD4_c14_Th17_CCR6'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD4') & df_metadata$meta.cluster == 'CD4_c15_Th17_IL23R'] <- 'CD4_c15_Th17_IL26'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD4') & df_metadata$meta.cluster == 'CD4_c17_TfhTh1_CXCL13'] <- 'CD4_c17_TfhTh1_IFNG'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD4') & df_metadata$meta.cluster == 'CD4_c18_Treg_RTKN2'] <- 'CD4_c18_Treg_TNFRSF9-'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD4') & df_metadata$meta.cluster == 'CD4_c21_Treg_OAS1'] <- 'CD4_c21_Treg_ISG'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD4') & df_metadata$meta.cluster == 'CD4_c22_ISG_IFIT1'] <- 'CD4_c22_Th_ISG'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD8') & df_metadata$meta.cluster == 'CD8_c01_Tn_MAL'] <- 'CD8_c01_Naive'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD8') & df_metadata$meta.cluster == 'CD8_c03_Tm_RPS12'] <- 'uncharacterized'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD8') & df_metadata$meta.cluster == 'CD8_c04_Tm_CD52'] <- 'CD8_c04_Tm_ZNF683'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD8') & df_metadata$meta.cluster == 'CD8_c05_Tem_CXCR5'] <- 'CD8_c05_Tem_Early'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD8') & df_metadata$meta.cluster == 'CD8_c08_Tk_TYROBP'] <- 'CD8_c08_NK-like_EOMES'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD8') & df_metadata$meta.cluster == 'CD8_c09_Tk_KIR2DL4'] <- 'CD8_c09_NK-like_TXK'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD8') & df_metadata$meta.cluster == 'CD8_c11_Tex_PDCD1'] <- 'CD8_c11_Tex_GZMK'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD8') & df_metadata$meta.cluster == 'CD8_c13_Tex_myl12a'] <- 'CD8_c13_Tex_OXPHOS-'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD8') & df_metadata$meta.cluster == 'CD8_c15_ISG_IFIT1'] <- 'CD8_c15_ISG'
df_metadata$meta.cluster[df_metadata$meta.cluster != 'uncharacterized'] <- str_replace(df_metadata$meta.cluster[df_metadata$meta.cluster != 'uncharacterized'], 'c\\d{2}\\_', '')
df_metadata <- filter(df_metadata, !meta.cluster %in% c('uncharacterized', 'CD4_Mix_NME1', 'CD4_Mix_NME2'))

ref_CD4 <- lapply(file_list[str_detect(file_list, 'CD4')], function(x){
  cancertype <- unique(str_split(x,'_',simplify = T)[,4])
  matrix_meta <- filter(df_metadata, cancerType == cancertype, str_detect(meta.cluster, 'CD4')) 
  matrix_count <- data.table::fread(x) %>% tibble::column_to_rownames(var = 'V1') %>% as.sparse()
  matrix_count <- matrix_count[, which(colnames(matrix_count) %in% rownames(matrix_meta))]
  sce <- SingleCellExperiment(assay = list(counts = matrix_count)) 
  sce$major <- 'CD4'
  return(sce)
})

ref_CD8 <- lapply(file_list[str_detect(file_list, 'CD8')], function(x){
  cancertype <- unique(str_split(x,'_',simplify = T)[,4])
  print(cancertype)
  matrix_meta <- filter(df_metadata, cancerType == cancertype, str_detect(meta.cluster, 'CD8')) 
  matrix_count <- data.table::fread(x) %>% tibble::column_to_rownames(var = 'V1') %>% as.sparse()
  matrix_count <- matrix_count[, which(colnames(matrix_count) %in% rownames(matrix_meta))]
  sce <- SingleCellExperiment(assay = list(counts = matrix_count)) 
  sce$major <- 'CD8'
  return(sce)
})
ref_list_T <- list()
for (i in 1:5){
  print(i)
  common_genes <- intersect(rownames(ref_CD4[[i]]), rownames(ref_CD8[[i]]))
  ref_T <- cbind(ref_CD4[[i]][common_genes,], ref_CD8[[i]][common_genes,])
  ref_list_T[[i]] <- ref_T
}

ref_T <- lapply(ref_list_T, function(sce){
  sce <- logNormCounts(sce) %>% aggregateReference(.$major)
  return(sce)
})
qsave(ref_T, "./PanCancer_ICI/data/Ref_SingleR/T_ref.qs")

# Fine level
ref_CD4 <- lapply(file_list[str_detect(file_list, 'CD4')], function(x){
  cancertype <- unique(str_split(x,'_',simplify = T)[,4])
  matrix_meta <- filter(df_metadata, cancerType == cancertype, str_detect(meta.cluster, 'CD4')) 
  matrix_count <- data.table::fread(x) %>% tibble::column_to_rownames(var = 'V1') %>% as.sparse()
  matrix_count <- matrix_count[, which(colnames(matrix_count) %in% rownames(matrix_meta))]
  sce <- SingleCellExperiment(assay = list(counts = matrix_count), colData = matrix_meta[colnames(matrix_count),]) %>% logNormCounts() %>% aggregateReference(.$meta.cluster)
  return(sce)
})
qsave(ref_CD4, "./PanCancer_ICI/data/Ref_SingleR/T_CD4_ref.qs")

ref_CD8 <- lapply(file_list[str_detect(file_list, 'CD8')], function(x){
  cancertype <- unique(str_split(x,'_',simplify = T)[,4])
  matrix_meta <- filter(df_metadata, cancerType == cancertype, str_detect(meta.cluster, 'CD8')) 
  matrix_count <- data.table::fread(x) %>% tibble::column_to_rownames(var = 'V1') %>% as.sparse()
  matrix_count <- matrix_count[, which(colnames(matrix_count) %in% rownames(matrix_meta))]
  sce <- SingleCellExperiment(assay = list(counts = matrix_count), colData = matrix_meta[colnames(matrix_count),]) %>% logNormCounts() %>% aggregateReference(.$meta.cluster)
  return(sce)
})
qsave(ref_CD8, "./PanCancer_ICI/data/Ref_SingleR/T_CD8_ref.qs")

# NK
data_dir <- './PanCancer_ICI/data/GSE212890(ref_NK)/'
list.files(data_dir) 
matrix_count <- readMM(paste0(data_dir, list.files(data_dir)[str_detect(list.files(data_dir) , 'counts')])) 
features <- read.csv(paste0(data_dir, list.files(data_dir)[str_detect(list.files(data_dir) , 'genes')]))
barcodes <- read.csv(paste0(data_dir, list.files(data_dir)[str_detect(list.files(data_dir) , 'barcodes')]))
metadata <- read.csv(paste0(data_dir, list.files(data_dir)[str_detect(list.files(data_dir) , 'metadata')]))
colnames(matrix_count) <- features$X0
rownames(matrix_count) <- barcodes$X0
matrix_count <- as(matrix_count, "CsparseMatrix")
# ref_NK <- SingleCellExperiment(assay = list(counts = t(matrix_count)))
# ref_NK <- ref_NK[, metadata$meta_tissue == 'Tumor']
# ref_NK$major <- 'NK'

ref_NK <- SingleCellExperiment(assay = list(counts = t(matrix_count)), colData = metadata) 
ref_NK <- ref_NK[, ref_NK$meta_tissue == 'Tumor']
cancertype <- unique(ref_NK$meta_histology)
ref_NK <- lapply(cancertype, function(x) {
  return(ref_NK[, ref_NK$meta_histology == x])
})
ref_main_NK <- lapply(ref_NK, function(x) {
  sce <- logNormCounts(x) %>% aggregateReference(.$Majortype)
  return(sce)
})
qsave(ref_main_NK, "./PanCancer_ICI/data/Ref_SingleR/NK_main_ref.qs")
ref_fine_NK <- lapply(ref_NK, function(x) {
  sce <- logNormCounts(x) %>% aggregateReference(.$celltype)
  return(sce)
})
qsave(ref_fine_NK, "./PanCancer_ICI/data/Ref_SingleR/NK_fine_ref.qs")

# CAFs
metadata <- read.table(gzfile('./PanCancer_ICI/data/GSE210347(ref_CAF)/GSE210347_meta.txt.gz'), row.names = 1, header = T) %>%
  filter(cluster %in% c('c1','c2','c4','c6','c7','c8'), group == 'Tumor')
metadata$cluster <- mapvalues(metadata$cluster,from = c('c1','c2','c4','c6','c7','c8'),
                              to = c('Myofibroblast','CAF_inflammatory','CAF_adipogenic','CAF_EndMT','CAF_PN','CAF_AP'))
matrix_count <- readRDS('./PanCancer_ICI/data/GSE210347(ref_CAF)/GSE210347_counts.Rds')[, rownames(metadata)]
caf_ref <- CreateSeuratObject(counts = matrix_count, meta.data = metadata)
sce <- SingleCellExperiment(assay = list(counts = matrix_count), colData = metadata)
sce$dataset <- str_split(sce$SampleID, "_", simplify = T)[,1]
dataset <- unique(sce$dataset)
ref_CAF <- lapply(dataset, function(x) {
  sce <- sce[, sce$dataset == x]
  sce <- logNormCounts(sce) %>% aggregateReference(.$cluster)
  return(sce)
})
qsave(ref_CAF, "./PanCancer_ICI/data/Ref_SingleR/CAF_ref.qs")

# CAF
metadata <- read.table(gzfile('./PanCancer_ICI/data/GSE210347(ref_CAF)/GSE210347_meta.txt.gz'), row.names = 1, header = T) %>%
  filter(cluster %in% c('c1','c2','c4','c6','c7','c8'), group == 'Tumor')
metadata$cluster <- mapvalues(metadata$cluster,from = c('c1','c2','c4','c6','c7','c8'),
                              to = c('Myofibroblast','CAF_inflammatory','CAF_adipogenic','CAF_EndMT','CAF_PN','CAF_AP'))
matrix_count <- readRDS('./PanCancer_ICI/data/GSE210347(ref_CAF)/GSE210347_counts.Rds')[, rownames(metadata)]
seu <- CreateSeuratObject(counts = matrix_count, meta.data = metadata)
seu$dataset <- str_split(seu$SampleID, "_", simplify = T)[,1]
sce <- seu |> as.SingleCellExperiment() |> logNormCounts() 
sce <- sce |> aggregateReference(sce$cluster)
qsave(sce, file = './PanCancer_ICI/data/Ref_SingleR/Fibro.qs')

# Pan-cancer Endothelial cells
seu <- schard::h5ad2seurat('.//PanCancer_ICI/data/pan-cancer_Endo/NSR.panE.exprs.all.h5ad.gz')
metadata <- read.csv('./PanCancer_ICI/data/pan-cancer_Endo/NSR.panE.obs.meta.csv', row.names = 1)
seu <- CreateSeuratObject(counts = seu@assays$RNA@counts, meta.data = metadata)
seu <- subset(seu, subset = tissue == 'T')
seu$major <- ifelse(seu$cellType.major %in% c('veins', 'hypoxia', 'tip cell', 'arteries', 'capillaries'), 'vascular', seu$cellType.major)
ref_list <- lapply(unique(seu$Dataset), function(dataset){
  print(dataset)
  sce <- subset(seu, subset = Dataset == dataset) |> 
    as.SingleCellExperiment() |> 
    logNormCounts() 
  sce <- sce |> aggregateReference(labels = sce$major)
  return(sce)
})
qsave(ref_list, file = './PanCancer_ICI/data/Ref_SingleR/Endo_major.qs')
# Endo_vascular
ref_list <- lapply(unique(seu$Dataset), function(dataset){
  print(dataset)
  sce <- subset(seu, subset = seu$cellType.major %in% c('veins', 'tip cell', 'arteries', 'capillaries') & Dataset == dataset) |> 
    as.SingleCellExperiment() |> 
    logNormCounts() 
  sce <- sce |> aggregateReference(labels = sce$cellType.major)
  return(sce)
})
qsave(ref_list, file = './PanCancer_ICI/data/Ref_SingleR/Endo_vascular.qs')
# Mural cells
endo_mc <- readRDS('./PanCancer_ICI/data/ref_Endo/panMC.rds')
endo_mc$dataset <- str_split(endo_mc$DonorID, '_', simplify = T)[,1]
sce <- endo_mc |> 
  as.SingleCellExperiment() |> 
  logNormCounts() 
sce <- sce |> aggregateReference(labels = sce$Anno_Tier1)
qsave(sce, file = './PanCancer_ICI/data/Ref_SingleR/Mural_ref.qs')

# B cells
# seu <- readRDS("./PanCancer_ICI/data/blueprint_B/scRNA_data/panB_Plasma_cell_selected_scRNA_processed_data.rds")
seu <- readRDS("./PanCancer_ICI/data/blueprint_B/scRNA_data/panB_scRNA_processed_data.rds")
genes_to_check = c('TCL1A','YBX3','FCER2',
                   'CD27','CD19','IFIT3', 'IFI44L', 'STAT1', 'ISG15',
                   'HSPA1A','DNAJB1','HSPA1B',
                   'MT1X','MT2A','SLC30A1',
                   'TNF','EGR1',
                   'ZNF331','VPS37B','NR4A2','CD69','CD86','IGHM','IGHD',
                   'CD83','CCR7','NFKBID',
                   'CRIP1','S100A10','S100A4', 'ITGB1',
                   'DUSP4','FCRL5','ZEB2','ITGAX',
                   'ENO1','PSME2','NME1',
                   'SUGCT','NEIL1','MEF2B',
                   'LMO2','GMDS','MARCKSL1',
                   'STMN1','HMGB2','MKI67','TOP2A',
                   'TXNDC5','MYDGF',
                   'XBP1','MZB1','JCHAIN','MS4A1')
DotPlot(seu, group.by = 'celltype', features = genes_to_check) + RotatedAxis() + scale_colour_gradient2(low = "#FFFFD9", mid = "#41B6C4", high = "#081D58")
seu <- CreateSeuratObject(counts = seu@assays$RNA@counts, meta.data = seu@meta.data)
seu$majortype <- ifelse(seu$celltype %in% c("B.14.Plasmablast", "B.15.Plasma cell"), 'Plasma',
                        ifelse(seu$celltype %in% c("B.12.LMO2+LZ_GCB", "B.11.SUGCT+DZ_GCB", "B.13.Cycling_GCB"), 'GCB', 'B'))
ref_B_major <- lapply(unique(seu$dataid), function(x) {
  print(x)
  seu_sub <- subset(seu, subset = dataid == x)
  sce <- SingleCellExperiment(assay = list(counts = seu_sub@assays$RNA$counts), colData = seu_sub@meta.data) |> 
    logNormCounts() 
  sce <- aggregateReference(sce,sce$majortype)
  return(sce)
})
qsave(ref_B_major, "./PanCancer_ICI/data/Ref_SingleR/Pan_B_major.qs")

ref_plasma <- lapply(unique(seu$dataid), function(x) {
  print(x)
  seu_sub <- subset(seu, subset = dataid == x & majortype == 'Plasma')
  sce <- SingleCellExperiment(assay = list(counts = seu_sub@assays$RNA$counts), colData = seu_sub@meta.data) |> 
    logNormCounts() 
  sce <- aggregateReference(sce,sce$celltype)
  return(sce)
})
qsave(ref_plasma, "./PanCancer_ICI/data/Ref_SingleR/Plasma.qs")

seu_gcb <- subset(seu, subset = majortype == 'GCB')
ref_GCB <- lapply(unique(seu_gcb$dataid)[table(seu_gcb$dataid)>50], function(x) {
  print(x)
  seu_sub <- subset(seu, subset = dataid == x & majortype == 'GCB')
  sce <- SingleCellExperiment(assay = list(counts = seu_sub@assays$RNA$counts), colData = seu_sub@meta.data) |> 
    logNormCounts() 
  sce <- aggregateReference(sce,sce$celltype)
  return(sce)
})
qsave(ref_GCB, "./PanCancer_ICI/data/Ref_SingleR/GCB.qs")

seu_b <- subset(seu, subset = majortype == 'B')
ref_B <- lapply(names(table(seu_b$dataid))[table(seu_b$dataid)>50], function(x) {
  print(x)
  seu_sub <- subset(seu, subset = dataid == x & majortype == 'B')
  sce <- SingleCellExperiment(assay = list(counts = seu_sub@assays$RNA$counts), colData = seu_sub@meta.data) |> 
    logNormCounts() 
  sce <- aggregateReference(sce,sce$celltype)
  return(sce)
})
qsave(ref_B, "./PanCancer_ICI/data/Ref_SingleR/B.qs")


