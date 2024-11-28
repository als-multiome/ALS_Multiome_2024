
  ### 0.0 Load libraries -------------------------------------------------------

source("~/ALS_Brain_Multiome.Rcfg") 

library(qs) 
library(Seurat)
library(Signac)
library(zellkonverter)




  ### 1.0 Load data ------------------------------------------------------------


      # M0_RNA Seurat single-cell data object ----------------------------------

M0_RNA <- qread(
  paste0(
    "../Data/SeuratObjects/", 
    "M0_RNA", 
    ".qrds"
  ), 
  nthr = nthr
)




  ### 2.0 Export M0_RNA as AnnData and metadata as CSV -------------------------

M0_RNA$Rand1_AllCases <- M0_RNA@misc$Sample_data$Rand1_AllCases
M0_RNA$Rand1_ALS <- M0_RNA@misc$Sample_data$Rand1_ALS 
M0_RNA$Rand1_ALSFTD <- M0_RNA@misc$Sample_data$Rand1_ALSFTD
M0_RNA$Rand1_nC9_ALSFTD <- M0_RNA@misc$Sample_data$Rand1_ALSFTD_NonC9
M0_RNA$Rand1_C9_ALSFTD <- M0_RNA@misc$Sample_data$Rand1_ALSFTD_C9
M0_RNA$Rand1_C9_ALSFTD <- M0_RNA@misc$Sample_data$Rand1_ALSFTD_C9
M0_RNA$Rand1_ALSFTD_C9_NonC9 <- M0_RNA@misc$Sample_data$Rand1_ALSFTD_C9vsALSFTD_NonC9

M0_RNA_sce <- as.SingleCellExperiment(
  M0_RNA, 
  assay="RNA"
)

writeH5AD(
  M0_RNA_sce, 
  paste0(
    "../Data/AnnDataObjects/", 
    "M0_RNA_sce", 
    ".h5ad"
  ), 
  X_name="counts"
)

write.csv(
  M0_RNA@meta.data, 
  paste0(
    "../Data/Input/", 
    "M0_RNA_metadata", 
    ".csv"
  ), 
  quote=FALSE, 
  row.names = TRUE
)






