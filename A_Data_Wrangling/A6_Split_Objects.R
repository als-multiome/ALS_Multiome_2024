

  ### 0.0 Load libraries -------------------------------------------------------

source("~/ALS_Brain_Multiome.Rcfg") 

library(qs) 
library(Seurat)
library(Signac)




  ### 1.0 Load data ------------------------------------------------------------


      # M0 Seurat single-cell data object --------------------------------------

M0 <- qread(
  paste0(
    "../Data/SeuratObjects/", 
    "M0", 
    ".qrds"
  ), 
  nthr = nthr
)




  ### 2.0 Save full object -----------------------------------------------------

qsave(
  M0, 
  paste0(
    "../Data/SeuratObjects/", 
    "M0_Full", 
    ".qrds"
  ), 
  nthr = nthr
) 




  ### 3.0 Subset Other Cells ---------------------------------------------------

M0 <- subset(M0, subset=WNN_L2=="Other_Cells", invert=TRUE)




  ### 4.0 Subset by Case -------------------------------------------------------

table(M0$Case)
table(M0$Case_Type)

M0_HC <- subset(M0, subset = Case == "HC") 
M0_HC 
qsave(
  M0_HC, 
  paste0(
    "../Data/SeuratObjects/", 
    "M0_HC", 
    ".qrds"
  ), 
  nthr = nthr
) 
rm(M0_HC)


M0_ALS <- subset(M0, subset = Case == "ALS") 
M0_ALS 
qsave(
  M0_ALS, 
  paste0(
    "../Data/SeuratObjects/", 
    "M0_ALS", 
    ".qrds"
  ), 
  nthr = nthr
) 
rm(M0_ALS)


M0_ALS_FTD_All <- subset(M0, subset = Case == "ALS_FTD") 
M0_ALS_FTD_All 
qsave(
  M0_ALS_FTD_All, 
  paste0(
    "../Data/SeuratObjects/", 
    "M0_ALS_FTD_All", 
    ".qrds"
  ), 
  nthr = nthr
) 
rm(M0_ALS_FTD_All)


M0_ALS_FTD <- subset(M0, subset = Case_Type == "ALS_FTD") 
M0_ALS_FTD 
qsave(
  M0_ALS_FTD, 
  paste0(
    "../Data/SeuratObjects/", 
    "M0_ALSFTD_nC9", 
    ".qrds"
  ), 
  nthr = nthr
) 
rm(M0_ALS_FTD) 


M0_C9_ALS_FTD <- subset(M0, subset = Case_Type == "C9_ALS_FTD") 
M0_C9_ALS_FTD 
qsave(
  M0_C9_ALS_FTD, 
  paste0(
    "../Data/SeuratObjects/", 
    "M0_ALSFTD_C9", 
    ".qrds"
  ), 
  nthr = nthr
) 
rm(M0_C9_ALS_FTD)



  ### 5.0 Subset by Assay ------------------------------------------------------

M0_RNA <- M0 
M0_RNA@assays$ATAC <- NULL 
M0_RNA@assays$RNATE <- NULL 
M0_RNA 
qsave(
  M0_RNA, 
  paste0(
    "../Data/SeuratObjects/", 
    "M0_RNA", 
    ".qrds"
  ), 
  nthr = nthr
) 
rm(M0_RNA) 

M0_ATAC <- M0 
M0_ATAC@assays$RNA <- NULL   
M0_ATAC@assays$RNATE <- NULL  
DefaultAssay(M0_ATAC) <- "ATAC"
M0_ATAC 
qsave(
  M0_ATAC, 
  paste0(
    "../Data/SeuratObjects/", 
    "M0_ATAC", 
    ".qrds"
  ), 
  nthr = nthr
) 
rm(M0_ATAC)


M0_RNATE <- M0 
M0_RNATE@assays$ATAC <- NULL  
M0_RNATE@assays$RNA <- NULL  
DefaultAssay(M0_RNATE) <- "RNATE"
M0_RNATE 
qsave(
  M0_RNATE, 
  paste0(
    "../Data/SeuratObjects/", 
    "M0_RNATE", 
    ".qrds"
  ), 
  nthr = nthr
) 
rm(M0_RNATE)




  ### 6.0 Save full object -----------------------------------------------------

M0@assays$RNATE <- NULL 
qsave(
  M0, 
  paste0(
    "../Data/SeuratObjects/", 
    "M0", 
    ".qrds"
  ), 
  nthr = nthr
)  


