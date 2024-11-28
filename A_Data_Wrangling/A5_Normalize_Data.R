

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




  ### 4.0 Normalize RNA Data with LogNormalize ---------------------------------

M0@assays$RNA <- NormalizeData(M0@assays$RNA, normalization.method = "LogNormalize")




  ### 5.0 Normalize ATAC Data with TFIDF ---------------------------------------

M0@assays$ATAC <- RunTFIDF(M0@assays$ATAC, method=1)

message("Done!")




  ### 6.0 Export Data ----------------------------------------------------------

qsave(
  M0, 
  paste0(
    "../Data/SeuratObjects/", 
    "M0", 
    ".qrds"
  ), 
  nthr = nthr
)



