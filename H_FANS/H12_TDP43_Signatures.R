# H12_TDP43_Signatures 



  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 

library(qs)  


library(Seurat)
library(SeuratObject) 
library(MAST)
library(zinbwave)
library(patchwork)
library(tidyverse)




qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)




  ### 1.0 Load data ------------------------------------------------------------

FANS_DE_TDP43_MAST <- qread(
  paste0(
    "../Data/DE/FANS/", 
    "FANS_DE_TDP43_MAST", 
    ".qrds"
  ), 
  nthr=nthr
)

F2 <- qread(
  paste0(
    "../Data/SeuratObjects/", 
    "F2", 
    ".qrds"
  ), 
  nthr=nthr
)

VlnPlot(F2, features="ATP1A3") 

TDP43_Signature <- FANS_DE_TDP43_MAST[
  FANS_DE_TDP43_MAST$fdr<0.1, 
]

