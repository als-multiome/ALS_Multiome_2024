
  ### 0.0 Load libraries -------------------------------------------------------


source("~/ALS_Brain_Multiome.Rcfg") 

library(qs)
library(tidyverse)
library(lme4)


qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)




  ### 1.0 Load data ------------------------------------------------------------

M0_RNA <- qread(
  paste0(
    "../Data/SeuratObjects/", 
    "M0_RNA", 
    ".qrds"
  ), 
    nthr=nthr
)


  

  ### 2.0 MASC Case ALS vs. HC -------------------------------------------------

dat.ALS <- M0_RNA@meta.data[
  M0_RNA@meta.data$Case %in% c("ALS", "HC"),
]
dat.ALS$Case <- factor(
  dat.ALS$Case, 
  levels = c("HC", "ALS")
)

dat.ALS <- dat.ALS %>% 
  group_by(Well) %>% 
    mutate(Well_CellLibrary_Size = n()) %>% 
      ungroup() 
  
dat.ALS$Well_CellLibrary_Size <- scale(dat.ALS$Well_CellLibrary_Size)

MASC_ALS_WNN_L25 <- MASC(
  dataset = dat.ALS, 
  cluster = "WNN_L25",  
  contrast = "Case",  
  random_effects = c("ID", "Well"), 
  fixed_effects = c("Sex", "Well_CellLibrary_Size"), 
  verbose = TRUE, 
  return_models = TRUE
)
rm(dat.ALS)




  ### 3.0 MASC Case ALS_FTD vs. HC ---------------------------------------------

dat.ALSFTD <- M0_RNA@meta.data[
  M0_RNA@meta.data$Case %in% c("ALS_FTD", "HC"),
]
dat.ALSFTD$Case <- factor(
  dat.ALSFTD$Case, 
  levels = c("HC", "ALS_FTD")
)

dat.ALSFTD <- dat.ALSFTD %>% 
  group_by(Well) %>% 
  mutate(Well_CellLibrary_Size = n()) %>% 
  ungroup() 

dat.ALSFTD$Well_CellLibrary_Size <- scale(dat.ALSFTD$Well_CellLibrary_Size)



MASC_ALSFTD_WNN_L25 <- MASC(
  dataset = dat.ALSFTD, 
  cluster = "WNN_L25",  
  contrast = "Case",  
  random_effects = c("ID"), 
  fixed_effects = c("Sex", "Well_CellLibrary_Size"), 
  verbose = TRUE, 
  return_models = TRUE
)
rm(dat.ALSFTD)




  ### 4.0 MASC Case C9_ALS_FTD vs. ALS_FTD -------------------------------------

dat.C9_ALSFTD <- M0_RNA@meta.data[
  M0_RNA@meta.data$Case_Type %in% c("C9_ALS_FTD", "ALS_FTD"),
]
dat.C9_ALSFTD$Case_Type <- factor(
  dat.C9_ALSFTD$Case_Type, 
  levels = c("ALS_FTD", "C9_ALS_FTD")
)

dat.C9_ALSFTD <- dat.C9_ALSFTD %>% 
  group_by(Well) %>% 
  mutate(Well_CellLibrary_Size = n()) %>% 
  ungroup() 

dat.C9_ALSFTD$Well_CellLibrary_Size <- scale(dat.C9_ALSFTD$Well_CellLibrary_Size)


MASC_C9_ALSFTD_WNN_L25 <- MASC(
  dataset = dat.C9_ALSFTD, 
  cluster = "WNN_L25",  
  contrast = "Case_Type",  
  random_effects = c("ID", "Well"), 
  fixed_effects = c("Sex", "Well_CellLibrary_Size"), 
  verbose = TRUE, 
  return_models = TRUE
)
rm(dat.C9_ALSFTD)




  ### 5.0 Save data ------------------------------------------------------------

qsave(
  MASC_ALS_WNN_L25, 
  paste0(
    "../Data/Coda/MASC/", 
    "MASC_ALS_WNN_L25", 
    ".qrds"
  ),  
  nthr = nthr
)

qsave(
  MASC_ALSFTD_WNN_L25, 
  paste0(
    "../Data/Coda/MASC/", 
    "MASC_ALSFTD_WNN_L25", 
    ".qrds"
  ),  
  nthr = nthr
)

qsave(
  MASC_C9_ALSFTD_WNN_L25, 
  paste0(
    "../Data/Coda/MASC/", 
    "MASC_C9_ALSFTD_WNN_L25", 
    ".qrds"
  ),  
  nthr = nthr
)
