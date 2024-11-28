

  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 
library(qs)

qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)




  ### 1.0 Load data ------------------------------------------------------------

M0_WNN_L25_HC_ModuleScores_Markers <- qread(
  paste0(
    "../Data/Annotations/Module_Scores/RNA/HC_Derived/", 
  "WNN_L25_HC_ModuleScores_Markers", 
  ".qrds"
  ), 
  nthr=nthr
)    




  ### 2.0 Write Table with N markers -------------------------------------------

write.table(
  M0_WNN_L25_HC_ModuleScores_Markers[,colnames(M0_WNN_L25_HC_ModuleScores_Markers) %in% c("CellType", "n", "fdr", "log2fc")], 
  file = paste0(
    "../Data/Visualization/Tables/SupplTab_CellType_Markers_N.txt" 
    
  ), 
  sep="\t", 
  quote=FALSE, 
  col.names = TRUE, 
  row.names = FALSE
)

