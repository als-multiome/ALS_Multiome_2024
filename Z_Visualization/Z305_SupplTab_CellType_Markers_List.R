
  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 

library(qs)
library(data.table) 
library(xlsx)
library(EnsDb.Hsapiens.v86) 


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


Markers_Extratelencephalic_Neurons <- fread(
  paste0(
    "../Data/Input/Markers_Extratelencephalic_Neurons", 
    ".txt"
  ), 
  header=FALSE, 
  data.table=FALSE
)   

Markers_Extratelencephalic_Neurons <- Markers_Extratelencephalic_Neurons$V1
Markers_Extratelencephalic_Neurons <- Markers_Extratelencephalic_Neurons[which(! Markers_Extratelencephalic_Neurons %in% c(
  "FAM19A1", 
  "C4orf22", 
  "FAM84B", 
  "KIAA0368"
))]
      # These markers were not found in the RNA Assay of the M0 data 




  ### 2.0 Generate a dictionary for ENSG IDs and HUGO Gene Symbols -------------

genes <- genes(EnsDb.Hsapiens.v86) 
features <- data.frame(
  SYMBOL = Markers_Extratelencephalic_Neurons
)

features$Gene_Name <- genes$gene_name[match(features$SYMBOL, genes$gene_name)]
features$ENSG <- genes$gene_id[match(features$SYMBOL, genes$gene_name)]

all(features$Gene_Name==features$Symbol, na.rm=TRUE)
features$Gene_biotype = genes$gene_biotype[match(features$ENSG, genes$gene_id)]

rm(genes)




  ### 3.0 Write XLSX Table with markers ----------------------------------------


for (CellType in sort(M0_WNN_L25_HC_ModuleScores_Markers$CellType)){
  tmp <- data.frame(
    ENSG = M0_WNN_L25_HC_ModuleScores_Markers$markers_ENSG[[which(M0_WNN_L25_HC_ModuleScores_Markers$CellType==CellType)]], 
    SYMBOL = M0_WNN_L25_HC_ModuleScores_Markers$markers_SYMB[[which(M0_WNN_L25_HC_ModuleScores_Markers$CellType==CellType)]] 
  )
  write.xlsx(
    x = tmp, 
    file = paste0(
      "../Data/Visualization/Tables/", 
      "SupplTab_CellType_Markers_List", 
      ".xlsx"
    ), 
    sheetName = str_replace_all(CellType, "_", " "), 
    col.names = TRUE, 
    row.names = FALSE, 
    append = TRUE
  )
}

write.xlsx(
  x = features[,colnames(features) %in% c("ENSG", "SYMBOL")][,c(2,1)], 
  file = paste0(
    "../Data/Visualization/Tables/", 
    "SupplTab_CellType_Markers_List", 
    ".xlsx"
  ), 
  sheetName = "Extratelencephalic Cell Markers", 
  col.names = TRUE, 
  row.names = FALSE, 
  append = TRUE
)
