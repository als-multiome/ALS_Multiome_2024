# A9_Merge_Seurat_Objects.R 

  ### 0.0 Load libraries ------------------------------------------------------- 

library(data.table)
library(qs2)
library(tidyverse)
library(Seurat)
  

  
  ### 1.0 Load data ------------------------------------------------------------
  
Samples <- qs_read(
  "../Data/FANS/Samples/Samples.qs2"
)

FANS_Seurats <- qs_read(
  "../Data/FANS/SeuratObjects/FANS_Single_SeuratObjects_Filtered_Singlets.qs2"
)




  ### 2.0 Merge Seurat Objects -------------------------------------------------

FANS <- merge(
  FANS_Seurats[[1]], 
  y = FANS_Seurats[c(2:length(FANS_Seurats))],
  add.cell.ids = names(FANS_Seurats)
)

FANS$ID <- FANS$orig.ident
FANS$Sample <- Samples$Name[match(FANS$ID, Samples$ID)]
FANS$Batch <- Samples$Batch[match(FANS$ID, Samples$ID)]
FANS$TDP43 <- Samples$TDP43[match(FANS$ID, Samples$ID)]

table(FANS$TDP43)




  ### 3.0 Save data ------------------------------------------------------------

qs_save(
  FANS, 
  "../Data/FANS/SeuratObjects/FANS_Merged.qs2"
)



