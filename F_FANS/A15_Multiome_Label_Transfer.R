# A10_Multiome_Label_Transfer.R 




  ### 0.0 Load libraries ------------------------------------------------------- 

library(qs)
library(qs2)
library(tidyverse)
library(Seurat)

qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)
  
  


  ### 1.0 Load data ------------------------------------------------------------

M0_RNA <- qread("../../../ALS_Brain_Multiome/Data/SeuratObjects/M0_RNA.qrds")
  
FANS <- qs_read(
  "../Data/FANS/SeuratObjects/FANS_Integrated.qs2"
)

Samples <- qs_read(
  "../Data/FANS/Samples/Samples.qs2"
)

FANS <- NormalizeData(FANS, assay = "RNA")
FANS <- JoinLayers(FANS, assay = "RNA")


FANS.anchors <- FindTransferAnchors(
  reference = M0_RNA, 
  reference.assay = "RNA",
  query.assay = "RNA", 
  query = FANS, 
  dims = 1:30,
  reference.reduction = NULL
) 

FANS_WNN_L1_Predictions <- TransferData(anchorset = FANS.anchors, refdata = M0_RNA$WNN_L1, query.assay = "RNA", dims = 1:30)
FANS_WNN_L15_Predictions <- TransferData(anchorset = FANS.anchors, refdata = M0_RNA$WNN_L15, query.assay = "RNA", dims = 1:30)
FANS_WNN_L2_Predictions <- TransferData(anchorset = FANS.anchors, refdata = M0_RNA$WNN_L2, query.assay = "RNA", dims = 1:30)
FANS_WNN_L25_Predictions <- TransferData(anchorset = FANS.anchors, refdata = M0_RNA$WNN_L25, query.assay = "RNA", dims = 1:30)
FANS_WNN_L3_Predictions <- TransferData(anchorset = FANS.anchors, refdata = M0_RNA$WNN_L3, query.assay = "RNA", dims = 1:30)
FANS_WNN_L4_Predictions <- TransferData(anchorset = FANS.anchors, refdata = M0_RNA$RNA_L4, query.assay = "RNA", dims = 1:30)

FANS <- AddMetaData(FANS, metadata = FANS_WNN_L1_Predictions$predicted.id, col.name = "ID_WNN_L1_Predicted")  
FANS <- AddMetaData(FANS, metadata = FANS_WNN_L15_Predictions$predicted.id, col.name = "ID_WNN_L15_Predicted")  
FANS <- AddMetaData(FANS, metadata = FANS_WNN_L2_Predictions$predicted.id, col.name = "ID_WNN_L2_Predicted")  
FANS <- AddMetaData(FANS, metadata = FANS_WNN_L25_Predictions$predicted.id, col.name = "ID_WNN_L25_Predicted")  
FANS <- AddMetaData(FANS, metadata = FANS_WNN_L3_Predictions$predicted.id, col.name = "ID_WNN_L3_Predicted")  
FANS <- AddMetaData(FANS, metadata = FANS_WNN_L4_Predictions$predicted.id, col.name = "ID_WNN_L4_Predicted")  

DimPlot(FANS, group.by = c("TDP43", "ID_WNN_L1_Predicted"))
DimPlot(FANS, group.by = c("TDP43", "ID_WNN_L15_Predicted"))
DimPlot(FANS, group.by = c("TDP43", "ID_WNN_L2_Predicted"))
DimPlot(FANS, group.by = c("TDP43", "ID_WNN_L25_Predicted"))
DimPlot(FANS, group.by = c("TDP43", "ID_WNN_L3_Predicted"))
DimPlot(FANS, group.by = c("TDP43", "ID_WNN_L4_Predicted"))




  ### 4.0 Export data ----------------------------------------------------------

qs_save(
  FANS, 
  "../Data/FANS/SeuratObjects/FANS_Integrated_LabelsTransferred.qs2"
)



