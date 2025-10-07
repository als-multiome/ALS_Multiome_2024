# A14_Integrate.R 

  ### 0.0 Load libraries ------------------------------------------------------- 

library(qs2)
library(tidyverse)
library(Seurat)
library(glmGamPoi)  

  

    
  ### 1.0 Load data ------------------------------------------------------------

FANS  <- qs_read(
  "../Data/FANS/SeuratObjects/FANS_Merged_QCFiltered_SexAdded_SampleDonorAdded.qs2"
)


Samples <- qs_read(
  "../Data/FANS/Samples/Samples.qs2"
)




  ### 2.0 Join Layers ----------------------------------------------------------

FANS
FANS <- SCTransform(FANS)
FANS <- RunPCA(FANS)

FANS <- IntegrateLayers(
  object = FANS,
  method = HarmonyIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca", 
  new.reduction = "harmony",
  verbose = T 
)




  ### 3.0 Cluster cells --------------------------------------------------------

FANS <- FindNeighbors(FANS, reduction = "harmony", dims = 1:30)
FANS <- FindClusters(FANS, cluster.name = "Seurat_Clusters")
FANS <- RunUMAP(FANS, reduction = "harmony", dims = 1:30, reduction.name = "Umap")

DimPlot(
  FANS,
  reduction = "Umap",
  group.by = c("TDP43", "ID", "Seurat_Clusters", "Sex", "Batch"),
  combine = TRUE, label.size = 2
)




  ### 4.0 Save data ------------------------------------------------------------

qs_save(
  FANS, 
  "../Data/FANS/SeuratObjects/FANS_Integrated.qs2"
)



