
  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 

library(qs) 
library(spatialLIBD)
library(ExperimentHub)
library(Seurat)
library(EnsDb.Hsapiens.v86) 


qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)




  ### 1.0 Load data ------------------------------------------------------------

spe <- fetch_data(type = "spe", eh = ExperimentHub())




  ### 2.0 Generate a SeuratObject for the spatial data -------------------------

vis_counts <- spe@assays@data$counts
all(colnames(vis_counts)==rownames(spe@colData))
colnames(vis_counts) <- spe@colData$key

vis_logcounts <- spe@assays@data$logcounts
all(colnames(vis_logcounts)==rownames(spe@colData))
colnames(vis_logcounts) <- spe@colData$key

metadata <- spe@colData
metadata$Original_Barcode <- rownames(metadata)
rownames(metadata) <- metadata$key


Spatial_Counts <- CreateAssayObject(counts=vis_counts)
Spatial_LogCounts <- CreateAssayObject(counts=vis_logcounts)
all(colnames(Spatial_Counts)==colnames(Spatial_LogCounts))
all(rownames(Spatial_Counts)==rownames(Spatial_LogCounts))
all(colnames(Spatial_Counts)==rownames(metadata))


Visium <- CreateSeuratObject(Spatial_Counts)
all(rownames(Visium@meta.data)==rownames(metadata))

Visium@meta.data <- cbind(Visium@meta.data, metadata)
Visium@assays$SpatialCounts <- Visium@assays$RNA
Visium@assays$SpatialLogCounts <- Spatial_LogCounts
DefaultAssay(Visium) <- "SpatialLogCounts"
Visium@assays$RNA <- NULL 

rm(vis_counts, vis_logcounts, Spatial_Counts, Spatial_LogCounts, metadata) 




  ### 3.0 Save Spatial Data Objects --------------------------------------------

qsave(
  spe, 
  paste0(
    "../Data/SpatialData/", 
    "LIBD_Spatial_Spe", 
    ".qrds"
  ), 
  nthr = nthr
)

qsave(
  Visium, 
  paste0(
    "../Data/SpatialData/", 
    "LIBD_Spatial_Seurat", 
    ".qrds"
  ), 
  nthr=nthr
)




