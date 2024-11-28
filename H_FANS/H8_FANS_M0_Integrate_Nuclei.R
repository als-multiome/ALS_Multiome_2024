
library(qs) 
library(Seurat)
library(tidyverse) 
library(SeuratObject) 
library(MAST)
library(zinbwave)

library(patchwork)
nthr=46



  ### 1.0 Load data ------------------------------------------------------------

M0_RNA <- qread(
  paste0(
    "../Data/SeuratObjects/", 
    "M0_RNA", 
    ".qrds"
  ), nthr=46
)
FC1 <- Read10X(
  "/home/veselin/Bioinfo_Data/Cellranger_Output/Human/scRNA/scRNA_FACS_Cortex_ALSFTS_Sarah/Cellranger_Count_Optimized/FC1/outs/filtered_feature_bc_matrix/"
)
FC1 <- CreateSeuratObject(FC1)
FC1$orig.ident <- "FC1"

FC2 <- Read10X(
  "/home/veselin/Bioinfo_Data/Cellranger_Output/Human/scRNA/scRNA_FACS_Cortex_ALSFTS_Sarah/Cellranger_Count_Optimized/FC2/outs/filtered_feature_bc_matrix/"
)
FC2 <- CreateSeuratObject(FC2)
FC2$orig.ident <- "FC2"

FC3 <- Read10X(
  "/home/veselin/Bioinfo_Data/Cellranger_Output/Human/scRNA/scRNA_FACS_Cortex_ALSFTS_Sarah/Cellranger_Count_Optimized/FC3/outs/filtered_feature_bc_matrix/"
)
FC3 <- CreateSeuratObject(FC3)
FC3$orig.ident <- "FC3"

FC4 <- Read10X(
  "/home/veselin/Bioinfo_Data/Cellranger_Output/Human/scRNA/scRNA_FACS_Cortex_ALSFTS_Sarah/Cellranger_Count_Optimized/FC4/outs/filtered_feature_bc_matrix/"
)
FC4 <- CreateSeuratObject(FC4)
FC4$orig.ident <- "FC4"

FC5 <- Read10X(
  "/home/veselin/Bioinfo_Data/Cellranger_Output/Human/scRNA/scRNA_FACS_Cortex_ALSFTS_Sarah/Cellranger_Count_Optimized/FC5/outs/filtered_feature_bc_matrix/"
)
FC5 <- CreateSeuratObject(FC5)
FC5$orig.ident <- "FC5"

FC6 <- Read10X(
  "/home/veselin/Bioinfo_Data/Cellranger_Output/Human/scRNA/scRNA_FACS_Cortex_ALSFTS_Sarah/Cellranger_Count_Optimized/FC6/outs/filtered_feature_bc_matrix/"
)
FC6 <- CreateSeuratObject(FC6)
FC6$orig.ident <- "FC6"

FC7 <- Read10X(
  "/home/veselin/Bioinfo_Data/Cellranger_Output/Human/scRNA/scRNA_FACS_Cortex_ALSFTS_Sarah/Cellranger_Count_Optimized/FC7/outs/filtered_feature_bc_matrix/"
)
FC7 <- CreateSeuratObject(FC7)
FC7$orig.ident <- "FC7"

FC8 <- Read10X(
  "/home/veselin/Bioinfo_Data/Cellranger_Output/Human/scRNA/scRNA_FACS_Cortex_ALSFTS_Sarah/Cellranger_Count_Optimized/FC8/outs/filtered_feature_bc_matrix/"
)
FC8 <- CreateSeuratObject(FC8)
FC8$orig.ident <- "FC8"




  ### 2.0 Create Merged Seurat Objects -----------------------------------------

FC12 <- merge(FC1, y=FC2)
FC12[["RNA"]] <- JoinLayers(FC12[["RNA"]])


FC34 <- merge(FC3, y=FC4)
FC34[["RNA"]] <- JoinLayers(FC34[["RNA"]])


FC78 <- merge(FC7, y=FC8)
FC78[["RNA"]] <- JoinLayers(FC78[["RNA"]])


F1 <- merge(FC12, y = c(FC34), add.cell.ids = c("FC12", "FC34"), project = "FANSSeq")
F1




  ### 3.0 Add PctMito and filter cells -----------------------------------------

F1[["PctMt"]] <- PercentageFeatureSet(F1, pattern = "^MT-")
summary(F1$PctMt)
table(F1$PctMt < 5)
table(F1$nFeature_RNA > 100, F1$orig.ident)
table(F1$nFeature_RNA < 2500, F1$orig.ident)
F2 <- subset(F1, subset = nFeature_RNA > 100 & nFeature_RNA < 2500 & PctMt < 5)




  ### 4.0 Normalize Data -------------------------------------------------------

F2 <- NormalizeData(F2) 




  ### 5.0 Add dimensionality reductions ----------------------------------------

F2 <- FindVariableFeatures(F2)
F2 <- ScaleData(F2)
F2 <- RunPCA(F2)

F2 <- RunUMAP(F2, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(F2, reduction = "umap.unintegrated", group.by = c("orig.ident")) 




  ### 6.0 Integrate Layers and join them ---------------------------------------

F2 <- IntegrateLayers(
  object = F2, 
  method = CCAIntegration, 
  orig.reduction = "pca", 
  new.reduction = "integrated.cca",
  verbose = TRUE, 
  dims=1:20, 
  k.weight=45
  )


F2[["RNA"]] <- JoinLayers(F2[["RNA"]])




  ### 7.0 Add dimentionality reduction of merged integrated object ------------- 

F2 <- FindNeighbors(F2, reduction = "integrated.cca", dims = 1:30)
F2 <- FindClusters(F2, resolution = 2)

F2 <- RunUMAP(F2, dims = 1:30, reduction = "integrated.cca", reduction.name = "UMAP_Integrated")




  ### 8.0 Add TDP43 FANS status data -------------------------------------------

F2$TDP43 <- "TDP43_Low"
F2$TDP43[F2$orig.ident %in% c("FC3", "FC4")] <- "TDP43_High"
DimPlot(F2, group.by="TDP43", reduction = "umap.unintegrated")

Idents(F2) <- F2$TDP43




  ### 9.0 Add Sex metadata -----------------------------------------------------

F2$Sex <- "M"
F2$Sex[F2$orig.ident=="FC2"] <- "F"




  ### 10.0 F2 Transfer Cell-Type Labels from M0 --------------------------------

M0_RNA.anchors <- FindTransferAnchors(reference = M0_RNA, query = F2, dims = 1:30,
                                      reference.reduction = NULL)

F2_M0_WNN_L25_Predictions <- TransferData(anchorset = M0_RNA.anchors, refdata = M0_RNA$WNN_L25, dims = 1:30)
F2_M0_WNN_L4_Predictions <- TransferData(anchorset = M0_RNA.anchors, refdata = M0_RNA$WNN_L4, dims = 1:30)
colnames(F2_M0_WNN_L25_Predictions)[1] <- "Predicted_ID_M0_WNN_L25"
colnames(F2_M0_WNN_L4_Predictions)[1] <- "Predicted_ID_M0_WNN_L4"


F2 <- AddMetaData(F2, metadata = F2_M0_WNN_L25_Predictions["Predicted_ID_M0_WNN_L25"])
F2 <- AddMetaData(F2, metadata = F2_M0_WNN_L4_Predictions["Predicted_ID_M0_WNN_L4"]) 




  ### 11.0 Export data ---------------------------------------------------------

qsave(
  F2, 
  paste0(
    "../Data/SeuratObjects/", 
    "F2", 
    ".qrds"
  ), 
  nthr=nthr
)


