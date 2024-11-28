

  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 

library(qs)
library(tidyverse)
library(Seurat)
library(EnsDb.Hsapiens.v86) 
library(tidyverse)
library(purrr) 


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

LIBD_Spatial_Spe <- qread(
  paste0(
    "../Data/SpatialData/", 
    "LIBD_Spatial_Spe", 
    ".qrds"
  ), 
  nthr=nthr
)

LIBD_Spatial_Seurat <- qread(
  paste0(
    "../Data/SpatialData/", 
    "LIBD_Spatial_Seurat", 
    ".qrds"
  ), 
  nthr=nthr
)

all(rownames(LIBD_Spatial_Spe)==rownames(LIBD_Spatial_Seurat))
all(colnames(LIBD_Spatial_Spe)==str_split(colnames(LIBD_Spatial_Seurat), "_", simplify=TRUE)[,2])


ALS_ALSFTD_WNN_L25_Signature <- qread(
  paste0(
    "../Data/Annotations/Signatures/RNA/", 
    "ALS_ALSFTD_WNN_L25_Signature", 
    ".qrds"
  ), 
  nthr=nthr
)




  ### 2.0 Generate a dictionary for ENSG IDs and HUGO Gene Symbols -------------

genes <- genes(EnsDb.Hsapiens.v86) 
features <- data.frame(
  ENSG = rownames(LIBD_Spatial_Seurat)
)

grep("\\.", features$ENS)
grep("LRG", features$ENS)

features$Gene_Name <- genes$gene_name[match(features$ENSG, genes$gene_id)]
features$Symbol <- genes$symbol[match(features$ENSG, genes$gene_id)]

all(features$Gene_Name==features$Symbol, na.rm=TRUE)
features$Gene_biotype = genes$gene_biotype[match(features$ENSG, genes$gene_id)]

features$Key <- features$ENSG
features$Key[!is.na(features$Symbol)] <- features$Symbol[!is.na(features$Symbol)]
length(grep("ENSG", features$Key, value=TRUE))==sum(is.na(features$Symbol))

grep("_", features$Key, value = TRUE)
grep("LRG", features$Key, value = TRUE)
features$Key[duplicated(features$Key)] |> sort()

rm(genes)




  ### 3.0 Generate M0_RNA ALS/ALSFTD WNN_L25 DEG signatures ModuleScores -------

M0_RNA <- AddModuleScore(
  M0_RNA, 
  features=list(
    ALS_ALSFTD_WNN_L25_Signature$SYMBOL[ALS_ALSFTD_WNN_L25_Signature$Signif_In=="ALS"]
  ), 
  ctrl=500, 
  name = "ALS_WNN_L25_DEGs_Score_"
)

M0_RNA <- AddModuleScore(
  M0_RNA, 
  features=list(
    ALS_ALSFTD_WNN_L25_Signature$SYMBOL[ALS_ALSFTD_WNN_L25_Signature$Signif_In=="ALSFTD"]
  ), 
  ctrl=500, 
  name = "ALSFTD_WNN_L25_DEGs_Score_"
)

M0_RNA <- AddModuleScore(
  M0_RNA, 
  features=list(
    ALS_ALSFTD_WNN_L25_Signature$SYMBOL[ALS_ALSFTD_WNN_L25_Signature$Signif_In=="Both"]
  ), 
  ctrl=500, 
  name = "ALS_ALSFTD_WNN_L25_DEGs_Score_"
)

M0_RNA <- AddModuleScore(
  M0_RNA, 
  features=list(
    ALS_ALSFTD_WNN_L25_Signature$SYMBOL
  ), 
  ctrl=500, 
  name = "ALS_ALSFTD_WNN_L25_DEGs_SignifInBoth_Score_"
)

M0_WNN_L25_DEG_Signatures_ModuleScores <- M0_RNA@meta.data

qsave(
  M0_WNN_L25_DEG_Signatures_ModuleScores, 
  paste0(
    "../Data/Annotations/Module_Scores/Signatures/", 
    "M0_WNN_L25_DEG_Signatures_ModuleScores", 
    ".qrds"
  ), 
  nthr=nthr
)




  ### 4.0 Generate LIBD Spatial Data WNN_L25 DEG signature Module Scores -------

fts <- ALS_ALSFTD_WNN_L25_Signature$Gene[ALS_ALSFTD_WNN_L25_Signature$Signif_ALS]
fts <- features$ENSG[match(fts, features$Symbol)]
table(is.na(fts))

LIBD_Spatial_Seurat <- AddModuleScore(
  LIBD_Spatial_Seurat, 
  features=list(fts), 
  ctrl=500, 
  name="ALS_WNN_L25_DEGs_Score_"
)

rm(fts)



fts <- ALS_ALSFTD_WNN_L25_Signature$Gene[ALS_ALSFTD_WNN_L25_Signature$Signif_ALSFTD]
fts <- features$ENSG[match(fts, features$Symbol)]
table(is.na(fts))

LIBD_Spatial_Seurat <- AddModuleScore(
  LIBD_Spatial_Seurat, 
  features=list(fts), 
  ctrl=500, 
  name="ALSFTD_WNN_L25_DEGs_Score_"
)

rm(fts)



fts <- ALS_ALSFTD_WNN_L25_Signature$Gene
fts <- features$ENSG[match(fts, features$Symbol)]
table(is.na(fts))

LIBD_Spatial_Seurat <- AddModuleScore(
  LIBD_Spatial_Seurat, 
  features=list(fts), 
  ctrl=500, 
  name="ALS_ALSFTD_WNN_L25_DEGs_Score_"
)

rm(fts) 



fts <- ALS_ALSFTD_WNN_L25_Signature$Gene[ALS_ALSFTD_WNN_L25_Signature$Signif_In=="Both"]
fts <- features$ENSG[match(fts, features$Symbol)]
table(is.na(fts))

LIBD_Spatial_Seurat <- AddModuleScore(
  LIBD_Spatial_Seurat, 
  features=list(fts), 
  ctrl=500, 
  name="ALS_ALSFTD_WNN_L25_DEGs_SignifInBoth_Score_"
)

rm(fts) 


LIBD_Spatial_Seurat_WNN_L25_DEG_Signatures_ModuleScores <- LIBD_Spatial_Seurat@meta.data

qsave(
  LIBD_Spatial_Seurat_WNN_L25_DEG_Signatures_ModuleScores, 
  paste0(
    "../Data/Annotations/Module_Scores/Signatures/", 
    "LIBD_Spatial_Seurat_WNN_L25_DEG_Signatures_ModuleScores", 
    ".qrds"
  ), 
  nthr=nthr
)
