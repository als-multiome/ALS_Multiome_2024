

  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 

library(qs)
library(tidyverse)
library(ExperimentHub)
library(Seurat)
library(EnsDb.Hsapiens.v86) 
library(tidyverse)


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



M0_RNA_HC <- subset(M0_RNA, subset=Case=="HC")

dim(M0_RNA_HC) 
table(M0_RNA_HC$Case) 




  ### 2.0 Generate a dictionary for ENSG IDs and HUGO Gene Symbols -------------

genes <- genes(EnsDb.Hsapiens.v86) 
features <- data.frame(
  ENSG = rownames(LIBD_Spatial_Seurat)
)

grep("\\.", features$ENS)

features$Gene_Name <- genes$gene_name[match(features$ENSG, genes$gene_id)]
features$Symbol <- genes$symbol[match(features$ENSG, genes$gene_id)]

all(features$Gene_Name==features$Symbol, na.rm=TRUE)
features$Gene_biotype = genes$gene_biotype[match(features$ENSG, genes$gene_id)]

features$Key <- features$ENSG
features$Key[!is.na(features$Symbol)] <- features$Symbol[!is.na(features$Symbol)]
length(grep("ENSG", features$Key, value=TRUE))==sum(is.na(features$Symbol))

grep("_", features$Key, value = TRUE)
features$Key[duplicated(features$Key)] |> sort()

rm(genes)




  ### 3.0 Generate CellType-specific signatures with the M0 HC dataset ---------



    ## 3.1 CellTypes Gene Modules WNN_L15 --------------------------------------

WNN_L15_Markers <- list()  
Idents(M0_RNA_HC) <- M0_RNA_HC$WNN_L15

for (CellType in unique(M0_RNA_HC$WNN_L15)){
  tryCatch({
    WNN_L15_Markers[[CellType]] <- FindMarkers(
      M0_RNA_HC, 
      ident.1 = CellType, 
      test.use="wilcox"
    ) 
  }, error=function(e){
    WNN_L15_Markers[[CellType]] <- "NA"
  }
  )
}
rm(CellType) 


WNN_L15_Markers <- lapply(
  WNN_L15_Markers, 
  FUN=function(x){
    tryCatch({
      x$SYMBOL <- rownames(x) 
      x$ENSG <- features$ENSG[match(x$SYMBOL, features$Symbol)] 
      x
    }, error=function(e){}
    )
  }
)

DefaultAssay(LIBD_Spatial_Seurat) <- "SpatialLogCounts"


for (CellType in names(WNN_L15_Markers)){
  
  message("##################################################") 
  message(
    paste0(
      Sys.time(), 
      ": ", 
      "Processing WNN_L15 CellType ", 
      CellType, 
      "..."
    )
  )
  
  tryCatch({
    
    fs <- WNN_L15_Markers[[CellType]] 
    fs <- fs$ENSG[fs$p_val_adj < 1e-10 & fs$avg_log2FC > 2.32]
    fs <- fs[!is.na(fs)]
    fs_s <- features$Symbol[match(fs, features$ENSG)] 
    stopifnot(!any(is.na(fs_s))) 
    
    
    M0_RNA_HC <- AddModuleScore(
      M0_RNA_HC, 
      features=list(fs_s), 
      name=paste0("WNN_L15_", CellType, "_ModuleScore_50fts_"), 
      ctrl=50 , 
      assay = "RNA"
    )
    
    LIBD_Spatial_Seurat <- AddModuleScore(
      LIBD_Spatial_Seurat, 
      features=list(fs), 
      name=paste0("WNN_L15_", CellType, "_ModuleScore_50fts_"), 
      ctrl=50 , 
      assay = "SpatialLogCounts"
    ) 
    
    LIBD_Spatial_Spe@colData[[paste0("WNN_L15_", CellType, "_ModuleScore_50fts_1")]] <- LIBD_Spatial_Seurat[[paste0("WNN_L15_", CellType, "_ModuleScore_50fts_1")]]
    
    df <- data.frame(
      SYMBOL = fs_s, 
      ENSEMBLG = fs 
    )
    
    xlsx::write.xlsx(
      df, 
      file = paste0(
        "../Data/Annotations/CellType_Markers/WNN_L15_Signatures_Wilcox_50fts_OnlyBothENSGandSYMBavail", 
        ".xlsx"
      ), 
      sheetName = CellType, 
      col.names = TRUE, 
      row.names = FALSE, 
      append = TRUE, 
      showNA = FALSE
    )
    
    rm(df, fs, fs_s)
    
  }, error=function(e){
    message(
      paste0(
        "Error processing cell type ", 
        CellType, 
        "... "
      )
    )
  })

}




    ## 3.2 CellTypes Gene Modules WNN_L25 --------------------------------------

WNN_L25_Markers <- list()  
Idents(M0_RNA_HC) <- M0_RNA_HC$WNN_L25

for (CellType in unique(M0_RNA_HC$WNN_L25)){
  tryCatch({
    WNN_L25_Markers[[CellType]] <- FindMarkers(
      M0_RNA_HC, 
      ident.1 = CellType, 
      test.use="wilcox"
    ) 
  }, error=function(e){
    WNN_L25_Markers[[CellType]] <- "NA"
  }
  )
}
rm(CellType) 


WNN_L25_Markers <- lapply(
  WNN_L25_Markers, 
  FUN=function(x){
    tryCatch({
      x$SYMBOL <- rownames(x) 
      x$ENSG <- features$ENSG[match(x$SYMBOL, features$Symbol)] 
      x
    }, error=function(e){}
    )
  }
)

DefaultAssay(LIBD_Spatial_Seurat) <- "SpatialLogCounts"

 
for (CellType in names(WNN_L25_Markers)){
  
  message("##################################################") 
  message(
    paste0(
      Sys.time(), 
      ": ", 
      "Processing WNN_L25 CellType ", 
      CellType, 
      "..."
    )
  )
  
  tryCatch({
    
    fs <- WNN_L25_Markers[[CellType]] 
    fs <- fs$ENSG[fs$p_val_adj < 1e-10 & fs$avg_log2FC > 2.32]
    fs <- fs[!is.na(fs)]
    fs_s <- features$Symbol[match(fs, features$ENSG)] 
    stopifnot(!any(is.na(fs_s))) 
    
    
    M0_RNA_HC <- AddModuleScore(
      M0_RNA_HC, 
      features=list(fs_s), 
      name=paste0("WNN_L25_", CellType, "_ModuleScore_50fts_"), 
      ctrl=50 , 
      assay = "RNA"
    )
    
    LIBD_Spatial_Seurat <- AddModuleScore(
      LIBD_Spatial_Seurat, 
      features=list(fs), 
      name=paste0("WNN_L25_", CellType, "_ModuleScore_50fts_"), 
      ctrl=50 , 
      assay = "SpatialLogCounts"
    ) 
    
    LIBD_Spatial_Spe@colData[[paste0("WNN_L25_", CellType, "_ModuleScore_50fts_1")]] <- LIBD_Spatial_Seurat[[paste0("WNN_L25_", CellType, "_ModuleScore_50fts_1")]]
    
    df <- data.frame(
      SYMBOL = fs_s, 
      ENSEMBLG = fs 
    )
    
    xlsx::write.xlsx(
      df, 
      file = paste0(
        "../Data/Annotations/CellType_Markers/WNN_L25_Signatures_Wilcox_50fts_OnlyBothENSGandSYMBavail", 
        ".xlsx"
      ), 
      sheetName = CellType, 
      col.names = TRUE, 
      row.names = FALSE, 
      append = TRUE, 
      showNA = FALSE
    )
    
    rm(df, fs, fs_s)
    
    
  }, error=function(e){
    message(
      paste0(
        "Error processing cell type ", 
        CellType, 
        "... "
      )
    )
  })
  
}
rm(CellType) 




    ## 3.3 CellTypes Gene Modules WNN_L4 --------------------------------------

WNN_L4_Markers <- list()  
Idents(M0_RNA_HC) <- M0_RNA_HC$WNN_L4

for (CellType in unique(M0_RNA_HC$WNN_L4)){
  tryCatch({
    WNN_L4_Markers[[CellType]] <- FindMarkers(
      M0_RNA_HC, 
      ident.1 = CellType, 
      test.use="wilcox"
    ) 
  }, error=function(e){
    WNN_L4_Markers[[CellType]] <- "NA"
  }
  )
}
rm(CellType) 


WNN_L4_Markers <- lapply(
  WNN_L4_Markers, 
  FUN=function(x){
    tryCatch({
      x$SYMBOL <- rownames(x) 
      x$ENSG <- features$ENSG[match(x$SYMBOL, features$Symbol)] 
      x
    }, error=function(e){})
  }
)

DefaultAssay(LIBD_Spatial_Seurat) <- "SpatialLogCounts"


for (CellType in names(WNN_L4_Markers)){
  
  message("##################################################") 
  message(
    paste0(
      Sys.time(), 
      ": ", 
      "Processing WNN_L4 CellType ", 
      CellType, 
      "..."
    )
  )
  
  tryCatch({
    
    fs <- WNN_L4_Markers[[CellType]] 
    fs <- fs$ENSG[fs$p_val_adj < 1e-10 & fs$avg_log2FC > 2.32]
    fs <- fs[!is.na(fs)]
    fs_s <- features$Symbol[match(fs, features$ENSG)] 
    stopifnot(!any(is.na(fs_s))) 
    
    
    M0_RNA_HC <- AddModuleScore(
      M0_RNA_HC, 
      features=list(fs_s), 
      name=paste0("WNN_L4_", CellType, "_ModuleScore_50fts_"), 
      ctrl=50 , 
      assay = "RNA"
    )
    
    LIBD_Spatial_Seurat <- AddModuleScore(
      LIBD_Spatial_Seurat, 
      features=list(fs), 
      name=paste0("WNN_L4_", CellType, "_ModuleScore_50fts_"), 
      ctrl=50 , 
      assay = "SpatialLogCounts"
    ) 
    
    LIBD_Spatial_Spe@colData[[paste0("WNN_L4_", CellType, "_ModuleScore_50fts_1")]] <- LIBD_Spatial_Seurat[[paste0("WNN_L4_", CellType, "_ModuleScore_50fts_1")]]
    
    df <- data.frame(
      SYMBOL = fs_s, 
      ENSEMBLG = fs 
    )
    
    xlsx::write.xlsx(
      df, 
      file = paste0(
        "../Data/Annotations/CellType_Markers/WNN_L4_Signatures_Wilcox_50fts_OnlyBothENSGandSYMBavail", 
        ".xlsx"
      ), 
      sheetName = CellType, 
      col.names = TRUE, 
      row.names = FALSE, 
      append = TRUE, 
      showNA = FALSE
    )
    
    rm(df, fs, fs_s)
    
  }, error=function(e){
    message(
      paste0(
        "Error processing cell type ", 
        CellType, 
        "... "
      )
    )
  })
}
rm(CellType) 




  ### 4.0 Save data ------------------------------------------------------------


qsave(
  M0_RNA_HC, 
  paste0(
    "../Data/SeuratObjects/", 
    "M0_RNA_HC_WithCellType_Signatures", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  LIBD_Spatial_Spe, 
  paste0(
    "../Data/SpatialData/", 
    "LIBD_Spatial_Spe", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  LIBD_Spatial_Seurat, 
  paste0(
    "../Data/SpatialData/", 
    "LIBD_Spatial_Seurat", 
    ".qrds"
  ), 
  nthr=nthr
)


qsave(
  WNN_L15_Markers, 
  paste0(
    "../Data/Annotations/CellType_Markers/", 
    "WNN_L25_Markers_Wilcox", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  WNN_L25_Markers, 
  paste0(
    "../Data/Annotations/CellType_Markers/", 
    "WNN_L15_Markers_Wilcox", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  WNN_L4_Markers, 
  paste0(
    "../Data/Annotations/CellType_Markers/", 
    "WNN_L15_Markers_Wilcox", 
    ".qrds"
  ), 
  nthr=nthr
)
