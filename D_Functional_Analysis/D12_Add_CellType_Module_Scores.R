

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




  ### 3.0 Generate CellType-specific signatures with HC samples ----------------



    ## 3.1 CellTypes Gene Modules WNN_L15 --------------------------------------

WNN_L15_Markers <- qread(
  paste0(
    "../Data/Annotations/CellType_Markers/", 
    "WNN_L15_Markers_HC_Wilcox", 
    ".qrds"
  ), 
  nthr=nthr
)

WNN_L15_HC_ModuleScores_Markers <- data.frame(
  CellTypeLevel = "WNN_L15",
  CellType = names(WNN_L15_Markers),
  n = NA, 
  fdr = NA, 
  log2fc = NA
)
WNN_L15_HC_ModuleScores_Markers$markers_ENSG = list(rep(NA, length(WNN_L15_Markers))) 
WNN_L15_HC_ModuleScores_Markers$markers_SYMB = list(rep(NA, length(WNN_L15_Markers)))

imap(
  WNN_L15_Markers, 
  .f=function(x, y){
   
    log = NULL
    fdr = NULL
    for (l in log2(c(5,4,3,2))){
      for (f in c(1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 5e-2)){
        if(sum(x$p_val_adj < f & x$avg_log2FC > l) >= 50){
          fdr <- f 
          break
        }
      }
      if(!is.null(fdr)){
        log <- l 
        break
      }
    } 
    
    if(!is.null(log) & !is.null(fdr)){
      
      ind <- which(x$avg_log2FC > log & x$p_val_adj < fdr)
      
      symb <- x$SYMBOL[ind] 
      ensg <- features$ENSG[match(symb, features$Symbol)]
      
      symb <- symb[!is.na(ensg)] 
      ensg <- ensg[!is.na(ensg)] 
      
      stopifnot(!any(is.na(symb))) 
      stopifnot(length(ensg)==length(symb))
      i = which(WNN_L15_HC_ModuleScores_Markers$CellType==y)
      WNN_L15_HC_ModuleScores_Markers$n[i] <<- length(ensg)
      WNN_L15_HC_ModuleScores_Markers$fdr[i] <<- fdr 
      WNN_L15_HC_ModuleScores_Markers$log2fc[i] <<- log 
      WNN_L15_HC_ModuleScores_Markers$markers_ENSG[[i]] <<- ensg
      WNN_L15_HC_ModuleScores_Markers$markers_SYMB[[i]] <<- symb      
    }
  }
)

qsave(
  WNN_L15_HC_ModuleScores_Markers, 
  paste0(
    "../Data/Annotations/Module_Scores/RNA/HC_Derived/", 
    "WNN_L15_HC_ModuleScores_Markers", 
    ".qrds"
  ), 
  nthr = nthr
)

rm(WNN_L15_Markers)

M0_WNN_L15_HC_ModuleScores <- data.frame(
  CellId=M0_RNA$CellId, 
  gex_id = M0_RNA$gex_barcode, 
  atac_barcode = M0_RNA$atac_barcode, 
  WNN_L1 = M0_RNA$WNN_L1, 
  WNN_L15 = M0_RNA$WNN_L15, 
  WNN_L2 = M0_RNA$WNN_L2, 
  WNN_L25 = M0_RNA$WNN_L25, 
  WNN_L3 = M0_RNA$WNN_L3, 
  WNN_L4 = M0_RNA$WNN_L4, 
  Case = M0_RNA$Case, 
  Case_Type = M0_RNA$Case_Type
)
rownames(M0_WNN_L15_HC_ModuleScores) <- colnames(M0_RNA)

LIBD_Spatial_Seurat_WNN_L15_HC_ModuleScores <- data.frame(
  orig.ident=LIBD_Spatial_Seurat$orig.ident, 
  sample_id = LIBD_Spatial_Seurat$sample_id, 
  subject = LIBD_Spatial_Seurat$subject, 
  position = LIBD_Spatial_Seurat$position, 
  replicate = LIBD_Spatial_Seurat$replicate, 
  key = LIBD_Spatial_Seurat$key, 
  cell_count = LIBD_Spatial_Seurat$cell_count, 
  Original_Barcode  = LIBD_Spatial_Seurat$Original_Barcode, 
  layer_guess = LIBD_Spatial_Seurat$layer_guess, 
  layer_guess_reordered = LIBD_Spatial_Seurat$layer_guess_reordered, 
  layer_guess_reordered_short = LIBD_Spatial_Seurat$layer_guess_reordered_short
)
rownames(LIBD_Spatial_Seurat_WNN_L15_HC_ModuleScores) <- colnames(LIBD_Spatial_Seurat)


for (CellType in WNN_L15_HC_ModuleScores_Markers$CellType){
  
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
    
    for (nfts in c(50, 100, 200, 500, 1234)){
      if(nfts==1234){nfts <- length(WNN_L15_HC_ModuleScores_Markers$markers_SYMB[[which(WNN_L15_HC_ModuleScores_Markers$CellType==CellType)]]) }
      message(paste0("Number of features: ", nfts))
      stopifnot(all(colnames(M0_RNA)==M0_WNN_L15_HC_ModuleScores$CellId))
      name = paste0("WNN_L15_", CellType, "_ModuleScore_", nfts, "fts_")
      fts = WNN_L15_HC_ModuleScores_Markers$markers_SYMB[which(WNN_L15_HC_ModuleScores_Markers$CellType==CellType)] 
      fts_ensg = WNN_L15_HC_ModuleScores_Markers$markers_ENSG[which(WNN_L15_HC_ModuleScores_Markers$CellType==CellType)]  
      
      M0_WNN_L15_HC_ModuleScores[[paste0(name, "1")]] <- AddModuleScore(
        M0_RNA, 
        features=fts, 
        name=name, 
        ctrl=nfts , 
        assay = "RNA"
      )[[paste0(name, "1")]][,1]
      
      LIBD_Spatial_Seurat_WNN_L15_HC_ModuleScores[[paste0(name, "1")]] <- AddModuleScore(
        LIBD_Spatial_Seurat, 
        features=fts_ensg, 
        name=name, 
        ctrl=nfts , 
        assay = "SpatialLogCounts"
      )[[paste0(name, "1")]][,1] 
      
      
      rm(name, fts, fts_ensg)
    }
  }, error=function(e){
  
  })
}    
    
qsave(
  M0_WNN_L15_HC_ModuleScores, 
  paste0(
    "../Data/Annotations/Module_Scores/RNA/HC_Derived/", 
    "M0_WNN_L15_HC_ModuleScores", 
    ".qrds"
  ), 
  nthr=nthr
)  
    
qsave(
  LIBD_Spatial_Seurat_WNN_L15_HC_ModuleScores, 
  paste0(
    "../Data/Annotations/Module_Scores/RNA/HC_Derived/", 
    "LIBD_Spatial_Seurat_WNN_L15_HC_ModuleScores", 
    ".qrds"
  ), 
  nthr=nthr
)    

#rm(M0_WNN_L15_HC_ModuleScores, LIBD_Spatial_Seurat_WNN_L15_HC_ModuleScores)  

    ## 3.2 CellTypes Gene Modules WNN_L25 --------------------------------------

WNN_L25_Markers <- qread(
  paste0(
    "../Data/Annotations/CellType_Markers/", 
    "WNN_L25_Markers_HC_Wilcox", 
    ".qrds"
  ), 
  nthr=nthr
)

WNN_L25_HC_ModuleScores_Markers <- data.frame(
  CellTypeLevel = "WNN_L25",
  CellType = names(WNN_L25_Markers),
  n = NA, 
  fdr = NA, 
  log2fc = NA
)
WNN_L25_HC_ModuleScores_Markers$markers_ENSG = list(rep(NA, length(WNN_L25_Markers))) 
WNN_L25_HC_ModuleScores_Markers$markers_SYMB = list(rep(NA, length(WNN_L25_Markers)))

imap(
  WNN_L25_Markers, 
  .f=function(x, y){
    
    log = NULL
    fdr = NULL
    for (l in log2(c(5,4,3,2))){
      for (f in c(1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 5e-2)){
        if(sum(x$p_val_adj < f & x$avg_log2FC > l) >= 50){
          fdr <- f 
          break
        }
      }
      if(!is.null(fdr)){
        log <- l 
        break
      }
    } 
    
    if(!is.null(log) & !is.null(fdr)){
      
      ind <- which(x$avg_log2FC > log & x$p_val_adj < fdr)
      
      symb <- x$SYMBOL[ind] 
      ensg <- features$ENSG[match(symb, features$Symbol)]
      
      symb <- symb[!is.na(ensg)] 
      ensg <- ensg[!is.na(ensg)] 
      
      stopifnot(!any(is.na(symb))) 
      stopifnot(length(ensg)==length(symb))
      i = which(WNN_L25_HC_ModuleScores_Markers$CellType==y)
      WNN_L25_HC_ModuleScores_Markers$n[i] <<- length(ensg)
      WNN_L25_HC_ModuleScores_Markers$fdr[i] <<- fdr 
      WNN_L25_HC_ModuleScores_Markers$log2fc[i] <<- log 
      WNN_L25_HC_ModuleScores_Markers$markers_ENSG[[i]] <<- ensg
      WNN_L25_HC_ModuleScores_Markers$markers_SYMB[[i]] <<- symb      
    }
  }
)

qsave(
  WNN_L25_HC_ModuleScores_Markers, 
  paste0(
    "../Data/Annotations/Module_Scores/RNA/HC_Derived/", 
    "WNN_L25_HC_ModuleScores_Markers", 
    ".qrds"
  ), 
  nthr = nthr
)

rm(WNN_L25_Markers)

M0_WNN_L25_HC_ModuleScores <- data.frame(
  CellId=M0_RNA$CellId, 
  gex_id = M0_RNA$gex_barcode, 
  atac_barcode = M0_RNA$atac_barcode, 
  WNN_L1 = M0_RNA$WNN_L1, 
  WNN_L15 = M0_RNA$WNN_L15, 
  WNN_L2 = M0_RNA$WNN_L2, 
  WNN_L25 = M0_RNA$WNN_L25, 
  WNN_L3 = M0_RNA$WNN_L3, 
  WNN_L4 = M0_RNA$WNN_L4, 
  Case = M0_RNA$Case, 
  Case_Type = M0_RNA$Case_Type
)
rownames(M0_WNN_L25_HC_ModuleScores) <- colnames(M0_RNA)

LIBD_Spatial_Seurat_WNN_L25_HC_ModuleScores <- data.frame(
  orig.ident=LIBD_Spatial_Seurat$orig.ident, 
  sample_id = LIBD_Spatial_Seurat$sample_id, 
  subject = LIBD_Spatial_Seurat$subject, 
  position = LIBD_Spatial_Seurat$position, 
  replicate = LIBD_Spatial_Seurat$replicate, 
  key = LIBD_Spatial_Seurat$key, 
  cell_count = LIBD_Spatial_Seurat$cell_count, 
  Original_Barcode  = LIBD_Spatial_Seurat$Original_Barcode, 
  layer_guess = LIBD_Spatial_Seurat$layer_guess, 
  layer_guess_reordered = LIBD_Spatial_Seurat$layer_guess_reordered, 
  layer_guess_reordered_short = LIBD_Spatial_Seurat$layer_guess_reordered_short
)
rownames(LIBD_Spatial_Seurat_WNN_L25_HC_ModuleScores) <- colnames(LIBD_Spatial_Seurat)


for (CellType in WNN_L25_HC_ModuleScores_Markers$CellType){
  
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
    
    for (nfts in c(50, 100, 200, 500, 1234)){
      if(nfts==1234){nfts <- length(WNN_L25_HC_ModuleScores_Markers$markers_SYMB[[which(WNN_L25_HC_ModuleScores_Markers$CellType==CellType)]]) }
      message(paste0("Number of features: ", nfts))
      stopifnot(all(colnames(M0_RNA)==M0_WNN_L25_HC_ModuleScores$CellId))
      name = paste0("WNN_L25_", CellType, "_ModuleScore_", nfts, "fts_")
      fts = WNN_L25_HC_ModuleScores_Markers$markers_SYMB[which(WNN_L25_HC_ModuleScores_Markers$CellType==CellType)] 
      fts_ensg = WNN_L25_HC_ModuleScores_Markers$markers_ENSG[which(WNN_L25_HC_ModuleScores_Markers$CellType==CellType)]  
      
      M0_WNN_L25_HC_ModuleScores[[paste0(name, "1")]] <- AddModuleScore(
        M0_RNA, 
        features=fts, 
        name=name, 
        ctrl=nfts , 
        assay = "RNA"
      )[[paste0(name, "1")]][,1]
      
      LIBD_Spatial_Seurat_WNN_L25_HC_ModuleScores[[paste0(name, "1")]] <- AddModuleScore(
        LIBD_Spatial_Seurat, 
        features=fts_ensg, 
        name=name, 
        ctrl=nfts , 
        assay = "SpatialLogCounts"
      )[[paste0(name, "1")]][,1] 
      
      
      rm(name, fts, fts_ensg)
    }
  }, error=function(e){
    
  })
}    

qsave(
  M0_WNN_L25_HC_ModuleScores, 
  paste0(
    "../Data/Annotations/Module_Scores/RNA/HC_Derived/", 
    "M0_WNN_L25_HC_ModuleScores", 
    ".qrds"
  ), 
  nthr=nthr
)  

qsave(
  LIBD_Spatial_Seurat_WNN_L25_HC_ModuleScores, 
  paste0(
    "../Data/Annotations/Module_Scores/RNA/HC_Derived/", 
    "LIBD_Spatial_Seurat_WNN_L25_HC_ModuleScores", 
    ".qrds"
  ), 
  nthr=nthr
)    

#rm(M0_WNN_L25_HC_ModuleScores, LIBD_Spatial_Seurat_WNN_L25_HC_ModuleScores)  

    ## 3.3 CellTypes Gene Modules WNN_L4 ---------------------------------------

WNN_L4_Markers <- qread(
  paste0(
    "../Data/Annotations/CellType_Markers/", 
    "WNN_L4_Markers_HC_Wilcox", 
    ".qrds"
  ), 
  nthr=nthr
)

WNN_L4_HC_ModuleScores_Markers <- data.frame(
  CellTypeLevel = "WNN_L4",
  CellType = names(WNN_L4_Markers),
  n = NA, 
  fdr = NA, 
  log2fc = NA
)
WNN_L4_HC_ModuleScores_Markers$markers_ENSG = list(rep(NA, length(WNN_L4_Markers))) 
WNN_L4_HC_ModuleScores_Markers$markers_SYMB = list(rep(NA, length(WNN_L4_Markers)))

imap(
  WNN_L4_Markers, 
  .f=function(x, y){
    
    log = NULL
    fdr = NULL
    for (l in log2(c(5,4,3,2))){
      for (f in c(1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 5e-2)){
        if(sum(x$p_val_adj < f & x$avg_log2FC > l) >= 50){
          fdr <- f 
          break
        }
      }
      if(!is.null(fdr)){
        log <- l 
        break
      }
    } 
    
    if(!is.null(log) & !is.null(fdr)){
      
      ind <- which(x$avg_log2FC > log & x$p_val_adj < fdr)
      
      symb <- x$SYMBOL[ind] 
      ensg <- features$ENSG[match(symb, features$Symbol)]
      
      symb <- symb[!is.na(ensg)] 
      ensg <- ensg[!is.na(ensg)] 
      
      stopifnot(!any(is.na(symb))) 
      stopifnot(length(ensg)==length(symb))
      i = which(WNN_L4_HC_ModuleScores_Markers$CellType==y)
      WNN_L4_HC_ModuleScores_Markers$n[i] <<- length(ensg)
      WNN_L4_HC_ModuleScores_Markers$fdr[i] <<- fdr 
      WNN_L4_HC_ModuleScores_Markers$log2fc[i] <<- log 
      WNN_L4_HC_ModuleScores_Markers$markers_ENSG[[i]] <<- ensg
      WNN_L4_HC_ModuleScores_Markers$markers_SYMB[[i]] <<- symb      
    }
  }
)

qsave(
  WNN_L4_HC_ModuleScores_Markers, 
  paste0(
    "../Data/Annotations/Module_Scores/RNA/HC_Derived/", 
    "WNN_L4_HC_ModuleScores_Markers", 
    ".qrds"
  ), 
  nthr = nthr
)

rm(WNN_L4_Markers)

M0_WNN_L4_HC_ModuleScores <- data.frame(
  CellId=M0_RNA$CellId, 
  gex_id = M0_RNA$gex_barcode, 
  atac_barcode = M0_RNA$atac_barcode, 
  WNN_L1 = M0_RNA$WNN_L1, 
  WNN_L15 = M0_RNA$WNN_L15, 
  WNN_L2 = M0_RNA$WNN_L2, 
  WNN_L25 = M0_RNA$WNN_L25, 
  WNN_L3 = M0_RNA$WNN_L3, 
  WNN_L4 = M0_RNA$WNN_L4, 
  Case = M0_RNA$Case, 
  Case_Type = M0_RNA$Case_Type
)
rownames(M0_WNN_L4_HC_ModuleScores) <- colnames(M0_RNA)

LIBD_Spatial_Seurat_WNN_L4_HC_ModuleScores <- data.frame(
  orig.ident=LIBD_Spatial_Seurat$orig.ident, 
  sample_id = LIBD_Spatial_Seurat$sample_id, 
  subject = LIBD_Spatial_Seurat$subject, 
  position = LIBD_Spatial_Seurat$position, 
  replicate = LIBD_Spatial_Seurat$replicate, 
  key = LIBD_Spatial_Seurat$key, 
  cell_count = LIBD_Spatial_Seurat$cell_count, 
  Original_Barcode  = LIBD_Spatial_Seurat$Original_Barcode, 
  layer_guess = LIBD_Spatial_Seurat$layer_guess, 
  layer_guess_reordered = LIBD_Spatial_Seurat$layer_guess_reordered, 
  layer_guess_reordered_short = LIBD_Spatial_Seurat$layer_guess_reordered_short
)
rownames(LIBD_Spatial_Seurat_WNN_L4_HC_ModuleScores) <- colnames(LIBD_Spatial_Seurat)


for (CellType in WNN_L4_HC_ModuleScores_Markers$CellType){
  
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
    
    for (nfts in c(50, 100, 200, 500, 1234)){
      if(nfts==1234){nfts <- length(WNN_L4_HC_ModuleScores_Markers$markers_SYMB[[which(WNN_L4_HC_ModuleScores_Markers$CellType==CellType)]]) }
      message(paste0("Number of features: ", nfts))
      stopifnot(all(colnames(M0_RNA)==M0_WNN_L4_HC_ModuleScores$CellId))
      name = paste0("WNN_L4_", CellType, "_ModuleScore_", nfts, "fts_")
      fts = WNN_L4_HC_ModuleScores_Markers$markers_SYMB[which(WNN_L4_HC_ModuleScores_Markers$CellType==CellType)] 
      fts_ensg = WNN_L4_HC_ModuleScores_Markers$markers_ENSG[which(WNN_L4_HC_ModuleScores_Markers$CellType==CellType)]  
      
      M0_WNN_L4_HC_ModuleScores[[paste0(name, "1")]] <- AddModuleScore(
        M0_RNA, 
        features=fts, 
        name=name, 
        ctrl=nfts , 
        assay = "RNA"
      )[[paste0(name, "1")]][,1]
      
      LIBD_Spatial_Seurat_WNN_L4_HC_ModuleScores[[paste0(name, "1")]] <- AddModuleScore(
        LIBD_Spatial_Seurat, 
        features=fts_ensg, 
        name=name, 
        ctrl=nfts , 
        assay = "SpatialLogCounts"
      )[[paste0(name, "1")]][,1] 
      
      
      rm(name, fts, fts_ensg)
    }
  }, error=function(e){
    
  })
}    

qsave(
  M0_WNN_L4_HC_ModuleScores, 
  paste0(
    "../Data/Annotations/Module_Scores/RNA/HC_Derived/", 
    "M0_WNN_L4_HC_ModuleScores", 
    ".qrds"
  ), 
  nthr=nthr
)  

qsave(
  LIBD_Spatial_Seurat_WNN_L4_HC_ModuleScores, 
  paste0(
    "../Data/Annotations/Module_Scores/RNA/HC_Derived/", 
    "LIBD_Spatial_Seurat_WNN_L4_HC_ModuleScores", 
    ".qrds"
  ), 
  nthr=nthr
)    

#rm(M0_WNN_L4_HC_ModuleScores, LIBD_Spatial_Seurat_WNN_L4_HC_ModuleScores)  



    ## 3.4 Extratelencephalic Cells Gene Score ---------------------------------


ETNC_Markers <- fread(
  paste0(
    "../Data/Input/", 
    "Markers_Extratelencephalic_Neurons", 
    ".txt"
  ), 
  header=FALSE
)
  

M0_ETNC_HC_ModuleScores <- data.frame(
  CellId=M0_RNA$CellId, 
  gex_id = M0_RNA$gex_barcode, 
  atac_barcode = M0_RNA$atac_barcode, 
  WNN_L1 = M0_RNA$WNN_L1, 
  WNN_L15 = M0_RNA$WNN_L15, 
  WNN_L2 = M0_RNA$WNN_L2, 
  WNN_L25 = M0_RNA$WNN_L25, 
  WNN_L3 = M0_RNA$WNN_L3, 
  WNN_L4 = M0_RNA$WNN_L4, 
  Case = M0_RNA$Case, 
  Case_Type = M0_RNA$Case_Type
)
rownames(M0_ETNC_HC_ModuleScores) <- colnames(M0_RNA)

LIBD_Spatial_Seurat_ETNC_HC_ModuleScores <- data.frame(
  orig.ident=LIBD_Spatial_Seurat$orig.ident, 
  sample_id = LIBD_Spatial_Seurat$sample_id, 
  subject = LIBD_Spatial_Seurat$subject, 
  position = LIBD_Spatial_Seurat$position, 
  replicate = LIBD_Spatial_Seurat$replicate, 
  key = LIBD_Spatial_Seurat$key, 
  cell_count = LIBD_Spatial_Seurat$cell_count, 
  Original_Barcode  = LIBD_Spatial_Seurat$Original_Barcode, 
  layer_guess = LIBD_Spatial_Seurat$layer_guess, 
  layer_guess_reordered = LIBD_Spatial_Seurat$layer_guess_reordered, 
  layer_guess_reordered_short = LIBD_Spatial_Seurat$layer_guess_reordered_short
)
rownames(LIBD_Spatial_Seurat_ETNC_HC_ModuleScores) <- colnames(LIBD_Spatial_Seurat)

for (nfts in c(50, 100, 200, 500, length(ETNC_Markers$V1))){
  
  message(paste0("Number of features: ", nfts))
  stopifnot(all(colnames(M0_RNA)==M0_ETNC_HC_ModuleScores$CellId))
  
  name = paste0("ETNC_ModuleScore_", nfts, "fts_")
  fts = ETNC_Markers$V1
  fts_ensg = features$ENSG[match(fts, features$Symbol)]  
  
  table(is.na(fts)) 
  table(is.na(fts_ensg))
  table(duplicated(fts))
  table(duplicated(fts_ensg)) 
  
  
  M0_ETNC_HC_ModuleScores[[paste0(name, "1")]] <- AddModuleScore(
    M0_RNA, 
    features=list(fts), 
    name=name, 
    ctrl=nfts , 
    assay = "RNA"
  )[[paste0(name, "1")]][,1]
  
  LIBD_Spatial_Seurat_ETNC_HC_ModuleScores[[paste0(name, "1")]] <- AddModuleScore(
    LIBD_Spatial_Seurat, 
    features=list(fts_ensg),
    name=name, 
    ctrl=nfts , 
    assay = "SpatialLogCounts"
  )[[paste0(name, "1")]][,1] 
  
  
  rm(name, fts, fts_ensg)
}


qsave(
  M0_ETNC_HC_ModuleScores, 
  paste0(
    "../Data/Annotations/Module_Scores/RNA/HC_Derived/", 
    "M0_ETNC_HC_ModuleScores", 
    ".qrds"
  ), 
  nthr=nthr
)  

qsave(
  LIBD_Spatial_Seurat_ETNC_HC_ModuleScores, 
  paste0(
    "../Data/Annotations/Module_Scores/RNA/HC_Derived/", 
    "LIBD_Spatial_Seurat_ETNC_HC_ModuleScores", 
    ".qrds"
  ), 
  nthr=nthr
)    

rm(M0_WNN_L4_HC_ModuleScores, LIBD_Spatial_Seurat_WNN_L4_HC_ModuleScores)  



  ### 4.0 Generate CellType-specific signatures with All samples ---------------



    ## 4.1 CellTypes Gene Modules WNN_L15 --------------------------------------

WNN_L15_Markers <- qread(
  paste0(
    "../Data/Annotations/CellType_Markers/", 
    "WNN_L15_Markers_HC_ALSFTD_Wilcox", 
    ".qrds"
  ), 
  nthr=nthr
)

WNN_L15_HC_ALSFTD_ModuleScores_Markers <- data.frame(
  CellTypeLevel = "WNN_L15",
  CellType = names(WNN_L15_Markers),
  n = NA, 
  fdr = NA, 
  log2fc = NA
)
WNN_L15_HC_ALSFTD_ModuleScores_Markers$markers_ENSG = list(rep(NA, length(WNN_L15_Markers))) 
WNN_L15_HC_ALSFTD_ModuleScores_Markers$markers_SYMB = list(rep(NA, length(WNN_L15_Markers)))

imap(
  WNN_L15_Markers, 
  .f=function(x, y){
    
    log = NULL
    fdr = NULL
    for (l in log2(c(5,4,3,2))){
      for (f in c(1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 5e-2)){
        if(sum(x$p_val_adj < f & x$avg_log2FC > l) >= 50){
          fdr <- f 
          break
        }
      }
      if(!is.null(fdr)){
        log <- l 
        break
      }
    } 

    if(!is.null(log) & !is.null(fdr)){
      
      ind <- which(x$avg_log2FC > log & x$p_val_adj < fdr)
      
      symb <- x$SYMBOL[ind] 
      ensg <- features$ENSG[match(symb, features$Symbol)]
      
      symb <- symb[!is.na(ensg)] 
      ensg <- ensg[!is.na(ensg)] 
      
      stopifnot(!any(is.na(symb))) 
      stopifnot(length(ensg)==length(symb))
      i = which(WNN_L15_HC_ALSFTD_ModuleScores_Markers$CellType==y)
      WNN_L15_HC_ALSFTD_ModuleScores_Markers$n[i] <<- length(ensg)
      WNN_L15_HC_ALSFTD_ModuleScores_Markers$fdr[i] <<- fdr 
      WNN_L15_HC_ALSFTD_ModuleScores_Markers$log2fc[i] <<- log 
      WNN_L15_HC_ALSFTD_ModuleScores_Markers$markers_ENSG[[i]] <<- ensg
      WNN_L15_HC_ALSFTD_ModuleScores_Markers$markers_SYMB[[i]] <<- symb      
    }
  }
)

qsave(
  WNN_L15_HC_ALSFTD_ModuleScores_Markers, 
  paste0(
    "../Data/Annotations/Module_Scores/RNA/HC_ALSFTD_Derived/", 
    "WNN_L15_HC_ALSFTD_ModuleScores_Markers", 
    ".qrds"
  ), 
  nthr = nthr
)

rm(WNN_L15_Markers)

M0_WNN_L15_HC_ALSFTD_ModuleScores <- data.frame(
  CellId=M0_RNA$CellId, 
  gex_id = M0_RNA$gex_barcode, 
  atac_barcode = M0_RNA$atac_barcode, 
  WNN_L1 = M0_RNA$WNN_L1, 
  WNN_L15 = M0_RNA$WNN_L15, 
  WNN_L2 = M0_RNA$WNN_L2, 
  WNN_L25 = M0_RNA$WNN_L25, 
  WNN_L3 = M0_RNA$WNN_L3, 
  WNN_L4 = M0_RNA$WNN_L4, 
  Case = M0_RNA$Case, 
  Case_Type = M0_RNA$Case_Type
)
rownames(M0_WNN_L15_HC_ALSFTD_ModuleScores) <- colnames(M0_RNA)

LIBD_Spatial_Seurat_WNN_L15_HC_ALSFTD_ModuleScores <- data.frame(
  orig.ident=LIBD_Spatial_Seurat$orig.ident, 
  sample_id = LIBD_Spatial_Seurat$sample_id, 
  subject = LIBD_Spatial_Seurat$subject, 
  position = LIBD_Spatial_Seurat$position, 
  replicate = LIBD_Spatial_Seurat$replicate, 
  key = LIBD_Spatial_Seurat$key, 
  cell_count = LIBD_Spatial_Seurat$cell_count, 
  Original_Barcode  = LIBD_Spatial_Seurat$Original_Barcode, 
  layer_guess = LIBD_Spatial_Seurat$layer_guess, 
  layer_guess_reordered = LIBD_Spatial_Seurat$layer_guess_reordered, 
  layer_guess_reordered_short = LIBD_Spatial_Seurat$layer_guess_reordered_short
)
rownames(LIBD_Spatial_Seurat_WNN_L15_HC_ALSFTD_ModuleScores) <- colnames(LIBD_Spatial_Seurat)


for (CellType in WNN_L15_HC_ALSFTD_ModuleScores_Markers$CellType){
  
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
    
    for (nfts in c(50, 100, 200, 500, 1234)){
      if(nfts==1234){nfts <- length(WNN_L15_HC_ALSFTD_ModuleScores_Markers$markers_SYMB[[which(WNN_L15_HC_ALSFTD_ModuleScores_Markers$CellType==CellType)]]) }
      message(paste0("Number of features: ", nfts))
      stopifnot(all(colnames(M0_RNA)==M0_WNN_L15_HC_ALSFTD_ModuleScores$CellId))
      name = paste0("WNN_L15_", CellType, "_ModuleScore_", nfts, "fts_")
      fts = WNN_L15_HC_ALSFTD_ModuleScores_Markers$markers_SYMB[which(WNN_L15_HC_ALSFTD_ModuleScores_Markers$CellType==CellType)] 
      fts_ensg = WNN_L15_HC_ALSFTD_ModuleScores_Markers$markers_ENSG[which(WNN_L15_HC_ALSFTD_ModuleScores_Markers$CellType==CellType)]  
      
      M0_WNN_L15_HC_ALSFTD_ModuleScores[[paste0(name, "1")]] <- AddModuleScore(
        M0_RNA, 
        features=fts, 
        name=name, 
        ctrl=nfts , 
        assay = "RNA"
      )[[paste0(name, "1")]][,1]
      
      LIBD_Spatial_Seurat_WNN_L15_HC_ALSFTD_ModuleScores[[paste0(name, "1")]] <- AddModuleScore(
        LIBD_Spatial_Seurat, 
        features=fts_ensg, 
        name=name, 
        ctrl=nfts , 
        assay = "SpatialLogCounts"
      )[[paste0(name, "1")]][,1] 
      
      
      rm(name, fts, fts_ensg)
    }
  }, error=function(e){
    
  })
}    

qsave(
  M0_WNN_L15_HC_ALSFTD_ModuleScores, 
  paste0(
    "../Data/Annotations/Module_Scores/RNA/HC_ALSFTD_Derived/", 
    "M0_WNN_L15_HC_ALSFTD_ModuleScores", 
    ".qrds"
  ), 
  nthr=nthr
)  

qsave(
  LIBD_Spatial_Seurat_WNN_L15_HC_ALSFTD_ModuleScores, 
  paste0(
    "../Data/Annotations/Module_Scores/RNA/HC_ALSFTD_Derived/", 
    "LIBD_Spatial_Seurat_WNN_L15_HC_ALSFTD_ModuleScores", 
    ".qrds"
  ), 
  nthr=nthr
)    

#rm(M0_WNN_L15_HC_ALSFTD_ModuleScores, LIBD_Spatial_Seurat_WNN_L15_HC_ALSFTD_ModuleScores)  

    ## 4.2 CellTypes Gene Modules WNN_L25 --------------------------------------

WNN_L25_Markers <- qread(
  paste0(
    "../Data/Annotations/CellType_Markers/", 
    "WNN_L25_Markers_HC_ALSFTD_Wilcox", 
    ".qrds"
  ), 
  nthr=nthr
)

WNN_L25_HC_ALSFTD_ModuleScores_Markers <- data.frame(
  CellTypeLevel = "WNN_L25",
  CellType = names(WNN_L25_Markers),
  n = NA, 
  fdr = NA, 
  log2fc = NA
)
WNN_L25_HC_ALSFTD_ModuleScores_Markers$markers_ENSG = list(rep(NA, length(WNN_L25_Markers))) 
WNN_L25_HC_ALSFTD_ModuleScores_Markers$markers_SYMB = list(rep(NA, length(WNN_L25_Markers)))

imap(
  WNN_L25_Markers, 
  .f=function(x, y){
    
    log = NULL
    fdr = NULL
    for (l in log2(c(5,4,3,2))){
      for (f in c(1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 5e-2)){
        if(sum(x$p_val_adj < f & x$avg_log2FC > l) >= 50){
          fdr <- f 
          break
        }
      }
      if(!is.null(fdr)){
        log <- l 
        break
      }
    } 
    
    if(!is.null(log) & !is.null(fdr)){
      
      ind <- which(x$avg_log2FC > log & x$p_val_adj < fdr)
      
      symb <- x$SYMBOL[ind] 
      ensg <- features$ENSG[match(symb, features$Symbol)]
      
      symb <- symb[!is.na(ensg)] 
      ensg <- ensg[!is.na(ensg)] 
      
      stopifnot(!any(is.na(symb))) 
      stopifnot(length(ensg)==length(symb))
      i = which(WNN_L25_HC_ALSFTD_ModuleScores_Markers$CellType==y)
      WNN_L25_HC_ALSFTD_ModuleScores_Markers$n[i] <<- length(ensg)
      WNN_L25_HC_ALSFTD_ModuleScores_Markers$fdr[i] <<- fdr 
      WNN_L25_HC_ALSFTD_ModuleScores_Markers$log2fc[i] <<- log 
      WNN_L25_HC_ALSFTD_ModuleScores_Markers$markers_ENSG[[i]] <<- ensg
      WNN_L25_HC_ALSFTD_ModuleScores_Markers$markers_SYMB[[i]] <<- symb      
    }
  }
)

qsave(
  WNN_L25_HC_ALSFTD_ModuleScores_Markers, 
  paste0(
    "../Data/Annotations/Module_Scores/RNA/HC_ALSFTD_Derived/", 
    "WNN_L25_HC_ALSFTD_ModuleScores_Markers", 
    ".qrds"
  ), 
  nthr = nthr
)

rm(WNN_L25_Markers)

M0_WNN_L25_HC_ALSFTD_ModuleScores <- data.frame(
  CellId=M0_RNA$CellId, 
  gex_id = M0_RNA$gex_barcode, 
  atac_barcode = M0_RNA$atac_barcode, 
  WNN_L1 = M0_RNA$WNN_L1, 
  WNN_L15 = M0_RNA$WNN_L15, 
  WNN_L2 = M0_RNA$WNN_L2, 
  WNN_L25 = M0_RNA$WNN_L25, 
  WNN_L3 = M0_RNA$WNN_L3, 
  WNN_L4 = M0_RNA$WNN_L4, 
  Case = M0_RNA$Case, 
  Case_Type = M0_RNA$Case_Type
)
rownames(M0_WNN_L25_HC_ALSFTD_ModuleScores) <- colnames(M0_RNA)

LIBD_Spatial_Seurat_WNN_L25_HC_ALSFTD_ModuleScores <- data.frame(
  orig.ident=LIBD_Spatial_Seurat$orig.ident, 
  sample_id = LIBD_Spatial_Seurat$sample_id, 
  subject = LIBD_Spatial_Seurat$subject, 
  position = LIBD_Spatial_Seurat$position, 
  replicate = LIBD_Spatial_Seurat$replicate, 
  key = LIBD_Spatial_Seurat$key, 
  cell_count = LIBD_Spatial_Seurat$cell_count, 
  Original_Barcode  = LIBD_Spatial_Seurat$Original_Barcode, 
  layer_guess = LIBD_Spatial_Seurat$layer_guess, 
  layer_guess_reordered = LIBD_Spatial_Seurat$layer_guess_reordered, 
  layer_guess_reordered_short = LIBD_Spatial_Seurat$layer_guess_reordered_short
)
rownames(LIBD_Spatial_Seurat_WNN_L25_HC_ALSFTD_ModuleScores) <- colnames(LIBD_Spatial_Seurat)


for (CellType in WNN_L25_HC_ALSFTD_ModuleScores_Markers$CellType){
  
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
    
    for (nfts in c(50, 100, 200, 500, 1234)){
      if(nfts==1234){nfts <- length(WNN_L25_HC_ALSFTD_ModuleScores_Markers$markers_SYMB[[which(WNN_L25_HC_ALSFTD_ModuleScores_Markers$CellType==CellType)]]) }
      message(paste0("Number of features: ", nfts))
      stopifnot(all(colnames(M0_RNA)==M0_WNN_L25_HC_ALSFTD_ModuleScores$CellId))
      name = paste0("WNN_L25_", CellType, "_ModuleScore_", nfts, "fts_")
      fts = WNN_L25_HC_ALSFTD_ModuleScores_Markers$markers_SYMB[which(WNN_L25_HC_ALSFTD_ModuleScores_Markers$CellType==CellType)] 
      fts_ensg = WNN_L25_HC_ALSFTD_ModuleScores_Markers$markers_ENSG[which(WNN_L25_HC_ALSFTD_ModuleScores_Markers$CellType==CellType)]  
      
      M0_WNN_L25_HC_ALSFTD_ModuleScores[[paste0(name, "1")]] <- AddModuleScore(
        M0_RNA, 
        features=fts, 
        name=name, 
        ctrl=nfts , 
        assay = "RNA"
      )[[paste0(name, "1")]][,1]
      
      LIBD_Spatial_Seurat_WNN_L25_HC_ALSFTD_ModuleScores[[paste0(name, "1")]] <- AddModuleScore(
        LIBD_Spatial_Seurat, 
        features=fts_ensg, 
        name=name, 
        ctrl=nfts , 
        assay = "SpatialLogCounts"
      )[[paste0(name, "1")]][,1] 
      
      
      rm(name, fts, fts_ensg)
    }
  }, error=function(e){
    
  })
}    

qsave(
  M0_WNN_L25_HC_ALSFTD_ModuleScores, 
  paste0(
    "../Data/Annotations/Module_Scores/RNA/HC_ALSFTD_Derived/", 
    "M0_WNN_L25_HC_ALSFTD_ModuleScores", 
    ".qrds"
  ), 
  nthr=nthr
)  

qsave(
  LIBD_Spatial_Seurat_WNN_L25_HC_ALSFTD_ModuleScores, 
  paste0(
    "../Data/Annotations/Module_Scores/RNA/HC_ALSFTD_Derived/", 
    "LIBD_Spatial_Seurat_WNN_L25_HC_ALSFTD_ModuleScores", 
    ".qrds"
  ), 
  nthr=nthr
)    

#rm(M0_WNN_L25_HC_ALSFTD_ModuleScores, LIBD_Spatial_Seurat_WNN_L25_HC_ALSFTD_ModuleScores)  

    ## 4.3 CellTypes Gene Modules WNN_L4 ---------------------------------------

WNN_L4_Markers <- qread(
  paste0(
    "../Data/Annotations/CellType_Markers/", 
    "WNN_L4_Markers_HC_ALSFTD_Wilcox", 
    ".qrds"
  ), 
  nthr=nthr
)

WNN_L4_HC_ALSFTD_ModuleScores_Markers <- data.frame(
  CellTypeLevel = "WNN_L4",
  CellType = names(WNN_L4_Markers),
  n = NA, 
  fdr = NA, 
  log2fc = NA
)
WNN_L4_HC_ALSFTD_ModuleScores_Markers$markers_ENSG = list(rep(NA, length(WNN_L4_Markers))) 
WNN_L4_HC_ALSFTD_ModuleScores_Markers$markers_SYMB = list(rep(NA, length(WNN_L4_Markers)))

imap(
  WNN_L4_Markers, 
  .f=function(x, y){
    
    log = NULL
    fdr = NULL
    for (l in log2(c(5,4,3,2))){
      for (f in c(1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 5e-2)){
        if(sum(x$p_val_adj < f & x$avg_log2FC > l) >= 50){
          fdr <- f 
          break
        }
      }
      if(!is.null(fdr)){
        log <- l 
        break
      }
    } 
    
    if(!is.null(log) & !is.null(fdr)){
      
      ind <- which(x$avg_log2FC > log & x$p_val_adj < fdr)
      
      symb <- x$SYMBOL[ind] 
      ensg <- features$ENSG[match(symb, features$Symbol)]
      
      symb <- symb[!is.na(ensg)] 
      ensg <- ensg[!is.na(ensg)] 
      
      stopifnot(!any(is.na(symb))) 
      stopifnot(length(ensg)==length(symb))
      i = which(WNN_L4_HC_ALSFTD_ModuleScores_Markers$CellType==y)
      WNN_L4_HC_ALSFTD_ModuleScores_Markers$n[i] <<- length(ensg)
      WNN_L4_HC_ALSFTD_ModuleScores_Markers$fdr[i] <<- fdr 
      WNN_L4_HC_ALSFTD_ModuleScores_Markers$log2fc[i] <<- log 
      WNN_L4_HC_ALSFTD_ModuleScores_Markers$markers_ENSG[[i]] <<- ensg
      WNN_L4_HC_ALSFTD_ModuleScores_Markers$markers_SYMB[[i]] <<- symb      
    }
  }
)

qsave(
  WNN_L4_HC_ALSFTD_ModuleScores_Markers, 
  paste0(
    "../Data/Annotations/Module_Scores/RNA/HC_ALSFTD_Derived/", 
    "WNN_L4_HC_ALSFTD_ModuleScores_Markers", 
    ".qrds"
  ), 
  nthr = nthr
)

rm(WNN_L4_Markers)

M0_WNN_L4_HC_ALSFTD_ModuleScores <- data.frame(
  CellId=M0_RNA$CellId, 
  gex_id = M0_RNA$gex_barcode, 
  atac_barcode = M0_RNA$atac_barcode, 
  WNN_L1 = M0_RNA$WNN_L1, 
  WNN_L15 = M0_RNA$WNN_L15, 
  WNN_L2 = M0_RNA$WNN_L2, 
  WNN_L25 = M0_RNA$WNN_L25, 
  WNN_L3 = M0_RNA$WNN_L3, 
  WNN_L4 = M0_RNA$WNN_L4, 
  Case = M0_RNA$Case, 
  Case_Type = M0_RNA$Case_Type
)
rownames(M0_WNN_L4_HC_ALSFTD_ModuleScores) <- colnames(M0_RNA)

LIBD_Spatial_Seurat_WNN_L4_HC_ALSFTD_ModuleScores <- data.frame(
  orig.ident=LIBD_Spatial_Seurat$orig.ident, 
  sample_id = LIBD_Spatial_Seurat$sample_id, 
  subject = LIBD_Spatial_Seurat$subject, 
  position = LIBD_Spatial_Seurat$position, 
  replicate = LIBD_Spatial_Seurat$replicate, 
  key = LIBD_Spatial_Seurat$key, 
  cell_count = LIBD_Spatial_Seurat$cell_count, 
  Original_Barcode  = LIBD_Spatial_Seurat$Original_Barcode, 
  layer_guess = LIBD_Spatial_Seurat$layer_guess, 
  layer_guess_reordered = LIBD_Spatial_Seurat$layer_guess_reordered, 
  layer_guess_reordered_short = LIBD_Spatial_Seurat$layer_guess_reordered_short
)
rownames(LIBD_Spatial_Seurat_WNN_L4_HC_ALSFTD_ModuleScores) <- colnames(LIBD_Spatial_Seurat)


for (CellType in WNN_L4_HC_ALSFTD_ModuleScores_Markers$CellType){
  
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
    
    for (nfts in c(50, 100, 200, 500, 1234)){
      if(nfts==1234){nfts <- length(WNN_L4_HC_ALSFTD_ModuleScores_Markers$markers_SYMB[[which(WNN_L4_HC_ALSFTD_ModuleScores_Markers$CellType==CellType)]]) }
      message(paste0("Number of features: ", nfts))
      stopifnot(all(colnames(M0_RNA)==M0_WNN_L4_HC_ALSFTD_ModuleScores$CellId))
      name = paste0("WNN_L4_", CellType, "_ModuleScore_", nfts, "fts_")
      fts = WNN_L4_HC_ALSFTD_ModuleScores_Markers$markers_SYMB[which(WNN_L4_HC_ALSFTD_ModuleScores_Markers$CellType==CellType)] 
      fts_ensg = WNN_L4_HC_ALSFTD_ModuleScores_Markers$markers_ENSG[which(WNN_L4_HC_ALSFTD_ModuleScores_Markers$CellType==CellType)]  
      
      M0_WNN_L4_HC_ALSFTD_ModuleScores[[paste0(name, "1")]] <- AddModuleScore(
        M0_RNA, 
        features=fts, 
        name=name, 
        ctrl=nfts , 
        assay = "RNA"
      )[[paste0(name, "1")]][,1]
      
      LIBD_Spatial_Seurat_WNN_L4_HC_ALSFTD_ModuleScores[[paste0(name, "1")]] <- AddModuleScore(
        LIBD_Spatial_Seurat, 
        features=fts_ensg, 
        name=name, 
        ctrl=nfts , 
        assay = "SpatialLogCounts"
      )[[paste0(name, "1")]][,1] 
      
      
      rm(name, fts, fts_ensg)
    }
  }, error=function(e){
    
  })
}    

qsave(
  M0_WNN_L4_HC_ALSFTD_ModuleScores, 
  paste0(
    "../Data/Annotations/Module_Scores/RNA/HC_ALSFTD_Derived/", 
    "M0_WNN_L4_HC_ALSFTD_ModuleScores", 
    ".qrds"
  ), 
  nthr=nthr
)  

qsave(
  LIBD_Spatial_Seurat_WNN_L4_HC_ALSFTD_ModuleScores, 
  paste0(
    "../Data/Annotations/Module_Scores/RNA/HC_ALSFTD_Derived/", 
    "LIBD_Spatial_Seurat_WNN_L4_HC_ALSFTD_ModuleScores", 
    ".qrds"
  ), 
  nthr=nthr
)    

#rm(M0_WNN_L4_HC_ALSFTD_ModuleScores, LIBD_Spatial_Seurat_WNN_L4_HC_ALSFTD_ModuleScores)  



