



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

all(rownames(LIBD_Spatial_Spe)==rownames(LIBD_Spatial_Seurat))
all(colnames(LIBD_Spatial_Spe)==str_split(colnames(LIBD_Spatial_Seurat), "_", simplify=TRUE)[,2])



M0_RNA_HC <- subset(M0_RNA, subset=Case=="HC")

dim(M0_RNA_HC) 
table(M0_RNA_HC$Case) 




  ### 2.0 Generate a dictionary for ENSG IDs and HUGO Gene Symbols -------------

genes <- genes(EnsDb.Hsapiens.v86) 
features <- data.frame(
  symbol = rownames(M0_RNA@assays$RNA)
)

features$Gene_Name <- genes$gene_name[match(features$symbol, genes$symbol)]
features$ENSGl <- genes$gene_id[match(features$symbol, genes$symbol)]

all(features$Gene_Name==features$Symbol, na.rm=TRUE)
features$Gene_biotype = genes$gene_biotype[match(features$ENSG, genes$gene_id)]
table(features$Gene_biotype)


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
      x$ENSG <- features$ENSG[match(x$SYMBOL, features$symbol)] 
      x
    }, error=function(e){}
    )
  }
)



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
      x$ENSG <- features$ENSG[match(x$SYMBOL, features$symbol)] 
      x
    }, error=function(e){}
    )
  }
)



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



    ## 3.4 Save Data -----------------------------------------------------------

qsave(
  WNN_L15_Markers, 
  paste0(
    "../Data/Annotations/CellType_Markers/", 
    "WNN_L15_Markers_HC_Wilcox", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  WNN_L25_Markers, 
  paste0(
    "../Data/Annotations/CellType_Markers/", 
    "WNN_L25_Markers_HC_Wilcox", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  WNN_L4_Markers, 
  paste0(
    "../Data/Annotations/CellType_Markers/", 
    "WNN_L4_Markers_HC_Wilcox", 
    ".qrds"
  ), 
  nthr=nthr
)

rm(WNN_L15_Markers, WNN_L25_Markers, WNN_L4_Markers)




  ### 4.0 Generate CellType-specific signatures with the M0 Full dataset -------



    ## 4.1 CellTypes Gene Modules WNN_L15 --------------------------------------
  
  WNN_L15_Markers <- list()  
  Idents(M0_RNA) <- M0_RNA$WNN_L15
  
  for (CellType in unique(M0_RNA$WNN_L15)){
    tryCatch({
      WNN_L15_Markers[[CellType]] <- FindMarkers(
        M0_RNA, 
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
        x$ENSG <- features$ENSG[match(x$SYMBOL, features$symbol)] 
        x
      }, error=function(e){}
      )
    }
  )
  

  
    ## 4.2 CellTypes Gene Modules WNN_L25 --------------------------------------
  
  WNN_L25_Markers <- list()  
  Idents(M0_RNA) <- M0_RNA$WNN_L25
  
  for (CellType in unique(M0_RNA$WNN_L25)){
    tryCatch({
      WNN_L25_Markers[[CellType]] <- FindMarkers(
        M0_RNA, 
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
        x$ENSG <- features$ENSG[match(x$SYMBOL, features$symbol)] 
        x
      }, error=function(e){}
      )
    }
  )

  
  
    ## 4.3 CellTypes Gene Modules WNN_L4 --------------------------------------
  
  WNN_L4_Markers <- list()  
  Idents(M0_RNA) <- M0_RNA$WNN_L4
  
  for (CellType in unique(M0_RNA$WNN_L4)){
    tryCatch({
      WNN_L4_Markers[[CellType]] <- FindMarkers(
        M0_RNA, 
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
        x$ENSG <- features$ENSG[match(x$SYMBOL, features$symbol)] 
        x
      }, error=function(e){})
    }
  )


    
    ## 4.4 Save Data -----------------------------------------------------------
  
  qsave(
    WNN_L15_Markers, 
    paste0(
      "../Data/Annotations/CellType_Markers/", 
      "WNN_L15_Markers_HC_ALSFTD_Wilcox", 
      ".qrds"
    ), 
    nthr=nthr
  )
  
  qsave(
    WNN_L25_Markers, 
    paste0(
      "../Data/Annotations/CellType_Markers/", 
      "WNN_L25_Markers_HC_ALSFTD_Wilcox", 
      ".qrds"
    ), 
    nthr=nthr
  )
  
  qsave(
    WNN_L4_Markers, 
    paste0(
      "../Data/Annotations/CellType_Markers/", 
      "WNN_L4_Markers_HC_ALSFTD_Wilcox", 
      ".qrds"
    ), 
    nthr=nthr
  )
  
  
  
  