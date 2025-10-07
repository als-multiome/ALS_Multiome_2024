# A5_QC_Filter_Cells_SingleObjects.R 




  ### 0.0 Load libraries ------------------------------------------------------- 

library(Seurat)
library(data.table)
library(qs2)
library(tidyverse)

  
  
  
  ### 1.0 Load data ------------------------------------------------------------
  
FANS_Seurats_Raw <- qs_read(
  "../Data/FANS/SeuratObjects/FANS_Single_SeuratObjects_Raw.qs2"
)

Samples <- qs_read(
  "../Data/FANS/Samples/Samples.qs2"
)




  ### 2.0 Define AUX functions -------------------------------------------------

listCols<-function(m){
  return(
    lapply(seq_len(ncol(m)), function(j) {
      idx_start <- m@p[j] + 1
      idx_end <- m@p[j + 1]
      if (idx_start <= idx_end) {
        m@x[idx_start:idx_end]
      } else {
        numeric(0)
      }
    })
  )
}

percentage_in_top_Pct <- function(Var, Pct=NULL, N=NULL){
  return(
    sum(
      sort(Var, decreasing = TRUE)[1:if(is.null(N)){ceiling(Pct*length(Var))}else{N}], na.rm = TRUE
    )/sum(
      Var, na.rm = TRUE
    )
  )
}

filter_by_MADs <- function(Var, nMads){
  
  MAD <- mad(Var, center = median(Var, na.rm = TRUE), constant = 1, na.rm = TRUE)
  return(
    Var < median(Var, na.rm = TRUE) - nMads*MAD | 
      Var > median(Var, na.rm = TRUE) + nMads*MAD
  ) 
}




  ### 3.0 Calculate QC metrics ------------------------------------------------- 

FANS_Seurats <- lapply(
  FANS_Seurats_Raw, 
  function(x){
    x$percent.mt <- PercentageFeatureSet(x, assay = "RNA", pattern = "^MT-")
    x$percent.ribo <- PercentageFeatureSet(x, assay = "RNA", pattern = "^RP[SL]")
    x$percent.soup <- 1-(colSums(x@assays$RNA$counts, na.rm = TRUE)/colSums(x@assays$RNA_Uncorrected$counts, na.rm = TRUE))
    x$percent.mt_withSoup <- PercentageFeatureSet(x, assay = "RNA_Uncorrected", pattern = "^MT-")
    x$percent.ribo_withSoup <- PercentageFeatureSet(x, assay = "RNA_Uncorrected", pattern = "^RP[SL]")
    x$Features_by_Counts <- x$nFeature_RNA/x$nCount_RNA 
    x$nCount_RNA_log1p <- log(x$nCount_RNA + 1)
    x$Pct_Counts_In_Top_20_Genes <- vapply(
      listCols(x@assays$RNA$counts), 
      percentage_in_top_Pct, 
      N=20, 
      FUN.VALUE=0.0
    ) 
    
    x$MAD_Filter_Flag <- filter_by_MADs(x$nCount_RNA_log1p, nMads = 5) | 
      x$percent.mt_withSoup >= 5 | 
      x$percent.ribo >= 5 
    return(x)
  }
)




  ### 4.0 Filter cells by QC --------------------------------------------------- 

FANS_Seurats <- lapply(
  FANS_Seurats, 
  function(x){
    return(
      subset(x, subset = CellBender_IsCell)
    )
  }
)
 
FANS_Seurats <- lapply(
  FANS_Seurats, 
  function(x){
    return(
      subset(x, subset = MAD_Filter_Flag, invert = TRUE)
    )
  }
)

lapply(
  FANS_Seurats_Raw, 
  dim
) %>% 
  data.frame() %>% 
  t() %>% 
  data.frame() %>% 
  rownames_to_column(var = "Sample") -> tmp 
  
Samples$nCells_CR <- tmp$X2[
  match(
    Samples$ID, 
    tmp$Sample
  )
]
rm(tmp)

lapply(
  FANS_Seurats, 
  dim
) %>% 
  data.frame() %>% 
  t() %>% 
  data.frame() %>% 
  rownames_to_column(var = "Sample") -> tmp 

Samples$nCells_QC_filtered <- tmp$X2[
  match(
    Samples$ID, 
    tmp$Sample
  )
]





  ### 5.0 Save data ------------------------------------------------------------ 

qs_save(
  FANS_Seurats, 
  "../Data/FANS/SeuratObjects/FANS_Single_SeuratObjects_Filtered.qs2"
)

qs_save(
  Samples, 
  "../Data/FANS/Samples/Samples.qs2"
)

