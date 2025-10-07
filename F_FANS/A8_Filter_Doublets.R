# A8_Filter_Doublets.R 




  ### 0.0 Load libraries -------------------------------------------------------

library(qs2)
library(Seurat)
library(tidyverse) 
library(data.table)




  ### 1.0 Load data ------------------------------------------------------------

Samples <- qs_read(
  "../Data/FANS/Samples/Samples.qs2"
)  

FANS <- qs_read(
  "../Data/FANS/SeuratObjects/FANS_Single_SeuratObjects_Filtered_scDblFinder_Marked_Genotypes_Assigned.qs2"
)




  ### 2.0 Add data to Samples data ---------------------------------------------

nDoublets_scDblFinder <- sapply(
  FANS, 
  function(x){
    return(
      sum(
        x$scDblFinder.class == "doublet"
      )
    )
  }
)

PctDoublets_scDblFinder <- sapply(
  FANS, 
  function(x){
    return(
      sum(
        x$scDblFinder.class == "doublet"
      )*100/nrow(x@meta.data)
    )
  }
)

nDoublets_Souporcell <- sapply(
  FANS, 
  function(x){
    return(
      sum(
        x$Souporcell_Status == "doublet"
      )
    )
  }
)

PctDoublets_Souporcell <- sapply(
  FANS, 
  function(x){
    return(
      sum(
        x$Souporcell_Status == "doublet"
      )*100/nrow(x@meta.data)
    )
  }
)

Samples$nDoublets_scDblFinder2 <- nDoublets_scDblFinder[match(Samples$ID, names(nDoublets_scDblFinder))]
Samples$PctDoublets_scDblFinder2 <- PctDoublets_scDblFinder[match(Samples$ID, names(PctDoublets_scDblFinder))]
all(Samples$PctDoublets_scDblFinder==Samples$PctDoublets_scDblFinder2)
Samples <- Samples %>% 
  select(-PctDoublets_scDblFinder2)

Samples$nDoublets_Souporcell <- nDoublets_Souporcell[match(Samples$ID, names(nDoublets_Souporcell))]
Samples$PctDoublets_Souporcell <- PctDoublets_Souporcell[match(Samples$ID, names(PctDoublets_Souporcell))]




  ### 3.0 Filter doublets ------------------------------------------------------

FANS_Seurats_Singlets <- lapply(
  FANS, 
  function(Seurat){
    return(
      subset(Seurat, subset = scDblFinder.class == "singlet")
    )
  }
)
sapply(FANS, function(x){dim(x)[2]}) |> sum()
sapply(FANS_Seurats_Singlets, function(x){dim(x)[2]}) |> sum()

FANS_Seurats_Singlets <- lapply(
  FANS, 
  function(Seurat){
    Seurat$Souporcell_Status[is.na(Seurat$Souporcell_Status)] <- "NA"
    return(
      subset(Seurat, subset = Souporcell_Status %in% c("NA", "singlet"))
    )
  }
)

sapply(FANS, function(x){dim(x)[2]}) |> sum()
sapply(FANS_Seurats_Singlets, function(x){dim(x)[2]}) |> sum()

lapply(
  FANS_Seurats_Singlets, 
  dim
) %>% 
  data.frame() %>% 
  t() %>% 
  data.frame() %>% 
  rownames_to_column(var = "Sample") -> tmp 


Samples$nCells_Filtered <- tmp$X2[
  match(
    Samples$ID, 
    tmp$Sample
  )
]




  ### 4.0 Save data ------------------------------------------------------------

qs_save(
  FANS_Seurats_Singlets, 
  "../Data/FANS/SeuratObjects/FANS_Single_SeuratObjects_Filtered_Singlets.qs2"
) 

qs_save(
  Samples, 
  "../Data/FANS/Samples/Samples.qs2"
) 




