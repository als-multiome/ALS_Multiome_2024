# A10_QC_Filter_Cells_Merged_Seurat_Objects.R 

  ### 0.0 Load libraries ------------------------------------------------------- 

library(qs2)
library(tidyverse)
library(Seurat)
  

  
  
  ### 1.0 Load data ------------------------------------------------------------
  
Samples <- qs_read(
  "../Data/FANS/Samples/Samples.qs2"
)

FANS  <- qs_read(
  "../Data/FANS/SeuratObjects/FANS_Merged.qs2"
)




  ### 2.0 Define AUX functions -------------------------------------------------

filter_by_MADs <- function(Var, nMads){
  
  MAD <- mad(Var, center = median(Var, na.rm = TRUE), constant = 1, na.rm = TRUE)
  return(
    Var < median(Var, na.rm = TRUE) - nMads*MAD | 
      Var > median(Var, na.rm = TRUE) + nMads*MAD
  ) 
}




  ### 3.0 QC Filter Cells ------------------------------------------------------

FANS$Filter <- !(
  filter_by_MADs(FANS$nCount_RNA_log1p, nMads = 3) | 
    filter_by_MADs(FANS$Features_by_Counts, nMads = 3) | 
    filter_by_MADs(FANS$Pct_Counts_In_Top_20_Genes, nMads = 5) | 
    FANS$percent.mt > 1 | 
    FANS$percent.ribo > 1 
)

table(FANS$Filter, useNA = "always")

FANS <- subset(
  FANS, 
  subset = Filter
)
FANS  

FANS@meta.data <- FANS@meta.data %>% 
  select(
    -c(
      CellBender_IsCell, 
      MAD_Filter_Flag, 
      scDblFinder.class, 
      Souporcell_Status, 
      Filter
    )
  )




  ### 4.0 Assay Uncorrected filter ---------------------------------------------

FANS@assays$RNA_Uncorrected <- NULL 
FANS




  ### 4.0 Save data ------------------------------------------------------------

qs_save(
  FANS, 
  "../Data/FANS/SeuratObjects/FANS_Merged_QCFiltered.qs2"
)






