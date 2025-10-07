# A19_Export_Barcodes.R




  ### 0.0 Load libraries ------------------------------------------------------- 

library(qs2)
library(tidyverse)
library(Seurat)

  


  ### 1.0 Load data ------------------------------------------------------------

FANS <- qs_read(
    "../Data/FANS/SeuratObjects/FANS.Unsubsetted.qs2"
)

Samples <- qs_read(
  "../Data/FANS/Samples/Samples.qs2"
)




  ### 2.0 Extract Barcodes by Sample -------------------------------------------

table(FANS$ID, FANS$Sample_donor, useNA="always")
table(FANS$ID, FANS$TDP43, useNA="always")
table(FANS$ID_WNN_L4_Predicted, useNA="always")



    ## 2.1 All Cells -----------------------------------------------------------

Barcodes_All <- lapply(
  setNames(
    Samples$ID, 
    nm = Samples$ID
  ), 
  function(o){
   FANS@meta.data %>% 
      filter(ID == o) %>% 
      pull(CB)
  }
)

sapply(
  names(Barcodes_All), 
  function(q){
    write.table(
      Barcodes_All[[q]], 
      paste0(
        "../Data/FANS/Bam/Barcodes/BySample/AllCells/", 
        q, 
        "_Barcodes.tsv"
      ), quote = FALSE, col.names = FALSE, row.names = FALSE
    )
  }
)

rm(Barcodes_All)




  ### 3.0 Extract Barcodes by Sample_CellType ----------------------------------

FANS$Sample_WNN_L15 <- paste0(
  FANS$ID, 
  "_", 
  "WNN_L15", 
  "_", 
  FANS$ID_WNN_L15_Predicted
)

FANS$Sample_WNN_L25 <- paste0(
  FANS$ID, 
  "_", 
  "WNN_L25", 
  "_", 
  FANS$ID_WNN_L25_Predicted
)

FANS$Sample_WNN_L4 <- paste0(
  FANS$ID, 
  "_", 
  "WNN_L4", 
  "_", 
  FANS$ID_WNN_L4_Predicted
)

table(FANS$Sample_WNN_L15) |> length()
table(FANS$Sample_WNN_L25) |> length() 
table(FANS$Sample_WNN_L4) |> length() 



FANS_Samples_Barcodes <- lapply(
  setNames(
    unique(FANS$ID), 
    nm = unique(FANS$ID)
  ), 
  function(Id){
    return(
      FANS@meta.data %>% 
        filter(ID==Id) %>% 
        mutate(
          V1=CB, 
          V2=paste0(
            Sample_WNN_L15,
            ",", 
            Sample_WNN_L25, 
            ",", 
            Sample_WNN_L4
          )
        ) %>% 
        select(V1,V2)
    )
  }
)




 ### 4.0 Extract Barcodes Per Cell ---------------------------------------------

table(FANS$ID)
FANS_Samples_Barcodes <- lapply(
  setNames(
    unique(FANS$ID), 
    nm = unique(FANS$ID)
  ), 
  function(Id){
    return(
      FANS@meta.data %>% 
        filter(ID==Id) %>% 
        mutate(V1=CB,V2=CellId) %>% 
        select(V1,V2)
    )
  }
)


lapply(
  names(FANS_Samples_Barcodes), 
  function(Id){
    df = FANS_Samples_Barcodes[[Id]]
    
      write.table(
        df, 
        paste0(
          "../Data/FANS/Bam/Barcodes/BySample/", 
          Id, 
          "_Sample_CellType_Barcodes.tsv"
        ), 
        col.names = FALSE, 
        row.names = FALSE, 
        quote = FALSE, 
        sep= '\t'
      )
    
  }
)

table(Samples$ID, Samples$TDP43)



