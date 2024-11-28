

  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 

library(qs)
library(tidyverse)
library(data.table)
library(xlsx)


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




  ### 2.0 Collect CellType percentages  ----------------------------------------

dfs <- list() 

for (CellTypeLevel in c("WNN_L1", "WNN_L15", "WNN_L25", "WNN_L4")){
  M0_RNA[[CellTypeLevel]] = as.character(M0_RNA[[CellTypeLevel]][,1])
  dfs[[CellTypeLevel]] = data.frame(
    CellType = names(table(M0_RNA[[CellTypeLevel]])), 
    N=as.numeric(table(M0_RNA[[CellTypeLevel]])), 
    Percent=if(CellTypeLevel=="WNN_L4") {
      as.numeric(round(prop.table(table(M0_RNA[[CellTypeLevel]]))*100,3))
    } else {
      as.numeric(round(prop.table(table(M0_RNA[[CellTypeLevel]]))*100,1))
    }
  )
}

dfs <- lapply(dfs, FUN=function(x){x[order(x$CellType),]})
dfs <- lapply(dfs, FUN=function(x){colnames(x)[1]<-"Cell Type"; return(x)})




  ### 3.0 Export Data ----------------------------------------------------------

CellTypeLevel_Names <- c(
  "WNN_L1"="Broad cell classes", 
  "WNN_L15"="Cell classes", 
  "WNN_L25"="Broad cell types", 
  "WNN_L4"="Cell types"
)

for (table in names(dfs)){
  write.xlsx(
    x = dfs[[table]], 
    file = paste0(
      "../Data/Visualization/Tables/", 
      "SupplTab_CellType_Percentages", 
      ".xlsx"
    ), 
    sheetName = CellTypeLevel_Names[table], 
    col.names = TRUE, 
    row.names = FALSE, 
    append = TRUE
  )
}


