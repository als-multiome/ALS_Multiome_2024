
  ### 0.0 Load libraries -------------------------------------------------------

source("~/ALS_Brain_Multiome.Rcfg")

library(qs) 
library(data.table)
library(Seurat)
library(Signac)
library(tidyverse) 
library(DESeq2)
library(magrittr)




  ### 1.0 Load data ------------------------------------------------------------


      # M0_RNA Seurat single-cell data object ----------------------------------

M0_RNA <- qread(
  paste0(
    "../Data/SeuratObjects/",  
    "M0_RNA",
    ".qrds"
  ),
  nthr=nthr 
)

M0_RNA 




  ### 2.0 Generate pseudobulk matrices -----------------------------------------



    ## 2.1 Check how many unique cell types across all levels there are --------

unique_levels <- sum(
  length(unique(M0_RNA$WNN_L1)), 
  length(unique(M0_RNA$WNN_L15)), 
  length(unique(M0_RNA$WNN_L2)), 
  length(unique(M0_RNA$WNN_L25)), 
  length(unique(M0_RNA$WNN_L3)), 
  length(unique(M0_RNA$WNN_L4))
)



    ## 2.2 Prepare List of pseudobulk objects and a data frame of their characteristics ---- 

RNA_Psdblk_WNN_Matrices <- list() 
RNA_Psdblk_WNN_Matrices_Index <-data.frame(
  Index = rep(NA, unique_levels+1), 
  CelltypeLevel = rep(NA, unique_levels+1), 
  CellType=rep(NA, unique_levels+1) 
)

for(i in 1:nrow(RNA_Psdblk_WNN_Matrices_Index)){
  RNA_Psdblk_WNN_Matrices[[i]] <- NA
}



    ## 2.3 Add Pseudobulk list first element: All cells ------------------------

Pseudobulk <- AggregateExpression(
  M0_RNA, 
  group.by = "ID", 
  assays = "RNA",
  slot="counts"
) 

RNA_Psdblk_WNN_Matrices[[1]] <- Pseudobulk$RNA  
RNA_Psdblk_WNN_Matrices_Index$Index[1] <- 1
RNA_Psdblk_WNN_Matrices_Index$CelltypeLevel[1] <- "AllCells"
RNA_Psdblk_WNN_Matrices_Index$CellType[1] <- "AllCells"
rm(Pseudobulk)



    ## 2.4 Add all other pseudobulk objects ------------------------------------

ind=2
for (CellTypeLevel in c("WNN_L1", "WNN_L15", "WNN_L2", "WNN_L25", "WNN_L3", "WNN_L4")){
  
  for(CellType in unique(M0_RNA@meta.data[CellTypeLevel][,1])){
    
    if(CellTypeLevel=="WNN_L1"){
      M1 <- subset(M0_RNA, subset = WNN_L1==CellType) 
    } 
    
    if(CellTypeLevel=="WNN_L15"){
      M1 <- subset(M0_RNA, subset = WNN_L15==CellType) 
    }
    if(CellTypeLevel=="WNN_L2"){
      M1 <- subset(M0_RNA, subset = WNN_L2==CellType) 
    } 
    if(CellTypeLevel=="WNN_L25"){
      M1 <- subset(M0_RNA, subset = WNN_L25==CellType) 
    } 
    if(CellTypeLevel=="WNN_L3"){
      M1 <- subset(M0_RNA, subset = WNN_L3==CellType) 
    }
    if(CellTypeLevel=="WNN_L4"){
      M1 <- subset(M0_RNA, subset = WNN_L4==CellType) 
    } 
    
    tryCatch({
      Pseudobulk <- AggregateExpression(
        M1, 
        group.by = "ID", 
        assays="RNA",
        slot="counts"
      )  
    
      RNA_Psdblk_WNN_Matrices[[ind]] <- Pseudobulk$RNA  
      RNA_Psdblk_WNN_Matrices_Index$Index[ind] <- ind
      RNA_Psdblk_WNN_Matrices_Index$CelltypeLevel[ind] <- CellTypeLevel
      RNA_Psdblk_WNN_Matrices_Index$CellType[ind] <- CellType
    }, error=function(e){})
    
    ind=ind+1 
  }
}

rm(M1, i, CellType, CellTypeLevel, ind, Pseudobulk)




  ### 3.0 Inspect Pseudobulk matrices ------------------------------------------



    ## 3.1 Quality  checks -----------------------------------------------------

RNA_Psdblk_WNN_Matrices_Index$MatrixGenerated <- !is.na(RNA_Psdblk_WNN_Matrices)
any(is.na(RNA_Psdblk_WNN_Matrices))

RNA_Psdblk_WNN_Matrices_Index$nSamples <- NA
RNA_Psdblk_WNN_Matrices_Index$nFeatures <- NA

for(i in 1:length(RNA_Psdblk_WNN_Matrices)){
  RNA_Psdblk_WNN_Matrices_Index$nFeatures[i] <- dim(RNA_Psdblk_WNN_Matrices[[i]])[1] 
  RNA_Psdblk_WNN_Matrices_Index$nSamples[i] <- dim(RNA_Psdblk_WNN_Matrices[[i]])[2]
}

all(RNA_Psdblk_WNN_Matrices_Index$nFeatures==dim(M0_RNA)[1])

table(RNA_Psdblk_WNN_Matrices_Index$nSamples[RNA_Psdblk_WNN_Matrices_Index$CelltypeLevel=="AllCells"])
table(RNA_Psdblk_WNN_Matrices_Index$nSamples[RNA_Psdblk_WNN_Matrices_Index$CelltypeLevel=="WNN_L1"])
table(RNA_Psdblk_WNN_Matrices_Index$nSamples[RNA_Psdblk_WNN_Matrices_Index$CelltypeLevel=="WNN_L15"])
table(RNA_Psdblk_WNN_Matrices_Index$nSamples[RNA_Psdblk_WNN_Matrices_Index$CelltypeLevel=="WNN_L2"])
table(RNA_Psdblk_WNN_Matrices_Index$nSamples[RNA_Psdblk_WNN_Matrices_Index$CelltypeLevel=="WNN_L25"])
table(RNA_Psdblk_WNN_Matrices_Index$nSamples[RNA_Psdblk_WNN_Matrices_Index$CelltypeLevel=="WNN_L3"])
table(RNA_Psdblk_WNN_Matrices_Index$nSamples[RNA_Psdblk_WNN_Matrices_Index$CelltypeLevel=="WNN_L4"])

Samples.tmp <- lapply(RNA_Psdblk_WNN_Matrices, function(x){return(colnames(x))}) 

SampData.tmp <-  M0_RNA@meta.data %>% 
  group_by(ID) %>% 
  summarise(Case=Case[1]) 



    ## 3.2 Count cases ---------------------------------------------------------

count.Cases <- function(x, Case){
  nSamples.tmp <- table(SampData.tmp$Case[match(x, SampData.tmp$ID)])
  return(as.numeric(nSamples.tmp[Case]))
}

for (Case.tmp in unique(SampData.tmp$Case)){
  RNA_Psdblk_WNN_Matrices_Index[[paste0("n", Case.tmp)]] <- sapply(
    Samples.tmp, 
    FUN = count.Cases, 
    Case=Case.tmp
  ) 
}

      # Double-Check 
all(
  RNA_Psdblk_WNN_Matrices_Index$nSamples==
    rowSums(RNA_Psdblk_WNN_Matrices_Index[c("nALS", "nHC", "nALS_FTD")], na.rm = TRUE)
)


rm(Samples.tmp, SampData.tmp, count.Cases, Case.tmp, i) 



    ## 3.3 Count cells ---------------------------------------------------------

RNA_Psdblk_WNN_Matrices_Index$nCells <- NA 
RNA_Psdblk_WNN_Matrices_Index$nCellsALS <- NA 
RNA_Psdblk_WNN_Matrices_Index$nCellsALSFTD <- NA 
RNA_Psdblk_WNN_Matrices_Index$nCellsHC <- NA 
RNA_Psdblk_WNN_Matrices_Index$nCellsALSandHC <- NA 


is.num.tmp <- function(x){return(if(length(x)==0) 0 else x)}

for (i in 1:nrow(RNA_Psdblk_WNN_Matrices_Index)){
  
  M0_RNA@meta.data %>% 
    group_by_at(RNA_Psdblk_WNN_Matrices_Index$CelltypeLevel[i]) %>% 
    summarise(nCells=length(CellId)) %>%
      filter(!!as.symbol(RNA_Psdblk_WNN_Matrices_Index$CelltypeLevel[i])==RNA_Psdblk_WNN_Matrices_Index$CellType[i]) %>% 
        select(nCells) |> as.numeric() |> is.num.tmp() -> 
          RNA_Psdblk_WNN_Matrices_Index$nCells[i]
  
  M0_RNA@meta.data %>% 
    group_by(!!as.symbol(RNA_Psdblk_WNN_Matrices_Index$CelltypeLevel[i]), Case) %>% 
      summarise(N=length(CellId)) %>% 
        filter(!!as.symbol(RNA_Psdblk_WNN_Matrices_Index$CelltypeLevel[i])==RNA_Psdblk_WNN_Matrices_Index$CellType[i] & Case=="ALS") %$% 
          N |> as.numeric() |> is.num.tmp() -> RNA_Psdblk_WNN_Matrices_Index$nCellsALS[i]
  
  M0_RNA@meta.data %>% 
    group_by(!!as.symbol(RNA_Psdblk_WNN_Matrices_Index$CelltypeLevel[i]), Case) %>% 
      summarise(N=length(CellId)) %>% 
        filter(!!as.symbol(RNA_Psdblk_WNN_Matrices_Index$CelltypeLevel[i])==RNA_Psdblk_WNN_Matrices_Index$CellType[i] & Case=="ALS_FTD") %$% 
          N |> as.numeric() |> is.num.tmp() -> RNA_Psdblk_WNN_Matrices_Index$nCellsALSFTD[i] 
  
  M0_RNA@meta.data %>% 
    group_by(!!as.symbol(RNA_Psdblk_WNN_Matrices_Index$CelltypeLevel[i]), Case) %>% 
      summarise(N=length(CellId)) %>% 
        filter(!!as.symbol(RNA_Psdblk_WNN_Matrices_Index$CelltypeLevel[i])==RNA_Psdblk_WNN_Matrices_Index$CellType[i] & Case=="HC") %$% 
          N |> as.numeric() |> is.num.tmp() -> RNA_Psdblk_WNN_Matrices_Index$nCellsHC[i]
  
}

rm(is.num.tmp, i)

RNA_Psdblk_WNN_Matrices_Index$nCellsALSandHC <- RNA_Psdblk_WNN_Matrices_Index$nCellsALS + RNA_Psdblk_WNN_Matrices_Index$nCellsHC


      # Double-check 

all(
  RNA_Psdblk_WNN_Matrices_Index$nCellsALS + 
    RNA_Psdblk_WNN_Matrices_Index$nCellsALSFTD + 
      RNA_Psdblk_WNN_Matrices_Index$nCellsHC == 
        RNA_Psdblk_WNN_Matrices_Index$nCells
)




  ### 4.0 Serialize RNA Pseudobulk data ----------------------------------------

qsave(
  RNA_Psdblk_WNN_Matrices, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/", 
    "RNA_Psdblk_WNN_Matrices", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  RNA_Psdblk_WNN_Matrices_Index, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/", 
    "RNA_Psdblk_WNN_Matrices_Index", 
    ".qrds"
  ), 
  nthr=nthr
)

