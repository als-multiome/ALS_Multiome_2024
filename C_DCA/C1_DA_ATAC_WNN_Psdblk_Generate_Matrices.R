
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


      # M0_ATAC Seurat single-cell data object ---------------------------------

M0_ATAC <- qread(
  paste0(
    "../Data/SeuratObjects/",  
    "M0_ATAC",
    ".qrds"
  ),
  nthr=nthr 
)




  ### 2.0 Generate pseudobulk matrices -----------------------------------------



    ## 2.1 Check how many unique cell types across all levels there are --------

unique_levels <- sum(
  length(unique(M0_ATAC$WNN_L1)), 
  length(unique(M0_ATAC$WNN_L15)), 
  length(unique(M0_ATAC$WNN_L2)), 
  length(unique(M0_ATAC$WNN_L25)), 
  length(unique(M0_ATAC$WNN_L3)), 
  length(unique(M0_ATAC$WNN_L4))
)



    ## 2.2 Prepare List of pseudobulk objects and a data frame of their characteristics ---- 

ATAC_Psdblk_WNN_Matrices <- list() 
ATAC_Psdblk_WNN_Matrices_Index <-data.frame(
  Index = rep(NA, unique_levels+1), 
  CelltypeLevel = rep(NA, unique_levels+1), 
  CellType=rep(NA, unique_levels+1) 
)

for(i in 1:nrow(ATAC_Psdblk_WNN_Matrices_Index)){
  ATAC_Psdblk_WNN_Matrices[[i]] <- NA
}



    ## 2.3 Add Pseudobulk list first element: All cells ------------------------

Pseudobulk <- AggregateExpression(
  M0_ATAC, 
  group.by = "ID", 
  assays = "ATAC",
  slot="counts"
) 

ATAC_Psdblk_WNN_Matrices[[1]] <- Pseudobulk$ATAC  
ATAC_Psdblk_WNN_Matrices_Index$Index[1] <- 1
ATAC_Psdblk_WNN_Matrices_Index$CelltypeLevel[1] <- "AllCells"
ATAC_Psdblk_WNN_Matrices_Index$CellType[1] <- "AllCells"
rm(Pseudobulk, i)



    ## 2.4 Add all other pseudobulk objects ------------------------------------

ind=2
for (CellTypeLevel in c("WNN_L1", "WNN_L15", "WNN_L2", "WNN_L25", "WNN_L3", "WNN_L4")){
  
  for(CellType in unique(M0_ATAC@meta.data[CellTypeLevel][,1])){
    
    if(CellTypeLevel=="WNN_L1"){
      M1 <- subset(M0_ATAC, subset = WNN_L1==CellType) 
    } 
    
    if(CellTypeLevel=="WNN_L15"){
      M1 <- subset(M0_ATAC, subset = WNN_L15==CellType) 
    }
    if(CellTypeLevel=="WNN_L2"){
      M1 <- subset(M0_ATAC, subset = WNN_L2==CellType) 
    } 
    if(CellTypeLevel=="WNN_L25"){
      M1 <- subset(M0_ATAC, subset = WNN_L25==CellType) 
    } 
    if(CellTypeLevel=="WNN_L3"){
      M1 <- subset(M0_ATAC, subset = WNN_L3==CellType) 
    }
    if(CellTypeLevel=="WNN_L4"){
      M1 <- subset(M0_ATAC, subset = WNN_L4==CellType) 
    } 
    
    tryCatch({
      Pseudobulk <- AggregateExpression(
        M1, 
        group.by = "ID", 
        assays="ATAC",
        slot="counts"
      )  
    
      ATAC_Psdblk_WNN_Matrices[[ind]] <- Pseudobulk$ATAC  
      ATAC_Psdblk_WNN_Matrices_Index$Index[ind] <- ind
      ATAC_Psdblk_WNN_Matrices_Index$CelltypeLevel[ind] <- CellTypeLevel
      ATAC_Psdblk_WNN_Matrices_Index$CellType[ind] <- CellType
    }, error=function(e){})
    
    message(paste0("ATAC Pseudobulk Matrix #", ind, " prepared... "))
    ind=ind+1 
  }
}

rm(M1, i, CellType, CellTypeLevel, ind, Pseudobulk)




  ### 3.0 Generate ATAC Link Peaks Pseudobulk matrices -------------------------

fts <- rownames(ATAC_Psdblk_WNN_Matrices[[1]])
all(sapply(ATAC_Psdblk_WNN_Matrices, FUN=function(x){return(all(rownames(x)==fts))}))


ind <- which(fts %in% ATAC_Links$peak)

ATAC_LinkPeaks_Psdblk_WNN_Matrices <- lapply(
  ATAC_Psdblk_WNN_Matrices, 
  FUN=function(x){
    return(
      x[ind,]
    )
  }
)

fts <- rownames(ATAC_LinkPeaks_Psdblk_WNN_Matrices[[1]]) 
all(sapply(ATAC_LinkPeaks_Psdblk_WNN_Matrices, FUN=function(x){return(all(rownames(x)==fts))}))

ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index <- ATAC_Psdblk_WNN_Matrices_Index




  ### 4.0 Inspect Pseudobulk matrices ------------------------------------------



    ## 4.1 Quality  checks -----------------------------------------------------

ATAC_Psdblk_WNN_Matrices_Index$MatrixGenerated <- !is.na(ATAC_Psdblk_WNN_Matrices)
any(is.na(ATAC_Psdblk_WNN_Matrices))

ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$MatrixGenerated <- !is.na(ATAC_LinkPeaks_Psdblk_WNN_Matrices)
any(is.na(ATAC_LinkPeaks_Psdblk_WNN_Matrices))


ATAC_Psdblk_WNN_Matrices_Index$nSamples <- NA
ATAC_Psdblk_WNN_Matrices_Index$nFeatures <- NA

ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$nSamples <- NA
ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$nFeatures <- NA

 
for(i in 1:length(ATAC_Psdblk_WNN_Matrices)){
  ATAC_Psdblk_WNN_Matrices_Index$nFeatures[i] <- dim(ATAC_Psdblk_WNN_Matrices[[i]])[1] 
  ATAC_Psdblk_WNN_Matrices_Index$nSamples[i] <- dim(ATAC_Psdblk_WNN_Matrices[[i]])[2]
}

for(i in 1:length(ATAC_LinkPeaks_Psdblk_WNN_Matrices)){
  ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$nFeatures[i] <- dim(ATAC_LinkPeaks_Psdblk_WNN_Matrices[[i]])[1] 
  ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$nSamples[i] <- dim(ATAC_LinkPeaks_Psdblk_WNN_Matrices[[i]])[2]
}

all(ATAC_Psdblk_WNN_Matrices_Index$nFeatures==dim(M0_ATAC)[1])
all(ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$nFeatures==length(fts))
rm(i, fts, ind)


any(
  sapply(
    ATAC_Psdblk_WNN_Matrices, 
    FUN=function(x){
      return(
        any(
          colSums(x)==0
        )
      )
    }
  )
)


any(
  sapply(
    ATAC_LinkPeaks_Psdblk_WNN_Matrices, 
    FUN=function(x){
      return(
        any(
          colSums(x)==0
        )
      )
    }
  )
)



table(ATAC_Psdblk_WNN_Matrices_Index$nSamples[ATAC_Psdblk_WNN_Matrices_Index$CelltypeLevel=="AllCells"])
table(ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$nSamples[ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$CelltypeLevel=="AllCells"])

table(ATAC_Psdblk_WNN_Matrices_Index$nSamples[ATAC_Psdblk_WNN_Matrices_Index$CelltypeLevel=="WNN_L1"])
table(ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$nSamples[ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$CelltypeLevel=="WNN_L1"])

table(ATAC_Psdblk_WNN_Matrices_Index$nSamples[ATAC_Psdblk_WNN_Matrices_Index$CelltypeLevel=="WNN_L15"])
table(ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$nSamples[ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$CelltypeLevel=="WNN_L15"])

table(ATAC_Psdblk_WNN_Matrices_Index$nSamples[ATAC_Psdblk_WNN_Matrices_Index$CelltypeLevel=="WNN_L2"])
table(ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$nSamples[ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$CelltypeLevel=="WNN_L2"])

table(ATAC_Psdblk_WNN_Matrices_Index$nSamples[ATAC_Psdblk_WNN_Matrices_Index$CelltypeLevel=="WNN_L25"])
table(ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$nSamples[ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$CelltypeLevel=="WNN_L25"])

table(ATAC_Psdblk_WNN_Matrices_Index$nSamples[ATAC_Psdblk_WNN_Matrices_Index$CelltypeLevel=="WNN_L3"])
table(ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$nSamples[ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$CelltypeLevel=="WNN_L3"])

table(ATAC_Psdblk_WNN_Matrices_Index$nSamples[ATAC_Psdblk_WNN_Matrices_Index$CelltypeLevel=="WNN_L4"])
table(ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$nSamples[ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$CelltypeLevel=="WNN_L4"])



Samples.tmp <- lapply(ATAC_Psdblk_WNN_Matrices, function(x){return(colnames(x))}) 

SampData.tmp <-  M0_ATAC@meta.data %>% 
  group_by(ID) %>% 
  summarise(Case=Case[1]) 


Samples_LinkPeaks.tmp <- lapply(ATAC_LinkPeaks_Psdblk_WNN_Matrices, function(x){return(colnames(x))}) 

SampData_LinkPeaks.tmp <-  M0_ATAC@meta.data %>% 
  group_by(ID) %>% 
  summarise(Case=Case[1]) 



    ## 4.2 Count cases ---------------------------------------------------------

count.Cases <- function(x, Case){
  nSamples.tmp <- table(SampData.tmp$Case[match(x, SampData.tmp$ID)])
  return(as.numeric(nSamples.tmp[Case]))
} 


for (Case.tmp in unique(SampData.tmp$Case)){
  ATAC_Psdblk_WNN_Matrices_Index[[paste0("n", Case.tmp)]] <- sapply(
    Samples.tmp, 
    FUN = count.Cases, 
    Case=Case.tmp
  ) 
}

for (Case.tmp in unique(SampData_LinkPeaks.tmp$Case)){
  ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index[[paste0("n", Case.tmp)]] <- sapply(
    Samples_LinkPeaks.tmp, 
    FUN = count.Cases, 
    Case=Case.tmp
  ) 
}


      # Double-Check 
all(
  ATAC_Psdblk_WNN_Matrices_Index$nSamples==
    rowSums(ATAC_Psdblk_WNN_Matrices_Index[c("nALS", "nHC", "nALS_FTD")], na.rm = TRUE)
)

all(
  ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$nSamples==
    rowSums(ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index[c("nALS", "nHC", "nALS_FTD")], na.rm = TRUE)
)

rm(Samples.tmp, Samples_LinkPeaks.tmp, SampData.tmp, SampData_LinkPeaks.tmp, count.Cases, Case.tmp) 



    ## 4.3 Count cells ---------------------------------------------------------

ATAC_Psdblk_WNN_Matrices_Index$nCells <- NA 
ATAC_Psdblk_WNN_Matrices_Index$nCellsALS <- NA 
ATAC_Psdblk_WNN_Matrices_Index$nCellsALSFTD <- NA 
ATAC_Psdblk_WNN_Matrices_Index$nCellsHC <- NA 
ATAC_Psdblk_WNN_Matrices_Index$nCellsALSandHC <- NA 

ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$nCells <- NA 
ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$nCellsALS <- NA 
ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$nCellsALSFTD <- NA 
ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$nCellsHC <- NA 
ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$nCellsALSandHC <- NA 


is.num.tmp <- function(x){return(if(length(x)==0) 0 else x)}

for (i in 1:nrow(ATAC_Psdblk_WNN_Matrices_Index)){
  
  M0_ATAC@meta.data %>% 
    group_by_at(ATAC_Psdblk_WNN_Matrices_Index$CelltypeLevel[i]) %>% 
    summarise(nCells=length(CellId)) %>%
      filter(!!as.symbol(ATAC_Psdblk_WNN_Matrices_Index$CelltypeLevel[i])==ATAC_Psdblk_WNN_Matrices_Index$CellType[i]) %>% 
        select(nCells) |> as.numeric() |> is.num.tmp() -> 
          ATAC_Psdblk_WNN_Matrices_Index$nCells[i]
  
  M0_ATAC@meta.data %>% 
    group_by(!!as.symbol(ATAC_Psdblk_WNN_Matrices_Index$CelltypeLevel[i]), Case) %>% 
      summarise(N=length(CellId)) %>% 
        filter(!!as.symbol(ATAC_Psdblk_WNN_Matrices_Index$CelltypeLevel[i])==ATAC_Psdblk_WNN_Matrices_Index$CellType[i] & Case=="ALS") %$% 
          N |> as.numeric() |> is.num.tmp() -> ATAC_Psdblk_WNN_Matrices_Index$nCellsALS[i]
  
  M0_ATAC@meta.data %>% 
    group_by(!!as.symbol(ATAC_Psdblk_WNN_Matrices_Index$CelltypeLevel[i]), Case) %>% 
      summarise(N=length(CellId)) %>% 
        filter(!!as.symbol(ATAC_Psdblk_WNN_Matrices_Index$CelltypeLevel[i])==ATAC_Psdblk_WNN_Matrices_Index$CellType[i] & Case=="ALS_FTD") %$% 
          N |> as.numeric() |> is.num.tmp() -> ATAC_Psdblk_WNN_Matrices_Index$nCellsALSFTD[i] 
  
  M0_ATAC@meta.data %>% 
    group_by(!!as.symbol(ATAC_Psdblk_WNN_Matrices_Index$CelltypeLevel[i]), Case) %>% 
      summarise(N=length(CellId)) %>% 
        filter(!!as.symbol(ATAC_Psdblk_WNN_Matrices_Index$CelltypeLevel[i])==ATAC_Psdblk_WNN_Matrices_Index$CellType[i] & Case=="HC") %$% 
          N |> as.numeric() |> is.num.tmp() -> ATAC_Psdblk_WNN_Matrices_Index$nCellsHC[i]
  
}


ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$nCellsALSandHC <- ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$nCellsALS + ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$nCellsHC


all(rownames(M0_ATAC)==rownames(ATAC_Psdblk_WNN_Matrices[[1]]))  
fts <- which(rownames(M0_ATAC) %in% ATAC_Links$peak)
ATAC_Links_Assay <- M0_ATAC@assays$ATAC$counts[fts,] 
all(rownames(M0_ATAC[fts])==rownames(ATAC_Links_Assay))
all(colnames(M0_ATAC)==colnames(ATAC_Links_Assay))

M0_ATAC$nCount_ATAC_LinkPeaks <- colSums(ATAC_Links_Assay)
any(
  M0_ATAC$nCount_ATAC_LinkPeaks==0
)

for (i in 1:nrow(ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index)){
  
  
  
  M0_ATAC@meta.data %>% 
    group_by_at(ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$CelltypeLevel[i]) %>% 
    summarise(nCells=length(CellId)) %>%
    filter(!!as.symbol(ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$CelltypeLevel[i])==ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$CellType[i]) %>% 
    select(nCells) |> as.numeric() |> is.num.tmp() -> 
    ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$nCells[i]
  
  M0_ATAC@meta.data %>% 
    group_by(!!as.symbol(ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$CelltypeLevel[i]), Case) %>% 
    summarise(N=length(CellId)) %>% 
    filter(!!as.symbol(ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$CelltypeLevel[i])==ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$CellType[i] & Case=="ALS") %$% 
    N |> as.numeric() |> is.num.tmp() -> ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$nCellsALS[i]
  
  M0_ATAC@meta.data %>% 
    group_by(!!as.symbol(ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$CelltypeLevel[i]), Case) %>% 
    summarise(N=length(CellId)) %>% 
    filter(!!as.symbol(ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$CelltypeLevel[i])==ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$CellType[i] & Case=="ALS_FTD") %$% 
    N |> as.numeric() |> is.num.tmp() -> ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$nCellsALSFTD[i] 
  
  M0_ATAC@meta.data %>% 
    group_by(!!as.symbol(ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$CelltypeLevel[i]), Case) %>% 
    summarise(N=length(CellId)) %>% 
    filter(!!as.symbol(ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$CelltypeLevel[i])==ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$CellType[i] & Case=="HC") %$% 
    N |> as.numeric() |> is.num.tmp() -> ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$nCellsHC[i]
  
}


ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$nCellsALSandHC <- ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$nCellsALS + ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$nCellsHC


      # Double-check 

all(
  ATAC_Psdblk_WNN_Matrices_Index$nCellsALS + 
    ATAC_Psdblk_WNN_Matrices_Index$nCellsALSFTD + 
      ATAC_Psdblk_WNN_Matrices_Index$nCellsHC == 
        ATAC_Psdblk_WNN_Matrices_Index$nCells
)

all(
  ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$nCellsALS + 
    ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$nCellsALSFTD + 
      ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$nCellsHC == 
        ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$nCells
)

cor.test(
  ATAC_Psdblk_WNN_Matrices_Index$nCells, 
  ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index$nCells
)




  ### 5.0 Serialize ATAC Pseudobulk data ----------------------------------------

qsave(
  ATAC_Psdblk_WNN_Matrices, 
  paste0(
    "../Data/DA/WNN/", 
    "ATAC_Psdblk_WNN_Matrices", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  ATAC_Psdblk_WNN_Matrices_Index, 
  paste0(
    "../Data/DA/WNN/", 
    "ATAC_Psdblk_WNN_Matrices_Index", 
    ".qrds"
  ), 
  nthr=nthr
)


qsave(
  ATAC_LinkPeaks_Psdblk_WNN_Matrices, 
  paste0(
    "../Data/DA/WNN/", 
    "ATAC_LinkPeaks_Psdblk_WNN_Matrices", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index, 
  paste0(
    "../Data/DA/WNN/", 
    "ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index", 
    ".qrds"
  ), 
  nthr=nthr
)