# A1_SoupX_Ambient_RNA_Correction.R 

 
### 0.0 Load libraries ------------------------------------------------------- 

library(SoupX)
library(Seurat)
library(glmGamPoi)
library(data.table)
library(qs2)
library(tidyverse)




  ### 1.0 Load data ------------------------------------------------------------

Samples <- data.table::fread(
  "../Data/Input/FANS_Seq_TDP43_Samples.txt", 
  data.table = FALSE
)
nrow(Samples) == length(unique(Samples$ID))


Data_Dirs <- read.table("../Data/cfg/Files_list.txt", header = FALSE)
Data_Dirs$V2 <- str_replace_all(Data_Dirs$V2, "\\$HOME", "~")
FANS_Cellranger_Output_Dir <- Data_Dirs$V2[Data_Dirs$V1=="FANS_Cellranger_Output_Dir"]




  ### 2.0 SoupX RNA Count matrix correction ------------------------------------

setequal(
  Samples$Name, 
  list.files(FANS_Cellranger_Output_Dir)  
) | all(
  Samples$Name %in% 
    list.files(FANS_Cellranger_Output_Dir)
)

SoupX_Objects <- list()


for (Sample in Samples$Name){
  
  cellranger_outs_raw <- file.path(FANS_Cellranger_Output_Dir, sprintf("%s/outs/raw_feature_bc_matrix/", Sample)) 
  cellranger_outs_filtered <- file.path(FANS_Cellranger_Output_Dir, sprintf("%s/outs/filtered_feature_bc_matrix/", Sample)) 
  
  raw_data <- Read10X(cellranger_outs_raw)  
  filtered_data <- Read10X(cellranger_outs_filtered) 
  
  clusters <- file.path(FANS_Cellranger_Output_Dir, sprintf("%s/outs/analysis/clustering/gene_expression_graphclust/clusters.csv", Sample))
  clusters <- fread(clusters, data.table = FALSE)
  
  sc = SoupChannel(raw_data, filtered_data)
  sc = setClusters(sc, clusters = setNames(clusters$Cluster, nm = clusters$Barcode))
  tryCatch(
    {
      sc = autoEstCont(sc, forceAccept = TRUE)
      out = adjustCounts(sc, roundToInt = TRUE) 
      SoupX_Objects[[Sample]] <- list(
        SoupChannel = sc, 
        AdjustedCounts = out, 
        RNA_SoupX_Rho = 1-(sum(out)/sum(filtered_data)) 
      )
    }, error = function(e){
      SoupX_Objects[[Sample]] <- NA
    }
  )
  

  rm(
    cellranger_outs_raw, 
    cellranger_outs_filtered, 
    raw_data, 
    filtered_data, 
    clusters, 
    sc, 
    out
  )
} 

all(names(SoupX_Objects) == Samples$Name)
SoupX_Rhos <- sapply(
  SoupX_Objects, 
  function(x){
    x$RNA_SoupX_Rho
  }
)




  ### 3.0 Samples SoupX Rhos ---------------------------------------------------

Samples$SoupX_Rho <- SoupX_Rhos[match(Samples$Name, names(SoupX_Rhos))]




  ### 4.0 Save data ------------------------------------------------------------

qs_save(
  SoupX_Objects, 
  "../Data/FANS/SoupX/SoupX_Objects.qs2"
)

qs_save(
  Samples, 
  "../Data/FANS/Samples/Samples.qs2"
)




