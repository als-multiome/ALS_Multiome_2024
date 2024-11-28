


  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 

library(qs)
library(spatialLIBD)
library(Seurat) 
library(Signac)
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
  )
)

M0_ATAC <- qread(
  paste0(
    "../Data/SeuratObjects/", 
    "M0_ATAC", 
    ".qrds"
  )
) 


LIBD_Spatial_Seurat <- qread( 
  paste0(
    "../Data/SpatialData/", 
    "LIBD_Spatial_Seurat", 
    ".qrds"
  ), 
  nthr = nthr
)


libd_layer_colors2 <- libd_layer_colors
libd_layer_colors2["WM"] <- "#A4A4A4"
libd_layer_colors3 <- libd_layer_colors2
libd_layer_colors3["WM"] <- "#515151" 




  ### 2.0 Calculate M0_RNA random score controls with HC samples ---------------

M0_HC <- subset(M0_RNA, subset=Case=="HC") 


dfs <- list() 

for (nfts in c(50, 100, 200, 500, 1000)){
  
  for (i in 1:10){
    
    Rand_fts <- sample(1:dim(M0_HC)[1], size = nfts, replace = FALSE)
    Rand_fts <- Features(M0_HC, assay = "RNA")[Rand_fts] 
    
    M0_HC <- AddModuleScore(
      M0_HC, 
      features = list(Rand_fts), 
      ctrl=500, 
      assay="RNA", 
      name=paste0(
        "Random_Score", 
        i, 
        "_"
      )
    )
  }
  
  df <- M0_HC@meta.data 
  df <- df[,colnames(df) %in% c("WNN_L25", paste0("Random_Score", 1:10, "_1"))]
  for (i in 2:11){
    df[,i] <- scale(df[,i])
  }
  df$ScoreScaledCombined <- rowMeans(df[,2:11])
  dfs[[paste0("NFTS_", nfts, "_df")]] <- df 
  rm(df)
} 
 


M0_HC <- subset(M0_RNA, subset=Case=="HC")  


nfts=500

  for (i in 1:100){
    
    print(Sys.time())
    message(paste0("Iteration ", i, "...")) 
    
    Rand_fts <- sample(1:dim(M0_HC)[1], size = nfts, replace = FALSE)
    Rand_fts <- Features(M0_HC, assay = "RNA")[Rand_fts] 
    
    M0_HC <- AddModuleScore(
      M0_HC, 
      features = list(Rand_fts), 
      ctrl=500, 
      assay="RNA", 
      name=paste0(
        "Random_Score", 
        i, 
        "_"
      )
    )
}

  
df <- M0_HC@meta.data 
df <- df[,colnames(df) %in% c("WNN_L25", paste0("Random_Score", 1:100, "_1"))]
 

for (i in 2:101){
  df[,i] <- scale(df[,i])
}

df$ScoreScaledCombined <- rowMeans(df[,2:101])


dfs[["NFTS_500_100x_df"]] <- df 
 
qsave(
  dfs, 
  paste0(
    "../Data/Annotations/Module_Scores/RNA/HC_Derived/", 
    "M0_HC_Random_ModuleScores", 
    ".qrds"
  ), 
  nthr=nthr
)
 
rm(M0_HC, df, dfs, i, nfts)




  ### 3.0 Calculate M0_RNA random score controls with ALS-FTD and HC samples ----

dfs <- list() 

for (nfts in c(50, 100, 200, 500, 1000)){
  
  for (i in 1:10){
    
    Rand_fts <- sample(1:dim(M0_RNA)[1], size = nfts, replace = FALSE)
    Rand_fts <- Features(M0_RNA, assay = "RNA")[Rand_fts] 
    
    M0_RNA <- AddModuleScore(
      M0_RNA, 
      features = list(Rand_fts), 
      ctrl=500, 
      assay="RNA", 
      name=paste0(
        "Random_Score", 
        i, 
        "_"
      )
    )
  }
  
  df <- M0_RNA@meta.data 
  df <- df[,colnames(df) %in% c("WNN_L25", paste0("Random_Score", 1:10, "_1"))]
  for (i in 2:11){
    df[,i] <- scale(df[,i])
  }
  df$ScoreScaledCombined <- rowMeans(df[,2:11])
  dfs[[paste0("NFTS_", nfts, "_df")]] <- df 
  rm(df)
} 


M0_RNA <- qread(
  paste0(
    "../Data/SeuratObjects/", 
    "M0_RNA", 
    ".qrds"
  )
)

nfts=500

for (i in 1:100){
  
  print(Sys.time())
  message(paste0("Iteration ", i, "...")) 
  
  Rand_fts <- sample(1:dim(M0_RNA)[1], size = nfts, replace = FALSE)
  Rand_fts <- Features(M0_RNA, assay = "RNA")[Rand_fts] 
  
  M0_RNA <- AddModuleScore(
    M0_RNA, 
    features = list(Rand_fts), 
    ctrl=500, 
    assay="RNA", 
    name=paste0(
      "Random_Score", 
      i, 
      "_"
    )
  )
}


df <- M0_RNA@meta.data 
df <- df[,colnames(df) %in% c("WNN_L25", paste0("Random_Score", 1:100, "_1"))]


for (i in 2:101){
  df[,i] <- scale(df[,i])
}

df$ScoreScaledCombined <- rowMeans(df[,2:101])


dfs[["NFTS_500_100x_df"]] <- df 

qsave(
  dfs, 
  paste0(
    "../Data/Annotations/Module_Scores/RNA/HC_ALSFTD_Derived/", 
    "M0_HC_ALSFTD_Random_ModuleScores", 
    ".qrds"
  ), 
  nthr=nthr
)

rm(M0_RNA, df, dfs, i, nfts)




  ### 4.0 Calculate Spatial random score controls with HC samples --------------

dfs <- list() 

for (nfts in c(50, 100, 200, 500, 1000)){
  
  for (i in 1:10){
    
    Rand_fts <- sample(1:dim(LIBD_Spatial_Seurat)[1], size = nfts, replace = FALSE)
    Rand_fts <- Features(LIBD_Spatial_Seurat, assay = "SpatialLogCounts")[Rand_fts] 
    
    LIBD_Spatial_Seurat <- AddModuleScore(
      LIBD_Spatial_Seurat, 
      features = list(Rand_fts), 
      ctrl=500, 
      assay="SpatialLogCounts", 
      name=paste0(
        "Random_Score", 
        i, 
        "_"
      )
    )
  }
  
  df <- LIBD_Spatial_Seurat@meta.data 
  df <- df[,colnames(df) %in% c("layer_guess_reordered", "key", paste0("Random_Score", 1:10, "_1"))]
  for (i in 3:12){
    df[,i] <- scale(df[,i])
  }
  df$ScoreScaledCombined <- rowMeans(df[,3:12])
  dfs[[paste0("NFTS_", nfts, "_df")]] <- df 
  rm(df)
} 


LIBD_Spatial_Seurat <- qread( 
  paste0(
    "../Data/SpatialData/", 
    "LIBD_Spatial_Seurat", 
    ".qrds"
  ), 
  nthr = nthr
)



nfts=500

for (i in 1:100){
  
  print(Sys.time())
  message(paste0("Iteration ", i, "...")) 
  
  Rand_fts <- sample(1:dim(LIBD_Spatial_Seurat)[1], size = nfts, replace = FALSE)
  Rand_fts <- Features(LIBD_Spatial_Seurat, assay = "SpatialLogCounts")[Rand_fts] 
  
  LIBD_Spatial_Seurat <- AddModuleScore(
    LIBD_Spatial_Seurat, 
    features = list(Rand_fts), 
    ctrl=500, 
    assay="SpatialLogCounts", 
    name=paste0(
      "Random_Score", 
      i, 
      "_"
    )
  )
}


df <- LIBD_Spatial_Seurat@meta.data 
df <- df[,colnames(df) %in% c("layer_guess_reordered", "key", paste0("Random_Score", 1:100, "_1"))]


for (i in 3:102){
  df[,i] <- scale(df[,i])
}

df$ScoreScaledCombined <- rowMeans(df[,3:102])


dfs[["NFTS_500_100x_df"]] <- df 

qsave(
  dfs, 
  paste0(
    "../Data/Annotations/Module_Scores/RNA/HC_Derived/", 
    "LIBD_Spatial_HC_Random_ModuleScores", 
    ".qrds"
  ), 
  nthr=nthr
)


  ### 5.0 Calculate M0_ATAC random score controls with HC samples -------------- 

M0_HC <- subset(M0_ATAC, subset=Case=="HC") 


dfs <- list() 

for (nfts in c(50, 100, 200, 500, 1000)){
  
  for (i in 1:10){
    
    Rand_fts <- sample(1:dim(M0_HC)[1], size = nfts, replace = FALSE)
    Rand_fts <- Features(M0_HC, assay = "ATAC")[Rand_fts] 
    
    M0_HC <- AddModuleScore(
      M0_HC, 
      features = list(Rand_fts), 
      ctrl=500, 
      assay="ATAC", 
      name=paste0(
        "Random_Score", 
        i, 
        "_"
      )
    )
  }
  
  df <- M0_HC@meta.data 
  df <- df[,colnames(df) %in% c("WNN_L25", paste0("Random_Score", 1:10, "_1"))]
  for (i in 2:11){
    df[,i] <- scale(df[,i])
  }
  df$ScoreScaledCombined <- rowMeans(df[,2:11])
  dfs[[paste0("NFTS_", nfts, "_df")]] <- df 
  rm(df)
} 



M0_HC <- subset(M0_ATAC, subset=Case=="HC") 


nfts=500

for (i in 1:100){
  
  print(Sys.time())
  message(paste0("Iteration ", i, "...")) 
  
  Rand_fts <- sample(1:dim(M0_HC)[1], size = nfts, replace = FALSE)
  Rand_fts <- Features(M0_HC, assay = "ATAC")[Rand_fts] 
  
  M0_HC <- AddModuleScore(
    M0_HC, 
    features = list(Rand_fts), 
    ctrl=500, 
    assay="ATAC", 
    name=paste0(
      "Random_Score", 
      i, 
      "_"
    )
  )
}


df <- M0_HC@meta.data 
df <- df[,colnames(df) %in% c("WNN_L25", paste0("Random_Score", 1:100, "_1"))]


for (i in 2:101){
  df[,i] <- scale(df[,i])
}

df$ScoreScaledCombined <- rowMeans(df[,2:101])


dfs[["NFTS_500_100x_df"]] <- df 

qsave(
  dfs, 
  paste0(
    "../Data/Annotations/Module_Scores/ATAC/HC_Derived/", 
    "M0_HC_Random_ModuleScores", 
    ".qrds"
  ), 
  nthr=nthr
)

rm(M0_HC, df, dfs, i, nfts)
