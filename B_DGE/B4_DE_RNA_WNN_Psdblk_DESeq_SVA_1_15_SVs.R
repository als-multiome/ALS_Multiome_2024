
  ### 0.0 Load libraries -------------------------------------------------------

source("~/ALS_Brain_Multiome.Rcfg")


library(qs) 
library(data.table)
library(Seurat)
library(tidyverse) 
library(DESeq2)
library(BiocParallel)
library(sva)


qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)




  ### 1.0 AllCase --------------------------------------------------------------



    ## 1.1 Load Data ----------------------------------------------------------- 

DDS_list <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/", 
    "DDS_DESEqed_list", 
    ".qrds"
  ), 
  nthr=nthr
)


DESeq_Results_Index <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/", 
    "DESeq_Results_Index", 
    ".qrds"
  ), 
  nthr=nthr
)    

DDS_list <- DDS_list[1:76] 
DESeq_Results_Index <- DESeq_Results_Index[1:76,1:9]



    ## 1.2 Find 15 SVs --------------------------------------------------------- 

SVAs_15SVs <- list() 
for (i in 1:nrow(DESeq_Results_Index)){
  SVAs_15SVs[[i]] <- NA
}

SVAs_15SVs_Index <- data.frame(
  Index=1:nrow(DESeq_Results_Index), 
  SVA=NA, 
  SVs = NA 
)

SVAs_15SVs_Index$CellTypeLevel <- DESeq_Results_Index$CellTypeLevel
SVAs_15SVs_Index$CellType <- DESeq_Results_Index$CellType
SVAs_15SVs_Index$Comparison <- DESeq_Results_Index$Comparison


SVAs_15SVs_Index$SVs <- NA 

for (i in 1:nrow(SVAs_15SVs_Index)){
  tryCatch({
    
    dat  <- counts(DDS_list[[i]], normalized = TRUE)
    idx  <- rowMeans(dat) > 1 
    dat  <- dat[idx, ] 
    if(DESeq_Results_Index$Comparison[i]=="All_Cases"){
      mod  <- model.matrix(~ Case, colData(DDS_list[[i]]))
    }else{
      mod  <- model.matrix(~ Rand, colData(DDS_list[[i]]))
    }
    mod0 <- model.matrix(~ 1, colData(DDS_list[[i]]))
    
    svseq <- svaseq(dat, mod, mod0, n.sv=15)
    SVAs_15SVs[[i]] <- svseq
    
    SVAs_15SVs_Index$SVA[i] <- TRUE 
    SVAs_15SVs_Index$SVs[i] <- 15 
    message(paste0("SVAs #", i, " generated! "))
  }, error=function(e){
    
    SVAs_15SVs_Index$SVA[i] <<- FALSE 
    SVAs_15SVs_Index$SVs[i] <<- 0 
    message(paste0("SVAs #", i, " could not be generated :("))
  })
}


qsave(
  SVAs_15SVs, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/1_15_SVs/", 
    "SVAs_15SVs", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  SVAs_15SVs_Index, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/1_15_SVs/", 
    "SVAs_15SVs_Index", 
    ".qrds"
  ), 
  nthr=nthr
)


all(SVAs_15SVs_Index$SVs[SVAs_15SVs_Index$SVA]==15)
rm(dat, mod, mod0, svseq, i, idx)



    ## 1.3 Add 1-15 SVs to DDS_list --------------------------------------------

DESeq_SVA_1_15_Results_Index <- DESeq_Results_Index
rm(DESeq_Results_Index)


for (i in 1:15){
  
  DDS_SVA_list <- list() 
  for (j in 1:length(DDS_list)){
    DDS_SVA_list[[j]] <- NA
  } 
  
  DESeq_SVA_1_15_Results_Index[[paste0("SVA",i)]]  <- NA
  
  
  for (k in 1:nrow(DESeq_SVA_1_15_Results_Index)){ 
    
    tryCatch({ 
      
      DDS.tmp <- DDS_list[[k]] 
      
      for (l in 1:i){
        
        DDS.tmp[[paste0("SV", l)]] <- SVAs_15SVs[[k]]$sv[,l]
      }
      
      
      if(DESeq_SVA_1_15_Results_Index$Comparison[k]=="All_Cases"){
        if(i==1){
          design(DDS.tmp) <- ~ SV1 + Case 
        }
        if(i==2){
          design(DDS.tmp) <- ~ SV1 + SV2 + Case 
        }
        if(i==3){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + Case 
        }
        if(i==4){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + Case 
        }
        if(i==5){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + Case 
        }
        if(i==6){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + Case 
        }
        if(i==7){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + Case 
        }
        if(i==8){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + Case 
        }
        if(i==9){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + Case 
        }
        if(i==10){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + Case 
        }
        if(i==11){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + Case 
        }
        if(i==12){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + Case 
        }
        if(i==13){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + Case 
        }
        if(i==14){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + SV14 + Case 
        }
        if(i==15){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + SV14 + SV15 + Case 
        }
        
      }else{
        if(i==1){
          design(DDS.tmp) <- ~ SV1 + Rand 
        }
        if(i==2){
          design(DDS.tmp) <- ~ SV1 + SV2 + Rand 
        }
        if(i==3){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + Rand 
        }
        if(i==4){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + Rand 
        }
        if(i==5){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + Rand 
        }
        if(i==6){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + Rand 
        }
        if(i==7){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + Rand 
        }
        if(i==8){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + Rand 
        }
        if(i==9){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + Rand 
        }
        if(i==10){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + Rand 
        }
        if(i==11){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + Rand 
        }
        if(i==12){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + Rand 
        }
        if(i==13){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + Rand 
        }
        if(i==14){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + SV14 + Rand 
        }
        if(i==15){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + SV14 + SV15 + Rand 
        }
        
      }
      
      DDS_SVA_list[[k]] <- DDS.tmp 
      rm(DDS.tmp)
      
      DESeq_SVA_1_15_Results_Index[[paste0("SVA",i)]][k] <- TRUE 
      message(k)
      message("########################################################")
      message(Sys.time())
    }, error=function(e){
      DESeq_SVA_1_15_Results_Index[[paste0("SVA",i)]][k] <<- FALSE 
      message("########################################################")
      message(Sys.time()) 
    })
  }  
  
  qsave(
    DDS_SVA_list, 
    paste0(
      "../Data/De/WNN/Pseudobulk/DESeq2/AllCase/SVA/1_15_SVs/",
      "DDS_", 
      i, 
      "SVs_", 
      "list", 
      ".qrds"
    ), 
    nthr=nthr
  )
  rm(DDS_SVA_list)
}


      # START HPC Cluster Outsourcing ------------------------------------------


      # Y_B4_1_DE_RNA_WNN_SVA_AllCase_1SVs.R -----------------------------------


      # Import HPC Output and summarize results --------------------------------

##$$$----


  ### 2.0 ALS ------------------------------------------------------------------



    ## 2.1 Load Data ----------------------------------------------------------- 

DDS_list <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALS/", 
    "DDS_DESEqed_list", 
    ".qrds"
  ), 
  nthr=nthr
)


DESeq_Results_Index <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALS/", 
    "DESeq_Results_Index", 
    ".qrds"
  ), 
  nthr=nthr
)    

DDS_list <- DDS_list[1:76] 
DESeq_Results_Index <- DESeq_Results_Index[1:76,1:9]



    ## 2.2 Find 15 SVs --------------------------------------------------------- 

SVAs_15SVs <- list() 
for (i in 1:nrow(DESeq_Results_Index)){
  SVAs_15SVs[[i]] <- NA
}

SVAs_15SVs_Index <- data.frame(
  Index=1:nrow(DESeq_Results_Index), 
  SVA=NA, 
  SVs = NA 
)

SVAs_15SVs_Index$CellTypeLevel <- DESeq_Results_Index$CellTypeLevel
SVAs_15SVs_Index$CellType <- DESeq_Results_Index$CellType
SVAs_15SVs_Index$Comparison <- DESeq_Results_Index$Comparison


SVAs_15SVs_Index$SVs <- NA 

for (i in 1:nrow(SVAs_15SVs_Index)){
  tryCatch({
    
    dat  <- counts(DDS_list[[i]], normalized = TRUE)
    idx  <- rowMeans(dat) > 1 
    dat  <- dat[idx, ] 
    if(DESeq_Results_Index$Comparison[i]=="ALSvsHC"){
      mod  <- model.matrix(~ Case, colData(DDS_list[[i]]))
    }else{
      mod  <- model.matrix(~ Rand, colData(DDS_list[[i]]))
    }
    mod0 <- model.matrix(~ 1, colData(DDS_list[[i]]))
    
    svseq <- svaseq(dat, mod, mod0, n.sv=15)
    SVAs_15SVs[[i]] <- svseq
    
    SVAs_15SVs_Index$SVA[i] <- TRUE 
    SVAs_15SVs_Index$SVs[i] <- 15 
    message(paste0("SVAs #", i, " generated! "))
  }, error=function(e){
    
    SVAs_15SVs_Index$SVA[i] <<- FALSE 
    SVAs_15SVs_Index$SVs[i] <<- 0 
    message(paste0("SVAs #", i, " could not be generated :("))
  })
}


qsave(
  SVAs_15SVs, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALS/SVA/1_15_SVs/", 
    "SVAs_15SVs", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  SVAs_15SVs_Index, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALS/SVA/1_15_SVs/", 
    "SVAs_15SVs_Index", 
    ".qrds"
  ), 
  nthr=nthr
)


all(SVAs_15SVs_Index$SVs[SVAs_15SVs_Index$SVA]==15)
rm(dat, mod, mod0, svseq, i, idx)



    ## 2.3 Add 1-15 SVs to DDS_list --------------------------------------------

DESeq_SVA_1_15_Results_Index <- DESeq_Results_Index
rm(DESeq_Results_Index)


for (i in 1:15){
  
  DDS_SVA_list <- list() 
  for (j in 1:length(DDS_list)){
    DDS_SVA_list[[j]] <- NA
  } 
  
  DESeq_SVA_1_15_Results_Index[[paste0("SVA",i)]]  <- NA
  
  
  for (k in 1:nrow(DESeq_SVA_1_15_Results_Index)){ 
    
    tryCatch({ 
      
      DDS.tmp <- DDS_list[[k]] 
      
      for (l in 1:i){
        
        DDS.tmp[[paste0("SV", l)]] <- SVAs_15SVs[[k]]$sv[,l]
      }
      
      
      if(DESeq_SVA_1_15_Results_Index$Comparison[k]=="ALSvsHC"){
        if(i==1){
          design(DDS.tmp) <- ~ SV1 + Case 
        }
        if(i==2){
          design(DDS.tmp) <- ~ SV1 + SV2 + Case 
        }
        if(i==3){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + Case 
        }
        if(i==4){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + Case 
        }
        if(i==5){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + Case 
        }
        if(i==6){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + Case 
        }
        if(i==7){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + Case 
        }
        if(i==8){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + Case 
        }
        if(i==9){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + Case 
        }
        if(i==10){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + Case 
        }
        if(i==11){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + Case 
        }
        if(i==12){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + Case 
        }
        if(i==13){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + Case 
        }
        if(i==14){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + SV14 + Case 
        }
        if(i==15){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + SV14 + SV15 + Case 
        }
        
      }else{
        if(i==1){
          design(DDS.tmp) <- ~ SV1 + Rand 
        }
        if(i==2){
          design(DDS.tmp) <- ~ SV1 + SV2 + Rand 
        }
        if(i==3){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + Rand 
        }
        if(i==4){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + Rand 
        }
        if(i==5){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + Rand 
        }
        if(i==6){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + Rand 
        }
        if(i==7){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + Rand 
        }
        if(i==8){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + Rand 
        }
        if(i==9){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + Rand 
        }
        if(i==10){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + Rand 
        }
        if(i==11){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + Rand 
        }
        if(i==12){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + Rand 
        }
        if(i==13){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + Rand 
        }
        if(i==14){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + SV14 + Rand 
        }
        if(i==15){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + SV14 + SV15 + Rand 
        }
        
      }
      
      DDS_SVA_list[[k]] <- DDS.tmp 
      rm(DDS.tmp)
      
      DESeq_SVA_1_15_Results_Index[[paste0("SVA",i)]][k] <- TRUE 
      message(k)
      message("########################################################")
      message(Sys.time())
    }, error=function(e){
      DESeq_SVA_1_15_Results_Index[[paste0("SVA",i)]][k] <<- FALSE 
      message("########################################################")
      message(Sys.time()) 
    })
  }  
  
  qsave(
    DDS_SVA_list, 
    paste0(
      "../Data/De/WNN/Pseudobulk/DESeq2/ALS/SVA/1_15_SVs/",
      "DDS_", 
      i, 
      "SVs_", 
      "list", 
      ".qrds"
    ), 
    nthr=nthr
  )
  rm(DDS_SVA_list)
}

##$$$----
      # START HPC Cluster Outsourcing ------------------------------------------


      # Y_B4_2_DE_RNA_WNN_SVA_ALS_1SVs.R -----------------------------------


      # Import HPC Output and summarize results --------------------------------





  ### 3.0 ALSFTD ---------------------------------------------------------------



    ## 3.1 Load Data ----------------------------------------------------------- 

DDS_list <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD/", 
    "DDS_DESEqed_list", 
    ".qrds"
  ), 
  nthr=nthr
)


DESeq_Results_Index <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD/", 
    "DESeq_Results_Index", 
    ".qrds"
  ), 
  nthr=nthr
)    

DDS_list <- DDS_list[1:76] 
DESeq_Results_Index <- DESeq_Results_Index[1:76,1:9]



    ## 3.2 Find 15 SVs --------------------------------------------------------- 

SVAs_15SVs <- list() 
for (i in 1:nrow(DESeq_Results_Index)){
  SVAs_15SVs[[i]] <- NA
}

SVAs_15SVs_Index <- data.frame(
  Index=1:nrow(DESeq_Results_Index), 
  SVA=NA, 
  SVs = NA 
)

SVAs_15SVs_Index$CellTypeLevel <- DESeq_Results_Index$CellTypeLevel
SVAs_15SVs_Index$CellType <- DESeq_Results_Index$CellType
SVAs_15SVs_Index$Comparison <- DESeq_Results_Index$Comparison


SVAs_15SVs_Index$SVs <- NA 

for (i in 1:nrow(SVAs_15SVs_Index)){
  tryCatch({
    
    dat  <- counts(DDS_list[[i]], normalized = TRUE)
    idx  <- rowMeans(dat) > 1 
    dat  <- dat[idx, ] 
    if(DESeq_Results_Index$Comparison[i]=="ALSFTDvsHC"){
      mod  <- model.matrix(~ Case, colData(DDS_list[[i]]))
    }else{
      mod  <- model.matrix(~ Rand, colData(DDS_list[[i]]))
    }
    mod0 <- model.matrix(~ 1, colData(DDS_list[[i]]))
    
    svseq <- svaseq(dat, mod, mod0, n.sv=15)
    SVAs_15SVs[[i]] <- svseq
    
    SVAs_15SVs_Index$SVA[i] <- TRUE 
    SVAs_15SVs_Index$SVs[i] <- 15 
    message(paste0("SVAs #", i, " generated! "))
  }, error=function(e){
    
    SVAs_15SVs_Index$SVA[i] <<- FALSE 
    SVAs_15SVs_Index$SVs[i] <<- 0 
    message(paste0("SVAs #", i, " could not be generated :("))
  })
}


qsave(
  SVAs_15SVs, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD/SVA/1_15_SVs/", 
    "SVAs_15SVs", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  SVAs_15SVs_Index, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD/SVA/1_15_SVs/", 
    "SVAs_15SVs_Index", 
    ".qrds"
  ), 
  nthr=nthr
)


all(SVAs_15SVs_Index$SVs[SVAs_15SVs_Index$SVA]==15)
rm(dat, mod, mod0, svseq, i, idx)



    ## 3.3 Add 1-15 SVs to DDS_list --------------------------------------------

DESeq_SVA_1_15_Results_Index <- DESeq_Results_Index
rm(DESeq_Results_Index)


for (i in 1:15){
  
  DDS_SVA_list <- list() 
  for (j in 1:length(DDS_list)){
    DDS_SVA_list[[j]] <- NA
  } 
  
  DESeq_SVA_1_15_Results_Index[[paste0("SVA",i)]]  <- NA
  
  
  for (k in 1:nrow(DESeq_SVA_1_15_Results_Index)){ 
    
    tryCatch({ 
      
      DDS.tmp <- DDS_list[[k]] 
      
      for (l in 1:i){
        
        DDS.tmp[[paste0("SV", l)]] <- SVAs_15SVs[[k]]$sv[,l]
      }
      
      
      if(DESeq_SVA_1_15_Results_Index$Comparison[k]=="ALSFTDvsHC"){
        if(i==1){
          design(DDS.tmp) <- ~ SV1 + Case 
        }
        if(i==2){
          design(DDS.tmp) <- ~ SV1 + SV2 + Case 
        }
        if(i==3){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + Case 
        }
        if(i==4){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + Case 
        }
        if(i==5){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + Case 
        }
        if(i==6){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + Case 
        }
        if(i==7){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + Case 
        }
        if(i==8){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + Case 
        }
        if(i==9){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + Case 
        }
        if(i==10){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + Case 
        }
        if(i==11){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + Case 
        }
        if(i==12){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + Case 
        }
        if(i==13){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + Case 
        }
        if(i==14){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + SV14 + Case 
        }
        if(i==15){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + SV14 + SV15 + Case 
        }
        
      }else{
        if(i==1){
          design(DDS.tmp) <- ~ SV1 + Rand 
        }
        if(i==2){
          design(DDS.tmp) <- ~ SV1 + SV2 + Rand 
        }
        if(i==3){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + Rand 
        }
        if(i==4){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + Rand 
        }
        if(i==5){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + Rand 
        }
        if(i==6){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + Rand 
        }
        if(i==7){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + Rand 
        }
        if(i==8){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + Rand 
        }
        if(i==9){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + Rand 
        }
        if(i==10){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + Rand 
        }
        if(i==11){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + Rand 
        }
        if(i==12){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + Rand 
        }
        if(i==13){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + Rand 
        }
        if(i==14){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + SV14 + Rand 
        }
        if(i==15){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + SV14 + SV15 + Rand 
        }
        
      }
      
      DDS_SVA_list[[k]] <- DDS.tmp 
      rm(DDS.tmp)
      
      DESeq_SVA_1_15_Results_Index[[paste0("SVA",i)]][k] <- TRUE 
      message(k)
      message("########################################################")
      message(Sys.time())
    }, error=function(e){
      DESeq_SVA_1_15_Results_Index[[paste0("SVA",i)]][k] <<- FALSE 
      message("########################################################")
      message(Sys.time()) 
    })
  }  
  
  qsave(
    DDS_SVA_list, 
    paste0(
      "../Data/De/WNN/Pseudobulk/DESeq2/ALSFTD/SVA/1_15_SVs/",
      "DDS_", 
      i, 
      "SVs_", 
      "list", 
      ".qrds"
    ), 
    nthr=nthr
  )
  rm(DDS_SVA_list)
}

##$$$----
      # START HPC Cluster Outsourcing ------------------------------------------


      # Y_B4_3_DE_RNA_WNN_SVA_AllCase_1SVs.R -----------------------------------


      # Import HPC Output and summarize results --------------------------------




  ### 4.0 ALSFTD_C9 ------------------------------------------------------------



    ## 4.1 Load Data ----------------------------------------------------------- 

DDS_list <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9/", 
    "DDS_DESEqed_list", 
    ".qrds"
  ), 
  nthr=nthr
)


DESeq_Results_Index <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9/", 
    "DESeq_Results_Index", 
    ".qrds"
  ), 
  nthr=nthr
)    

DDS_list <- DDS_list[1:76] 
DESeq_Results_Index <- DESeq_Results_Index[1:76,1:9]



    ## 2.2 Find 15 SVs --------------------------------------------------------- 

SVAs_15SVs <- list() 
for (i in 1:nrow(DESeq_Results_Index)){
  SVAs_15SVs[[i]] <- NA
}

SVAs_15SVs_Index <- data.frame(
  Index=1:nrow(DESeq_Results_Index), 
  SVA=NA, 
  SVs = NA 
)

SVAs_15SVs_Index$CellTypeLevel <- DESeq_Results_Index$CellTypeLevel
SVAs_15SVs_Index$CellType <- DESeq_Results_Index$CellType
SVAs_15SVs_Index$Comparison <- DESeq_Results_Index$Comparison


SVAs_15SVs_Index$SVs <- NA 

for (i in 1:nrow(SVAs_15SVs_Index)){
  tryCatch({
    
    dat  <- counts(DDS_list[[i]], normalized = TRUE)
    idx  <- rowMeans(dat) > 1 
    dat  <- dat[idx, ] 
    if(DESeq_Results_Index$Comparison[i]=="C9_ALSFTDvsHC"){
      mod  <- model.matrix(~ Case, colData(DDS_list[[i]]))
    }else{
      mod  <- model.matrix(~ Rand, colData(DDS_list[[i]]))
    }
    mod0 <- model.matrix(~ 1, colData(DDS_list[[i]]))
    
    svseq <- svaseq(dat, mod, mod0, n.sv=15)
    SVAs_15SVs[[i]] <- svseq
    
    SVAs_15SVs_Index$SVA[i] <- TRUE 
    SVAs_15SVs_Index$SVs[i] <- 15 
    message(paste0("SVAs #", i, " generated! "))
  }, error=function(e){
    
    SVAs_15SVs_Index$SVA[i] <<- FALSE 
    SVAs_15SVs_Index$SVs[i] <<- 0 
    message(paste0("SVAs #", i, " could not be generated :("))
  })
}


qsave(
  SVAs_15SVs, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9/SVA/1_15_SVs/", 
    "SVAs_15SVs", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  SVAs_15SVs_Index, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9/SVA/1_15_SVs/", 
    "SVAs_15SVs_Index", 
    ".qrds"
  ), 
  nthr=nthr
)


all(SVAs_15SVs_Index$SVs[SVAs_15SVs_Index$SVA]==15)
rm(dat, mod, mod0, svseq, i, idx)



    ## 4.3 Add 1-15 SVs to DDS_list --------------------------------------------

DESeq_SVA_1_15_Results_Index <- DESeq_Results_Index
rm(DESeq_Results_Index)


for (i in 1:15){
  
  DDS_SVA_list <- list() 
  for (j in 1:length(DDS_list)){
    DDS_SVA_list[[j]] <- NA
  } 
  
  DESeq_SVA_1_15_Results_Index[[paste0("SVA",i)]]  <- NA
  
  
  for (k in 1:nrow(DESeq_SVA_1_15_Results_Index)){ 
    
    tryCatch({ 
      
      DDS.tmp <- DDS_list[[k]] 
      
      for (l in 1:i){
        
        DDS.tmp[[paste0("SV", l)]] <- SVAs_15SVs[[k]]$sv[,l]
      }
      
      
      if(DESeq_SVA_1_15_Results_Index$Comparison[k]=="C9_ALSFTDvsHC"){
        if(i==1){
          design(DDS.tmp) <- ~ SV1 + Case 
        }
        if(i==2){
          design(DDS.tmp) <- ~ SV1 + SV2 + Case 
        }
        if(i==3){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + Case 
        }
        if(i==4){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + Case 
        }
        if(i==5){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + Case 
        }
        if(i==6){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + Case 
        }
        if(i==7){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + Case 
        }
        if(i==8){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + Case 
        }
        if(i==9){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + Case 
        }
        if(i==10){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + Case 
        }
        if(i==11){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + Case 
        }
        if(i==12){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + Case 
        }
        if(i==13){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + Case 
        }
        if(i==14){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + SV14 + Case 
        }
        if(i==15){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + SV14 + SV15 + Case 
        }
        
      }else{
        if(i==1){
          design(DDS.tmp) <- ~ SV1 + Rand 
        }
        if(i==2){
          design(DDS.tmp) <- ~ SV1 + SV2 + Rand 
        }
        if(i==3){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + Rand 
        }
        if(i==4){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + Rand 
        }
        if(i==5){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + Rand 
        }
        if(i==6){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + Rand 
        }
        if(i==7){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + Rand 
        }
        if(i==8){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + Rand 
        }
        if(i==9){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + Rand 
        }
        if(i==10){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + Rand 
        }
        if(i==11){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + Rand 
        }
        if(i==12){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + Rand 
        }
        if(i==13){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + Rand 
        }
        if(i==14){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + SV14 + Rand 
        }
        if(i==15){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + SV14 + SV15 + Rand 
        }
        
      }
      
      DDS_SVA_list[[k]] <- DDS.tmp 
      rm(DDS.tmp)
      
      DESeq_SVA_1_15_Results_Index[[paste0("SVA",i)]][k] <- TRUE 
      message(k)
      message("########################################################")
      message(Sys.time())
    }, error=function(e){
      DESeq_SVA_1_15_Results_Index[[paste0("SVA",i)]][k] <<- FALSE 
      message("########################################################")
      message(Sys.time()) 
    })
  }  
  
  qsave(
    DDS_SVA_list, 
    paste0(
      "../Data/De/WNN/Pseudobulk/DESeq2/ALSFTD_C9/SVA/1_15_SVs/",
      "DDS_", 
      i, 
      "SVs_", 
      "list", 
      ".qrds"
    ), 
    nthr=nthr
  )
  rm(DDS_SVA_list)
}

##$$$----
      # START HPC Cluster Outsourcing ------------------------------------------


      # Y_B4_4_DE_RNA_WNN_SVA_AllCase_1SVs.R -----------------------------------


      # Import HPC Output and summarize results --------------------------------





  ### 5.0 ALSFTD_NonC9 ---------------------------------------------------------



    ## 5.1 Load Data ----------------------------------------------------------- 

DDS_list <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_NonC9/", 
    "DDS_DESEqed_list", 
    ".qrds"
  ), 
  nthr=nthr
)


DESeq_Results_Index <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_NonC9/", 
    "DESeq_Results_Index", 
    ".qrds"
  ), 
  nthr=nthr
)    

DDS_list <- DDS_list[1:76] 
DESeq_Results_Index <- DESeq_Results_Index[1:76,1:9]



    ## 5.2 Find 15 SVs ---------------------------------------------------------



SVAs_15SVs <- list() 
for (i in 1:nrow(DESeq_Results_Index)){
  SVAs_15SVs[[i]] <- NA
}

SVAs_15SVs_Index <- data.frame(
  Index=1:nrow(DESeq_Results_Index), 
  SVA=NA, 
  SVs = NA 
)

SVAs_15SVs_Index$CellTypeLevel <- DESeq_Results_Index$CellTypeLevel
SVAs_15SVs_Index$CellType <- DESeq_Results_Index$CellType
SVAs_15SVs_Index$Comparison <- DESeq_Results_Index$Comparison


SVAs_15SVs_Index$SVs <- NA 

for (i in 1:nrow(SVAs_15SVs_Index)){
  tryCatch({
    
    dat  <- counts(DDS_list[[i]], normalized = TRUE)
    idx  <- rowMeans(dat) > 1 
    dat  <- dat[idx, ] 
    if(DESeq_Results_Index$Comparison[i]=="NonC9_ALSFTDvsHC"){
      mod  <- model.matrix(~ Case, colData(DDS_list[[i]]))
    }else{
      mod  <- model.matrix(~ Rand, colData(DDS_list[[i]]))
    }
    mod0 <- model.matrix(~ 1, colData(DDS_list[[i]]))
    
    svseq <- svaseq(dat, mod, mod0, n.sv=15)
    SVAs_15SVs[[i]] <- svseq
    
    SVAs_15SVs_Index$SVA[i] <- TRUE 
    SVAs_15SVs_Index$SVs[i] <- 15 
    message(paste0("SVAs #", i, " generated! "))
  }, error=function(e){
    
    SVAs_15SVs_Index$SVA[i] <<- FALSE 
    SVAs_15SVs_Index$SVs[i] <<- 0 
    message(paste0("SVAs #", i, " could not be generated :("))
  })
}


qsave(
  SVAs_15SVs, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_NonC9/SVA/1_15_SVs/", 
    "SVAs_15SVs", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  SVAs_15SVs_Index, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_NonC9/SVA/1_15_SVs/", 
    "SVAs_15SVs_Index", 
    ".qrds"
  ), 
  nthr=nthr
)


all(SVAs_15SVs_Index$SVs[SVAs_15SVs_Index$SVA]==15)
rm(dat, mod, mod0, svseq, i, idx)



    ## 5.3 Add 1-15 SVs to DDS_list --------------------------------------------

DESeq_SVA_1_15_Results_Index <- DESeq_Results_Index
rm(DESeq_Results_Index)


for (i in 1:15){
  
  DDS_SVA_list <- list() 
  for (j in 1:length(DDS_list)){
    DDS_SVA_list[[j]] <- NA
  } 
  
  DESeq_SVA_1_15_Results_Index[[paste0("SVA",i)]]  <- NA
  
  
  for (k in 1:nrow(DESeq_SVA_1_15_Results_Index)){ 
    
    tryCatch({ 
      
      DDS.tmp <- DDS_list[[k]] 
      
      for (l in 1:i){
        
        DDS.tmp[[paste0("SV", l)]] <- SVAs_15SVs[[k]]$sv[,l]
      }
      
      
      if(DESeq_SVA_1_15_Results_Index$Comparison[k]=="NonC9_ALSFTDvsHC"){
        if(i==1){
          design(DDS.tmp) <- ~ SV1 + Case 
        }
        if(i==2){
          design(DDS.tmp) <- ~ SV1 + SV2 + Case 
        }
        if(i==3){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + Case 
        }
        if(i==4){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + Case 
        }
        if(i==5){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + Case 
        }
        if(i==6){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + Case 
        }
        if(i==7){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + Case 
        }
        if(i==8){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + Case 
        }
        if(i==9){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + Case 
        }
        if(i==10){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + Case 
        }
        if(i==11){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + Case 
        }
        if(i==12){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + Case 
        }
        if(i==13){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + Case 
        }
        if(i==14){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + SV14 + Case 
        }
        if(i==15){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + SV14 + SV15 + Case 
        }
        
      }else{
        if(i==1){
          design(DDS.tmp) <- ~ SV1 + Rand 
        }
        if(i==2){
          design(DDS.tmp) <- ~ SV1 + SV2 + Rand 
        }
        if(i==3){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + Rand 
        }
        if(i==4){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + Rand 
        }
        if(i==5){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + Rand 
        }
        if(i==6){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + Rand 
        }
        if(i==7){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + Rand 
        }
        if(i==8){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + Rand 
        }
        if(i==9){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + Rand 
        }
        if(i==10){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + Rand 
        }
        if(i==11){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + Rand 
        }
        if(i==12){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + Rand 
        }
        if(i==13){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + Rand 
        }
        if(i==14){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + SV14 + Rand 
        }
        if(i==15){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + SV14 + SV15 + Rand 
        }
        
      }
      
      DDS_SVA_list[[k]] <- DDS.tmp 
      rm(DDS.tmp)
      
      DESeq_SVA_1_15_Results_Index[[paste0("SVA",i)]][k] <- TRUE 
      message(k)
      message("########################################################")
      message(Sys.time())
    }, error=function(e){
      DESeq_SVA_1_15_Results_Index[[paste0("SVA",i)]][k] <<- FALSE 
      message("########################################################")
      message(Sys.time()) 
    })
  }  
  
  qsave(
    DDS_SVA_list, 
    paste0(
      "../Data/De/WNN/Pseudobulk/DESeq2/ALSFTD_NonC9/SVA/1_15_SVs/",
      "DDS_", 
      i, 
      "SVs_", 
      "list", 
      ".qrds"
    ), 
    nthr=nthr
  )
  rm(DDS_SVA_list)
}

##$$$----
      # START HPC Cluster Outsourcing ------------------------------------------


      # Y_B4_5_DE_RNA_WNN_SVA_AllCase_1SVs.R -----------------------------------


      # Import HPC Output and summarize results --------------------------------




  ### 6.0 ALSFTD_C9_ALSFTD_NonC9 -----------------------------------------------



    ## 6.1 Load Data ----------------------------------------------------------- 

DDS_list <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9_NonC9/", 
    "DDS_DESEqed_list", 
    ".qrds"
  ), 
  nthr=nthr
)


DESeq_Results_Index <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9_NonC9/", 
    "DESeq_Results_Index", 
    ".qrds"
  ), 
  nthr=nthr
)    

DDS_list <- DDS_list[1:76] 
DESeq_Results_Index <- DESeq_Results_Index[1:76,1:9]



    ## 6.2 Find 14 SVs ---------------------------------------------------------

SVAs_14SVs <- list() 
for (i in 1:nrow(DESeq_Results_Index)){
  SVAs_14SVs[[i]] <- NA
}

SVAs_14SVs_Index <- data.frame(
  Index=1:nrow(DESeq_Results_Index), 
  SVA=NA, 
  SVs = NA 
)

SVAs_14SVs_Index$CellTypeLevel <- DESeq_Results_Index$CellTypeLevel
SVAs_14SVs_Index$CellType <- DESeq_Results_Index$CellType
SVAs_14SVs_Index$Comparison <- DESeq_Results_Index$Comparison


SVAs_14SVs_Index$SVs <- NA 

for (i in 1:nrow(SVAs_14SVs_Index)){
  tryCatch({
    
    dat  <- counts(DDS_list[[i]], normalized = TRUE)

    idx  <- rowMeans(dat) > 1 
    dat  <- dat[idx, ] 
    if(DESeq_Results_Index$Comparison[i]=="ALSFTD_C9vsALSFTD_NonC9"){
      colData(DDS_list[[i]])$Case_Type <- as.factor(colData(DDS_list[[i]])$Case_Type)
      mod  <- model.matrix(~ Case_Type, colData(DDS_list[[i]]))
    }else{
      mod  <- model.matrix(~ Rand, colData(DDS_list[[i]]))
    }
    mod0 <- model.matrix(~ 1, colData(DDS_list[[i]]))
    
    svseq <- svaseq(dat, mod, mod0, n.sv=14)
    SVAs_14SVs[[i]] <- svseq
    
    SVAs_14SVs_Index$SVA[i] <- TRUE 
    SVAs_14SVs_Index$SVs[i] <- 14 
    message(paste0("SVAs #", i, " generated! "))
  }, error=function(e){
    
    SVAs_14SVs_Index$SVA[i] <<- FALSE 
    SVAs_14SVs_Index$SVs[i] <<- 0 
    message(paste0("SVAs #", i, " could not be generated :("))
  })
}


qsave(
  SVAs_14SVs, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9_NonC9/SVA/1_14_SVs/", 
    "SVAs_14SVs", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  SVAs_14SVs_Index, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9_NonC9/SVA/1_14_SVs/", 
    "SVAs_14SVs_Index", 
    ".qrds"
  ), 
  nthr=nthr
)


all(SVAs_14SVs_Index$SVs[SVAs_14SVs_Index$SVA]==14)
rm(dat, mod, mod0, svseq, i, idx)



    ## 6.3 Add 1-14 SVs to DDS_list --------------------------------------------

DESeq_SVA_1_14_Results_Index <- DESeq_Results_Index
rm(DESeq_Results_Index)


for (i in 1:14){
  
  DDS_SVA_list <- list() 
  for (j in 1:length(DDS_list)){
    DDS_SVA_list[[j]] <- NA
  } 
  
  DESeq_SVA_1_14_Results_Index[[paste0("SVA",i)]]  <- NA
  
  
  for (k in 1:nrow(DESeq_SVA_1_14_Results_Index)){ 
    
    tryCatch({ 
      
      DDS.tmp <- DDS_list[[k]] 
      
      for (l in 1:i){
        
        DDS.tmp[[paste0("SV", l)]] <- SVAs_14SVs[[k]]$sv[,l]
      }
      
      
      if(DESeq_SVA_1_14_Results_Index$Comparison[k]=="ALSFTD_C9vsALSFTD_NonC9"){
        if(i==1){
          design(DDS.tmp) <- ~ SV1 + Case_Type 
        }
        if(i==2){
          design(DDS.tmp) <- ~ SV1 + SV2 + Case_Type 
        }
        if(i==3){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + Case_Type 
        }
        if(i==4){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + Case_Type 
        }
        if(i==5){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + Case_Type 
        }
        if(i==6){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + Case_Type 
        }
        if(i==7){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + Case_Type 
        }
        if(i==8){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + Case_Type 
        }
        if(i==9){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + Case_Type 
        }
        if(i==10){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + Case_Type 
        }
        if(i==11){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + Case_Type 
        }
        if(i==12){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + Case_Type 
        }
        if(i==13){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + Case_Type 
        }
        if(i==14){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + SV14 + Case_Type 
        }
        
      }else{
        if(i==1){
          design(DDS.tmp) <- ~ SV1 + Rand 
        }
        if(i==2){
          design(DDS.tmp) <- ~ SV1 + SV2 + Rand 
        }
        if(i==3){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + Rand 
        }
        if(i==4){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + Rand 
        }
        if(i==5){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + Rand 
        }
        if(i==6){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + Rand 
        }
        if(i==7){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + Rand 
        }
        if(i==8){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + Rand 
        }
        if(i==9){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + Rand 
        }
        if(i==10){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + Rand 
        }
        if(i==11){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + Rand 
        }
        if(i==12){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + Rand 
        }
        if(i==13){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + Rand 
        }
        if(i==14){
          design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + SV13 + SV14 + Rand 
        }
        
      }
      
      DDS_SVA_list[[k]] <- DDS.tmp 
      rm(DDS.tmp)
      
      DESeq_SVA_1_14_Results_Index[[paste0("SVA",i)]][k] <- TRUE 
      message(k)
      message("########################################################")
      message(Sys.time())
    }, error=function(e){
      DESeq_SVA_1_14_Results_Index[[paste0("SVA",i)]][k] <<- FALSE 
      message("########################################################")
      message(Sys.time()) 
    })
  }  
  
  qsave(
    DDS_SVA_list, 
    paste0(
      "../Data/De/WNN/Pseudobulk/DESeq2/ALSFTD_C9_NonC9/SVA/1_14_SVs/",
      "DDS_", 
      i, 
      "SVs_", 
      "list", 
      ".qrds"
    ), 
    nthr=nthr
  )
  rm(DDS_SVA_list)
}

##$$$----
      # START HPC Cluster Outsourcing ------------------------------------------


      # Y_B4_6_DE_RNA_WNN_SVA_ALSFTD_C9_NonC9_14SVs.R ---


      # Import HPC Output and summarize results --------------------------------








