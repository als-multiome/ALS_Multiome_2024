
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



    ## 1.2 Find 12 SVs --------------------------------------------------------- 

SVAs_12SVs <- list() 
for (i in 1:nrow(DESeq_Results_Index)){
  SVAs_12SVs[[i]] <- NA
}

SVAs_12SVs_Index <- data.frame(
  Index=1:nrow(DESeq_Results_Index), 
  SVA=NA, 
  SVs = 12 
)

SVAs_12SVs_Index$CellTypeLevel <- DESeq_Results_Index$CellTypeLevel
SVAs_12SVs_Index$CellType <- DESeq_Results_Index$CellType
SVAs_12SVs_Index$Comparison <- DESeq_Results_Index$Comparison


SVAs_12SVs_Index$SVs <- NA 

for (i in 1:nrow(SVAs_12SVs_Index)){
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
    
    svseq <- svaseq(dat, mod, mod0, n.sv=12)
    SVAs_12SVs[[i]] <- svseq
    
    SVAs_12SVs_Index$SVA[i] <- TRUE 
    SVAs_12SVs_Index$SVs[i] <- 12 
    message(paste0("SVAs #", i, " generated! "))
  }, error=function(e){
    
    SVAs_12SVs_Index$SVA[i] <<- FALSE 
    SVAs_12SVs_Index$SVs[i] <<- 0 
    message(paste0("SVAs #", i, " could not be generated :("))
  })
}


qsave(
  SVAs_12SVs, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/12_SVs/", 
    "SVAs_12SVs", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  SVAs_12SVs_Index, 
  paste0(
    "../Data/DE/WNN/AllCase/SVA/12_SVs/", 
    "SVAs_12SVs_Index", 
    ".qrds"
  ), 
  nthr=nthr
)

all(SVAs_12SVs_Index$SVs[SVAs_12SVs_Index$SVA]==12)
rm(dat, mod, mod0, svseq, i, idx)



    ## 1.3 DESeq2 with 12 SVs --------------------------------------------------

register(MulticoreParam(nthr)) 

DDS_SVA_12sv_list <- list() 
for (i in 1:length(DDS_list)){
  DDS_SVA_12sv_list[[i]] <- NA
}

DESeq_SVA_12sv_Results_Index <- DESeq_Results_Index[,c(1:9)] 
rm(DESeq_Results_Index)
message(
  paste0(
    "Length of DDS_SVA_12sv_list same as Nrow DESeq_SVA_12sv_Results_Index: ", 
    length(DDS_SVA_12sv_list)==nrow(DESeq_SVA_12sv_Results_Index)
  )
)
DESeq_SVA_12sv_Results_Index$SVA <- NA 

for (i in 1:nrow(DESeq_SVA_12sv_Results_Index)){ 
  tryCatch({
    DDS.tmp <- DDS_list[[i]] 
    
    DDS.tmp$SV1 <- SVAs_12SVs[[i]]$sv[,1]
    DDS.tmp$SV2 <- SVAs_12SVs[[i]]$sv[,2]
    DDS.tmp$SV3 <- SVAs_12SVs[[i]]$sv[,3]
    DDS.tmp$SV4 <- SVAs_12SVs[[i]]$sv[,4] 
    DDS.tmp$SV5 <- SVAs_12SVs[[i]]$sv[,5]
    DDS.tmp$SV6 <- SVAs_12SVs[[i]]$sv[,6]
    DDS.tmp$SV7 <- SVAs_12SVs[[i]]$sv[,7]
    DDS.tmp$SV8 <- SVAs_12SVs[[i]]$sv[,8]
    DDS.tmp$SV9 <- SVAs_12SVs[[i]]$sv[,9] 
    DDS.tmp$SV10 <- SVAs_12SVs[[i]]$sv[,10]
    DDS.tmp$SV11 <- SVAs_12SVs[[i]]$sv[,11]
    DDS.tmp$SV12 <- SVAs_12SVs[[i]]$sv[,12]
    
    if(DESeq_SVA_12sv_Results_Index$Comparison[i]=="All_Cases"){
      design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + Case
    }else{
      design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + Rand 
    }
    DDS_SVA_12sv_list[[i]] <- DESeq(DDS.tmp, test="Wald", quiet=FALSE, parallel=TRUE) 
    DESeq_SVA_12sv_Results_Index$SVA[i] <- TRUE 
    message("########################################################")
    message(Sys.time())
    message(paste0("DESeq analysis #", i, " sucessfull!...")) 
    message("########################################################")
  }, error=function(e){
    DESeq_SVA_12sv_Results_Index$SVA[i] <<- FALSE 
    message("########################################################")
    message(Sys.time()) 
    message(paste0("DESeq analysis #", i, " failed :( "))
    message("########################################################")
  })
  tryCatch({rm(DDS.tmp)}, error=function(e){})
  tryCatch({
    if(!exists("it.tmp")) it.tmp <- 0 
    if(!exists("Systimes.tmp")) Systimes.tmp <- as.character(Sys.time())
    it.tmp <- c(it.tmp, i) 
    Systimes.tmp <- as.character(c(Systimes.tmp, as.character(Sys.time())))
    write.csv(data.frame(It=it.tmp, Time=Systimes.tmp), "../Data/Benchmarks/Loop_State_B3_AllCase.csv", row.names = FALSE)}, error=function(e){})
}

qsave(
  DDS_SVA_12sv_list, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/12_SVs/", 
    "DDS_SVA_12sv_list", 
    ".qrds"
  ), 
  nthr=nthr
)

write.csv(resultsNames(DDS_SVA_12sv_list[[1]]), "./Loop done.csv")



    ## 1.4 DESeq SVA results summary -------------------------------------------

resultsNames(DDS_SVA_12sv_list[[1]])
resultsNames(DDS_SVA_12sv_list[[2]])

DESeq_Results_12SVs_ALS <- list()  
DESeq_Results_12SVs_ALSFTD <- list()  


DESeq_Results_12SVs_ALS <- lapply(
  DDS_SVA_12sv_list, 
  FUN=function(x){
    return(
      if(is(x, "DESeqDataSet")){
          if("Case_ALS_vs_HC" %in% resultsNames(x)){
            DESeq2::results(x, alpha=0.05, lfcThreshold = 0, name = "Case_ALS_vs_HC")
          }else{
            DESeq2::results(x, alpha=0.05, lfcThreshold = 0, name = "Rand_R2_vs_R1")
          }
      }else{
        NA
      }
    )
  }
)


DESeq_Results_12SVs_ALSFTD <- lapply(
  DDS_SVA_12sv_list, 
  FUN=function(x){
    return(
      if(is(x, "DESeqDataSet")){
        if("Case_ALS_FTD_vs_HC" %in% resultsNames(x)){
          DESeq2::results(x, alpha=0.05, lfcThreshold = 0, name = "Case_ALS_FTD_vs_HC")
        }else{
          DESeq2::results(x, alpha=0.05, lfcThreshold = 0, name = "Rand_R3_vs_R1")
        }
      }else{
        NA
      }
    )
  }
)

DESeq_Results_12SVs_Index$All_q_0.05_ALS <- NA 
DESeq_Results_12SVs_Index$All_q_0.05_ALSFTD <- NA 

DESeq_Results_12SVs_Index$Up_q_0.05_ALS <- NA   
DESeq_Results_12SVs_Index$Up_q_0.05_ALSFTD <- NA  

DESeq_Results_12SVs_Index$Down_q_0.05_ALS <- NA 
DESeq_Results_12SVs_Index$Down_q_0.05_FTD <- NA 

DESeq_Results_12SVs_Index$nSamples <- NA 
DESeq_Results_12SVs_Index$nSamples_HC <- NA  
DESeq_Results_12SVs_Index$nSamples_ALS <- NA 
DESeq_Results_12SVs_Index$nSamples_ALSFTD <- NA 

DESeq_Results_12SVs_Index$TotalReads <- NA 

DESeq_Results_12SVs_Index$MeanReads <- NA 

DESeq_Results_12SVs_Index$NonzeroCounts_ALS <- NA 
DESeq_Results_12SVs_Index$NonzeroCounts_ALSFTD <- NA 

DESeq_Results_12SVs_Index$LowCounts_ALS <- NA 
DESeq_Results_12SVs_Index$LowCounts_ALSFTD <- NA 

sample_data <- DDS_SVA_12sv_list[[1]]@colData
DESeq_Results_12SVs_Index$Covariates <- NA 

for (i in 1:length(DDS_SVA_12sv_list)){
  if(!DESeq_Results_12SVs_Index$SVA[i]){
    
  }else{
    DESeq_Results_12SVs_Index$Covariates[i] <- "12_SVs"
    DESeq_Results_12SVs_Index$NonzeroCounts_ALS[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_12SVs_ALS[[i]]))[2], fixed("of "))[[1]][2], " with")[[1]][1]) 
    DESeq_Results_12SVs_Index$NonzeroCounts_ALSFTD[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_12SVs_ALSFTD[[i]]))[2], fixed("of "))[[1]][2], " with")[[1]][1]) 
    
    DESeq_Results_12SVs_Index$Up_q_0.05_ALS[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_12SVs_ALS[[i]]))[4], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_12SVs_Index$Up_q_0.05_ALSFTD[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_12SVs_ALSFTD[[i]]))[4], fixed(": "))[[1]][2], ", ")[[1]][1]) 
    
    DESeq_Results_12SVs_Index$Down_q_0.05_ALS[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_12SVs_ALS[[i]]))[5], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_12SVs_Index$Down_q_0.05_ALSFTD[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_12SVs_ALSFTD[[i]]))[5], fixed(": "))[[1]][2], ", ")[[1]][1])
    
    DESeq_Results_12SVs_Index$LowCounts_ALS[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_12SVs_ALS[[i]]))[7], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_12SVs_Index$LowCounts_ALSFTD[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_12SVs_ALSFTD[[i]]))[7], fixed(": "))[[1]][2], ", ")[[1]][1])
    
    DESeq_Results_12SVs_Index$nSamples[i] <- if(is.na(DDS_SVA_12sv_list[i])){NA} else{dim(DDS_SVA_12sv_list[[i]]@assays@data$counts)[2]}
    
    
    DESeq_Results_12SVs_Index$nSamples_HC[i] <- if(is.na(DDS_SVA_12sv_list[i])){NA} else{table(sample_data$Case[sample_data$ID %in% colnames(DDS_SVA_12sv_list[[i]]@assays@data$counts)])[1]}
    DESeq_Results_12SVs_Index$nSamples_ALS[i] <- if(is.na(DDS_SVA_12sv_list[i])){NA} else{table(sample_data$Case[sample_data$ID %in% colnames(DDS_SVA_12sv_list[[i]]@assays@data$counts)])[2]}
    DESeq_Results_12SVs_Index$nSamples_ALSFTD[i] <- if(is.na(DDS_SVA_12sv_list[i])){NA} else{table(sample_data$Case[sample_data$ID %in% colnames(DDS_SVA_12sv_list[[i]]@assays@data$counts)])[3]}
    
    DESeq_Results_12SVs_Index$TotalReads[i] <- if(is.na(DDS_SVA_12sv_list[i])){NA} else{sum(colSums(DDS_SVA_12sv_list[[i]]@assays@data$counts))}
    DESeq_Results_12SVs_Index$MeanReads[i] <- if(is.na(DDS_SVA_12sv_list[i])){NA} else{sum(colMeans(DDS_SVA_12sv_list[[i]]@assays@data$counts))}
    DESeq_Results_12SVs_Index$nFeatures[i] <- if(is.na(DDS_SVA_12sv_list[i])){NA} else{sum(rowMeans(DDS_SVA_12sv_list[[i]]@assays@data$counts)>0)}
  }
  
}

DESeq_Results_12SVs_Index$All_q_0.05_ALS <- DESeq_Results_12SVs_Index$Up_q_0.05_ALS + DESeq_Results_12SVs_Index$Down_q_0.05_ALS 
DESeq_Results_12SVs_Index$All_q_0.05_ALSFTD <- DESeq_Results_12SVs_Index$Up_q_0.05_ALSFTD + DESeq_Results_12SVs_Index$Down_q_0.05_ALSFTD 



    ## 1.5 Save Data -----------------------------------------------------------

qsave(
  DESeq_Results_12SVs_ALS, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/12_SVs/", 
    "DESeq_Results_12SVs_ALS", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_12SVs_ALSFTD, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/12_SVs/", 
    "DESeq_Results_12SVs_ALSFTD", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_12SVs_Index, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/12_SVs/", 
    "DESeq_Results_12SVs_Index", 
    ".qrds"
  ), 
  nthr=nthr
)

rm(DDS_list, DDS_SVA_12sv_list, DESeq_Results_12SVs_ALS, DESeq_Results_12SVs_ALSFTD, DESeq_Results_12SVs_Index, sample_data, i, SVAs_12SVs, SVAs_12SVs_Index, it.tmp, Systimes.tmp) 




  ### 2.0 ALSvsHC -------------------------------------------------------------- 



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



    ## 2.2 Find 12 SVs --------------------------------------------------------- 

SVAs_12SVs <- list() 
for (i in 1:nrow(DESeq_Results_Index)){
  SVAs_12SVs[[i]] <- NA
}


SVAs_12SVs_Index <- data.frame(
  Index=1:nrow(DESeq_Results_Index), 
  SVA=NA, 
  SVs = 12 
)

SVAs_12SVs_Index$CellTypeLevel <- DESeq_Results_Index$CellTypeLevel
SVAs_12SVs_Index$CellType <- DESeq_Results_Index$CellType
SVAs_12SVs_Index$Comparison <- DESeq_Results_Index$Comparison

SVAs_12SVs_Index$SVs <- NA  

for (i in 1:nrow(SVAs_12SVs_Index)){
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
    
    svseq <- svaseq(dat, mod, mod0, n.sv=12)
    SVAs_12SVs[[i]] <- svseq
    
    SVAs_12SVs_Index$SVA[i] <- TRUE 
    SVAs_12SVs_Index$SVs[i] <- 12 
    message(paste0("SVAs #", i, " generated! "))
  }, error=function(e){
    
    SVAs_12SVs_Index$SVA[i] <<- FALSE 
    SVAs_12SVs_Index$SVs[i] <<- 0 
    message(paste0("SVAs #", i, " could not be generated :("))
  })
}


qsave(
  SVAs_12SVs, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALS/SVA/12_SVs/", 
    "SVAs_12SVs", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  SVAs_12SVs_Index, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALS/SVA/12_SVs/", 
    "SVAs_12SVs_Index", 
    ".qrds"
  ), 
  nthr=nthr
)

all(SVAs_12SVs_Index$SVs[SVAs_12SVs_Index$SVA]==12)
rm(dat, mod, mod0, svseq, i, idx)



    ## 2.3 DESeq2 with 12 SVs --------------------------------------------------

register(MulticoreParam(nthr)) 

DDS_SVA_12sv_list <- list() 
for (i in 1:length(DDS_list)){
  DDS_SVA_12sv_list[[i]] <- NA
}

DESeq_SVA_12sv_Results_Index <- DESeq_Results_Index[,c(1:9)] 
rm(DESeq_Results_Index)

message(
  paste0(
    "Length of DDS_SVA_12sv_list same as Nrow DESeq_SVA_12sv_Results_Index: ", 
    length(DDS_SVA_12sv_list)==nrow(DESeq_SVA_12sv_Results_Index)
  )
)
DESeq_SVA_12sv_Results_Index$SVA <- NA 

for (i in 1:nrow(DESeq_SVA_12sv_Results_Index)){ 
  tryCatch({
    DDS.tmp <- DDS_list[[i]] 
    
    DDS.tmp$SV1 <- SVAs_12SVs[[i]]$sv[,1]
    DDS.tmp$SV2 <- SVAs_12SVs[[i]]$sv[,2]
    DDS.tmp$SV3 <- SVAs_12SVs[[i]]$sv[,3]
    DDS.tmp$SV4 <- SVAs_12SVs[[i]]$sv[,4] 
    DDS.tmp$SV5 <- SVAs_12SVs[[i]]$sv[,5]
    DDS.tmp$SV6 <- SVAs_12SVs[[i]]$sv[,6]
    DDS.tmp$SV7 <- SVAs_12SVs[[i]]$sv[,7]
    DDS.tmp$SV8 <- SVAs_12SVs[[i]]$sv[,8]
    DDS.tmp$SV9 <- SVAs_12SVs[[i]]$sv[,9] 
    DDS.tmp$SV10 <- SVAs_12SVs[[i]]$sv[,10]
    DDS.tmp$SV11 <- SVAs_12SVs[[i]]$sv[,11]
    DDS.tmp$SV12 <- SVAs_12SVs[[i]]$sv[,12]
    
    if(DESeq_SVA_12sv_Results_Index$Comparison[i]=="ALSvsHC"){
      design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + Case
    }else{
      design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + Rand 
    }
    DDS_SVA_12sv_list[[i]] <- DESeq(DDS.tmp, test="Wald", quiet=FALSE, parallel=TRUE) 
    DESeq_SVA_12sv_Results_Index$SVA[i] <- TRUE 
    message("########################################################")
    message(Sys.time())
    message(paste0("DESeq analysis #", i, " sucessfull!...")) 
    message("########################################################")
  }, error=function(e){
    DESeq_SVA_12sv_Results_Index$SVA[i] <<- FALSE 
    message("########################################################")
    message(Sys.time()) 
    message(paste0("DESeq analysis #", i, " failed :( "))
    message("########################################################")
  })
  tryCatch({rm(DDS.tmp)}, error=function(e){})
  tryCatch({
    if(!exists("it.tmp")) it.tmp <- 0 
    if(!exists("Systimes.tmp")) Systimes.tmp <- as.character(Sys.time())
    it.tmp <- c(it.tmp, i) 
    Systimes.tmp <- as.character(c(Systimes.tmp, as.character(Sys.time())))
    write.csv(data.frame(It=it.tmp, Time=Systimes.tmp), "../Data/Benchmarks/Loop_State_B3_ALSvsHC.csv", row.names = FALSE)}, error=function(e){})
}

qsave(
  DDS_SVA_12sv_list, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALS/SVA/12_SVs/", 
    "DDS_SVA_12sv_list", 
    ".qrds"
  ), 
  nthr=nthr
)



    ## 2.4 DESeq SVA results summary -------------------------------------------

resultsNames(DDS_SVA_12sv_list[[1]])
resultsNames(DDS_SVA_12sv_list[[2]])

DESeq_Results_12SVs <- list() 
DESeq_Results_12SVs <- lapply(
  DDS_SVA_12sv_list, 
  FUN=function(x){
    return(
      if(is(x, "DESeqDataSet")){
        results(x, alpha=0.05, lfcThreshold = 0)
      }else{
        NA
      }
    )
  }
)


DESeq_Results_12SVs_Index <- DESeq_SVA_12sv_Results_Index 
rm(DESeq_SVA_12sv_Results_Index)

DESeq_Results_12SVs_Index$Covariates <- NA
DESeq_Results_12SVs_Index$Covariates[DESeq_Results_12SVs_Index$SVA] <- "12_SVs"
DESeq_Results_12SVs_Index$SVA[is.na(DESeq_Results_12SVs_Index$SVA)] <- FALSE
DESeq_Results_12SVs_Index$All_q_0.05 <- NA 
DESeq_Results_12SVs_Index$Up_q_0.05 <- NA  
DESeq_Results_12SVs_Index$Down_q_0.05 <- NA 
DESeq_Results_12SVs_Index$nSamples <- NA 
DESeq_Results_12SVs_Index$nSamples_HC <- NA  
DESeq_Results_12SVs_Index$nSamples_ALS <- NA 
DESeq_Results_12SVs_Index$TotalReads <- NA 
DESeq_Results_12SVs_Index$MeanReads <- NA 
DESeq_Results_12SVs_Index$NonzeroCounts <- NA 
DESeq_Results_12SVs_Index$nFeatures <- NA 
DESeq_Results_12SVs_Index$LowCounts <- NA 

sample_data <- colData(DDS_SVA_12sv_list[[1]])

for (i in 1:length(DDS_SVA_12sv_list)){
  if(!DESeq_Results_12SVs_Index$SVA[i]){
    
  }else{
    DESeq_Results_12SVs_Index$NonzeroCounts[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_12SVs[[i]]))[2], fixed("of "))[[1]][2], " with")[[1]][1])
    DESeq_Results_12SVs_Index$Up_q_0.05[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_12SVs[[i]]))[4], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_12SVs_Index$Down_q_0.05[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_12SVs[[i]]))[5], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_12SVs_Index$LowCounts[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_12SVs[[i]]))[7], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_12SVs_Index$nSamples[i] <- dim(DDS_SVA_12sv_list[[i]]@assays@data$counts)[2] 
    DESeq_Results_12SVs_Index$nSamples_HC[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_SVA_12sv_list[[i]]@assays@data$counts)])[1]
    DESeq_Results_12SVs_Index$nSamples_ALS[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_SVA_12sv_list[[i]]@assays@data$counts)])[2]
    DESeq_Results_12SVs_Index$TotalReads[i] <- sum(colSums(DDS_SVA_12sv_list[[i]]@assays@data$counts)) 
    DESeq_Results_12SVs_Index$MeanReads[i] <- sum(colMeans(DDS_SVA_12sv_list[[i]]@assays@data$counts)) 
    DESeq_Results_12SVs_Index$nFeatures[i] <- sum(rowMeans(DDS_SVA_12sv_list[[i]]@assays@data$counts)>0) 
  }
  
}

DESeq_Results_12SVs_Index$All_q_0.05 <- DESeq_Results_12SVs_Index$Up_q_0.05 + DESeq_Results_12SVs_Index$Down_q_0.05 


    ## 2.5 Save Data -----------------------------------------------------------

qsave(
  DESeq_Results_12SVs, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALS/SVA/12_SVs/", 
    "DESeq_Results_12SVs", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_12SVs_Index, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALS/SVA/12_SVs/", 
    "DESeq_Results_12SVs_Index", 
    ".qrds"
  ), 
  nthr=nthr
)

rm(DDS_list, DDS_SVA_12sv_list, DESeq_Results_12SVs, DESeq_Results_12SVs_Index, sample_data, i) 




  ### 3.0 ALSFTDvsHC -----------------------------------------------------------



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



    ## 3.2 Find 12 SVs --------------------------------------------------------- 

SVAs_12SVs <- list() 
for (i in 1:nrow(DESeq_Results_Index)){
  SVAs_12SVs[[i]] <- NA
}

SVAs_12SVs_Index <- data.frame(
  Index=1:nrow(DESeq_Results_Index), 
  SVA=NA, 
  SVs = 12 
)

SVAs_12SVs_Index$CellTypeLevel <- DESeq_Results_Index$CellTypeLevel
SVAs_12SVs_Index$CellType <- DESeq_Results_Index$CellType
SVAs_12SVs_Index$Comparison <- DESeq_Results_Index$Comparison


DESeq_Results_Index$SVs <- NA 

for (i in 1:nrow(DESeq_Results_Index)){
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
    
    svseq <- svaseq(dat, mod, mod0, n.sv=12)
    DDS_list[[i]]$SV1 <- svseq$sv[,1]
    DDS_list[[i]]$SV2 <- svseq$sv[,2]
    DDS_list[[i]]$SV3 <- svseq$sv[,3]
    DDS_list[[i]]$SV4 <- svseq$sv[,4] 
    DDS_list[[i]]$SV5 <- svseq$sv[,5]
    DDS_list[[i]]$SV6 <- svseq$sv[,6]
    DDS_list[[i]]$SV7 <- svseq$sv[,7]
    DDS_list[[i]]$SV8 <- svseq$sv[,8]
    DDS_list[[i]]$SV9 <- svseq$sv[,9] 
    DDS_list[[i]]$SV10 <- svseq$sv[,10]
    DDS_list[[i]]$SV11 <- svseq$sv[,11]
    DDS_list[[i]]$SV12 <- svseq$sv[,12]
    DESeq_Results_Index$SVs[i] <- 12 
    SVAs_12SVs_Index$SVA[i] <- TRUE 
    SVAs_12SVs_Index$SVs[i] <- 12 
    message(paste0("SVAs #", i, " generated! "))
  }, error=function(e){
    DESeq_Results_Index$SVs[i] <<- 0 
    SVAs_12SVs_Index$SVA[i] <<- FALSE 
    SVAs_12SVs_Index$SVs[i] <<- 0 
    message(paste0("SVAs #", i, " could not be generated :("))
  })
}


qsave(
  SVAs_12SVs, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD/SVA/12_SVs/", 
    "SVAs_12SVs", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  SVAs_12SVs_Index, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD/SVA/12_SVs/", 
    "SVAs_12SVs_Index", 
    ".qrds"
  ), 
  nthr=nthr
)

all(SVAs_12SVs_Index$SVs[SVAs_12SVs_Index$SVA]==12)
rm(dat, mod, mod0, svseq, i, idx)



    ## 3.3 DESeq2 with 12 SVs --------------------------------------------------


      # Start HPC Cluster Outsourcing ------------------------------------------


      # Export environment -----------------------------------------------------

do.call(
  qsavem, 
  c(
    lapply(
      list(
        "DDS_list", 
        "DESeq_Results_Index"
      ), 
      as.symbol
    ), 
    file=paste0(
      "../Data/tmp/", 
      "B5_3_DE_12SVA_ALSFTD_Pre", 
      ".qrdata"
    ), 
    nthr=nthr)
) 

      # Y_B5_3_DE_12SVA_ALSFTD_Slurm.R ----------------------------------------- 


      # Import HPC environment with results ------------------------------------

rm(DDS_list, DESeq_Results_Index, SVAs_12SVs, SVAs_12SVs_Index)

qreadm(
  paste0(
    "../Data/tmp/",
    "B5_3_DE_12SVA_ALSFTD_Post", 
    ".qrdata"
  ),
  nthr=nthr
)


      # END HPC Cluster Outsourcing --------------------------------------------



    ## 3.4 DESeq SVA results summary -------------------------------------------

qsave(
  DDS_SVA_12sv_list, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD/SVA/12_SVs/", 
    "DDS_SVA_12sv_list", 
    ".qrds"
  ), 
  nthr=nthr
)


DESeq_Results_12SVs <- list() 
DESeq_Results_12SVs <- lapply(
  DDS_SVA_12sv_list, 
  FUN=function(x){
    return(
      if(is(x, "DESeqDataSet")){
        results(x, alpha=0.05, lfcThreshold = 0)
      }else{
        NA
      }
    )
  }
)


DESeq_Results_12SVs_Index <- DESeq_SVA_12sv_Results_Index 
rm(DESeq_SVA_12sv_Results_Index)

DESeq_Results_12SVs_Index$Covariates <- "12_SVs"
DESeq_Results_12SVs_Index$SVA[is.na(DESeq_Results_12SVs_Index$SVA)] <- FALSE
DESeq_Results_12SVs_Index$All_q_0.05 <- NA 
DESeq_Results_12SVs_Index$Up_q_0.05 <- NA  
DESeq_Results_12SVs_Index$Down_q_0.05 <- NA 
DESeq_Results_12SVs_Index$nSamples <- NA 
DESeq_Results_12SVs_Index$nSamples_HC <- NA  
DESeq_Results_12SVs_Index$nSamples_ALSFTD <- NA 
DESeq_Results_12SVs_Index$TotalReads <- NA 
DESeq_Results_12SVs_Index$MeanReads <- NA 
DESeq_Results_12SVs_Index$NonzeroCounts <- NA 
DESeq_Results_12SVs_Index$nFeatures <- NA 
DESeq_Results_12SVs_Index$LowCounts <- NA 

sample_data <- colData(DDS_SVA_12sv_list[[1]])

for (i in 1:length(DDS_SVA_12sv_list)){
  if(!is(DESeq_Results_12SVs[[i]], "DESeqResults")){
    
  }else{
    DESeq_Results_12SVs_Index$NonzeroCounts[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_12SVs[[i]]))[2], fixed("of "))[[1]][2], " with")[[1]][1])
    DESeq_Results_12SVs_Index$Up_q_0.05[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_12SVs[[i]]))[4], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_12SVs_Index$Down_q_0.05[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_12SVs[[i]]))[5], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_12SVs_Index$LowCounts[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_12SVs[[i]]))[7], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_12SVs_Index$nSamples[i] <- dim(DDS_SVA_12sv_list[[i]]@assays@data$counts)[2] 
    DESeq_Results_12SVs_Index$nSamples_HC[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_SVA_12sv_list[[i]]@assays@data$counts)])[1]
    DESeq_Results_12SVs_Index$nSamples_ALSFTD[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_SVA_12sv_list[[i]]@assays@data$counts)])[2]
    DESeq_Results_12SVs_Index$TotalReads[i] <- sum(colSums(DDS_SVA_12sv_list[[i]]@assays@data$counts)) 
    DESeq_Results_12SVs_Index$MeanReads[i] <- sum(colMeans(DDS_SVA_12sv_list[[i]]@assays@data$counts)) 
    DESeq_Results_12SVs_Index$nFeatures[i] <- sum(rowMeans(DDS_SVA_12sv_list[[i]]@assays@data$counts)>0) 
  }
  
}

DESeq_Results_12SVs_Index$All_q_0.05 <- DESeq_Results_12SVs_Index$Up_q_0.05 + DESeq_Results_12SVs_Index$Down_q_0.05 



    ## 3.5 Save Data -----------------------------------------------------------

qsave(
  DESeq_Results_12SVs, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD/SVA/12_SVs", 
    "DESeq_Results_12SVs", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_12SVs_Index, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD/SVA/12_SVs", 
    "DESeq_Results_12SVs_Index", 
    ".qrds"
  ), 
  nthr=nthr
)

rm(DDS_list, DDS_SVA_12sv_list, DESeq_Results_12SVs, DESeq_Results_12SVs_Index, sample_data, i) 




  ### 4.0 C9_ALSFTDvsHC --------------------------------------------------------



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



    ## 4.2 Find 12 SVs --------------------------------------------------------- 

SVAs_12SVs <- list() 
for (i in 1:nrow(DESeq_Results_Index)){
  SVAs_12SVs[[i]] <- NA
}

SVAs_12SVs_Index <- data.frame(
  Index=1:nrow(DESeq_Results_Index), 
  SVA=NA, 
  SVs = 12 
)

SVAs_12SVs_Index$CellTypeLevel <- DESeq_Results_Index$CellTypeLevel
SVAs_12SVs_Index$CellType <- DESeq_Results_Index$CellType
SVAs_12SVs_Index$Comparison <- DESeq_Results_Index$Comparison


DESeq_Results_Index$SVs <- NA 

for (i in 1:nrow(DESeq_Results_Index)){
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
    
    svseq <- svaseq(dat, mod, mod0, n.sv=12)
    DDS_list[[i]]$SV1 <- svseq$sv[,1]
    DDS_list[[i]]$SV2 <- svseq$sv[,2]
    DDS_list[[i]]$SV3 <- svseq$sv[,3]
    DDS_list[[i]]$SV4 <- svseq$sv[,4] 
    DDS_list[[i]]$SV5 <- svseq$sv[,5]
    DDS_list[[i]]$SV6 <- svseq$sv[,6]
    DDS_list[[i]]$SV7 <- svseq$sv[,7]
    DDS_list[[i]]$SV8 <- svseq$sv[,8]
    DDS_list[[i]]$SV9 <- svseq$sv[,9] 
    DDS_list[[i]]$SV10 <- svseq$sv[,10]
    DDS_list[[i]]$SV11 <- svseq$sv[,11]
    DDS_list[[i]]$SV12 <- svseq$sv[,12]
    DESeq_Results_Index$SVs[i] <- 12 
    SVAs_12SVs_Index$SVA[i] <- TRUE 
    SVAs_12SVs_Index$SVs[i] <- 12 
    message(paste0("SVAs #", i, " generated! "))
  }, error=function(e){
    DESeq_Results_Index$SVs[i] <<- 0 
    SVAs_12SVs_Index$SVA[i] <<- FALSE 
    SVAs_12SVs_Index$SVs[i] <<- 0 
    message(paste0("SVAs #", i, " could not be generated :("))
  })
}


qsave(
  SVAs_12SVs, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9/SVA/12_SVs/", 
    "SVAs_12SVs", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  SVAs_12SVs_Index, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9/SVA/12_SVs/", 
    "SVAs_12SVs_Index", 
    ".qrds"
  ), 
  nthr=nthr
)

all(SVAs_12SVs_Index$SVs[SVAs_12SVs_Index$SVA]==12)
rm(dat, mod, mod0, svseq, i, idx)



    ## 4.3 1507NotDOneOn DESeq2 with 12 SVs --------------------------------------------------


      # Start HPC Cluster Outsourcing ------------------------------------------


      # Export environment -----------------------------------------------------

do.call(
  qsavem, 
  c(
    lapply(
      list(
        "DDS_list", 
        "DESeq_Results_Index"
      ), 
      as.symbol
    ), 
    file=paste0(
      "../Data/tmp/", 
      "B5_4_DE_12SVA_ALSFTD_C9_Pre", 
      ".qrdata"
    ), 
    nthr=nthr)
) 


      # Y_B5_4_DE_12SVA_C9_ALSFTD_Slurm.R -------------------------------------- 


      # Import HPC environment with results ------------------------------------

rm(DDS_list, DESeq_Results_Index, SVAs_12SVs, SVAs_12SVs_Index)

qreadm(
  paste0(
    "../Data/tmp/",
    "B5_4_DE_12SVA_ALSFTD_C9_Post", 
    ".qrdata"
  ),
  nthr=nthr
)


      # END HPC Cluster Outsourcing --------------------------------------------



    ## 4.4 DESeq SVA results summary -------------------------------------------

qsave(
  DDS_SVA_12sv_list, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9/SVA/12_SVs/", 
    "DDS_SVA_12sv_list", 
    ".qrds"
  ), 
  nthr=nthr
)


DESeq_Results_12SVs <- list() 
DESeq_Results_12SVs <- lapply(
  DDS_SVA_12sv_list, 
  FUN=function(x){
    return(
      if(is(x, "DESeqDataSet")){
        results(x, alpha=0.05, lfcThreshold = 0)
      }else{
        NA
      }
    )
  }
)


DESeq_Results_12SVs_Index <- DESeq_SVA_12sv_Results_Index 
rm(DESeq_SVA_12sv_Results_Index)

DESeq_Results_12SVs_Index$Covariates <- "12_SVs"
DESeq_Results_12SVs_Index$SVA[is.na(DESeq_Results_12SVs_Index$SVA)] <- FALSE
DESeq_Results_12SVs_Index$All_q_0.05 <- NA 
DESeq_Results_12SVs_Index$Up_q_0.05 <- NA  
DESeq_Results_12SVs_Index$Down_q_0.05 <- NA 
DESeq_Results_12SVs_Index$nSamples <- NA 
DESeq_Results_12SVs_Index$nSamples_HC <- NA  
DESeq_Results_12SVs_Index$nSamples_ALSFTD <- NA 
DESeq_Results_12SVs_Index$TotalReads <- NA 
DESeq_Results_12SVs_Index$MeanReads <- NA 
DESeq_Results_12SVs_Index$NonzeroCounts <- NA 
DESeq_Results_12SVs_Index$nFeatures <- NA 
DESeq_Results_12SVs_Index$LowCounts <- NA 

sample_data <- colData(DDS_SVA_12sv_list[[1]])

for (i in 1:length(DDS_SVA_12sv_list)){
  if(!is(DESeq_Results_12SVs[[i]], "DESeqResults")){
    
  }else{
    DESeq_Results_12SVs_Index$NonzeroCounts[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_12SVs[[i]]))[2], fixed("of "))[[1]][2], " with")[[1]][1])
    DESeq_Results_12SVs_Index$Up_q_0.05[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_12SVs[[i]]))[4], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_12SVs_Index$Down_q_0.05[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_12SVs[[i]]))[5], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_12SVs_Index$LowCounts[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_12SVs[[i]]))[7], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_12SVs_Index$nSamples[i] <- dim(DDS_SVA_12sv_list[[i]]@assays@data$counts)[2] 
    DESeq_Results_12SVs_Index$nSamples_HC[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_SVA_12sv_list[[i]]@assays@data$counts)])[1]
    DESeq_Results_12SVs_Index$nSamples_ALSFTD[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_SVA_12sv_list[[i]]@assays@data$counts)])[2]
    DESeq_Results_12SVs_Index$TotalReads[i] <- sum(colSums(DDS_SVA_12sv_list[[i]]@assays@data$counts)) 
    DESeq_Results_12SVs_Index$MeanReads[i] <- sum(colMeans(DDS_SVA_12sv_list[[i]]@assays@data$counts)) 
    DESeq_Results_12SVs_Index$nFeatures[i] <- sum(rowMeans(DDS_SVA_12sv_list[[i]]@assays@data$counts)>0) 
  }
  
}

DESeq_Results_12SVs_Index$All_q_0.05 <- DESeq_Results_12SVs_Index$Up_q_0.05 + DESeq_Results_12SVs_Index$Down_q_0.05 



    ## 4.5 Save Data -----------------------------------------------------------

qsave(
  DESeq_Results_12SVs, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9/SVA/12_SVs", 
    "DESeq_Results_12SVs", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_12SVs_Index, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9/SVA/12_SVs", 
    "DESeq_Results_12SVs_Index", 
    ".qrds"
  ), 
  nthr=nthr
)

rm(DDS_list, DDS_SVA_12sv_list, DESeq_Results_12SVs, DESeq_Results_12SVs_Index, sample_data, i) 




  ### 5.0 NonC9_ALSFTDvsHC -----------------------------------------------------



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



    ## 5.2 Find 12 SVs --------------------------------------------------------- 

SVAs_12SVs <- list() 
for (i in 1:nrow(DESeq_Results_Index)){
  SVAs_12SVs[[i]] <- NA
}

SVAs_12SVs_Index <- data.frame(
  Index=1:nrow(DESeq_Results_Index), 
  SVA=NA, 
  SVs = 12 
)

SVAs_12SVs_Index$CellTypeLevel <- DESeq_Results_Index$CellTypeLevel
SVAs_12SVs_Index$CellType <- DESeq_Results_Index$CellType
SVAs_12SVs_Index$Comparison <- DESeq_Results_Index$Comparison


DESeq_Results_Index$SVs <- NA 

for (i in 1:nrow(DESeq_Results_Index)){
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
    
    svseq <- svaseq(dat, mod, mod0, n.sv=12)
    DDS_list[[i]]$SV1 <- svseq$sv[,1]
    DDS_list[[i]]$SV2 <- svseq$sv[,2]
    DDS_list[[i]]$SV3 <- svseq$sv[,3]
    DDS_list[[i]]$SV4 <- svseq$sv[,4] 
    DDS_list[[i]]$SV5 <- svseq$sv[,5]
    DDS_list[[i]]$SV6 <- svseq$sv[,6]
    DDS_list[[i]]$SV7 <- svseq$sv[,7]
    DDS_list[[i]]$SV8 <- svseq$sv[,8]
    DDS_list[[i]]$SV9 <- svseq$sv[,9] 
    DDS_list[[i]]$SV10 <- svseq$sv[,10]
    DDS_list[[i]]$SV11 <- svseq$sv[,11]
    DDS_list[[i]]$SV12 <- svseq$sv[,12]
    DESeq_Results_Index$SVs[i] <- 12 
    SVAs_12SVs_Index$SVA[i] <- TRUE 
    SVAs_12SVs_Index$SVs[i] <- 12 
    message(paste0("SVAs #", i, " generated! "))
  }, error=function(e){
    DESeq_Results_Index$SVs[i] <<- 0 
    SVAs_12SVs_Index$SVA[i] <<- FALSE 
    SVAs_12SVs_Index$SVs[i] <<- 0 
    message(paste0("SVAs #", i, " could not be generated :("))
  })
}


qsave(
  SVAs_12SVs, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_NonC9/SVA/12_SVs/", 
    "SVAs_12SVs", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  SVAs_12SVs_Index, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_NonC9/SVA/12_SVs/", 
    "SVAs_12SVs_Index", 
    ".qrds"
  ), 
  nthr=nthr
)

all(SVAs_12SVs_Index$SVs[SVAs_12SVs_Index$SVA]==12)
rm(dat, mod, mod0, svseq, i, idx)



    ## 5.3 DESeq2 with 12 SVs --------------------------------------------------


      # Start HPC Cluster Outsourcing ------------------------------------------


      # Export environment -----------------------------------------------------

do.call(
  qsavem, 
  c(
    lapply(
      list(
        "DDS_list", 
        "DESeq_Results_Index"
      ), 
      as.symbol
    ), 
    file=paste0(
      "../Data/tmp/", 
      "B5_5_DE_12SVA_ALSFTD_NonC9_Pre", 
      ".qrdata"
    ), 
    nthr=nthr)
) 

      # Y_B3_5_DE_12SVA_NonC9_ALSFTD_Slurm.R -----------------------------------


      # Import HPC environment with results ------------------------------------

rm(DDS_list, DESeq_Results_Index, SVAs_12SVs, SVAs_12SVs_Index)

qreadm(
  paste0(
    "../Data/tmp/",
    "B5_5_DE_12SVA_ALSFTD_NonC9_Post", 
    ".qrdata"
  ),
  nthr=nthr
)


      # END HPC Cluster Outsourcing --------------------------------------------


    ## 5.4 DESeq SVA results summary -------------------------------------------

qsave(
  DDS_SVA_12sv_list, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_NonC9/SVA/12_SVs/", 
    "DDS_SVA_12sv_list", 
    ".qrds"
  ), 
  nthr=nthr
)


DESeq_Results_12SVs <- list() 
DESeq_Results_12SVs <- lapply(
  DDS_SVA_12sv_list, 
  FUN=function(x){
    return(
      if(is(x, "DESeqDataSet")){
        results(x, alpha=0.05, lfcThreshold = 0)
      }else{
        NA
      }
    )
  }
)


DESeq_Results_12SVs_Index <- DESeq_SVA_12sv_Results_Index 
rm(DESeq_SVA_12sv_Results_Index)

DESeq_Results_12SVs_Index$Covariates <- "12_SVs"
DESeq_Results_12SVs_Index$SVA[is.na(DESeq_Results_12SVs_Index$SVA)] <- FALSE
DESeq_Results_12SVs_Index$All_q_0.05 <- NA 
DESeq_Results_12SVs_Index$Up_q_0.05 <- NA  
DESeq_Results_12SVs_Index$Down_q_0.05 <- NA 
DESeq_Results_12SVs_Index$nSamples <- NA 
DESeq_Results_12SVs_Index$nSamples_HC <- NA  
DESeq_Results_12SVs_Index$nSamples_ALSFTD <- NA 
DESeq_Results_12SVs_Index$TotalReads <- NA 
DESeq_Results_12SVs_Index$MeanReads <- NA 
DESeq_Results_12SVs_Index$NonzeroCounts <- NA 
DESeq_Results_12SVs_Index$nFeatures <- NA 
DESeq_Results_12SVs_Index$LowCounts <- NA 

sample_data <- colData(DDS_SVA_12sv_list[[1]])

for (i in 1:length(DDS_SVA_12sv_list)){
  if(!is(DESeq_Results_12SVs[[i]], "DESeqResults")){
    
  }else{
    DESeq_Results_12SVs_Index$NonzeroCounts[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_12SVs[[i]]))[2], fixed("of "))[[1]][2], " with")[[1]][1])
    DESeq_Results_12SVs_Index$Up_q_0.05[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_12SVs[[i]]))[4], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_12SVs_Index$Down_q_0.05[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_12SVs[[i]]))[5], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_12SVs_Index$LowCounts[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_12SVs[[i]]))[7], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_12SVs_Index$nSamples[i] <- dim(DDS_SVA_12sv_list[[i]]@assays@data$counts)[2] 
    DESeq_Results_12SVs_Index$nSamples_HC[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_SVA_12sv_list[[i]]@assays@data$counts)])[1]
    DESeq_Results_12SVs_Index$nSamples_ALSFTD[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_SVA_12sv_list[[i]]@assays@data$counts)])[2]
    DESeq_Results_12SVs_Index$TotalReads[i] <- sum(colSums(DDS_SVA_12sv_list[[i]]@assays@data$counts)) 
    DESeq_Results_12SVs_Index$MeanReads[i] <- sum(colMeans(DDS_SVA_12sv_list[[i]]@assays@data$counts)) 
    DESeq_Results_12SVs_Index$nFeatures[i] <- sum(rowMeans(DDS_SVA_12sv_list[[i]]@assays@data$counts)>0) 
  }
  
}

DESeq_Results_12SVs_Index$All_q_0.05 <- DESeq_Results_12SVs_Index$Up_q_0.05 + DESeq_Results_12SVs_Index$Down_q_0.05 



    ## 5.5 Save Data -----------------------------------------------------------

qsave(
  DESeq_Results_12SVs, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_NonC9/SVA/12_SVs", 
    "DESeq_Results_12SVs", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_12SVs_Index, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_NonC9/SVA/12_SVs", 
    "DESeq_Results_12SVs_Index", 
    ".qrds"
  ), 
  nthr=nthr
)

rm(DDS_list, DDS_SVA_12sv_list, DESeq_Results_12SVs, DESeq_Results_12SVs_Index, sample_data, i) 




  ### 6.0 C9_ALSFTDvsNonC9_ALSFTD ----------------------------------------------



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



    ## 6.2 Find 12 SVs --------------------------------------------------------- 

SVAs_12SVs <- list() 
for (i in 1:nrow(DESeq_Results_Index)){
  SVAs_12SVs[[i]] <- NA
}

SVAs_12SVs_Index <- data.frame(
  Index=1:nrow(DESeq_Results_Index), 
  SVA=NA, 
  SVs = 12 
)

SVAs_12SVs_Index$CellTypeLevel <- DESeq_Results_Index$CellTypeLevel
SVAs_12SVs_Index$CellType <- DESeq_Results_Index$CellType
SVAs_12SVs_Index$Comparison <- DESeq_Results_Index$Comparison


DESeq_Results_Index$SVs <- NA 

for (i in 1:nrow(DESeq_Results_Index)){
  tryCatch({
    
    dat  <- counts(DDS_list[[i]], normalized = TRUE)
    idx  <- rowMeans(dat) > 1 
    dat  <- dat[idx, ] 
    if(DESeq_Results_Index$Comparison[i]=="ALSFTD_C9vsALSFTD_NonC9"){
      mod  <- model.matrix(~ Case_Type, colData(DDS_list[[i]]))
    }else{
      mod  <- model.matrix(~ Rand, colData(DDS_list[[i]]))
    }
    mod0 <- model.matrix(~ 1, colData(DDS_list[[i]]))
    
    svseq <- svaseq(dat, mod, mod0, n.sv=12)
    DDS_list[[i]]$SV1 <- svseq$sv[,1]
    DDS_list[[i]]$SV2 <- svseq$sv[,2]
    DDS_list[[i]]$SV3 <- svseq$sv[,3]
    DDS_list[[i]]$SV4 <- svseq$sv[,4] 
    DDS_list[[i]]$SV5 <- svseq$sv[,5]
    DDS_list[[i]]$SV6 <- svseq$sv[,6]
    DDS_list[[i]]$SV7 <- svseq$sv[,7]
    DDS_list[[i]]$SV8 <- svseq$sv[,8]
    DDS_list[[i]]$SV9 <- svseq$sv[,9] 
    DDS_list[[i]]$SV10 <- svseq$sv[,10]
    DDS_list[[i]]$SV11 <- svseq$sv[,11]
    DDS_list[[i]]$SV12 <- svseq$sv[,12]
    DESeq_Results_Index$SVs[i] <- 12 
    SVAs_12SVs_Index$SVA[i] <- TRUE 
    SVAs_12SVs_Index$SVs[i] <- 12 
    message(paste0("SVAs #", i, " generated! "))
  }, error=function(e){
    DESeq_Results_Index$SVs[i] <<- 0 
    SVAs_12SVs_Index$SVA[i] <<- FALSE 
    SVAs_12SVs_Index$SVs[i] <<- 0 
    message(paste0("SVAs #", i, " could not be generated :("))
  })
}


qsave(
  SVAs_12SVs, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9_NonC9/SVA/12_SVs/", 
    "SVAs_12SVs", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  SVAs_12SVs_Index, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9_NonC9/SVA/12_SVs/", 
    "SVAs_12SVs_Index", 
    ".qrds"
  ), 
  nthr=nthr
)

all(SVAs_12SVs_Index$SVs[SVAs_12SVs_Index$SVA]==12)
rm(dat, mod, mod0, svseq, i, idx)



    ## 6.3 DESeq2 with 12 SVs --------------------------------------------------


      # Start HPC Cluster Outsourcing ------------------------------------------


      # Export environment -----------------------------------------------------

do.call(
  qsavem, 
  c(
    lapply(
      list(
        "DDS_list", 
        "DESeq_Results_Index"
      ), 
      as.symbol
    ), 
    file=paste0(
      "../Data/tmp/", 
      "B5_6_DE_12SVA_ALSFTD_C9_NonC9_Pre", 
      ".qrdata"
    ), 
    nthr=nthr)
) 

      # Y_B3_6_DE_12SVA_ALSFTD_C9_NonC9_Slurm.R --------------------------------




      # Not Done 15.07 Import HPC environment with results ---------------------
rm(DDS_list, DESeq_Results_Index, SVAs_12SVs, SVAs_12SVs_Index)

qreadm(
  paste0(
    "../Data/tmp/",
    "B5_6_DE_12SVA_ALSFTD_C9_NonC9_Post", 
    ".qrdata"
  ),
  nthr=nthr
)


      # END HPC Cluster Outsourcing --------------------------------------------



    ## 6.4 DESeq SVA results summary -------------------------------------------

qsave(
  DDS_SVA_12sv_list, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9_NonC9/SVA/12_SVs/", 
    "DDS_SVA_12sv_list", 
    ".qrds"
  ), 
  nthr=nthr
)


DESeq_Results_12SVs <- list() 
DESeq_Results_12SVs <- lapply(
  DDS_SVA_12sv_list, 
  FUN=function(x){
    return(
      if(is(x, "DESeqDataSet")){
        results(x, alpha=0.05, lfcThreshold = 0)
      }else{
        NA
      }
    )
  }
)


DESeq_Results_12SVs_Index <- DESeq_SVA_12sv_Results_Index 
rm(DESeq_SVA_12sv_Results_Index)

DESeq_Results_12SVs_Index$Covariates <- "12_SVs"
DESeq_Results_12SVs_Index$SVA[is.na(DESeq_Results_12SVs_Index$SVA)] <- FALSE
DESeq_Results_12SVs_Index$All_q_0.05 <- NA 
DESeq_Results_12SVs_Index$Up_q_0.05 <- NA  
DESeq_Results_12SVs_Index$Down_q_0.05 <- NA 
DESeq_Results_12SVs_Index$nSamples <- NA 
DESeq_Results_12SVs_Index$nSamples_ALSFTD_C9 <- NA  
DESeq_Results_12SVs_Index$nSamples_ALSFTD_NonC9 <- NA 
DESeq_Results_12SVs_Index$TotalReads <- NA 
DESeq_Results_12SVs_Index$MeanReads <- NA 
DESeq_Results_12SVs_Index$NonzeroCounts <- NA 
DESeq_Results_12SVs_Index$nFeatures <- NA 
DESeq_Results_12SVs_Index$LowCounts <- NA 

sample_data <- colData(DDS_SVA_12sv_list[[1]])

for (i in 1:length(DDS_SVA_12sv_list)){
  if(!is(DESeq_Results_12SVs[[i]], "DESeqResults")){
    
  }else{
    DESeq_Results_12SVs_Index$NonzeroCounts[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_12SVs[[i]]))[2], fixed("of "))[[1]][2], " with")[[1]][1])
    DESeq_Results_12SVs_Index$Up_q_0.05[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_12SVs[[i]]))[4], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_12SVs_Index$Down_q_0.05[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_12SVs[[i]]))[5], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_12SVs_Index$LowCounts[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_12SVs[[i]]))[7], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_12SVs_Index$nSamples[i] <- dim(DDS_SVA_12sv_list[[i]]@assays@data$counts)[2] 
    DESeq_Results_12SVs_Index$nSamples_ALSFTD_C9[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_SVA_12sv_list[[i]]@assays@data$counts)])[1]
    DESeq_Results_12SVs_Index$nSamples_ALSFTD_NonC9[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_SVA_12sv_list[[i]]@assays@data$counts)])[2]
    DESeq_Results_12SVs_Index$TotalReads[i] <- sum(colSums(DDS_SVA_12sv_list[[i]]@assays@data$counts)) 
    DESeq_Results_12SVs_Index$MeanReads[i] <- sum(colMeans(DDS_SVA_12sv_list[[i]]@assays@data$counts)) 
    DESeq_Results_12SVs_Index$nFeatures[i] <- sum(rowMeans(DDS_SVA_12sv_list[[i]]@assays@data$counts)>0) 
  }
  
}

DESeq_Results_12SVs_Index$All_q_0.05 <- DESeq_Results_12SVs_Index$Up_q_0.05 + DESeq_Results_12SVs_Index$Down_q_0.05 



    ## 6.5 Save Data -----------------------------------------------------------

qsave(
  DESeq_Results_12SVs, 
  paste0(
    "../Data/DE/WNN/APseudobulk/DESeq2/ALSFTD_C9_NonC9/SVA/12_SVs", 
    "DESeq_Results_12SVs", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_12SVs_Index, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9_NonC9/SVA/12_SVs", 
    "DESeq_Results_12SVs_Index", 
    ".qrds"
  ), 
  nthr=nthr
)

rm(DDS_list, DDS_SVA_12sv_list, DESeq_Results_12SVs, DESeq_Results_12SVs_Index, sample_data, i) 





