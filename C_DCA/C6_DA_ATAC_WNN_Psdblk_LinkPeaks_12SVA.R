

  ### 0.0 Load libraries -------------------------------------------------------

source("~/ALS_Brain_Multiome.Rcfg")


library(qs) 
library(data.table)
library(Seurat)
library(tidyverse) 
library(DESeq2)
library(sva)



qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)




  ### 1.0 AllCase --------------------------------------------------------------



    ## 1.1 Load Data ----------------------------------------------------------- 

DDS_list <- qread(
  paste0(
    "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/", 
    "DDS_DESeqed_list", 
    ".qrds"
  ), 
  nthr=nthr
)


DESeq_Results_Index <- qread(
  paste0(
    "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/", 
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
    "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/SVA/12_SVs/", 
    "SVAs_12SVs", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  SVAs_12SVs_Index, 
  paste0(
    "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/SVA/12_SVs/", 
    "SVAs_12SVs_Index", 
    ".qrds"
  ), 
  nthr=nthr
)


all(SVAs_12SVs_Index$SVs[SVAs_12SVs_Index$SVA]==12)
rm(dat, mod, mod0, svseq, i, idx)



    ## 1.3 Add 12 SVs to DDS_list ----------------------------------------------

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
    
    DDS_SVA_12sv_list[[i]] <- DDS.tmp 
    rm(DDS.tmp)
    
    DESeq_SVA_12sv_Results_Index$SVA[i] <- TRUE 
    message(i)
    message("########################################################")
    message(Sys.time())
    }, error=function(e){
    DESeq_SVA_12sv_Results_Index$SVA[i] <<- FALSE 
    message("########################################################")
    message(Sys.time()) 
  })
}  

rm(i, DDS.tmp)

colnames(colData(DDS_SVA_12sv_list[[201]]))
design(DDS_SVA_12sv_list[[201]])

colnames(colData(DDS_SVA_12sv_list[[202]]))
design(DDS_SVA_12sv_list[[202]])


      # START HPC Cluster Outsourcing ------------------------------------------


      # Export Environment -----------------------------------------------------

do.call(
  qsavem,
  c(
    lapply(
      c(
        "DDS_SVA_12sv_list", 
        "DESeq_SVA_12sv_Results_Index"
      ), 
      as.symbol
    ), 
    file=paste0(
      "../Data/tmp/", 
      "C6_1_DA_ATAC_WNN_SVA_AllCase_LinkPeaks_Pre", 
      ".qrdata"
    ), 
    nthr=nthr
  )
)         

rm(DDS_list, DDS_SVA_12sv_list, DESeq_SVA_12sv_Results_Index, SVAs_12SVs, SVAs_12SVs_Index)


      # Y_C6_1_DA_ATAC_WNN_SVA_AllCase_LinkPeaks_Slurm. R ----------------------



      # Import HPC Environment  ------------------------------------------------

qreadm(
  paste0(
    "../Data/tmp/",
    "C6_1_DA_ATAC_WNN_SVA_AllCase_LinkPeaks_Post", 
    ".qrdata"
  ), 
  nthr=nthr
)


      # END HPC Cluster Outsourcing --------------------------------------------



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

DESeq_Results_12SVs_Index <- DESeq_SVA_12sv_Results_Index 
rm(DESeq_SVA_12sv_Results_Index)

DESeq_Results_12SVs_Index$All_q_0.05_ALS <- NA 
DESeq_Results_12SVs_Index$All_q_0.05_ALSFTD <- NA 

DESeq_Results_12SVs_Index$Up_q_0.05_ALS <- NA   
DESeq_Results_12SVs_Index$Up_q_0.05_ALSFTD <- NA  

DESeq_Results_12SVs_Index$Down_q_0.05_ALS <- NA 
DESeq_Results_12SVs_Index$Down_q_0.05_ALSFTD <- NA 

DESeq_Results_12SVs_Index$nSamples <- NA 
DESeq_Results_12SVs_Index$nSamples_HC <- NA  
DESeq_Results_12SVs_Index$nSamples_ALS <- NA 
DESeq_Results_12SVs_Index$nSamples_ALSFTD <- NA 

DESeq_Results_12SVs_Index$TotalReads <- NA 

DESeq_Results_12SVs_Index$MeanReads <- NA 

DESeq_Results_12SVs_Index$nFeatures <- NA 

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
    
    
    if(is(object = DDS_SVA_12sv_list[[i]], class="DESeqDataSet")){
      
      DESeq_Results_12SVs_Index$nSamples[i] <- dim(DDS_SVA_12sv_list[[i]]@assays@data$counts)[2] 
      
      DESeq_Results_12SVs_Index$nSamples_HC[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_SVA_12sv_list[[i]]@assays@data$counts)])[1]
      DESeq_Results_12SVs_Index$nSamples_ALS[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_SVA_12sv_list[[i]]@assays@data$counts)])[2]
      DESeq_Results_12SVs_Index$nSamples_ALSFTD[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_SVA_12sv_list[[i]]@assays@data$counts)])[3]
      
      DESeq_Results_12SVs_Index$TotalReads[i] <- sum(colSums(DDS_SVA_12sv_list[[i]]@assays@data$counts)) 
      DESeq_Results_12SVs_Index$MeanReads[i] <- sum(colMeans(DDS_SVA_12sv_list[[i]]@assays@data$counts)) 
      DESeq_Results_12SVs_Index$nFeatures[i] <- sum(rowMeans(DDS_SVA_12sv_list[[i]]@assays@data$counts)>0) 
    }else{
      DESeq_Results_12SVs_Index$nSamples[i] <- NA
      
      DESeq_Results_12SVs_Index$nSamples_HC[i] <- NA
      DESeq_Results_12SVs_Index$nSamples_ALS[i] <- NA
      DESeq_Results_12SVs_Index$nSamples_ALSFTD[i] <- NA
      
      DESeq_Results_12SVs_Index$TotalReads[i] <- NA
      DESeq_Results_12SVs_Index$MeanReads[i] <- NA
      DESeq_Results_12SVs_Index$nFeatures[i] <- NA
    }
    
  }
  
}

DESeq_Results_12SVs_Index$All_q_0.05_ALS <- DESeq_Results_12SVs_Index$Up_q_0.05_ALS + DESeq_Results_12SVs_Index$Down_q_0.05_ALS 
DESeq_Results_12SVs_Index$All_q_0.05_ALSFTD <- DESeq_Results_12SVs_Index$Up_q_0.05_ALSFTD + DESeq_Results_12SVs_Index$Down_q_0.05_ALSFTD 



    ## 1.5 Save Data -----------------------------------------------------------

qsave(
  DDS_SVA_12sv_list, 
  paste0(
    "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/SVA/12_SVs/", 
    "DDS_SVA_12sv_list", 
    ".qrds"
  ), 
  nthr=nthr
)


qsave(
  DESeq_Results_12SVs_ALS, 
  paste0(
    "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/SVA/12_SVs/", 
    "DESeq_Results_12SVs_ALS", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_12SVs_ALSFTD, 
  paste0(
    "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/SVA/12_SVs/", 
    "DESeq_Results_12SVs_ALSFTD", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_12SVs_Index, 
  paste0(
    "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/SVA/12_SVs/", 
    "DESeq_Results_12SVs_Index", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  sample_data, 
  paste0(
    "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/SVA/12_SVs/", 
    "sample_data", 
    ".qrds"
  ), 
  nthr=nthr
)

rm(
  DDS_SVA_12sv_list, 
  DESeq_Results_12SVs_ALS, 
  DESeq_Results_12SVs_ALSFTD, 
  DESeq_Results_12SVs_Index, 
  sample_data, 
  i
) 



