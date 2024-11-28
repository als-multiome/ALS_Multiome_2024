

  ### 0.0 Load libraries -------------------------------------------------------

source("~/ALS_Brain_Multiome.Rcfg")


library(qs) 
library(data.table)
library(Seurat)
library(tidyverse) 
library(DESeq2)
library(sva)




  ### 1.0 AllCase --------------------------------------------------------------



    ## 1.1 Load Data ----------------------------------------------------------- 

DDS_list <- qread(
  paste0(
    "../Data/DA/WNN/AllCase/", 
    "DDS_DESEqed_list", 
    ".qrds"
  ), 
  nthr=nthr
)


DESeq_Results_Index <- qread(
  paste0(
    "../Data/DE/WNN/AllCase/", 
    "DESeq_Results_Index", 
    ".qrds"
  ), 
  nthr=nthr
)    

DDS_list <- DDS_list[1:76] 
DESeq_Results_Index <- DESeq_Results_Index[1:76,]



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
    "../Data/DA/WNN/AllCase/SVA/1_15_SVs/", 
    "SVAs_15SVs", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  SVAs_15SVs_Index, 
  paste0(
    "../Data/DA/WNN/AllCase/SVA/1_15_SVs/", 
    "SVAs_15SVs_Index", 
    ".qrds"
  ), 
  nthr=nthr
)


all(SVAs_15SVs_Index$SVs[SVAs_15SVs_Index$SVA]==15)
rm(dat, mod, mod0, svseq, i, idx)



    ## 1.3 Add 1-15 SVs to DDS_list --------------------------------------

DESeq_SVA_1_15_Results_Index <- DESeq_Results_Index[,c(1:9)] 
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
      "../Data/DA/WNN/AllCase/SVA/1_15_SVs/",
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
      "C3_1_DA_ATAC_WNN_SVA_AllCase_Pre", 
      ".qrdata"
    ), 
    nthr=nthr
  )
)         

rm(list=ls())


      # Y_C4_1_DA_ATAC_WNN_SVA_AllCase_Slurm.R, Files 1-15 ---------------------


      # Import HPC Output and summarize results --------------------------------

for (i in 1:15){
  
  DDS_List <- qread(
    file=paste0(
      "../Data/DA/WNN/AllCase/SVA/1_15_SVs/DDS_", 
      i, 
      "SVs_list", 
      ".qrds"
    ), 
    nthr=nthr
  )
  
  
  Results_Index <- qread(
    paste0(
      "../Data/DA/WNN/AllCase/SVA/1_15_SVs/",
      "SVAs_15SVs_Index", 
      ".qrds"
    ), 
    nthr=nthr
  )[,c(1,4,5,6)]
  
  Results_Index$SVs <- i 
  Results_Index$SVA <- unlist(lapply(DDS_List, FUN=function(x){is(x, "DESeqDataSet")}))
  
  
  Results_ALS <- list()  
  Results_ALSFTD <- list()  
  
  
  Results_ALS <- lapply(
    DDS_List, 
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
  
  
  Results_ALSFTD <- lapply(
    DDS_List, 
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
  
  
  Results_Index$All_q_0.05_ALS <- NA 
  Results_Index$All_q_0.05_ALSFTD <- NA 
  
  Results_Index$Up_q_0.05_ALS <- NA   
  Results_Index$Up_q_0.05_ALSFTD <- NA  
  
  Results_Index$Down_q_0.05_ALS <- NA 
  Results_Index$Down_q_0.05_ALSFTD <- NA 
  
  Results_Index$nSamples <- NA 
  Results_Index$nSamples_HC <- NA  
  Results_Index$nSamples_ALS <- NA 
  Results_Index$nSamples_ALSFTD <- NA 
  
  Results_Index$TotalReads <- NA 
  
  Results_Index$MeanReads <- NA 
  
  Results_Index$nFeatures <- NA 
  
  Results_Index$NonzeroCounts_ALS <- NA 
  Results_Index$NonzeroCounts_ALSFTD <- NA 
  
  Results_Index$LowCounts_ALS <- NA 
  Results_Index$LowCounts_ALSFTD <- NA 
  
  sample_data <- DDS_List[[i]]@colData
  Results_Index$Covariates <- NA 
  
  for (j in 1:length(DDS_List)){
    if(!Results_Index$SVA[j]){
      
    }else{
      Results_Index$Covariates[j] <- paste0(i, "_SVs")
      Results_Index$NonzeroCounts_ALS[j] <- as.numeric(str_split(str_split(capture.output(summary(Results_ALS[[j]]))[2], fixed("of "))[[1]][2], " with")[[1]][1]) 
      Results_Index$NonzeroCounts_ALSFTD[j] <- as.numeric(str_split(str_split(capture.output(summary(Results_ALSFTD[[j]]))[2], fixed("of "))[[1]][2], " with")[[1]][1]) 
      
      Results_Index$Up_q_0.05_ALS[j] <- as.numeric(str_split(str_split(capture.output(summary(Results_ALS[[j]]))[4], fixed(": "))[[1]][2], ", ")[[1]][1])
      Results_Index$Up_q_0.05_ALSFTD[j] <- as.numeric(str_split(str_split(capture.output(summary(Results_ALSFTD[[j]]))[4], fixed(": "))[[1]][2], ", ")[[1]][1]) 
      
      Results_Index$Down_q_0.05_ALS[j] <- as.numeric(str_split(str_split(capture.output(summary(Results_ALS[[j]]))[5], fixed(": "))[[1]][2], ", ")[[1]][1])
      Results_Index$Down_q_0.05_ALSFTD[j] <- as.numeric(str_split(str_split(capture.output(summary(Results_ALSFTD[[j]]))[5], fixed(": "))[[1]][2], ", ")[[1]][1])
      
      Results_Index$LowCounts_ALS[j] <- as.numeric(str_split(str_split(capture.output(summary(Results_ALS[[j]]))[7], fixed(": "))[[1]][2], ", ")[[1]][1])
      Results_Index$LowCounts_ALSFTD[j] <- as.numeric(str_split(str_split(capture.output(summary(Results_ALSFTD[[j]]))[7], fixed(": "))[[1]][2], ", ")[[1]][1])
      
      
      if(is(object = DDS_List[[j]], class="DESeqDataSet")){
        
        Results_Index$nSamples[j] <- dim(DDS_List[[j]]@assays@data$counts)[2] 
        
        Results_Index$nSamples_HC[j] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_List[[j]]@assays@data$counts)])[1]
        Results_Index$nSamples_ALS[j] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_List[[j]]@assays@data$counts)])[2]
        Results_Index$nSamples_ALSFTD[j] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_List[[j]]@assays@data$counts)])[3]
        
        Results_Index$TotalReads[j] <- sum(colSums(DDS_List[[j]]@assays@data$counts)) 
        Results_Index$MeanReads[j] <- sum(colMeans(DDS_List[[j]]@assays@data$counts)) 
        Results_Index$nFeatures[j] <- sum(rowMeans(DDS_List[[j]]@assays@data$counts)>0) 
        
      }else{
        
        Results_12SVs_Index$nSamples[j] <- NA
        
        Results_Index$nSamples_HC[j] <- NA
        Results_Index$nSamples_ALS[j] <- NA
        Results_Index$nSamples_ALSFTD[j] <- NA
        
        Results_Index$TotalReads[j] <- NA
        Results_Index$MeanReads[j] <- NA
        Results_Index$nFeatures[j] <- NA
      }
      
    }
    
  }
  
  Results_Index$All_q_0.05_ALS <- Results_Index$Up_q_0.05_ALS + Results_Index$Down_q_0.05_ALS 
  Results_Index$All_q_0.05_ALSFTD <- Results_Index$Up_q_0.05_ALSFTD + Results_Index$Down_q_0.05_ALSFTD 
  
  
  qsave(
    DDS_List, 
    paste0(
      "../Data/DA/WNN/AllCase/SVA/1_15_SVs/", 
      "DDS_SVA_",
      i, 
      "sv_list", 
      ".qrds"
    ), 
    nthr=nthr
  )
  
  qsave(
    Results_ALS, 
    paste0(
      "../Data/DA/WNN/AllCase/SVA/1_15_SVs/", 
      "DESeq_Results_", 
      i, 
      "SVs_ALS", 
      ".qrds"
    ), 
    nthr=nthr
  )
  
  qsave(
    Results_ALSFTD, 
    paste0(
      "../Data/DA/WNN/AllCase/SVA/1_15_SVs/", 
      "DESeq_Results_",
      i, 
      "SVs_ALSFTD", 
      ".qrds"
    ), 
    nthr=nthr
  )
  
  qsave(
    Results_Index, 
    paste0(
      "../Data/DA/WNN/AllCase/SVA/1_15_SVs/", 
      "DESeq_Results_", 
      i, 
      "SVs_Index", 
      ".qrds"
    ), 
    nthr=nthr
  )
  
  qsave(
    sample_data, 
    paste0(
      "../Data/DA/WNN/AllCase/SVA/1_15_SVs/", 
      "Sample_data_", 
      i,
      "SVs", 
      ".qrds"
    ), 
    nthr=nthr
  )
  
  rm(
    DDS_List, 
    Results_Index, 
    Results_ALS, 
    Results_ALSFTD, 
    sample_data, 
    j
  ) 
  
  message("######################################################")
  print(Sys.time())
  message(
    paste0(
      "Processed and saved data with ", 
      i, 
      "SVs!"
    )
  )
  message("######################################################")
  
}
  
rm(i) 



