
  ### 0.0 Load libraries -------------------------------------------------------

source("~/ALS_Brain_Multiome.Rcfg")

library(qs) 
library(data.table)
library(Seurat)
library(tidyverse) 
library(DESeq2)

qs::set_trust_promises(TRUE)




  ### 1.0 AllCase ------------------------------------------------------------------



    ## 1.1 Load data -----------------------------------------------------------

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



    ### 1.2 DE Analysis --------------------------------------------------------

DESeq_Sex_PctMito_Results_Index <- DESeq_Results_Index[,c(1:9)] 


DDS_list_Covs <- list() 
for (i in 1:nrow(DESeq_Sex_PctMito_Results_Index)){
  DDS_list_Covs[[i]] <- NA
}

rm(DESeq_Results_Index, i) 

for (i in 1:nrow(DESeq_Sex_PctMito_Results_Index)){ 
  tryCatch({
    DDS.tmp <- DDS_list[[i]] 
    DDS.tmp$Sex <- factor(DDS.tmp$Sex, levels = c("f", "m"))
    if(DESeq_Sex_PctMito_Results_Index$Comparison[i]=="All_Cases"){
      design(DDS.tmp) <-  ~ Sex + pctMito_RNA_scaled + Case
    }else{
      design(DDS.tmp) <- ~ Sex + pctMito_RNA_scaled + Rand 
    } 
    DDS_list_Covs[[i]] <- DDS.tmp 
    rm(DDS.tmp)
  }, error=function(e){
      message(
        paste0(
          "Design of Object # ", 
          i, 
          " could not be set!"
        )
      )
    }
  ) 
}


      # START HPC Cluster Outsourcing ------------------------------------------


      # Export Environment -----------------------------------------------------

do.call(
  qsavem,
  c(
    lapply(
      c(
        "DDS_list_Covs", 
        "DESeq_Sex_PctMito_Results_Index"
      ), 
      as.symbol
    ), 
    file=paste0(
      "../Data/tmp/", 
      "B3_1_DE_RNA_WNN_Pre", 
      ".qrdata"
    ), 
    nthr=nthr
  )
)         


      # Y_B3_1_RNA_DE_Sex_PctMito_AllCases_Slurm.R -----------------------------


      # Import HPC Environment  ------------------------------------------------

rm(
  i, 
  DDS_list, 
  DDS_list_Covs, 
  DESeq_Sex_PctMito_Results_Index
)


qreadm(
  paste0(
    "../Data/tmp/", 
    "B3_1_DE_RNA_WNN_Post", 
    ".qrdata"
  ), 
  nthr=nthr
)



    ## 1.3 Collect DE stats ----------------------------------------------------

resultsNames(DDS_List[[1]])
resultsNames(DDS_List[[2]])

DESeq_Results_ALS <- list() 
DESeq_Results_ALSFTD <- list()  

for (i in 1:length(DDS_List)){
  
  if(i%%2==1){
    tryCatch({
      DESeq_Results_ALS[[i]] <- results(
        DDS_List[[i]], 
        name="Case_ALS_vs_HC", 
        alpha=0.05, 
        lfcThreshold = 0
      )
    },error=function(e){
      message(paste0("Results #", i, " , ALS, could not be generated")) 
      DESeq_Results_ALS[[i]] <<- NA  
    }) 
    
    tryCatch({
      DESeq_Results_ALSFTD[[i]] <- results(
        DDS_List[[i]], 
        name="Case_ALS_FTD_vs_HC", 
        alpha=0.05, 
        lfcThreshold = 0
      )
    },error=function(e){
      message(paste0("Results #", i, " , ALSFTD, could not be generated")) 
      DESeq_Results_ALS[[i]] <<- NA  
    })
  }else{
    tryCatch({
      DESeq_Results_ALS[[i]] <- results(
        DDS_List[[i]], 
        name="Rand_R2_vs_R1", 
        alpha=0.05, 
        lfcThreshold = 0
      )
    },error=function(e){
      message(paste0("Results #", i, " , ALS/Rand, could not be generated")) 
      DESeq_Results_ALS[[i]] <<- NA  
    }) 
    
    tryCatch({
      DESeq_Results_ALSFTD[[i]] <- results(
        DDS_List[[i]], 
        name="Rand_R3_vs_R1", 
        alpha=0.05, 
        lfcThreshold = 0
      )
    },error=function(e){
      message(paste0("Results #", i, " , ALSFTD/Rand, could not be generated")) 
      DESeq_Results_ALS[[i]] <<- NA  
    })
  }
  
}

DESeq_Results_Index <- DESeq_Sex_PctMito_Results_Index
rm(DESeq_Sex_PctMito_Results_Index) 

DESeq_Results_Index$All_q_0.05_ALS <- NA 
DESeq_Results_Index$All_q_0.05_ALSFTD <- NA 

DESeq_Results_Index$Up_q_0.05_ALS <- NA   
DESeq_Results_Index$Up_q_0.05_ALSFTD <- NA  

DESeq_Results_Index$Down_q_0.05_ALS <- NA 
DESeq_Results_Index$Down_q_0.05_ALSFTD <- NA 

DESeq_Results_Index$nSamples <- NA 
DESeq_Results_Index$nSamples_HC <- NA  
DESeq_Results_Index$nSamples_ALS <- NA 
DESeq_Results_Index$nSamples_ALSFTD <- NA 

DESeq_Results_Index$TotalReads <- NA 

DESeq_Results_Index$MeanReads <- NA 

DESeq_Results_Index$NonzeroCounts_ALS <- NA 
DESeq_Results_Index$NonzeroCounts_ALSFTD <- NA 

DESeq_Results_Index$LowCounts_ALS <- NA 
DESeq_Results_Index$LowCounts_ALSFTD <- NA 


for (i in 1:length(DDS_List)){
  if(!is.na(DESeq_Results_Index$Remarks[i])){
    
  }else{
    
    sample_data = colData(DDS_List[[i]])
    DESeq_Results_Index$NonzeroCounts_ALS[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_ALS[[i]]))[2], fixed("of "))[[1]][2], " with")[[1]][1]) 
    DESeq_Results_Index$NonzeroCounts_ALSFTD[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_ALSFTD[[i]]))[2], fixed("of "))[[1]][2], " with")[[1]][1]) 
    
    DESeq_Results_Index$Up_q_0.05_ALS[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_ALS[[i]]))[4], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_Index$Up_q_0.05_ALSFTD[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_ALSFTD[[i]]))[4], fixed(": "))[[1]][2], ", ")[[1]][1]) 
    
    DESeq_Results_Index$Down_q_0.05_ALS[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_ALS[[i]]))[5], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_Index$Down_q_0.05_ALSFTD[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_ALSFTD[[i]]))[5], fixed(": "))[[1]][2], ", ")[[1]][1])
    
    DESeq_Results_Index$LowCounts_ALS[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_ALS[[i]]))[7], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_Index$LowCounts_ALSFTD[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_ALSFTD[[i]]))[7], fixed(": "))[[1]][2], ", ")[[1]][1])
    
    DESeq_Results_Index$nSamples[i] <- dim(DDS_List[[i]]@assays@data$counts)[2] 
    
    
    DESeq_Results_Index$nSamples_HC[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_List[[i]]@assays@data$counts)])[1]
    DESeq_Results_Index$nSamples_ALS[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_List[[i]]@assays@data$counts)])[2]
    DESeq_Results_Index$nSamples_ALSFTD[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_List[[i]]@assays@data$counts)])[3]
    
    DESeq_Results_Index$TotalReads[i] <- sum(colSums(DDS_List[[i]]@assays@data$counts)) 
    DESeq_Results_Index$MeanReads[i] <- sum(colMeans(DDS_List[[i]]@assays@data$counts)) 
    DESeq_Results_Index$nFeatures[i] <- sum(rowMeans(DDS_List[[i]]@assays@data$counts)>0) 
  }
  
}

DESeq_Results_Index$All_q_0.05_ALS <- DESeq_Results_Index$Up_q_0.05_ALS + DESeq_Results_Index$Down_q_0.05_ALS 
DESeq_Results_Index$All_q_0.05_ALSFTD <- DESeq_Results_Index$Up_q_0.05_ALSFTD + DESeq_Results_Index$Down_q_0.05_ALSFTD 




    ## 1.4 Export data ---------------------------------------------------------

qsave(
  DDS_List, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/Sex_PctMitoScaled/", 
    "DDS_DESEqed_list", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_ALS, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/Sex_PctMitoScaled/", 
    "DESeq_Results_ALS", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_ALSFTD, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/Sex_PctMitoScaled/", 
    "DESeq_Results_ALSFTD", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_Index, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/Sex_PctMitoScaled/", 
    "DESeq_Results_Index", 
    ".qrds"
  ), 
  nthr=nthr
) 

qsave(
  sample_data, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/Sex_PctMitoScaled/", 
    "sample_data", 
    ".qrds"
  ), 
  nthr=nthr
)

rm(list=ls())

source("~/ALS_Brain_Multiome.Rcfg")




  ### 2.0 ALSvsHC --------------------------------------------------------------



    ### 2.1 Load data ----------------------------------------------------------

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



    ## 2.2 DE Analysis ---------------------------------------------------------

DESeq_Sex_PctMito_Results_Index <- DESeq_Results_Index[,c(1:9)] 


DDS_list_Covs <- list() 
for (i in 1:nrow(DESeq_Sex_PctMito_Results_Index)){
  DDS_list_Covs[[i]] <- NA
}

rm(DESeq_Results_Index, i) 

for (i in 1:nrow(DESeq_Sex_PctMito_Results_Index)){ 
  tryCatch({
    DDS.tmp <- DDS_list[[i]] 
    DDS.tmp$Sex <- factor(DDS.tmp$Sex, levels = c("f", "m"))
    if(DESeq_Sex_PctMito_Results_Index$Comparison[i]=="ALSvsHC"){
      design(DDS.tmp) <-  ~ Sex + pctMito_RNA_scaled + Case
    }else{
      design(DDS.tmp) <- ~ Sex + pctMito_RNA_scaled + Rand 
    } 
    DDS_list_Covs[[i]] <- DDS.tmp 
    rm(DDS.tmp)
  }, error=function(e){
    message(
      paste0(
        "Design of Object # ", 
        i, 
        " could not be set!"
      )
    )
  }
  ) 
}


      # START HPC Cluster Outsourcing ------------------------------------------


      # Export Environment -----------------------------------------------------

do.call(
  qsavem,
  c(
    lapply(
      c(
        "DDS_list_Covs", 
        "DESeq_Sex_PctMito_Results_Index"
      ), 
      as.symbol
    ), 
    file=paste0(
      "../Data/tmp/", 
      "B3_2_DE_RNA_WNN_Pre", 
      ".qrdata"
    ), 
    nthr=nthr
  )
)         


      # Y_B3_2_RNA_DE_Sex_PctMito_ALS_Slurm.R ----------------------------------


      # Import HPC Environment  ------------------------------------------------

rm(
  i, 
  DDS_list, 
  DDS_list_Covs, 
  DESeq_Sex_PctMito_Results_Index
)


qreadm(
  paste0(
    "../Data/tmp/", 
    "B3_2_DE_RNA_WNN_Post", 
    ".qrdata"
  ), 
  nthr=nthr
)


      # END HPC Cluster Outsourcing --------------------------------------------



    ## 2.3 Collect DE stats ----------------------------------------------------


DESeq_Results <- list() 

for (i in 1:length(DDS_List)){
  
  if(i%%2==1){
    tryCatch({
      DESeq_Results[[i]] <- results(
        DDS_List[[i]], 
        name="Case_ALS_vs_HC", 
        alpha=0.05, 
        lfcThreshold = 0
      )
    },error=function(e){
      message(paste0("Results #", i, " , ALS, could not be generated")) 
      DESeq_Results[[i]] <<- NA  
    }) 
    
    
  }else{
    tryCatch({
      DESeq_Results[[i]] <- results(
        DDS_List[[i]], 
        name="Rand_R2_vs_R1", 
        alpha=0.05, 
        lfcThreshold = 0
      )
    },error=function(e){
      message(paste0("Results #", i, " , ALS/Rand, could not be generated")) 
      DESeq_Results[[i]] <<- NA  
    }) 
    
  }
  
}


DESeq_Results_Index <- DESeq_Sex_PctMito_Results_Index
rm(DESeq_Sex_PctMito_Results_Index)

DESeq_Results_Index$All_q_0.05 <- NA 
DESeq_Results_Index$Up_q_0.05 <- NA  
DESeq_Results_Index$Down_q_0.05 <- NA 
DESeq_Results_Index$nSamples <- NA 
DESeq_Results_Index$nSamples_HC <- NA  
DESeq_Results_Index$nSamples_ALS <- NA 
DESeq_Results_Index$TotalReads <- NA 
DESeq_Results_Index$MeanReads <- NA 
DESeq_Results_Index$NonzeroCounts <- NA 
DESeq_Results_Index$nFeatures <- NA 
DESeq_Results_Index$LowCounts <- NA 

sample_data <- colData(DDS_List[[1]])

for (i in 1:length(DDS_List)){
  if(!is.na(DESeq_Results_Index$Remarks[i])){
    
  }else{
    DESeq_Results_Index$NonzeroCounts[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[2], fixed("of "))[[1]][2], " with")[[1]][1])
    DESeq_Results_Index$Up_q_0.05[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[4], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_Index$Down_q_0.05[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[5], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_Index$LowCounts[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[7], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_Index$nSamples[i] <- dim(DDS_List[[i]]@assays@data$counts)[2] 
    DESeq_Results_Index$nSamples_HC[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_List[[i]]@assays@data$counts)])[1]
    DESeq_Results_Index$nSamples_ALS[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_List[[i]]@assays@data$counts)])[2]
    DESeq_Results_Index$TotalReads[i] <- sum(colSums(DDS_List[[i]]@assays@data$counts)) 
    DESeq_Results_Index$MeanReads[i] <- sum(colMeans(DDS_List[[i]]@assays@data$counts)) 
    DESeq_Results_Index$nFeatures[i] <- sum(rowMeans(DDS_List[[i]]@assays@data$counts)>0) 
  }
  
}

DESeq_Results_Index$All_q_0.05 <- DESeq_Results_Index$Up_q_0.05 + DESeq_Results_Index$Down_q_0.05 



    ## 2.4 Export data ---------------------------------------------------------

qsave(
  DDS_List, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALS/Sex_PctMitoScaled/", 
    "DDS_DESEqed_list", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALS/Sex_PctMitoScaled/", 
    "DESeq_Results", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_Index, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALS/Sex_PctMitoScaled/", 
    "DESeq_Results_Index", 
    ".qrds"
  ), 
  nthr=nthr
) 

qsave(
  sample_data, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALS/Sex_PctMitoScaled/", 
    "sample_data", 
    ".qrds"
  ), 
  nthr=nthr
)


rm(list=ls())
source("~/ALS_Brain_Multiome.Rcfg")
  



  ### 3.0 ALSFTDvsHC -----------------------------------------------------------



    ## 3.1 Load data -----------------------------------------------------------

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



    ## 3.2 DE Analysis ---------------------------------------------------------

DESeq_Sex_PctMito_Results_Index <- DESeq_Results_Index[,c(1:9)] 


DDS_list_Covs <- list() 
for (i in 1:nrow(DESeq_Sex_PctMito_Results_Index)){
  DDS_list_Covs[[i]] <- NA
}

rm(DESeq_Results_Index, i) 

for (i in 1:nrow(DESeq_Sex_PctMito_Results_Index)){ 
  tryCatch({
    DDS.tmp <- DDS_list[[i]] 
    DDS.tmp$Sex <- factor(DDS.tmp$Sex, levels = c("f", "m"))
    if(DESeq_Sex_PctMito_Results_Index$Comparison[i]=="ALSFTDvsHC"){
      design(DDS.tmp) <-  ~ Sex + pctMito_RNA_scaled + Case
    }else{
      design(DDS.tmp) <- ~ Sex + pctMito_RNA_scaled + Rand 
    } 
    DDS_list_Covs[[i]] <- DDS.tmp 
    rm(DDS.tmp)
  }, error=function(e){
    message(
      paste0(
        "Design of Object # ", 
        i, 
        " could not be set!"
      )
    )
  }
  ) 
}


      # START HPC Cluster Outsourcing ------------------------------------------


      # Export Environment -----------------------------------------------------

do.call(
  qsavem,
  c(
    lapply(
      c(
        "DDS_list_Covs", 
        "DESeq_Sex_PctMito_Results_Index"
      ), 
      as.symbol
    ), 
    file=paste0(
      "../Data/tmp/", 
      "B3_3_DE_RNA_WNN_Pre", 
      ".qrdata"
    ), 
    nthr=nthr
  )
)         


      # Y_B3_3_RNA_DE_Sex_PctMito_ALSFTD_Slurm.R -------------------------------


      # Import HPC Environment  ------------------------------------------------

rm(
  i, 
  DDS_list, 
  DDS_list_Covs, 
  DESeq_Sex_PctMito_Results_Index
)


qreadm(
  paste0(
    "../Data/tmp/", 
    "B3_3_DE_RNA_WNN_Post", 
    ".qrdata"
  ), 
  nthr=nthr
)


      # END HPC Cluster Outsourcing --------------------------------------------



    ## 3.3 Collect DE stats ----------------------------------------------------


DESeq_Results <- list() 

for (i in 1:length(DDS_List)){
  
  if(i%%2==1){
    tryCatch({
      DESeq_Results[[i]] <- results(
        DDS_List[[i]], 
        name="Case_ALS_FTD_vs_HC", 
        alpha=0.05, 
        lfcThreshold = 0
      )
    },error=function(e){
      message(paste0("Results #", i, " , ALSFTD, could not be generated")) 
      DESeq_Results[[i]] <<- NA  
    }) 
    
    
  }else{
    tryCatch({
      DESeq_Results[[i]] <- results(
        DDS_List[[i]], 
        name="Rand_R2_vs_R1", 
        alpha=0.05, 
        lfcThreshold = 0
      )
    },error=function(e){
      message(paste0("Results #", i, " , ALSFTD/Rand, could not be generated")) 
      DESeq_Results[[i]] <<- NA  
    }) 
    
  }
  
}


DESeq_Results_Index <- DESeq_Sex_PctMito_Results_Index
rm(DESeq_Sex_PctMito_Results_Index)

DESeq_Results_Index$All_q_0.05 <- NA 
DESeq_Results_Index$Up_q_0.05 <- NA  
DESeq_Results_Index$Down_q_0.05 <- NA 
DESeq_Results_Index$nSamples <- NA 
DESeq_Results_Index$nSamples_HC <- NA  
DESeq_Results_Index$nSamples_ALSFTD <- NA 
DESeq_Results_Index$TotalReads <- NA 
DESeq_Results_Index$MeanReads <- NA 
DESeq_Results_Index$NonzeroCounts <- NA 
DESeq_Results_Index$nFeatures <- NA 
DESeq_Results_Index$LowCounts <- NA 

sample_data <- colData(DDS_List[[1]])

for (i in 1:length(DDS_List)){
  if(!is.na(DESeq_Results_Index$Remarks[i])){
    
  }else{
    DESeq_Results_Index$NonzeroCounts[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[2], fixed("of "))[[1]][2], " with")[[1]][1])
    DESeq_Results_Index$Up_q_0.05[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[4], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_Index$Down_q_0.05[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[5], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_Index$LowCounts[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[7], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_Index$nSamples[i] <- dim(DDS_List[[i]]@assays@data$counts)[2] 
    DESeq_Results_Index$nSamples_HC[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_List[[i]]@assays@data$counts)])[1]
    DESeq_Results_Index$nSamples_ALSFTD[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_List[[i]]@assays@data$counts)])[2]
    DESeq_Results_Index$TotalReads[i] <- sum(colSums(DDS_List[[i]]@assays@data$counts)) 
    DESeq_Results_Index$MeanReads[i] <- sum(colMeans(DDS_List[[i]]@assays@data$counts)) 
    DESeq_Results_Index$nFeatures[i] <- sum(rowMeans(DDS_List[[i]]@assays@data$counts)>0) 
  }
  
}

DESeq_Results_Index$All_q_0.05 <- DESeq_Results_Index$Up_q_0.05 + DESeq_Results_Index$Down_q_0.05 



    ## 3.4 Export data ---------------------------------------------------------

qsave(
  DDS_List, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD/Sex_PctMitoScaled/", 
    "DDS_DESEqed_list", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD/Sex_PctMitoScaled/", 
    "DESeq_Results", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_Index, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD/Sex_PctMitoScaled/", 
    "DESeq_Results_Index", 
    ".qrds"
  ), 
  nthr=nthr
) 

qsave(
  sample_data, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD/Sex_PctMitoScaled/", 
    "sample_data", 
    ".qrds"
  ), 
  nthr=nthr
)


rm(list=ls())
source("~/ALS_Brain_Multiome.Rcfg")


  ### 4.0 ALSFTD_C9vsHC --------------------------------------------------------



    ## 4.1 Load data -----------------------------------------------------------

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



    ## 4.2 DE Analysis ---------------------------------------------------------

DESeq_Sex_PctMito_Results_Index <- DESeq_Results_Index[,c(1:9)] 


DDS_list_Covs <- list() 
for (i in 1:nrow(DESeq_Sex_PctMito_Results_Index)){
  DDS_list_Covs[[i]] <- NA
}

rm(DESeq_Results_Index, i) 

for (i in 1:nrow(DESeq_Sex_PctMito_Results_Index)){ 
  tryCatch({
    DDS.tmp <- DDS_list[[i]] 
    DDS.tmp$Sex <- factor(DDS.tmp$Sex, levels = c("f", "m"))
    if(DESeq_Sex_PctMito_Results_Index$Comparison[i]=="C9_ALSFTDvsHC"){
      design(DDS.tmp) <-  ~ Sex + pctMito_RNA_scaled + Case
    }else{
      design(DDS.tmp) <- ~ Sex + pctMito_RNA_scaled + Rand 
    } 
    DDS_list_Covs[[i]] <- DDS.tmp 
    rm(DDS.tmp)
  }, error=function(e){
    message(
      paste0(
        "Design of Object # ", 
        i, 
        " could not be set!"
      )
    )
  }
  ) 
}


      # START HPC Cluster Outsourcing ------------------------------------------


      # Export Environment -----------------------------------------------------

do.call(
  qsavem,
  c(
    lapply(
      c(
        "DDS_list_Covs", 
        "DESeq_Sex_PctMito_Results_Index"
      ), 
      as.symbol
    ), 
    file=paste0(
      "../Data/tmp/", 
      "B3_4_DE_RNA_WNN_Pre", 
      ".qrdata"
    ), 
    nthr=nthr
  )
)         


      # Y_B3_4_RNA_DE_Sex_PctMito_ALSFTD_Slurm.R -------------------------------


      # Import HPC Environment  ------------------------------------------------

rm(
  i, 
  DDS_list, 
  DDS_list_Covs, 
  DESeq_Sex_PctMito_Results_Index
)


qreadm(
  paste0(
    "../Data/tmp/", 
    "B3_4_DE_RNA_WNN_Post", 
    ".qrdata"
  ), 
  nthr=nthr
)


      # END HPC Cluster Outsourcing --------------------------------------------


    ## 4.3 Collect DE stats ----------------------------------------------------


DESeq_Results <- list() 

for (i in 1:length(DDS_List)){
  
  if(i%%2==1){
    tryCatch({
      DESeq_Results[[i]] <- results(
        DDS_List[[i]], 
        name="Case_ALS_FTD_vs_HC", 
        alpha=0.05, 
        lfcThreshold = 0
      )
    },error=function(e){
      message(paste0("Results #", i, " , ALSFTD, could not be generated")) 
      DESeq_Results[[i]] <<- NA  
    }) 
    
    
  }else{
    tryCatch({
      DESeq_Results[[i]] <- results(
        DDS_List[[i]], 
        name="Rand_R2_vs_R1", 
        alpha=0.05, 
        lfcThreshold = 0
      )
    },error=function(e){
      message(paste0("Results #", i, " , ALSFTD/Rand, could not be generated")) 
      DESeq_Results[[i]] <<- NA  
    }) 
    
  }
  
}


DESeq_Results_Index <- DESeq_Sex_PctMito_Results_Index
rm(DESeq_Sex_PctMito_Results_Index)

DESeq_Results_Index$All_q_0.05 <- NA 
DESeq_Results_Index$Up_q_0.05 <- NA  
DESeq_Results_Index$Down_q_0.05 <- NA 
DESeq_Results_Index$nSamples <- NA 
DESeq_Results_Index$nSamples_HC <- NA  
DESeq_Results_Index$nSamples_ALSFTD <- NA 
DESeq_Results_Index$TotalReads <- NA 
DESeq_Results_Index$MeanReads <- NA 
DESeq_Results_Index$NonzeroCounts <- NA 
DESeq_Results_Index$nFeatures <- NA 
DESeq_Results_Index$LowCounts <- NA 

sample_data <- colData(DDS_List[[1]])

for (i in 1:length(DDS_List)){
  if(!is.na(DESeq_Results_Index$Remarks[i])){
    
  }else{
    DESeq_Results_Index$NonzeroCounts[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[2], fixed("of "))[[1]][2], " with")[[1]][1])
    DESeq_Results_Index$Up_q_0.05[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[4], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_Index$Down_q_0.05[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[5], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_Index$LowCounts[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[7], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_Index$nSamples[i] <- dim(DDS_List[[i]]@assays@data$counts)[2] 
    DESeq_Results_Index$nSamples_HC[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_List[[i]]@assays@data$counts)])[1]
    DESeq_Results_Index$nSamples_ALSFTD[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_List[[i]]@assays@data$counts)])[2]
    DESeq_Results_Index$TotalReads[i] <- sum(colSums(DDS_List[[i]]@assays@data$counts)) 
    DESeq_Results_Index$MeanReads[i] <- sum(colMeans(DDS_List[[i]]@assays@data$counts)) 
    DESeq_Results_Index$nFeatures[i] <- sum(rowMeans(DDS_List[[i]]@assays@data$counts)>0) 
  }
  
}

DESeq_Results_Index$All_q_0.05 <- DESeq_Results_Index$Up_q_0.05 + DESeq_Results_Index$Down_q_0.05 



    ## 4.4 Export data ---------------------------------------------------------

qsave(
  DDS_List, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9/Sex_PctMitoScaled/", 
    "DDS_DESEqed_list", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9/Sex_PctMitoScaled/", 
    "DESeq_Results", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_Index, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9/Sex_PctMitoScaled/", 
    "DESeq_Results_Index", 
    ".qrds"
  ), 
  nthr=nthr
) 

qsave(
  sample_data, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9/Sex_PctMitoScaled/", 
    "sample_data", 
    ".qrds"
  ), 
  nthr=nthr
)


rm(list=ls())
source("~/ALS_Brain_Multiome.Rcfg")



  ### 5.0 ALSFTD_NonC9vsHC -----------------------------------------------------



    ## 5.1 Load data -----------------------------------------------------------

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



    ## 5.2 DE Analysis ---------------------------------------------------------

DESeq_Sex_PctMito_Results_Index <- DESeq_Results_Index[,c(1:9)] 


DDS_list_Covs <- list() 
for (i in 1:nrow(DESeq_Sex_PctMito_Results_Index)){
  DDS_list_Covs[[i]] <- NA
}

rm(DESeq_Results_Index, i) 

for (i in 1:nrow(DESeq_Sex_PctMito_Results_Index)){ 
  tryCatch({
    DDS.tmp <- DDS_list[[i]] 
    DDS.tmp$Sex <- factor(DDS.tmp$Sex, levels = c("f", "m"))
    if(DESeq_Sex_PctMito_Results_Index$Comparison[i]=="NonC9_ALSFTDvsHC"){
      design(DDS.tmp) <-  ~ Sex + pctMito_RNA_scaled + Case
    }else{
      design(DDS.tmp) <- ~ Sex + pctMito_RNA_scaled + Rand 
    } 
    DDS_list_Covs[[i]] <- DDS.tmp 
    rm(DDS.tmp)
  }, error=function(e){
    message(
      paste0(
        "Design of Object # ", 
        i, 
        " could not be set!"
      )
    )
  }
  ) 
}


      # START HPC Cluster Outsourcing ------------------------------------------


      # Export Environment -----------------------------------------------------

do.call(
  qsavem,
  c(
    lapply(
      c(
        "DDS_list_Covs", 
        "DESeq_Sex_PctMito_Results_Index"
      ), 
      as.symbol
    ), 
    file=paste0(
      "../Data/tmp/", 
      "B3_5_DE_RNA_WNN_Pre", 
      ".qrdata"
    ), 
    nthr=nthr
  )
)         


      # Y_B3_5_RNA_DE_Sex_PctMito_ALSFTD_NonC9_Slurm.R -------------------------


      # Import HPC Environment  ------------------------------------------------

rm(
  i, 
  DDS_list, 
  DDS_list_Covs, 
  DESeq_Sex_PctMito_Results_Index
)


qreadm(
  paste0(
    "../Data/tmp/", 
    "B3_5_DE_RNA_WNN_Post", 
    ".qrdata"
  ), 
  nthr=nthr
)




      # END HPC Cluster Outsourcing --------------------------------------------


    ## 5.3 Collect DE stats ----------------------------------------------------


DESeq_Results <- list() 

for (i in 1:length(DDS_List)){
  
  if(i%%2==1){
    tryCatch({
      DESeq_Results[[i]] <- results(
        DDS_List[[i]], 
        name="Case_ALS_FTD_vs_HC", 
        alpha=0.05, 
        lfcThreshold = 0
      )
    },error=function(e){
      message(paste0("Results #", i, " , ALSFTD, could not be generated")) 
      DESeq_Results[[i]] <<- NA  
    }) 
    
    
  }else{
    tryCatch({
      DESeq_Results[[i]] <- results(
        DDS_List[[i]], 
        name="Rand_R2_vs_R1", 
        alpha=0.05, 
        lfcThreshold = 0
      )
    },error=function(e){
      message(paste0("Results #", i, " , ALSFTD/Rand, could not be generated")) 
      DESeq_Results[[i]] <<- NA  
    }) 
    
  }
  
}


DESeq_Results_Index <- DESeq_Sex_PctMito_Results_Index
rm(DESeq_Sex_PctMito_Results_Index)

DESeq_Results_Index$All_q_0.05 <- NA 
DESeq_Results_Index$Up_q_0.05 <- NA  
DESeq_Results_Index$Down_q_0.05 <- NA 
DESeq_Results_Index$nSamples <- NA 
DESeq_Results_Index$nSamples_HC <- NA  
DESeq_Results_Index$nSamples_ALSFTD <- NA 
DESeq_Results_Index$TotalReads <- NA 
DESeq_Results_Index$MeanReads <- NA 
DESeq_Results_Index$NonzeroCounts <- NA 
DESeq_Results_Index$nFeatures <- NA 
DESeq_Results_Index$LowCounts <- NA 

sample_data <- colData(DDS_List[[1]])

for (i in 1:length(DDS_List)){
  if(!is.na(DESeq_Results_Index$Remarks[i])){
    
  }else{
    DESeq_Results_Index$NonzeroCounts[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[2], fixed("of "))[[1]][2], " with")[[1]][1])
    DESeq_Results_Index$Up_q_0.05[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[4], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_Index$Down_q_0.05[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[5], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_Index$LowCounts[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[7], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_Index$nSamples[i] <- dim(DDS_List[[i]]@assays@data$counts)[2] 
    DESeq_Results_Index$nSamples_HC[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_List[[i]]@assays@data$counts)])[1]
    DESeq_Results_Index$nSamples_ALSFTD[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_List[[i]]@assays@data$counts)])[2]
    DESeq_Results_Index$TotalReads[i] <- sum(colSums(DDS_List[[i]]@assays@data$counts)) 
    DESeq_Results_Index$MeanReads[i] <- sum(colMeans(DDS_List[[i]]@assays@data$counts)) 
    DESeq_Results_Index$nFeatures[i] <- sum(rowMeans(DDS_List[[i]]@assays@data$counts)>0) 
  }
  
}

DESeq_Results_Index$All_q_0.05 <- DESeq_Results_Index$Up_q_0.05 + DESeq_Results_Index$Down_q_0.05 



    ## 5.4 Export data ---------------------------------------------------------

qsave(
  DDS_List, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_NonC9/Sex_PctMitoScaled/", 
    "DDS_DESEqed_list", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_NonC9/Sex_PctMitoScaled/", 
    "DESeq_Results", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_Index, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_NonC9/Sex_PctMitoScaled/", 
    "DESeq_Results_Index", 
    ".qrds"
  ), 
  nthr=nthr
) 

qsave(
  sample_data, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_NonC9/Sex_PctMitoScaled/", 
    "sample_data", 
    ".qrds"
  ), 
  nthr=nthr
)


rm(list=ls())
source("~/ALS_Brain_Multiome.Rcfg")




  ### 6.0 ALSFTD_C9vsALSFTD_NonC9 ----------------------------------------------



    ## 6.1 Load data -----------------------------------------------------------

DDS_list <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9vsALSFTD_NonC9/", 
    "DDS_DESEqed_list", 
    ".qrds"
  ), 
  nthr=nthr
)


DESeq_Results_Index <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9vsALSFTD_NonC9/",  
    "DESeq_Results_Index", 
    ".qrds"
  ), 
  nthr=nthr
)    



    ## 6.2 DE Analysis ---------------------------------------------------------

DESeq_Sex_PctMito_Results_Index <- DESeq_Results_Index[,c(1:9)] 


DDS_list_Covs <- list() 
for (i in 1:nrow(DESeq_Sex_PctMito_Results_Index)){
  DDS_list_Covs[[i]] <- NA
}

rm(DESeq_Results_Index, i) 

for (i in 1:nrow(DESeq_Sex_PctMito_Results_Index)){ 
  tryCatch({
    DDS.tmp <- DDS_list[[i]] 
    DDS.tmp$Sex <- factor(DDS.tmp$Sex, levels = c("f", "m"))
    if(DESeq_Sex_PctMito_Results_Index$Comparison[i]=="ALSFTD_C9vsALSFTD_NonC9"){
      design(DDS.tmp) <-  ~ Sex + pctMito_RNA_scaled + Case_Type
    }else{
      design(DDS.tmp) <- ~ Sex + pctMito_RNA_scaled + Rand 
    } 
    DDS_list_Covs[[i]] <- DDS.tmp 
    rm(DDS.tmp)
  }, error=function(e){
    message(
      paste0(
        "Design of Object # ", 
        i, 
        " could not be set!"
      )
    )
  }
  ) 
}


      # START HPC Cluster Outsourcing ------------------------------------------


      # Export Environment -----------------------------------------------------

do.call(
  qsavem,
  c(
    lapply(
      c(
        "DDS_list_Covs", 
        "DESeq_Sex_PctMito_Results_Index"
      ), 
      as.symbol
    ), 
    file=paste0(
      "../Data/tmp/", 
      "B3_6_DE_RNA_WNN_Pre", 
      ".qrdata"
    ), 
    nthr=nthr
  )
)         


      # Y_B3_6_RNA_DE_Sex_PctMito_ALSFTD_C9vsALSFTD_NonC9_Slurm.R --------------


      # Import HPC Environment  ------------------------------------------------

rm(
  i, 
  DDS_list, 
  DDS_list_Covs, 
  DESeq_Sex_PctMito_Results_Index
)


qreadm(
  paste0(
    "../Data/tmp/", 
    "B3_6_DE_RNA_WNN_Post", 
    ".qrdata"
  ), 
  nthr=nthr
)





      # END HPC Cluster Outsourcing --------------------------------------------


      ## 6.3 Collect DE stats ----------------------------------------------------

DESeq_Results <- list() 

for (i in 1:length(DDS_List)){
  
  if(i%%2==1){
    tryCatch({
      DESeq_Results[[i]] <- results(
        DDS_List[[i]], 
        name="Case_Type_ALS_FTD_vs_C9_ALS_FTD", 
        alpha=0.05, 
        lfcThreshold = 0
      )
    },error=function(e){
      message(paste0("Results #", i, " , ALSFTD, could not be generated")) 
      DESeq_Results[[i]] <<- NA  
    }) 
    
    
  }else{
    tryCatch({
      DESeq_Results[[i]] <- results(
        DDS_List[[i]], 
        name="Rand_R2_vs_R1", 
        alpha=0.05, 
        lfcThreshold = 0
      )
    },error=function(e){
      message(paste0("Results #", i, " , ALSFTD/Rand, could not be generated")) 
      DESeq_Results[[i]] <<- NA  
    }) 
    
  }
  
}


DESeq_Results_Index <- DESeq_Sex_PctMito_Results_Index
rm(DESeq_Sex_PctMito_Results_Index)

DESeq_Results_Index$All_q_0.05 <- NA 
DESeq_Results_Index$Up_q_0.05 <- NA  
DESeq_Results_Index$Down_q_0.05 <- NA 
DESeq_Results_Index$nSamples <- NA 
DESeq_Results_Index$nSamples_NonC9_ALSFTD <- NA  
DESeq_Results_Index$nSamples_C9_ALSFTD<- NA 
DESeq_Results_Index$TotalReads <- NA 
DESeq_Results_Index$MeanReads <- NA 
DESeq_Results_Index$NonzeroCounts <- NA 
DESeq_Results_Index$nFeatures <- NA 
DESeq_Results_Index$LowCounts <- NA 

sample_data <- colData(DDS_List[[1]])

for (i in 1:length(DDS_List)){
  if(!is.na(DESeq_Results_Index$Remarks[i])){
    
  }else{
    DESeq_Results_Index$NonzeroCounts[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[2], fixed("of "))[[1]][2], " with")[[1]][1])
    DESeq_Results_Index$Up_q_0.05[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[4], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_Index$Down_q_0.05[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[5], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_Index$LowCounts[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[7], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_Index$nSamples[i] <- dim(DDS_List[[i]]@assays@data$counts)[2] 
    DESeq_Results_Index$nSamples_NonC9_ALSFTD[i] <- table(sample_data$Case_Type[sample_data$ID %in% colnames(DDS_List[[i]]@assays@data$counts)])[2]
    DESeq_Results_Index$nSamples_C9_ALSFTD[i] <- table(sample_data$Case_Type[sample_data$ID %in% colnames(DDS_List[[i]]@assays@data$counts)])[1]
    DESeq_Results_Index$TotalReads[i] <- sum(colSums(DDS_List[[i]]@assays@data$counts)) 
    DESeq_Results_Index$MeanReads[i] <- sum(colMeans(DDS_List[[i]]@assays@data$counts)) 
    DESeq_Results_Index$nFeatures[i] <- sum(rowMeans(DDS_List[[i]]@assays@data$counts)>0) 
  }
  
}

DESeq_Results_Index$All_q_0.05 <- DESeq_Results_Index$Up_q_0.05 + DESeq_Results_Index$Down_q_0.05 



      ## 6.4 Export data ---------------------------------------------------------

qsave(
  DDS_List, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9vsALSFTD_NonC9/Sex_PctMitoScaled/", 
    "DDS_DESEqed_list", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9vsALSFTD_NonC9/Sex_PctMitoScaled/", 
    "DESeq_Results", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_Index, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9vsALSFTD_NonC9/Sex_PctMitoScaled/", 
    "DESeq_Results_Index", 
    ".qrds"
  ), 
  nthr=nthr
) 

qsave(
  sample_data, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9vsALSFTD_NonC9/Sex_PctMitoScaled/", 
    "sample_data", 
    ".qrds"
  ), 
  nthr=nthr
)


rm(list=ls())



