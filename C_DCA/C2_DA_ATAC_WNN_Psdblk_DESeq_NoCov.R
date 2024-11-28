
  
  ### 0.0 Load libraries -------------------------------------------------------

source("~/ALS_Brain_Multiome.Rcfg")

library(qs) 
library(data.table)
library(Seurat)
library(tidyverse) 
library(DESeq2)



qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)




  ### 1.0 Load data ------------------------------------------------------------


      # M0_ATAC Seurat single-cell data object  --------------------------------

M0_ATAC <- qread(
  paste0(
    "../Data/SeuratObjects/",  
    "M0_ATAC",
    ".qrds"
  ),
  nthr=nthr 
)

M0_ATAC 



      # ATAC WNN Pseudobulk matrices -------------------------------------------

ATAC_Psdblk_WNN_Matrices <- qread(
  paste0(
    "../Data/DA/WNN/Pseudobulk/", 
    "ATAC_Psdblk_WNN_Matrices", 
    ".qrds"
  ), 
  nthr=nthr
)

ATAC_Psdblk_WNN_Matrices_Index <- qread(
  paste0(
    "../Data/DA/WNN/Pseudobulk/", 
    "ATAC_Psdblk_WNN_Matrices_Index", 
    ".qrds"
  ), 
  nthr=nthr
)




  ### 2.0 All Cases DESeq2 -----------------------------------------------------



    ## 2.1 Extract sample data -------------------------------------------------

sample_data <- M0_ATAC@misc$Sample_metrics 
sample_data$Sex <- M0_ATAC$Sex[match(sample_data$ID, M0_ATAC$ID)]
setequal(
  unique(paste0(M0_ATAC@misc$Sample_data$ID, "_", M0_ATAC@misc$Sample_data$Sex)), 
  paste0(sample_data$ID, "_", sample_data$Sex)
)
sample_data <- sample_data[,c(1,17:19, 2:16)]

sample_data$Rand <- M0_ATAC@misc$Sample_data$Rand1_AllCases[match(sample_data$ID, M0_ATAC@misc$Sample_data$ID)] 
table(is.na(sample_data$Rand))
table(sample_data$Case, sample_data$Rand)



    ## 2.3 Format grouping variables -------------------------------------------

      # Define Case_Type and Rand as factors and order levels (for DESeq2) -----

sample_data$Case <- factor(sample_data$Case, levels=c("HC", "ALS", "ALS_FTD"))
sample_data$Case_Type <- factor(sample_data$Case_Type, levels=c("HC", "ALS", "C9_ALS_FTD", "ALS_FTD")) 
sample_data$Rand <- factor(sample_data$Rand, levels=c("R1", "R2", "R3"))



    ## 2.4 DA Analysis  --------------------------------------------------------

      # Prepare data containers and counters -----------------------------------

dds_list <- list() 
DDS_list <- list() 
DESeq_Results_ALS <- list()   
DESeq_Results_ALSFTD <- list()   

nR <- nrow(ATAC_Psdblk_WNN_Matrices_Index)

DESeq_Results_Index <- data.frame(
  Index=rep(NA, nR*2), 
  CellTypeLevel=rep(NA, nR*2), 
  CellType=rep(NA, nR*2), 
  Assay=rep("ATAC", nR*2), 
  Test=rep("WALD", nR*2), 
  Covariates=rep("NoCov", nR*2), 
  SizeFactors=rep("Poscounts", nR*2),
  Comparison=rep(NA, nR*2), 
  Remarks=rep(NA, nR*2)
) 


for(i in 1:2*nrow(ATAC_Psdblk_WNN_Matrices_Index)){
  dds_list[[i]] <- NA 
  DDS_list[[i]] <- NA
  DESeq_Results_ALS[[i]] <- NA 
  DESeq_Results_ALSFTD[[i]] <- NA 
}

rm(nR, i) 


      # Main DESeq Loop  -------------------------------------------------------


      # START HPC Cluster Outsourcing ------------------------------------------


      # Export Environment -----------------------------------------------------

do.call(
  qsavem,
  c(
    lapply(
      c(
        "dds_list", 
        "DDS_list", 
        "DESeq_Results_ALS", 
        "DESeq_Results_ALSFTD", 
        "DESeq_Results_Index", 
        "ATAC_Psdblk_WNN_Matrices", 
        "ATAC_Psdblk_WNN_Matrices_Index", 
        "sample_data"
      ), 
      as.symbol
    ), 
    file=paste0(
      "../Data/tmp/", 
      "C2_1_DA_ATAC_WNN_Pre", 
      ".qrdata"
    ), 
    nthr=nthr
  )
)         


      # Y_C2_1_ATAC_DA_Slurm.R -------------------------------------------------


      # Import HPC Environment  ------------------------------------------------

rm(
  "dds_list", 
  "DDS_list", 
  "DESeq_Results_ALS", 
  "DESeq_Results_ALSFTD", 
  "DESeq_Results_Index", 
  "ATAC_Psdblk_WNN_Matrices", 
  "ATAC_Psdblk_WNN_Matrices_Index", 
  "sample_data"   
)



qreadm(
  paste0(
    "../Data/tmp/",
    "C2_1_DA_ATAC_WNN_Post", 
    ".qrdata"
  ), 
  nthr=nthr
)


      # END HPC Cluster Outsourcing --------------------------------------------



    ## 2.5 Collect DA stats ----------------------------------------------------
 
resultsNames(DDS_list[[1]])
resultsNames(DDS_list[[2]])


DESeq_Results_ALS <- list() 
DESeq_Results_ALSFTD <- list()  

for (i in 1:length(DDS_list)){
  
  if(i%%2==1){
    tryCatch({
      DESeq_Results_ALS[[i]] <- results(
        DDS_list[[i]], 
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
        DDS_list[[i]], 
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
        DDS_list[[i]], 
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
        DDS_list[[i]], 
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


for (i in 1:nrow(DESeq_Results_Index)){
  print(i)
  if(!is.na(DESeq_Results_Index$Remarks[i])){
    
  }else{
    DESeq_Results_Index$NonzeroCounts_ALS[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_ALS[[i]]))[2], fixed("of "))[[1]][2], " with")[[1]][1]) 
    DESeq_Results_Index$NonzeroCounts_ALSFTD[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_ALSFTD[[i]]))[2], fixed("of "))[[1]][2], " with")[[1]][1]) 
    
    DESeq_Results_Index$Up_q_0.05_ALS[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_ALS[[i]]))[4], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_Index$Up_q_0.05_ALSFTD[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_ALSFTD[[i]]))[4], fixed(": "))[[1]][2], ", ")[[1]][1]) 
    
    DESeq_Results_Index$Down_q_0.05_ALS[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_ALS[[i]]))[5], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_Index$Down_q_0.05_ALSFTD[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_ALSFTD[[i]]))[5], fixed(": "))[[1]][2], ", ")[[1]][1])
    
    DESeq_Results_Index$LowCounts_ALS[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_ALS[[i]]))[7], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_Index$LowCounts_ALSFTD[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results_ALSFTD[[i]]))[7], fixed(": "))[[1]][2], ", ")[[1]][1])
    
    DESeq_Results_Index$nSamples[i] <- dim(DDS_list[[i]]@assays@data$counts)[2] 
    
    
    DESeq_Results_Index$nSamples_HC[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_list[[i]]@assays@data$counts)])[1]
    DESeq_Results_Index$nSamples_ALS[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_list[[i]]@assays@data$counts)])[2]
    DESeq_Results_Index$nSamples_ALSFTD[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_list[[i]]@assays@data$counts)])[3]
    
    DESeq_Results_Index$TotalReads[i] <- sum(colSums(DDS_list[[i]]@assays@data$counts)) 
    DESeq_Results_Index$MeanReads[i] <- sum(colMeans(DDS_list[[i]]@assays@data$counts)) 
    DESeq_Results_Index$nFeatures[i] <- sum(rowMeans(DDS_list[[i]]@assays@data$counts)>0) 
  }
  
}

DESeq_Results_Index$All_q_0.05_ALS <- DESeq_Results_Index$Up_q_0.05_ALS + DESeq_Results_Index$Down_q_0.05_ALS 
DESeq_Results_Index$All_q_0.05_ALSFTD <- DESeq_Results_Index$Up_q_0.05_ALSFTD + DESeq_Results_Index$Down_q_0.05_ALSFTD 



    ## 2.6 Export data ---------------------------------------------------------

qsave(
  dds_list, 
  paste0(
    "../Data/DA/WNN/AllCase/", 
    "dds_list", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DDS_list, 
  paste0(
    "../Data/DA/WNN/AllCase/", 
    "DDS_DESEqed_list", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_ALS, 
  paste0(
    "../Data/DA/WNN/AllCase/", 
    "DESeq_Results_ALS", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_ALSFTD, 
  paste0(
    "../Data/DA/WNN/AllCase/", 
    "DESeq_Results_ALSFTD", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_Index, 
  paste0(
    "../Data/DA/WNN/AllCase/", 
    "DESeq_Results_Index", 
    ".qrds"
  ), 
  nthr=nthr
) 

qsave(
  sample_data, 
  paste0(
    "../Data/DA/WNN/AllCase/", 
    "sample_data", 
    ".qrds"
  ), 
  nthr=nthr
)


rm(dds_list, DDS_list, DESeq_Results_ALS, DESeq_Results_ALSFTD, DESeq_Results_Index, sample_data, i)
gc() 




  ### 3.0 ALSvsHC ATAC DA Analysis with DESeq2 ----------------------------------



    ## 3.1 Extract sample data and subset ALS and HC persons -------------------

sample_data <- M0_RNA@misc$Sample_metrics 
sample_data$Sex <- M0_RNA$Sex[match(sample_data$ID, M0_RNA$ID)]
setequal(
  unique(paste0(M0_RNA@misc$Sample_data$ID, "_", M0_RNA@misc$Sample_data$Sex)), 
  paste0(sample_data$ID, "_", sample_data$Sex)
)
sample_data <- sample_data[,c(1,17:19, 2:16)]

sample_data <- sample_data[sample_data$Case %in% c("ALS", "HC"),]

sample_data$Rand <- M0_RNA@misc$Sample_data$Rand1_ALS[match(sample_data$ID, M0_RNA@misc$Sample_data$ID)]

table(sample_data$Rand)
table(sample_data$Case, sample_data$Rand) 
table(sample_data$Sex, sample_data$Case)
table(sample_data$Sex, sample_data$Rand) 



    ## 3.3 Format grouping variables -------------------------------------------

        # Define Case and Sex as factors and order levels (for DESeq2) 

sample_data$Case <- factor(sample_data$Case, levels=c("HC", "ALS")) 
sample_data$Sex <- factor(sample_data$Sex, levels=c("f", "m"))
sample_data$Rand <- factor(sample_data$Rand, levels=c("R1", "R2"))



    ## 3.4 DA Analysis  --------------------------------------------------------

      # Prepare data containers and counters -----------------------------------

dds_list <- list() 
DDS_list <- list() 
DESeq_Results <- list()   

nR <- nrow(ATAC_Psdblk_WNN_Matrices_Index)

DESeq_Results_Index <- data.frame(
  Index=rep(NA, nR*2), 
  CellTypeLevel=rep(NA, nR*2), 
  CellType=rep(NA, nR*2), 
  Assay=rep("ATAC", nR*2), 
  Test=rep("WALD", nR*2), 
  Covariates=rep("NoCov", nR*2), 
  SizeFactors=rep("Poscounts", nR*2),
  Comparison=rep(NA, nR*2), 
  Remarks=rep(NA, nR*2)
) 


for(i in 1:2*nrow(ATAC_Psdblk_WNN_Matrices_Index)){
  dds_list[[i]] <- NA 
  DDS_list[[i]] <- NA
  DESeq_Results[[i]] <- NA 
}

rm(nR, i) 


      # Main DESeq Loop  -------------------------------------------------------

      # START HPC Cluster Outsourcing ------------------------------------------

      # Export Environment -----------------------------------------------------

do.call(
  qsavem,
  c(
    lapply(
      c(
        "dds_list", 
        "DDS_list", 
        "DESeq_Results", 
        "DESeq_Results_Index", 
        "ATAC_Psdblk_WNN_Matrices", 
        "ATAC_Psdblk_WNN_Matrices_Index", 
        "sample_data"
      ), 
      as.symbol
    ), 
    file=paste0(
      "../Data/tmp/", 
      "C2_2_DA_ATAC_WNN_Pre", 
      ".qrdata"
    ), 
    nthr=nthr
  )
)     

      # Y_C2_2_ATAC_DA_Slurm.R --------------------------------------------------

      # Import HPC Environment  ------------------------------------------------

rm(
  "dds_list", 
  "DDS_list", 
  "DESeq_Results", 
  "DESeq_Results_Index", 
  "nthr", 
  "ATAC_Psdblk_WNN_Matrices", 
  "ATAC_Psdblk_WNN_Matrices_Index", 
  "sample_data"   
)

source("~/ALS_Brain_Multiome.Rcfg") 

qreadm(
  paste0(
    "../Data/tmp/", 
    "C2_2_DA_ATAC_WNN_Post", 
    ".qrdata"
  ), 
  nthr=nthr
)


      # END HPC Cluster Outsourcing --------------------------------------------



    ## 3.5 Collect DA stats ----------------------------------------------------

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

for (i in 1:length(dds_list)){
  if(!is.na(DESeq_Results_Index$Remarks[i])){
    
  }else{
    DESeq_Results_Index$NonzeroCounts[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[2], fixed("of "))[[1]][2], " with")[[1]][1])
    DESeq_Results_Index$Up_q_0.05[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[4], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_Index$Down_q_0.05[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[5], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_Index$LowCounts[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[7], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_Index$nSamples[i] <- dim(DDS_list[[i]]@assays@data$counts)[2] 
    DESeq_Results_Index$nSamples_HC[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_list[[i]]@assays@data$counts)])[1]
    DESeq_Results_Index$nSamples_ALS[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_list[[i]]@assays@data$counts)])[2]
    DESeq_Results_Index$TotalReads[i] <- sum(colSums(DDS_list[[i]]@assays@data$counts)) 
    DESeq_Results_Index$MeanReads[i] <- sum(colMeans(DDS_list[[i]]@assays@data$counts)) 
    DESeq_Results_Index$nFeatures[i] <- sum(rowMeans(DDS_list[[i]]@assays@data$counts)>0) 
  }
  
}

DESeq_Results_Index$All_q_0.05 <- DESeq_Results_Index$Up_q_0.05 + DESeq_Results_Index$Down_q_0.05 



    ## 3.6 Export data ---------------------------------------------------------

qsave(
  dds_list, 
  paste0(
    "../Data/DA/WNN/ALS/", 
    "dds_list", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DDS_list, 
  paste0(
    "../Data/DA/WNN/ALS/", 
    "DDS_DESEqed_list", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results, 
  paste0(
    "../Data/DA/WNN/ALS/", 
    "DESeq_Results", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_Index, 
  paste0(
    "../Data/DA/WNN/ALS/", 
    "DESeq_Results_Index", 
    ".qrds"
  ), 
  nthr=nthr
) 

qsave(
  sample_data, 
  paste0(
    "../Data/DA/WNN/ALS/", 
    "sample_data", 
    ".qrds"
  ), 
  nthr=nthr
)


rm(list=ls()[-which(ls() %in% c("M0_RNA", "nthr", "ATAC_Psdblk_WNN_Matrices", "ATAC_Psdblk_WNN_Matrices_Index"))])




  ### 4.0 ALSFTDvsHC ATAC DA Analysis with DESeq2 -------------------------------



    ## 4.1 Extract sample data and subset ALSFTD and HC persons ----------------

sample_data <- M0_RNA@misc$Sample_metrics 
sample_data$Sex <- M0_RNA$Sex[match(sample_data$ID, M0_RNA$ID)]
setequal(
  unique(paste0(M0_RNA@misc$Sample_data$ID, "_", M0_RNA@misc$Sample_data$Sex)), 
  paste0(sample_data$ID, "_", sample_data$Sex)
)
sample_data <- sample_data[,c(1,17:19, 2:16)]

sample_data <- sample_data[sample_data$Case %in% c("ALS_FTD", "HC"),]

sample_data$Rand <- M0_RNA@misc$Sample_data$Rand1_ALSFTD[match(sample_data$ID, M0_RNA@misc$Sample_data$ID)]
    

table(sample_data$Rand)
table(sample_data$Case, sample_data$Rand) 
table(sample_data$Sex, sample_data$Case)
table(sample_data$Sex, sample_data$Rand) 



    ## 4.3 Format grouping variables -------------------------------------------

      # Define Case and Sex as factors and order levels (for DESeq2)

sample_data$Case <- factor(sample_data$Case, levels=c("HC", "ALS_FTD")) 
sample_data$Sex <- factor(sample_data$Sex, levels=c("f", "m"))
sample_data$Rand <- factor(sample_data$Rand, levels=c("R1", "R2"))



    ## 4.4 DA Analysis  --------------------------------------------------------

      # Prepare data containers and counters -----------------------------------

dds_list <- list() 
DDS_list <- list() 
DESeq_Results <- list()   

nR <- nrow(ATAC_Psdblk_WNN_Matrices_Index)

DESeq_Results_Index <- data.frame(
  Index=rep(NA, nR*2), 
  CellTypeLevel=rep(NA, nR*2), 
  CellType=rep(NA, nR*2), 
  Assay=rep("ATAC", nR*2), 
  Test=rep("WALD", nR*2), 
  Covariates=rep("NoCov", nR*2), 
  SizeFactors=rep("Poscounts", nR*2),
  Comparison=rep(NA, nR*2), 
  Remarks=rep(NA, nR*2)
) 


for(i in 1:2*nrow(ATAC_Psdblk_WNN_Matrices_Index)){
  dds_list[[i]] <- NA 
  DDS_list[[i]] <- NA
  DESeq_Results[[i]] <- NA 
}

rm(nR, i) 


      # Main DESeq Loop  -------------------------------------------------------


      # START HPC Cluster Outsourcing ------------------------------------------


      # Export Environment -----------------------------------------------------

do.call(
  qsavem,
  c(
    lapply(
      c(
        "dds_list", 
        "DDS_list", 
        "DESeq_Results", 
        "DESeq_Results_Index", 
        "ATAC_Psdblk_WNN_Matrices", 
        "ATAC_Psdblk_WNN_Matrices_Index", 
        "sample_data"
      ), 
      as.symbol
    ), 
    file=paste0(
      "../Data/tmp/", 
      "C2_3_DA_ATAC_WNN_Pre", 
      ".qrdata"
    ), 
    nthr=nthr
  )
)     


      # Y_C2_3_ATAC_DA_Slurm.R --------------------------------------------------


      # Import HPC Environment  ------------------------------------------------

rm(
  "dds_list", 
  "DDS_list", 
  "DESeq_Results", 
  "DESeq_Results_Index", 
  "nthr", 
  "ATAC_Psdblk_WNN_Matrices", 
  "ATAC_Psdblk_WNN_Matrices_Index", 
  "sample_data"   
)

source("~/ALS_Brain_Multiome.Rcfg")  

qreadm(
  paste0(
    "../Data/tmp/", 
    "",
    "C2_3_DA_ATAC_WNN_Post", 
    ".qrdata"
  ), 
  nthr=nthr 
)



  #### END HPC Cluster Outsourcing ---------------------------------------------

    ## 4.5 Collect DA stats ----------------------------------------------------

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

for (i in 1:length(dds_list)){
  if(!is.na(DESeq_Results_Index$Remarks[i])){
    
  }else{
    DESeq_Results_Index$NonzeroCounts[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[2], fixed("of "))[[1]][2], " with")[[1]][1])
    DESeq_Results_Index$Up_q_0.05[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[4], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_Index$Down_q_0.05[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[5], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_Index$LowCounts[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[7], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_Index$nSamples[i] <- dim(DDS_list[[i]]@assays@data$counts)[2] 
    DESeq_Results_Index$nSamples_HC[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_list[[i]]@assays@data$counts)])[1]
    DESeq_Results_Index$nSamples_ALS[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_list[[i]]@assays@data$counts)])[2]
    DESeq_Results_Index$TotalReads[i] <- sum(colSums(DDS_list[[i]]@assays@data$counts)) 
    DESeq_Results_Index$MeanReads[i] <- sum(colMeans(DDS_list[[i]]@assays@data$counts)) 
    DESeq_Results_Index$nFeatures[i] <- sum(rowMeans(DDS_list[[i]]@assays@data$counts)>0) 
  }
  
}

DESeq_Results_Index$All_q_0.05 <- DESeq_Results_Index$Up_q_0.05 + DESeq_Results_Index$Down_q_0.05 
rm(i)



    ## 4.6 Export data ---------------------------------------------------------

qsave(
  dds_list, 
  paste0(
    "../Data/DA/WNN/ALSFTD/", 
    "dds_list", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DDS_list, 
  paste0(
    "../Data/DA/WNN/ALSFTD/", 
    "DDS_DESEqed_list", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results, 
  paste0(
    "../Data/DA/WNN/ALSFTD/", 
    "DESeq_Results", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_Index, 
  paste0(
    "../Data/DA/WNN/ALSFTD/", 
    "DESeq_Results_Index", 
    ".qrds"
  ), 
  nthr=nthr
) 

qsave(
  sample_data, 
  paste0(
    "../Data/DA/WNN/ALSFTD/", 
    "sample_data", 
    ".qrds"
  ), 
  nthr=nthr
)


rm(list=ls()[-which(ls() %in% c("M0_RNA", "nthr", "ATAC_Psdblk_WNN_Matrices", "ATAC_Psdblk_WNN_Matrices_Index"), )])




  ### 5.0 ALSFTD_C9vsHC ATAC DA Analysis with DESeq2 ----------------------------



    ## 5.1 Extract sample data and subset C9_ALSFTD and HC persons -------------

sample_data <- M0_RNA@misc$Sample_metrics 
sample_data$Sex <- M0_RNA$Sex[match(sample_data$ID, M0_RNA$ID)]
setequal(
  unique(paste0(M0_RNA@misc$Sample_data$ID, "_", M0_RNA@misc$Sample_data$Sex)), 
  paste0(sample_data$ID, "_", sample_data$Sex)
)
sample_data <- sample_data[,c(1,17:19, 2:16)]


sample_data <- sample_data[sample_data$Case_Type %in% c("C9_ALS_FTD", "HC"),]


sample_data$Rand <- M0_RNA@misc$Sample_data$Rand1_ALSFTD_C9[match(sample_data$ID, M0_RNA@misc$Sample_data$ID)]

table(sample_data$Rand)
table(sample_data$Case, sample_data$Rand) 
table(sample_data$Sex, sample_data$Case)
table(sample_data$Sex, sample_data$Rand) 



    ## 5.3 Format grouping variables -------------------------------------------

# Define Case and Sex as factors and order levels (for DESeq2)

sample_data$Case <- factor(sample_data$Case, levels=c("HC", "ALS_FTD")) 
sample_data$Sex <- factor(sample_data$Sex, levels=c("f", "m"))
sample_data$Rand <- factor(sample_data$Rand, levels=c("R1", "R2"))



    ## 5.4 DA Analysis  --------------------------------------------------------

      # Prepare data containers and counters -----------------------------------

dds_list <- list() 
DDS_list <- list() 
DESeq_Results <- list()   

nR <- nrow(ATAC_Psdblk_WNN_Matrices_Index)

DESeq_Results_Index <- data.frame(
  Index=rep(NA, nR*2), 
  CellTypeLevel=rep(NA, nR*2), 
  CellType=rep(NA, nR*2), 
  Assay=rep("ATAC", nR*2), 
  Test=rep("WALD", nR*2), 
  Covariates=rep("NoCov", nR*2), 
  SizeFactors=rep("Poscounts", nR*2),
  Comparison=rep(NA, nR*2), 
  Remarks=rep(NA, nR*2)
) 


for(i in 1:2*nrow(ATAC_Psdblk_WNN_Matrices_Index)){
  dds_list[[i]] <- NA 
  DDS_list[[i]] <- NA
  DESeq_Results[[i]] <- NA 
}

rm(nR, i) 


      # Main DESeq Loop  -------------------------------------------------------


      # START HPC Cluster Outsourcing ---------------------------------------


      # Export Environment -----------------------------------------------------

do.call(
  qsavem,
  c(
    lapply(
      c(
        "dds_list", 
        "DDS_list", 
        "DESeq_Results", 
        "DESeq_Results_Index", 
        "ATAC_Psdblk_WNN_Matrices", 
        "ATAC_Psdblk_WNN_Matrices_Index", 
        "sample_data"
      ), 
      as.symbol
    ), 
    file=paste0(
      "../Data/tmp/", 
      "C2_4_DA_ATAC_WNN_Pre", 
      ".qrdata"
    ), 
    nthr=nthr
  )
)    


      # Y_C2_4_ATAC_DA_Slurm.R --------------------------------------------------

      # Import HPC Environment  ------------------------------------------------

rm(
  "dds_list", 
  "DDS_list", 
  "DESeq_Results", 
  "DESeq_Results_Index", 
  "nthr", 
  "ATAC_Psdblk_WNN_Matrices", 
  "ATAC_Psdblk_WNN_Matrices_Index", 
  "sample_data"   
)

source("~/ALS_Brain_Multiome.Rcfg")  

qreadm(
  paste0(
    "../Data/tmp/", 
    "",
    "C2_4_DA_ATAC_WNN_Post", 
    ".qrdata"
  ), 
  nthr=nthr
)



      #### END HPC Cluster Outsourcing -----------------------------------------


    ## 5.5 Collect DA stats ----------------------------------------------------

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

for (i in 1:length(dds_list)){
  if(!is.na(DESeq_Results_Index$Remarks[i])){
    
  }else{
    DESeq_Results_Index$NonzeroCounts[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[2], fixed("of "))[[1]][2], " with")[[1]][1])
    DESeq_Results_Index$Up_q_0.05[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[4], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_Index$Down_q_0.05[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[5], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_Index$LowCounts[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[7], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_Index$nSamples[i] <- dim(DDS_list[[i]]@assays@data$counts)[2] 
    DESeq_Results_Index$nSamples_HC[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_list[[i]]@assays@data$counts)])[1]
    DESeq_Results_Index$nSamples_ALS[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_list[[i]]@assays@data$counts)])[2]
    DESeq_Results_Index$TotalReads[i] <- sum(colSums(DDS_list[[i]]@assays@data$counts)) 
    DESeq_Results_Index$MeanReads[i] <- sum(colMeans(DDS_list[[i]]@assays@data$counts)) 
    DESeq_Results_Index$nFeatures[i] <- sum(rowMeans(DDS_list[[i]]@assays@data$counts)>0) 
  }
  
}

DESeq_Results_Index$All_q_0.05 <- DESeq_Results_Index$Up_q_0.05 + DESeq_Results_Index$Down_q_0.05 
rm(i)



    ## 5.6 Export data ---------------------------------------------------------

qsave(
  dds_list, 
  paste0(
    "../Data/DA/WNN/ALSFTD_C9/", 
    "dds_list", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DDS_list, 
  paste0(
    "../Data/DA/WNN/ALSFTD_C9/", 
    "DDS_DESEqed_list", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results, 
  paste0(
    "../Data/DA/WNN/ALSFTD_C9/", 
    "DESeq_Results", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_Index, 
  paste0(
    "../Data/DA/WNN/ALSFTD_C9/", 
    "DESeq_Results_Index", 
    ".qrds"
  ), 
  nthr=nthr
) 

qsave(
  sample_data, 
  paste0(
    "../Data/DA/WNN/ALSFTD_C9/", 
    "sample_data", 
    ".qrds"
  ), 
  nthr=nthr
)



rm(list=ls()[-which(ls() %in% c("M0_RNA", "nthr", "ATAC_Psdblk_WNN_Matrices", "ATAC_Psdblk_WNN_Matrices_Index"), )])


  ### 6.0 ALSFTD_nonC9vsHC ATAC DA Analysis with DESeq2 ----------------------------



    ## 6.1 Extract sample data and subset C9_ALSFTD and HC persons -------------

sample_data <- M0_RNA@misc$Sample_metrics 
sample_data$Sex <- M0_RNA$Sex[match(sample_data$ID, M0_RNA$ID)]
setequal(
  unique(paste0(M0_RNA@misc$Sample_data$ID, "_", M0_RNA@misc$Sample_data$Sex)), 
  paste0(sample_data$ID, "_", sample_data$Sex)
)
sample_data <- sample_data[,c(1,17:19, 2:16)]

table(sample_data$Case_Type)
sample_data <- sample_data[sample_data$Case_Type %in% c("ALS_FTD", "HC"),]


sample_data$Rand <- M0_RNA@misc$Sample_data$Rand1_ALSFTD_NonC9[match(sample_data$ID, M0_RNA@misc$Sample_data$ID)]

table(sample_data$Rand)
table(sample_data$Case, sample_data$Rand) 
table(sample_data$Sex, sample_data$Case)
table(sample_data$Sex, sample_data$Rand) 



    ## 6.3 Format grouping variables -------------------------------------------

      # Define Case and Sex as factors and order levels (for DESeq2)

sample_data$Case <- factor(sample_data$Case, levels=c("HC", "ALS_FTD")) 
sample_data$Sex <- factor(sample_data$Sex, levels=c("f", "m"))
sample_data$Rand <- factor(sample_data$Rand, levels=c("R1", "R2"))



    ## 6.4 DA Analysis  --------------------------------------------------------

      # Prepare data containers and counters -----------------------------------

dds_list <- list() 
DDS_list <- list() 
DESeq_Results <- list()   

nR <- nrow(ATAC_Psdblk_WNN_Matrices_Index)

DESeq_Results_Index <- data.frame(
  Index=rep(NA, nR*2), 
  CellTypeLevel=rep(NA, nR*2), 
  CellType=rep(NA, nR*2), 
  Assay=rep("ATAC", nR*2), 
  Test=rep("WALD", nR*2), 
  Covariates=rep("NoCov", nR*2), 
  SizeFactors=rep("Poscounts", nR*2),
  Comparison=rep(NA, nR*2), 
  Remarks=rep(NA, nR*2)
) 


for(i in 1:2*nrow(ATAC_Psdblk_WNN_Matrices_Index)){
  dds_list[[i]] <- NA 
  DDS_list[[i]] <- NA
  DESeq_Results[[i]] <- NA 
}

rm(nR, i) 


      # Main DESeq Loop  -------------------------------------------------------


      # START HPC Cluster Outsourcing ------------------------------------------


      # Export Environment -----------------------------------------------------

do.call(
  qsavem,
  c(
    lapply(
      c(
        "dds_list", 
        "DDS_list", 
        "DESeq_Results", 
        "DESeq_Results_Index", 
        "ATAC_Psdblk_WNN_Matrices", 
        "ATAC_Psdblk_WNN_Matrices_Index", 
        "sample_data"
      ), 
      as.symbol
    ), 
    file=paste0(
      "../Data/tmp/", 
      "C2_5_DA_ATAC_WNN_Pre", 
      ".qrdata"
    ), 
    nthr=nthr
  )
)    


      # Y_C2_5_ATAC_DA_Slurm.R --------------------------------------------------

      # Import HPC Environment  ------------------------------------------------

rm(
  "dds_list", 
  "DDS_list", 
  "DESeq_Results", 
  "DESeq_Results_Index", 
  "nthr", 
  "ATAC_Psdblk_WNN_Matrices", 
  "ATAC_Psdblk_WNN_Matrices_Index", 
  "sample_data"   
)

source("~/ALS_Brain_Multiome.Rcfg")  

qreadm(
  paste0(
    "../Data/tmp/", 
    "",
    "C2_5_DA_ATAC_WNN_Post", 
    ".qrdata"
  ), 
  nthr=36
)


      # END HPC Cluster Outsourcing --------------------------------------------



    ## 6.5 Collect DA stats ----------------------------------------------------

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

for (i in 1:length(dds_list)){
  if(!is.na(DESeq_Results_Index$Remarks[i])){
    
  }else{
    DESeq_Results_Index$NonzeroCounts[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[2], fixed("of "))[[1]][2], " with")[[1]][1])
    DESeq_Results_Index$Up_q_0.05[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[4], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_Index$Down_q_0.05[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[5], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_Index$LowCounts[i] <- as.numeric(str_split(str_split(capture.output(summary(DESeq_Results[[i]]))[7], fixed(": "))[[1]][2], ", ")[[1]][1])
    DESeq_Results_Index$nSamples[i] <- dim(DDS_list[[i]]@assays@data$counts)[2] 
    DESeq_Results_Index$nSamples_HC[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_list[[i]]@assays@data$counts)])[1]
    DESeq_Results_Index$nSamples_ALS[i] <- table(sample_data$Case[sample_data$ID %in% colnames(DDS_list[[i]]@assays@data$counts)])[2]
    DESeq_Results_Index$TotalReads[i] <- sum(colSums(DDS_list[[i]]@assays@data$counts)) 
    DESeq_Results_Index$MeanReads[i] <- sum(colMeans(DDS_list[[i]]@assays@data$counts)) 
    DESeq_Results_Index$nFeatures[i] <- sum(rowMeans(DDS_list[[i]]@assays@data$counts)>0) 
  }
  
}

DESeq_Results_Index$All_q_0.05 <- DESeq_Results_Index$Up_q_0.05 + DESeq_Results_Index$Down_q_0.05 
rm(i)



    ## 6.6 Export data ---------------------------------------------------------

qsave(
  dds_list, 
  paste0(
    "../Data/DA/WNN/ALSFTD_NonC9/", 
    "dds_list", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DDS_list, 
  paste0(
    "../Data/DA/WNN/ALSFTD_NonC9/", 
    "DDS_DESEqed_list", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results, 
  paste0(
    "../Data/DA/WNN/ALSFTD_NonC9/", 
    "DESeq_Results", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_Index, 
  paste0(
    "../Data/DA/WNN/ALSFTD_NonC9/", 
    "DESeq_Results_Index", 
    ".qrds"
  ), 
  nthr=nthr
) 

qsave(
  sample_data, 
  paste0(
    "../Data/DA/WNN/ALSFTD_NonC9/", 
    "sample_data", 
    ".qrds"
  ), 
  nthr=nthr
)



rm(list=ls())
sessionInfo()
 

