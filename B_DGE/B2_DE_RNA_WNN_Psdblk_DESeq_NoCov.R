
  ### 0.0 Load libraries -------------------------------------------------------

source("~/ALS_Brain_Multiome.Rcfg")

library(qs) 
library(data.table)
library(Seurat)
library(tidyverse) 
library(DESeq2)


qs::set_trust_promises(TRUE) 




  ### 1.0 Load data ------------------------------------------------------------


      # M0_RNA Seurat single-cell data object  ---------------------------------

M0_RNA <- qread(
  paste0(
    "../Data/SeuratObjects/",  
    "M0_RNA",
    ".qrds"
  ),
  nthr=nthr 
)

M0_RNA 


      # RNA WNN Pseudobulk matrices --------------------------------------------


RNA_Psdblk_WNN_Matrices <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/", 
    "RNA_Psdblk_WNN_Matrices", 
    ".qrds"
  ), 
  nthr=nthr
)

RNA_Psdblk_WNN_Matrices_Index <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/", 
    "RNA_Psdblk_WNN_Matrices_Index", 
    ".qrds"
  ), 
  nthr=nthr
)




  ### 2.0 All Cases DESeq2 -----------------------------------------------------



    ## 2.1 Extract sample data -------------------------------------------------

sample_data <- M0_RNA@misc$Sample_metrics 
sample_data$Sex <- M0_RNA$Sex[match(sample_data$ID, M0_RNA$ID)]
setequal(
  unique(paste0(M0_RNA@misc$Sample_data$ID, "_", M0_RNA@misc$Sample_data$Sex)), 
  paste0(sample_data$ID, "_", sample_data$Sex)
)
sample_data <- sample_data[,c(1,17:19, 2:16)]

sample_data$Rand <- M0_RNA@misc$Sample_data$Rand1_AllCases[match(sample_data$ID, M0_RNA@misc$Sample_data$ID)] 
table(is.na(sample_data$Rand))
table(sample_data$Case, sample_data$Rand)



    ## 2.2 Format grouping variables -------------------------------------------

      # Define Case_Type and Rand as factors and order levels (for DESeq2) -----

sample_data$Case <- factor(sample_data$Case, levels=c("HC", "ALS", "ALS_FTD"))
sample_data$Case_Type <- factor(sample_data$Case_Type, levels=c("HC", "ALS", "C9_ALS_FTD", "ALS_FTD")) 
sample_data$Rand <- factor(sample_data$Rand, levels=c("R1", "R2", "R3"))



    ## 2.3 DE Analysis  --------------------------------------------------------

      # Prepare data containers and counters -----------------------------------

dds_list <- list() 
DDS_list <- list() 
DESeq_Results <- list()   

nR <- nrow(RNA_Psdblk_WNN_Matrices_Index)

DESeq_Results_Index <- data.frame(
  Index=rep(NA, nR*2), 
  CellTypeLevel=rep(NA, nR*2), 
  CellType=rep(NA, nR*2), 
  Assay=rep("RNA", nR*2), 
  Test=rep("WALD", nR*2), 
  Covariates=rep("NoCov", nR*2), 
  SizeFactors=rep("Ratio", nR*2),
  Comparison=rep(NA, nR*2), 
  Remarks=rep(NA, nR*2)
) 


for(i in 1:2*nrow(RNA_Psdblk_WNN_Matrices_Index)){
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
        "RNA_Psdblk_WNN_Matrices", 
        "RNA_Psdblk_WNN_Matrices_Index", 
        "sample_data"
      ), 
      as.symbol
    ), 
    file=paste0(
      "../Data/tmp/", 
      "B2_1_DE_RNA_WNN_Pre", 
      ".qrdata"
    ), 
    nthr=nthr
  )
)         


      # Y_B2_1_RNA_DE_Slurm.R --------------------------------------------------


      # Import HPC Environment  ------------------------------------------------

rm(
  "dds_list", 
  "DDS_list", 
  "DESeq_Results", 
  "DESeq_Results_Index", 
  "RNA_Psdblk_WNN_Matrices", 
  "RNA_Psdblk_WNN_Matrices_Index", 
  "sample_data"   
)

source("~/ALS_Brain_Multiome.Rcfg") 

qreadm(
  paste0(
    "../Data/tmp/", 
    "B2_1_DE_RNA_WNN_Post", 
    ".qrdata"
  ), 
  nthr=nthr
)



      # END HPC Cluster Outsourcing --------------------------------------------



    ## 2.4 Collect DE stats ----------------------------------------------------

resultsNames(DDS_list[[1]])
resultsNames(DDS_list[[2]])

rm(DESeq_Results)

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
      message(paste0("Results #", i, " , ALS, could not be genarated")) 
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
      message(paste0("Results #", i, " , ALSFTD, could not be genarated")) 
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
      message(paste0("Results #", i, " , ALS/Rand, could not be genarated")) 
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
      message(paste0("Results #", i, " , ALSFTD/Rand, could not be genarated")) 
      DESeq_Results_ALS[[i]] <<- NA  
    })
  }
  
}


DESeq_Results_Index$All_q_0.05_ALS <- NA 
DESeq_Results_Index$All_q_0.05_ALSFTD <- NA 

DESeq_Results_Index$Up_q_0.05_ALS <- NA   
DESeq_Results_Index$Up_q_0.05_ALSFTD <- NA  

DESeq_Results_Index$Down_q_0.05_ALS <- NA 
DESeq_Results_Index$Down_q_0.05_FTD <- NA 

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


for (i in 1:length(DDS_list)){
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



    ## 2.5 Export data ---------------------------------------------------------

qsave(
  dds_list, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/", 
    "dds_list", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DDS_list, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/", 
    "DDS_DESEqed_list", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_ALS, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/", 
    "DESeq_Results_ALS", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_ALSFTD, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/", 
    "DESeq_Results_ALSFTD", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_Index, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/", 
    "DESeq_Results_Index", 
    ".qrds"
  ), 
  nthr=nthr
) 

qsave(
  sample_data, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/", 
    "sample_data", 
    ".qrds"
  ), 
  nthr=nthr
)

{
  rm(list=ls())

  source("~/ALS_Brain_Multiome.Rcfg")
  
  RNA_Psdblk_WNN_Matrices <- qread(
    paste0(
      "../Data/DE/WNN/Pseudobulk/DESeq2/", 
      "RNA_Psdblk_WNN_Matrices", 
      ".qrds"
    ), 
    nthr=nthr
  )
  
  RNA_Psdblk_WNN_Matrices_Index <- qread(
    paste0(
      "../Data/DE/WNN/Pseudobulk/DESeq2/", 
      "RNA_Psdblk_WNN_Matrices_Index", 
      ".qrds"
    ), 
    nthr=nthr
  )
}




  ### 3.0 ALSvsHC RNA DE Analysis with DESeq2 ----------------------------------



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



    ## 3.2 Format grouping variables -------------------------------------------

# Define Case and Sex as factors and order levels (for DESeq2) 

sample_data$Case <- factor(sample_data$Case, levels=c("HC", "ALS")) 
sample_data$Sex <- factor(sample_data$Sex, levels=c("f", "m"))
sample_data$Rand <- factor(sample_data$Rand, levels=c("R1", "R2"))



    ## 3.3 DE Analysis  --------------------------------------------------------

      # Prepare data containers and counters -----------------------------------

dds_list <- list() 
DDS_list <- list() 
DESeq_Results <- list()   

nR <- nrow(RNA_Psdblk_WNN_Matrices_Index)

DESeq_Results_Index <- data.frame(
  Index=rep(NA, nR*2), 
  CellTypeLevel=rep(NA, nR*2), 
  CellType=rep(NA, nR*2), 
  Assay=rep("RNA", nR*2), 
  Test=rep("WALD", nR*2), 
  Covariates=rep("NoCov", nR*2), 
  SizeFactors=rep("Ratio", nR*2),
  Comparison=rep(NA, nR*2), 
  Remarks=rep(NA, nR*2)
) 


for(i in 1:2*nrow(RNA_Psdblk_WNN_Matrices_Index)){
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
        "RNA_Psdblk_WNN_Matrices", 
        "RNA_Psdblk_WNN_Matrices_Index", 
        "sample_data"
      ), 
      as.symbol
    ), 
    file=paste0(
      "../Data/tmp/", 
      "B2_2_DE_RNA_WNN_Pre", 
      ".qrdata"
    ), 
    nthr=nthr
  )
)     

      # Y_B2_2_RNA_DE_Slurm.R --------------------------------------------------

      # Import HPC Environment  ------------------------------------------------

rm(
  "dds_list", 
  "DDS_list", 
  "DESeq_Results", 
  "DESeq_Results_Index", 
  "nthr", 
  "RNA_Psdblk_WNN_Matrices", 
  "RNA_Psdblk_WNN_Matrices_Index", 
  "sample_data"   
)

source("~/ALS_Brain_Multiome.Rcfg") 

qreadm(
  paste0(
    "../Data/tmp/", 
    "B2_2_DE_RNA_WNN_Post", 
    ".qrdata"
  ), 
  nthr=nthr
)



      # END HPC Cluster Outsourcing --------------------------------------------



    ## 3.4 Collect DE stats ----------------------------------------------------

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



    ## 3.5 Export data ---------------------------------------------------------

qsave(
  dds_list, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALS/", 
    "dds_list", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DDS_list, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALS/", 
    "DDS_DESEqed_list", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALS/", 
    "DESeq_Results", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_Index, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALS/", 
    "DESeq_Results_Index", 
    ".qrds"
  ), 
  nthr=nthr
) 

qsave(
  sample_data, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALS/", 
    "sample_data", 
    ".qrds"
  ), 
  nthr=nthr
)


{
  rm(list=ls())
  source("~/ALS_Brain_Multiome.Rcfg")
  RNA_Psdblk_WNN_Matrices <- qread(
    paste0(
      "../Data/DE/WNN/Pseudobulk/", 
      "RNA_Psdblk_WNN_Matrices", 
      ".qrds"
    ), 
    nthr=nthr
  )
  
  RNA_Psdblk_WNN_Matrices_Index <- qread(
    paste0(
      "../Data/DE/WNN/Pseudobulk/", 
      "RNA_Psdblk_WNN_Matrices_Index", 
      ".qrds"
    ), 
    nthr=nthr
  )
}




  ### 4.0 ALSFTDvsHC RNA DE Analysis with DESeq2 -------------------------------



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



    ## 4.2 Format grouping variables -------------------------------------------

# Define Case and Sex as factors and order levels (for DESeq2)

sample_data$Case <- factor(sample_data$Case, levels=c("HC", "ALS_FTD")) 
sample_data$Sex <- factor(sample_data$Sex, levels=c("f", "m"))
sample_data$Rand <- factor(sample_data$Rand, levels=c("R1", "R2"))



    ## 4.3 DE Analysis  --------------------------------------------------------

      # Prepare data containers and counters -----------------------------------

dds_list <- list() 
DDS_list <- list() 
DESeq_Results <- list()   

nR <- nrow(RNA_Psdblk_WNN_Matrices_Index)

DESeq_Results_Index <- data.frame(
  Index=rep(NA, nR*2), 
  CellTypeLevel=rep(NA, nR*2), 
  CellType=rep(NA, nR*2), 
  Assay=rep("RNA", nR*2), 
  Test=rep("WALD", nR*2), 
  Covariates=rep("NoCov", nR*2), 
  SizeFactors=rep("Ratio", nR*2),
  Comparison=rep(NA, nR*2), 
  Remarks=rep(NA, nR*2)
) 


for(i in 1:2*nrow(RNA_Psdblk_WNN_Matrices_Index)){
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
        "RNA_Psdblk_WNN_Matrices", 
        "RNA_Psdblk_WNN_Matrices_Index", 
        "sample_data"
      ), 
      as.symbol
    ), 
    file=paste0(
      "../Data/tmp/", 
      "B2_3_DE_RNA_WNN_Pre", 
      ".qrdata"
    ), 
    nthr=nthr
  )
)     


      # Y_B2_3_RNA_DE_Slurm.R --------------------------------------------------


      # Import HPC Environment  ------------------------------------------------

rm(
  "dds_list", 
  "DDS_list", 
  "DESeq_Results", 
  "DESeq_Results_Index", 
  "nthr", 
  "RNA_Psdblk_WNN_Matrices", 
  "RNA_Psdblk_WNN_Matrices_Index", 
  "sample_data"   
)

source("~/ALS_Brain_Multiome.Rcfg")  

qreadm(
  paste0(
    "../Data/tmp/", 
    "",
    "B2_3_DE_RNA_WNN_Post", 
    ".qrdata"
  ), 
  nthr=nthr 
)


      # END HPC Cluster Outsourcing --------------------------------------------

    ## 4.4 Collect DE stats ----------------------------------------------------

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



    ## 4.5 Export data ---------------------------------------------------------

qsave(
  dds_list, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD/", 
    "dds_list", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DDS_list, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD/", 
    "DDS_DESEqed_list", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD/", 
    "DESeq_Results", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_Index, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD/", 
    "DESeq_Results_Index", 
    ".qrds"
  ), 
  nthr=nthr
) 

qsave(
  sample_data, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD/", 
    "sample_data", 
    ".qrds"
  ), 
  nthr=nthr
)


{
  rm(list=ls())
  source("~/ALS_Brain_Multiome.Rcfg")
  RNA_Psdblk_WNN_Matrices <- qread(
    paste0(
      "../Data/DE/WNN/Pseudobulk/", 
      "RNA_Psdblk_WNN_Matrices", 
      ".qrds"
    ), 
    nthr=nthr
  )
  
  RNA_Psdblk_WNN_Matrices_Index <- qread(
    paste0(
      "../Data/DE/WNN/Pseudobulk/", 
      "RNA_Psdblk_WNN_Matrices_Index", 
      ".qrds"
    ), 
    nthr=nthr
  )
}




  ### 5.0 ALSFTD_C9vsHC RNA DE Analysis with DESeq2 ----------------------------



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



    ## 5.2 Format grouping variables -------------------------------------------

# Define Case and Sex as factors and order levels (for DESeq2)

sample_data$Case <- factor(sample_data$Case, levels=c("HC", "ALS_FTD")) 
sample_data$Sex <- factor(sample_data$Sex, levels=c("f", "m"))
sample_data$Rand <- factor(sample_data$Rand, levels=c("R1", "R2"))



    ## 5.3 DE Analysis  --------------------------------------------------------

      # Prepare data containers and counters -----------------------------------

dds_list <- list() 
DDS_list <- list() 
DESeq_Results <- list()   

nR <- nrow(RNA_Psdblk_WNN_Matrices_Index)

DESeq_Results_Index <- data.frame(
  Index=rep(NA, nR*2), 
  CellTypeLevel=rep(NA, nR*2), 
  CellType=rep(NA, nR*2), 
  Assay=rep("RNA", nR*2), 
  Test=rep("WALD", nR*2), 
  Covariates=rep("NoCov", nR*2), 
  SizeFactors=rep("Ratio", nR*2),
  Comparison=rep(NA, nR*2), 
  Remarks=rep(NA, nR*2)
) 


for(i in 1:2*nrow(RNA_Psdblk_WNN_Matrices_Index)){
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
        "RNA_Psdblk_WNN_Matrices", 
        "RNA_Psdblk_WNN_Matrices_Index", 
        "sample_data"
      ), 
      as.symbol
    ), 
    file=paste0(
      "../Data/tmp/", 
      "B2_4_DE_RNA_WNN_Pre", 
      ".qrdata"
    ), 
    nthr=nthr
  )
)    


      # Y_B2_4_RNA_DE_Slurm.R --------------------------------------------------

      # Import HPC Environment  ------------------------------------------------

rm(
  "dds_list", 
  "DDS_list", 
  "DESeq_Results", 
  "DESeq_Results_Index", 
  "nthr", 
  "RNA_Psdblk_WNN_Matrices", 
  "RNA_Psdblk_WNN_Matrices_Index", 
  "sample_data"   
)

source("~/ALS_Brain_Multiome.Rcfg")  

qreadm(
  paste0(
    "../Data/tmp/", 
    "",
    "B2_4_DE_RNA_WNN_Post", 
    ".qrdata"
  ), 
  nthr=nthr
)


      # END HPC Cluster Outsourcing --------------------------------------------


    ## 5.4 Collect DE stats ----------------------------------------------------

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



    ## 5.5 Export data ---------------------------------------------------------

qsave(
  dds_list, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9/", 
    "dds_list", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DDS_list, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9/", 
    "DDS_DESEqed_list", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9/", 
    "DESeq_Results", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_Index, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9/", 
    "DESeq_Results_Index", 
    ".qrds"
  ), 
  nthr=nthr
) 

qsave(
  sample_data, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9/", 
    "sample_data", 
    ".qrds"
  ), 
  nthr=nthr
)



{
  rm(list=ls())
  source("~/ALS_Brain_Multiome.Rcfg")
  RNA_Psdblk_WNN_Matrices <- qread(
    paste0(
      "../Data/DE/WNN/Pseudobulk/", 
      "RNA_Psdblk_WNN_Matrices", 
      ".qrds"
    ), 
    nthr=nthr
  )
  
  RNA_Psdblk_WNN_Matrices_Index <- qread(
    paste0(
      "../Data/DE/WNN/Pseudobulk/", 
      "RNA_Psdblk_WNN_Matrices_Index", 
      ".qrds"
    ), 
    nthr=nthr
  )
}




  ### 6.0 ALSFTD_NonC9vsHC RNA DE Analysis with DESeq2 ----------------------------



    ## 6.1 Extract sample data and subset Non_C9_ALSFTD and HC persons ---------

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



    ## 6.2 Format grouping variables -------------------------------------------

# Define Case and Sex as factors and order levels (for DESeq2)

sample_data$Case <- factor(sample_data$Case, levels=c("HC", "ALS_FTD")) 
sample_data$Sex <- factor(sample_data$Sex, levels=c("f", "m"))
sample_data$Rand <- factor(sample_data$Rand, levels=c("R1", "R2"))



    ## 6.3 DE Analysis  --------------------------------------------------------

      # Prepare data containers and counters -----------------------------------

dds_list <- list() 
DDS_list <- list() 
DESeq_Results <- list()   

nR <- nrow(RNA_Psdblk_WNN_Matrices_Index)

DESeq_Results_Index <- data.frame(
  Index=rep(NA, nR*2), 
  CellTypeLevel=rep(NA, nR*2), 
  CellType=rep(NA, nR*2), 
  Assay=rep("RNA", nR*2), 
  Test=rep("WALD", nR*2), 
  Covariates=rep("NoCov", nR*2), 
  SizeFactors=rep("Ratio", nR*2),
  Comparison=rep(NA, nR*2), 
  Remarks=rep(NA, nR*2)
) 


for(i in 1:2*nrow(RNA_Psdblk_WNN_Matrices_Index)){
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
        "RNA_Psdblk_WNN_Matrices", 
        "RNA_Psdblk_WNN_Matrices_Index", 
        "sample_data"
      ), 
      as.symbol
    ), 
    file=paste0(
      "../Data/tmp/", 
      "B2_5_DE_RNA_WNN_Pre", 
      ".qrdata"
    ), 
    nthr=nthr
  )
)    


      # Y_B2_5_RNA_DE_Slurm.R --------------------------------------------------

      # Import HPC Environment  ------------------------------------------------

rm(
  "dds_list", 
  "DDS_list", 
  "DESeq_Results", 
  "DESeq_Results_Index", 
  "nthr", 
  "RNA_Psdblk_WNN_Matrices", 
  "RNA_Psdblk_WNN_Matrices_Index", 
  "sample_data"   
)

source("~/ALS_Brain_Multiome.Rcfg")  

qreadm(
  paste0(
    "../Data/tmp/", 
    "",
    "B2_5_DE_RNA_WNN_Post", 
    ".qrdata"
  ), 
  nthr=36
)


      # END HPC Cluster Outsourcing --------------------------------------------



    ## 6.4 Collect DE stats ----------------------------------------------------

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



    ## 6.5 Export data ---------------------------------------------------------

qsave(
  dds_list, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_NonC9/", 
    "dds_list", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DDS_list, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_NonC9/", 
    "DDS_DESEqed_list", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_NonC9/", 
    "DESeq_Results", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_Index, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_NonC9/", 
    "DESeq_Results_Index", 
    ".qrds"
  ), 
  nthr=nthr
) 

qsave(
  sample_data, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_NonC9/", 
    "sample_data", 
    ".qrds"
  ), 
  nthr=nthr
)


{
  rm(list=ls())
  source("~/ALS_Brain_Multiome.Rcfg")
  RNA_Psdblk_WNN_Matrices <- qread(
    paste0(
      "../Data/DE/WNN/Pseudobulk/", 
      "RNA_Psdblk_WNN_Matrices", 
      ".qrds"
    ), 
    nthr=nthr
  )
  
  RNA_Psdblk_WNN_Matrices_Index <- qread(
    paste0(
      "../Data/DE/WNN/Pseudobulk/", 
      "RNA_Psdblk_WNN_Matrices_Index", 
      ".qrds"
    ), 
    nthr=nthr
  )
}




  ### 7.0 ALSFTD_C9vsALSFTD_NonC9 RNA DE Analysis with DESeq2 ------------------



    ## 7.1 Extract sample data and subset C9_ALSFTD and HC persons -------------

sample_data <- M0_RNA@misc$Sample_metrics 
sample_data$Sex <- M0_RNA$Sex[match(sample_data$ID, M0_RNA$ID)]
setequal(
  unique(paste0(M0_RNA@misc$Sample_data$ID, "_", M0_RNA@misc$Sample_data$Sex)), 
  paste0(sample_data$ID, "_", sample_data$Sex)
)
sample_data <- sample_data[,c(1,17:19, 2:16)]

table(sample_data$Case_Type)
sample_data <- sample_data[sample_data$Case_Type %in% c("C9_ALS_FTD", "ALS_FTD"),]



sample_data$Rand <- M0_RNA@misc$Sample_data$Rand1_ALSFTD_C9vsALSFTD_NonC9[match(sample_data$ID, M0_RNA@misc$Sample_data$ID)]

table(sample_data$Rand)
table(sample_data$Case_Type, sample_data$Rand) 
table(sample_data$Sex, sample_data$Case_Type)
table(sample_data$Sex, sample_data$Rand) 



    ## 7.2 Format grouping variables -------------------------------------------

# Define Case and Sex as factors and order levels (for DESeq2)

sample_data$Case_Type <- factor(sample_data$Case_Type, levels=c("C9_ALS_FTD", "ALS_FTD")) 
sample_data$Sex <- factor(sample_data$Sex, levels=c("f", "m"))
sample_data$Rand <- factor(sample_data$Rand, levels=c("R1", "R2"))




    ## 7.3 DE Analysis  --------------------------------------------------------


      # Prepare data containers and counters -----------------------------------

dds_list <- list() 
DDS_list <- list() 
DESeq_Results <- list()   

nR <- nrow(RNA_Psdblk_WNN_Matrices_Index)

DESeq_Results_Index <- data.frame(
  Index=rep(NA, nR*2), 
  CellTypeLevel=rep(NA, nR*2), 
  CellType=rep(NA, nR*2), 
  Assay=rep("RNA", nR*2), 
  Test=rep("WALD", nR*2), 
  Covariates=rep("NoCov", nR*2), 
  SizeFactors=rep("Ratio", nR*2),
  Comparison=rep(NA, nR*2), 
  Remarks=rep(NA, nR*2)
) 


for(i in 1:2*nrow(RNA_Psdblk_WNN_Matrices_Index)){
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
        "RNA_Psdblk_WNN_Matrices", 
        "RNA_Psdblk_WNN_Matrices_Index", 
        "sample_data"
      ), 
      as.symbol
    ), 
    file=paste0(
      "../Data/tmp/", 
      "B2_6_DE_RNA_WNN_Pre", 
      ".qrdata"
    ), 
    nthr=nthr
  )
)    


      # Y_B2_6_RNA_DE_Slurm.R --------------------------------------------------


      # Import HPC Environment  ------------------------------------------------

rm(
  "dds_list", 
  "DDS_list", 
  "DESeq_Results", 
  "DESeq_Results_Index", 
  "RNA_Psdblk_WNN_Matrices", 
  "RNA_Psdblk_WNN_Matrices_Index", 
  "sample_data"   
)

qreadm(
  paste0(
    "../Data/tmp/", 
    "",
    "B2_6_DE_RNA_WNN_Post",  
    ".qrdata"
  ), 
  nthr=36
)


      # END HPC Cluster Outsourcing --------------------------------------------



    ## 7.4 Collect DE stats ----------------------------------------------------

DESeq_Results_Index$All_q_0.05 <- NA 
DESeq_Results_Index$Up_q_0.05 <- NA  
DESeq_Results_Index$Down_q_0.05 <- NA 
DESeq_Results_Index$nSamples <- NA 
DESeq_Results_Index$nSamples_ALSFTD_C9 <- NA  
DESeq_Results_Index$nSamples_ALSFTD_NonC9 <- NA 
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
    DESeq_Results_Index$nSamples_ALSFTD_C9[i] <- table(sample_data$Case_Type[sample_data$ID %in% colnames(DDS_list[[i]]@assays@data$counts)])[1]
    DESeq_Results_Index$nSamples_ALSFTD_NonC9[i] <- table(sample_data$Case_Type[sample_data$ID %in% colnames(DDS_list[[i]]@assays@data$counts)])[2]
    DESeq_Results_Index$TotalReads[i] <- sum(colSums(DDS_list[[i]]@assays@data$counts)) 
    DESeq_Results_Index$MeanReads[i] <- sum(colMeans(DDS_list[[i]]@assays@data$counts)) 
    DESeq_Results_Index$nFeatures[i] <- sum(rowMeans(DDS_list[[i]]@assays@data$counts)>0) 
  }
  
}

DESeq_Results_Index$All_q_0.05 <- DESeq_Results_Index$Up_q_0.05 + DESeq_Results_Index$Down_q_0.05 
rm(i)



    ## 7.5 Export data ---------------------------------------------------------

qsave(
  dds_list, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9vsALSFTD_NonC9/", 
    "dds_list", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DDS_list, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9vsALSFTD_NonC9/", 
    "DDS_DESEqed_list", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9vsALSFTD_NonC9/", 
    "DESeq_Results", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_Index, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9vsALSFTD_NonC9/", 
    "DESeq_Results_Index", 
    ".qrds"
  ), 
  nthr=nthr
) 

qsave(
  sample_data, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9vsALSFTD_NonC9/", 
    "sample_data", 
    ".qrds"
  ), 
  nthr=nthr
)



