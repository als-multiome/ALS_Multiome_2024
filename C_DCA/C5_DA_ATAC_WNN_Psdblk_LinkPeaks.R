

  ### 0.0 Load libraries -------------------------------------------------------

source("~/ALS_Brain_Multiome.Rcfg")

library(qs) 
library(data.table)
library(Seurat)
library(tidyverse) 
library(DESeq2)
library(BiocParallel)



qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)




  ### 1.0 Load data ------------------------------------------------------------

M0_ATAC <- qread(
  paste0(
    "../Data/SeuratObjects/",  
    "M0_ATAC",
    ".qrds"
  ),
  nthr=nthr 
)




      # ATAC LinkPeaks WNN Pseudobulk matrices ---------------------------------

ATAC_Psdblk_WNN_Matrices <- qread(
  paste0(
    "../Data/DA/WNN/Pseudobulk/", 
    "ATAC_LinkPeaks_Psdblk_WNN_Matrices", 
    ".qrds"
  ), 
  nthr=nthr
)

ATAC_Psdblk_WNN_Matrices_Index <- qread(
  paste0(
    "../Data/DA/WNN/Pseudobulk/", 
    "ATAC_LinkPeaks_Psdblk_WNN_Matrices_Index", 
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


register(MulticoreParam(34)) 


index=1
for (pseudobulk_index in 1:nrow(ATAC_Psdblk_WNN_Matrices_Index)){
  
  ######### Loop into new celltype ################################ 
  
  message("    ")
  message("#################################################################")
  
  message(paste0(
    "Working on ", 
    ATAC_Psdblk_WNN_Matrices_Index$CelltypeLevel[pseudobulk_index], 
    " with cell type ", 
    ATAC_Psdblk_WNN_Matrices_Index$CellType[pseudobulk_index], 
    "...")
  )
  message("#################################################################")  
  message("    ")
  
  coldata <- sample_data
  counts <- ATAC_Psdblk_WNN_Matrices[[pseudobulk_index]]  
  counts <- counts[,which(colnames(counts) %in% coldata$ID)] 
  coldata <- coldata[which(coldata$ID %in% colnames(counts)),]
  
  coldata=coldata[match(colnames(counts), coldata$ID),] 
  rownames(coldata) <- coldata$ID 
  
  message(
    paste0(
      "Pseudobulk columns and Sample_data rows match: ", 
      all(colnames(counts)==rownames(coldata)))  
  )
  
  
  ############ Case ################## ----
  
  
  message(paste0("Comparing HC, ALS and ALS-FTD and asssigning the results under index nr. ", index, "... ")) 
  
  DESeq_Results_Index$Index[index] <- index 
  DESeq_Results_Index$CellTypeLevel[index] <- ATAC_Psdblk_WNN_Matrices_Index$CelltypeLevel[pseudobulk_index]
  DESeq_Results_Index$CellType[index] <- ATAC_Psdblk_WNN_Matrices_Index$CellType[pseudobulk_index]
  DESeq_Results_Index$Comparison[index] <- "All_Cases" 
  
  
  tryCatch({
    
    
    
    dds <- DESeqDataSetFromMatrix(
      countData = counts, 
      colData = coldata, 
      design = ~ Case 
    )
    dds_list[[index]] <- dds             
    
    DDS <- estimateSizeFactors(dds, type="poscounts")
    DDS <- DESeq(DDS, test="Wald", quiet = FALSE,  parallel=TRUE)
    DDS_list[[index]] <- DDS
    
    
    tryCatch({
      DESeq_Results_ALS[[index]] <- results(DDS, alpha=0.05, lfcThreshold = 0, name="Case_ALS_vs_HC") 
      DESeq_Results_ALSFTD[[index]] <- results(DDS, alpha=0.05, lfcThreshold = 0, name="Case_ALS_FTD_vs_HC") 
    }, error=function(e){})
    
    
    
  }, error=function(e){
    
    message(paste0("Result #", index, " could not be generated "))
    
    DESeq_Results_Index$Remarks[index] <<- "Results could not be calculated"    
    
    dds_list[[index]] <<- "NA"            
    DDS_list[[index]] <<- "NA" 
    DESeq_Results_ALS[[index]] <<- "NA"
    DESeq_Results_ALSFTD[[index]] <<- "NA"
  }
  ) 
  
  suppressWarnings(rm(dds, DDS))
  
  index=index+1                               
  
  
  
  ############ Rand ################## ----
  
  
  message(paste0("Comparing Rand1, Rand2 and Rand3, asssigning the results under index nr. ", index, "... ")) 
  
  DESeq_Results_Index$Index[index] <- index 
  DESeq_Results_Index$CellTypeLevel[index] <- ATAC_Psdblk_WNN_Matrices_Index$CelltypeLevel[pseudobulk_index]
  DESeq_Results_Index$CellType[index] <- ATAC_Psdblk_WNN_Matrices_Index$CellType[pseudobulk_index]
  DESeq_Results_Index$Comparison[index] <- "Rand" 
  
  
  tryCatch({
    
    
    
    dds.r <- DESeqDataSetFromMatrix(
      countData = counts, 
      colData = coldata, 
      design = ~ Rand 
    )
    dds_list[[index]] <- dds.r             
    
    
    DDS.r <- estimateSizeFactors(dds.r, type="poscounts")
    DDS.r <- DESeq(DDS.r, test="Wald", quiet = FALSE, parallel=TRUE)
    DDS_list[[index]] <- DDS.r
    
    
    tryCatch({
      DESeq_Results_ALS[[index]] <- results(DDS.r, alpha=0.05, lfcThreshold = 0, name="Rand_R2_vs_R1") 
      DESeq_Results_ALSFTD[[index]] <- results(DDS.r, alpha=0.05, lfcThreshold = 0, name="Rand_R3_vs_R1") 
    }, error=function(e){})
    
    
    
  }, error=function(e){
    
    message(paste0("Result #", index, " could not be generated "))
    
    DESeq_Results_Index$Remarks[index] <<- "Results could not be calculated"    
    
    dds_list[[index]] <<- "NA"            
    DDS_list[[index]] <<- "NA" 
    DESeq_Results_ALS[[index]] <<- "NA" 
    DESeq_Results_ALSFTD[[index]] <<- "NA"
  }
  ) 
  
  suppressWarnings(rm(dds, dds.r, DDS, DDS.r))
  
  index=index+1  
  
  ############ Print Info ############ ----  
  
  message(
    paste0(
      "Done with CellType ",
      ATAC_Psdblk_WNN_Matrices_Index$CellType[pseudobulk_index], 
      " on cell type level ", 
      ATAC_Psdblk_WNN_Matrices_Index$CelltypeLevel[pseudobulk_index]
    )
  ) 
  
  if(index%%5==0){
    qsave(
      DDS_list, 
      paste0(
        "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/", 
        "DDS_DESeqed_list", 
        ".qrds"
      ), 
      nthr=nthr
    )
    
  }
  
  print(Sys.time())
  
}
suppressWarnings(rm(counts, coldata, index, pseudobulk_index)) 



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
    "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/", 
    "dds_list", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DDS_list, 
  paste0(
    "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/", 
    "DDS_DESEqed_list", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_ALS, 
  paste0(
    "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/", 
    "DESeq_Results_ALS", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_ALSFTD, 
  paste0(
    "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/", 
    "DESeq_Results_ALSFTD", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  DESeq_Results_Index, 
  paste0(
    "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/", 
    "DESeq_Results_Index", 
    ".qrds"
  ), 
  nthr=nthr
) 

qsave(
  sample_data, 
  paste0(
    "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/", 
    "sample_data", 
    ".qrds"
  ), 
  nthr=nthr
) 

qsave(
  DESeq_Results_Index, 
  paste0(
    "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/", 
    "DESeq_Results_Index", 
    ".qrds"
  ), 
  nthr=nthr
) 

qsave(
  sample_data, 
  paste0(
    "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/", 
    "sample_data", 
    ".qrds"
  ), 
  nthr=nthr
)


rm(dds_list, DDS_list, DESeq_Results_ALS, DESeq_Results_ALSFTD, DESeq_Results_Index, sample_data, i)


