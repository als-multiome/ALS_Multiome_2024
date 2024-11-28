.libPaths(c("/home/genomics/R/x86_64-pc-linux-gnu-library/4.4/", .libPaths()))

library(qs) 
library(DESeq2) 
library(BiocParallel) 

qreadm("./B2_6_DE_RNA_WNN_Pre.qrdata",nthr=40)

register(MulticoreParam(40)) 

# Loop ---- 

index=1
for (pseudobulk_index in 1:nrow(RNA_Psdblk_WNN_Matrices_Index)){
  
  ######### Loop into new celltype ################################ 
  
  message("    ")
  message("#################################################################")
  
  message(paste0(
    "Working on ", 
    RNA_Psdblk_WNN_Matrices_Index$CelltypeLevel[pseudobulk_index], 
    " with cell type ", 
    RNA_Psdblk_WNN_Matrices_Index$CellType[pseudobulk_index], 
    "...")
  )
  message("#################################################################")  
  message("    ")
  
  coldata <- sample_data
  counts <- RNA_Psdblk_WNN_Matrices[[pseudobulk_index]]  
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
  
  
  message(paste0("Comparing ALSFTD_C9 and ALSFTD_NonC9, asssigning the results under index nr. ", index, "... ")) 
  
  DESeq_Results_Index$Index[index] <- index 
  DESeq_Results_Index$CellTypeLevel[index] <- RNA_Psdblk_WNN_Matrices_Index$CelltypeLevel[pseudobulk_index]
  DESeq_Results_Index$CellType[index] <- RNA_Psdblk_WNN_Matrices_Index$CellType[pseudobulk_index]
  DESeq_Results_Index$Comparison[index] <- "ALSFTD_C9vsALSFTD_NonC9"  
  
  
  tryCatch({
    
    
    
    dds <- DESeqDataSetFromMatrix(
      countData = counts, 
      colData = coldata, 
      design = ~ Case_Type
    )
    dds_list[[index]] <- dds             
    
    DDS <- estimateSizeFactors(dds, type="ratio")
    DDS <- DESeq(DDS, test="Wald", quiet = TRUE,  parallel=TRUE)
    DDS_list[[index]] <- DDS
    
    
    DESeq_Results[[index]] <- results(DDS, alpha=0.05, lfcThreshold = 0) 
    
    
    
    
  }, error=function(e){
    
    message(paste0("Result #", index, " could not be generated "))
    
    DESeq_Results_Index$Remarks[index] <<- "Results could not be calculated"    
    
    dds_list[[index]] <<- "NA"            
    DDS_list[[index]] <<- "NA" 
    DESeq_Results[[index]] <<- "NA"
  }
  ) 
  
  suppressWarnings(rm(dds, DDS))
  
  index=index+1                               
  
  
  
  ############ Rand ################## ----
  
  
  message(paste0("Comparing Rand1 and Rand2, asssigning the results under index nr. ", index, "... ")) 
  
  DESeq_Results_Index$Index[index] <- index 
  DESeq_Results_Index$CellTypeLevel[index] <- RNA_Psdblk_WNN_Matrices_Index$CelltypeLevel[pseudobulk_index]
  DESeq_Results_Index$CellType[index] <- RNA_Psdblk_WNN_Matrices_Index$CellType[pseudobulk_index]
  DESeq_Results_Index$Comparison[index] <- "Rand2vsRand1" 
  
  
  tryCatch({
    
    
    
    dds.r <- DESeqDataSetFromMatrix(
      countData = counts, 
      colData = coldata, 
      design = ~ Rand 
    )
    dds_list[[index]] <- dds.r             
    
    
    DDS.r <- estimateSizeFactors(dds.r, type="ratio")
    DDS.r <- DESeq(DDS.r, test="Wald", quiet = TRUE, parallel=TRUE)
    DDS_list[[index]] <- DDS.r
    
    
    DESeq_Results[[index]] <- results(DDS.r, alpha=0.05, lfcThreshold = 0) 
    
    
    
    
  }, error=function(e){
    
    message(paste0("Result #", index, " could not be generated "))
    
    DESeq_Results_Index$Remarks[index] <<- "Results could not be calculated"    
    
    dds_list[[index]] <<- "NA"            
    DDS_list[[index]] <<- "NA" 
    DESeq_Results[[index]] <<- "NA"
  }
  ) 
  
  suppressWarnings(rm(dds, dds.r, DDS, DDS.r))
  
  index=index+1  
  
  ############ Print Info ############ ----  
  
  message(
    paste0(
      "Done with CellType ",
      RNA_Psdblk_WNN_Matrices_Index$CellType[pseudobulk_index], 
      " on cell type level ", 
      RNA_Psdblk_WNN_Matrices_Index$CelltypeLevel[pseudobulk_index]
    )
  ) 
  
  print(Sys.time())
  
}
suppressWarnings(rm(counts, coldata, index, pseudobulk_index))




do.call(qsavem,c(lapply(ls(), as.symbol), file="./B2_6_DE_RNA_WNN_Post.qrdata", nthr=40))
print("Job done!...")
