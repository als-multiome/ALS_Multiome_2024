
  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 

library(qs)
library(tidyverse)
library(spatialLIBD)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)


qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)




  ### 1.0 Load data ------------------------------------------------------------

SVA_DDS_List_Case <- list() 
SVA_DDS_List_Rand <- list() 

SVAs_15SVs_Index <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/1_15_SVs/SVAs_15SVs_Index.qrds"
  )
)

for (i in 12:12){
  tmp <- qread(
    paste0(
      "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/1_15_SVs/DDS_", 
      i, 
      "SVs_list_DESeqed", 
      ".qrds"
    ),
    nthr=nthr
  )
  SVA_DDS_List_Case[[paste0(i, "SVs")]] <- tmp[[1]]
  SVA_DDS_List_Rand[[paste0(i, "SVs")]] <- tmp[[2]] 
  rm(tmp)
}

DDS_12SVs_Case <- SVA_DDS_List_Case[["12SVs"]]
DDS_12SVs_Rand <- SVA_DDS_List_Rand[["12SVs"]] 


GEX_Features <- qread(
  paste0(
    "../Data/Annotations/Genes/", 
    "GEX_Features", 
    ".qrds"
  ), 
  nthr=nthr
)


rm(SVA_DDS_List_Case, SVA_DDS_List_Rand, SVAs_15SVs_Index)
  



  ### 2.0 Define AUX functions -------------------------------------------------

get_sign_genes_df <- function(res, fdr=0.05){
  res <- data.frame(res)
  res <- res[!is.na(res$padj),]
  res <- res[res$padj<fdr,]
  return(res)
}



  ### 3.0 Collect DEGs ---------------------------------------------------------



    ## 3.1 ALS DESeq -----------------------------------------------------------

ALS_DESeq <- as.data.frame(
  results(DDS_12SVs_Case, name="Case_ALS_vs_HC", alpha=0.05, lfcThreshold = 0)
)

ALS_DESeq$SYMBOL <- GEX_Features$SYMBOL[match(rownames(ALS_DESeq), GEX_Features$ID_10X)]
ALS_DESeq$SYMBOL[is.na(ALS_DESeq$SYMBOL)] <- rownames(ALS_DESeq)[is.na(ALS_DESeq$SYMBOL)]
ALS_DESeq <- ALS_DESeq[order(ALS_DESeq$padj),] 



    ## 3.2 ALSFTD DESeq -----------------------------------------------------------

ALSFTD_DESeq <- as.data.frame(
  results(DDS_12SVs_Case, name="Case_ALS_FTD_vs_HC", alpha=0.05, lfcThreshold = 0)
)

ALSFTD_DESeq$SYMBOL <- GEX_Features$SYMBOL[match(rownames(ALSFTD_DESeq), GEX_Features$ID_10X)]
ALSFTD_DESeq$SYMBOL[is.na(ALSFTD_DESeq$SYMBOL)] <- rownames(ALSFTD_DESeq)[is.na(ALSFTD_DESeq$SYMBOL)]
ALSFTD_DESeq <- ALSFTD_DESeq[order(ALSFTD_DESeq$padj),] 



    ## 3.3 ALS & ALSFTD DESeq --------------------------------------------------

ALS_ALSFTD_DESeq <- unique(
  c(
    rownames(ALS_DESeq)[which(ALS_DESeq$padj<0.05)], 
    rownames(ALSFTD_DESeq)[which(ALSFTD_DESeq$padj<0.05)]
  )
)
ALS_ALSFTD_DESeq <- data.frame(
  ID_CellRanger = ALS_ALSFTD_DESeq, 
  SYMBOL = GEX_Features$SYMBOL[match(ALS_ALSFTD_DESeq, GEX_Features$ID_10X)], 
  GENETYPE = GEX_Features$GENETYPE[match(ALS_ALSFTD_DESeq, GEX_Features$ID_10X)]
)

ALS_ALSFTD_DESeq$DEG_In <- "NA"
ALS_ALSFTD_DESeq$DEG_In[
  ALS_ALSFTD_DESeq$ID_CellRanger %in% rownames(ALS_DESeq)[which(ALS_DESeq$padj<0.05)]
] <- "ALS" 

ALS_ALSFTD_DESeq$DEG_In[
  ALS_ALSFTD_DESeq$ID_CellRanger %in% rownames(ALSFTD_DESeq)[which(ALSFTD_DESeq$padj<0.05)]
] <- "ALS-FTD"

ALS_ALSFTD_DESeq$DEG_In[
  ALS_ALSFTD_DESeq$ID_CellRanger %in% rownames(ALS_DESeq)[which(ALS_DESeq$padj<0.05)] & 
  ALS_ALSFTD_DESeq$ID_CellRanger %in% rownames(ALSFTD_DESeq)[which(ALSFTD_DESeq$padj<0.05)]
] <- "ALS & ALS-FTD"


table(ALS_ALSFTD_DESeq$DEG_In)




  ### 4.0 Write data -----------------------------------------------------------

write.csv(
  ALS_DESeq, 
  paste0(
    "../Data/Visualization/Tables/",
    "SupplTab_RNA_Complete_Pseudobulk_DEGs_ALS", 
    ".csv"
  ), 
  row.names = FALSE, 
  quote=FALSE
) 

write.csv(
  ALSFTD_DESeq, 
  paste0(
    "../Data/Visualization/Tables/",
    "SupplTab_RNA_Complete_Pseudobulk_DEGs_ALSFTD", 
    ".csv"
  ), 
  row.names = FALSE, 
  quote=FALSE
)

write.csv(
  ALS_ALSFTD_DESeq, 
  paste0(
    "../Data/Visualization/Tables/",
    "SupplTab_RNA_Complete_Pseudobulk_DEGs_ALS_ALSFTD", 
    ".csv"
  ), 
  row.names = FALSE, 
  quote=FALSE
)
