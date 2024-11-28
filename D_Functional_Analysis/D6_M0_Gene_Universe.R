

  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 

library(qs)
library(tidyverse)
library(DESeq2)



qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)




  ### 1.0 Load data ------------------------------------------------------------

DESeq_Results_12SVs_Disease_AllCells_list <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/Disease/SVA/12_SVs/", 
    "DESeq_Results_12SVs_Disease_AllCells_list", 
    ".qrds"
  ), 
  nthr=nthr
)

DESeq_Results_12SVs_Disease_AllCells_Index <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/Disease/SVA/12_SVs/", 
    "DESeq_Results_12SVs_Disease_AllCells_Index", 
    ".qrds"
  ), 
  nthr=nthr
)

DESeq_Results_12SVs_Disease_WNN_L25_list <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/Disease/SVA/12_SVs/", 
    "DESeq_Results_12SVs_Disease_WNN_L25_list", 
    ".qrds"
  ), 
  nthr=nthr
)

DESeq_Results_12SVs_Disease_WNN_L25_Index <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/Disease/SVA/12_SVs/", 
    "DESeq_Results_12SVs_Disease_WNN_L25_Index", 
    ".qrds"
  ), 
  nthr=nthr
)


DESeq_Results_12SVs_AllCase_AllCells_list <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/1_15_SVs/", 
    "DDS_12SVs_list_DESeqed", 
    ".qrds"
  ), 
  nthr=nthr
)

DESeq_Results_12SVs_Disease_AllCells_Index <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/1_15_SVs/", 
    "SVAs_15SVs_Index", 
    ".qrds"
  ), 
  nthr=nthr
)


DESeq_Results_12SVs_AllCase_WNN_L25_list <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/12_SVs/", 
    "DDS_SVA_12sv_list", 
    ".qrds"
  ), 
  nthr=nthr
)

DESeq_Results_12SVs_AllCase_WNN_L25_Index <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/12_SVs/", 
    "DESeq_Results_12SVs_Index", 
    ".qrds"
  ), 
  nthr=nthr
)


GEX_Features <- qread(
  paste0(
    "../Data/Annotations/Genes/", 
    "GEX_Features", 
    ".qrds"
  ), 
  nthr=nthr
)




  ### 2.0 Define AUX functions -------------------------------------------------

get_sign_genes_df <- function(res, fdr=0.05){
  res <- data.frame(res)
  res <- res[!is.na(res$padj),]
  res <- res[res$padj<fdr,]
  return(res)
}

get_nDEGs <- function(x, name){
  return(nrow(get_sign_genes_df(results(x, alpha=0.05, lfcThreshold = 0, name = name))))
}




  ### 3.0 Generate Gene Universes for 'Disease' comparisons --------------------

Universe_Disease_List <- list()




    ## 3.1 AllCells comparisons ------------------------------------------------

tmp <- as.data.frame(DESeq_Results_12SVs_Disease_AllCells_list[[1]])
tmp <- tmp[tmp$baseMean >= min(tmp$baseMean[tmp$padj<0.05], na.rm = TRUE),]

tmp <- data.frame(
  ID_10X = rownames(tmp)
)
tmp$SYMBOL <- GEX_Features$SYMBOL[match(tmp$ID_10X, GEX_Features$ID_10X)] 
Universe_Disease_List[["AllCells"]] <- tmp 
rm(tmp)



    ## 3.2 WNN_L25 comparisons -------------------------------------------------

ind <- which(
  DESeq_Results_12SVs_Disease_WNN_L25_Index$CellTypeLevel=="WNN_L25" & 
    DESeq_Results_12SVs_Disease_WNN_L25_Index$Comparison=="Disease"
)

for(i in ind){
  
  tmp <- as.data.frame(DESeq_Results_12SVs_Disease_WNN_L25_list[i]) 
  tmp <- tmp[tmp$baseMean >= min(tmp$baseMean[tmp$padj<0.05], na.rm = TRUE),]
  
  tmp <- data.frame(
    ID_10X = rownames(tmp)
  )
  tmp$SYMBOL <- GEX_Features$SYMBOL[match(tmp$ID_10X, GEX_Features$ID_10X)] 
  Universe_Disease_List[[DESeq_Results_12SVs_Disease_WNN_L25_Index$CellType[i]]] <- tmp 
  rm(tmp)
  
}



    ## 3.3 Save data -----------------------------------------------------------

qsave(
  Universe_Disease_List, 
  paste0(
    "../Data/Annotations/Genes/", 
    "Universe_Disease_List", 
    ".qrds"
  ), 
  nthr=nthr
)

rm(
  ind, 
  DESeq_Results_12SVs_Disease_AllCells_list, 
  DESeq_Results_12SVs_Disease_AllCells_Index, 
  DESeq_Results_12SVs_Disease_WNN_L25_list, 
  DESeq_Results_12SVs_Disease_WNN_L25_Index, 
  i, 
  Universe_Disease_List
)




  ### 4.0 Generate Gene Universes for 'AllCase' comparisons --------------------

Universe_AllCase_ALS_List <- list()
Universe_AllCase_ALSFTD_List <- list() 
Universe_AllCase_ALS_ALSFTD_List <- list()



    ## 4.1 AllCells comparisons ------------------------------------------------

ind <- which(
  DESeq_Results_12SVs_AllCase_WNN_L25_Index$CellTypeLevel=="AllCells" & 
    DESeq_Results_12SVs_AllCase_WNN_L25_Index$Comparison=="All_Cases"
)


tmp_ALS <- as.data.frame(
  results(
    DESeq_Results_12SVs_AllCase_AllCells_list[[ind]], 
    name="Case_ALS_vs_HC", 
    alpha=0.05, 
    lfcThreshold = 0
  )
)

tmp_ALS <- tmp_ALS[tmp_ALS$baseMean >= min(tmp_ALS$baseMean[tmp_ALS$padj<0.05], na.rm = TRUE),]

tmp_ALS <- data.frame(
  ID_10X = rownames(tmp_ALS)
)
tmp_ALS$SYMBOL <- GEX_Features$SYMBOL[match(tmp_ALS$ID_10X, GEX_Features$ID_10X)] 
Universe_AllCase_ALS_List[["AllCells"]] <- tmp_ALS 



tmp_ALSFTD <- as.data.frame(
  results(
    DESeq_Results_12SVs_AllCase_AllCells_list[[ind]], 
    name="Case_ALS_FTD_vs_HC", 
    alpha=0.05, 
    lfcThreshold = 0
  )
)

tmp_ALSFTD <- tmp_ALSFTD[tmp_ALSFTD$baseMean >= min(tmp_ALSFTD$baseMean[tmp_ALSFTD$padj<0.05], na.rm = TRUE),]

tmp_ALSFTD <- data.frame(
  ID_10X = rownames(tmp_ALSFTD)
)
tmp_ALSFTD$SYMBOL <- GEX_Features$SYMBOL[match(tmp_ALSFTD$ID_10X, GEX_Features$ID_10X)] 
Universe_AllCase_ALSFTD_List[["AllCells"]] <- tmp_ALSFTD 


tmp_ALS_ALSFTD <- data.frame(
  ID_10X = unique(
    c(
      tmp_ALS$ID_10X,
      tmp_ALSFTD$ID_10X
    )
  )
)
tmp_ALS_ALSFTD$SYMBOL <- GEX_Features$SYMBOL[match(tmp_ALS_ALSFTD$ID_10X, GEX_Features$ID_10X)] 
Universe_AllCase_ALS_ALSFTD_List[["AllCells"]] <- tmp_ALS_ALSFTD 


rm(tmp_ALS, tmp_ALSFTD, tmp_ALS_ALSFTD, ind)



    ## 4.2 WNN_L25 comparisons -------------------------------------------------

ind <- which(
  DESeq_Results_12SVs_AllCase_WNN_L25_Index$CellTypeLevel=="WNN_L25" & 
    DESeq_Results_12SVs_AllCase_WNN_L25_Index$Comparison=="All_Cases"
)

for(i in ind){
  
  tmp_ALS <- as.data.frame(
    results(
      DESeq_Results_12SVs_AllCase_WNN_L25_list[[i]], 
      name="Case_ALS_vs_HC", 
      alpha=0.05, 
      lfcThreshold = 0
    )
  )
  
  
  if(any(tmp_ALS$padj<0.05, na.rm=TRUE)){
    
    tmp_ALS <- tmp_ALS[tmp_ALS$baseMean >= min(tmp_ALS$baseMean[tmp_ALS$padj<0.05], na.rm = TRUE),]
    
    tmp_ALS <- data.frame(
      ID_10X = rownames(tmp_ALS)
    )
    tmp_ALS$SYMBOL <- GEX_Features$SYMBOL[match(tmp_ALS$ID_10X, GEX_Features$ID_10X)] 
    
  } else {
    
    tmp_ALS <- NULL
    
  }
  
  if(!is.null(tmp_ALS)) {
    Universe_AllCase_ALS_List[[DESeq_Results_12SVs_AllCase_WNN_L25_Index$CellType[i]]] <- tmp_ALS 
  } else {
    Universe_AllCase_ALS_List[[DESeq_Results_12SVs_AllCase_WNN_L25_Index$CellType[i]]] <- "NA"
  }
  
  
  
  tmp_ALSFTD <- as.data.frame(
    results(
      DESeq_Results_12SVs_AllCase_WNN_L25_list[[i]], 
      name="Case_ALS_FTD_vs_HC", 
      alpha=0.05, 
      lfcThreshold = 0
    )
  )
  
  if(any(tmp_ALSFTD$padj<0.05)){
    
    tmp_ALSFTD <- tmp_ALSFTD[tmp_ALSFTD$baseMean >= min(tmp_ALSFTD$baseMean[tmp_ALSFTD$padj<0.05], na.rm = TRUE),]
    
    tmp_ALSFTD <- data.frame(
      ID_10X = rownames(tmp_ALSFTD)
    )
    tmp_ALSFTD$SYMBOL <- GEX_Features$SYMBOL[match(tmp_ALSFTD$ID_10X, GEX_Features$ID_10X)] 
    
  } else {
    
    tmp_ALSFTD <- NULL
    
  }
  
  if(!is.null(tmp_ALSFTD)) {
    Universe_AllCase_ALSFTD_List[[DESeq_Results_12SVs_AllCase_WNN_L25_Index$CellType[i]]] <- tmp_ALSFTD 
  } else {
    Universe_AllCase_ALSFTD_List[[DESeq_Results_12SVs_AllCase_WNN_L25_Index$CellType[i]]] <- "NA" 
  }
  
  
  if(all(is.null(c(tmp_ALS, tmp_ALSFTD)))){
    
    tmp_ALS_ALSFTD <- "NA"
    
  } else {
    
    
    tmp_ALS_ALSFTD <- data.frame(
      ID_10X = unique(
        c(
          tmp_ALS$ID_10X,
          tmp_ALSFTD$ID_10X
        )
      )
    )
    tmp_ALS_ALSFTD$SYMBOL <- GEX_Features$SYMBOL[match(tmp_ALS_ALSFTD$ID_10X, GEX_Features$ID_10X)] 
    
  }
    
  Universe_AllCase_ALS_ALSFTD_List[[DESeq_Results_12SVs_AllCase_WNN_L25_Index$CellType[i]]] <- tmp_ALS_ALSFTD 

  
  rm(tmp_ALS, tmp_ALSFTD, tmp_ALS_ALSFTD)
  
}
rm(ind, i)



    ## 4.3 Save data -----------------------------------------------------------

qsave(
  Universe_AllCase_ALS_List, 
  paste0(
    "../Data/Annotations/Genes/", 
    "Universe_AllCase_ALS_List",
    ".qrds"
  ), 
  nthr=nthr
) 

qsave(
  Universe_AllCase_ALSFTD_List, 
  paste0(
    "../Data/Annotations/Genes/", 
    "Universe_AllCase_ALSFTD_List",
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  Universe_AllCase_ALS_ALSFTD_List, 
  paste0(
    "../Data/Annotations/Genes/", 
    "Universe_AllCase_ALS_ALSFTD_List",
    ".qrds"
  ), 
  nthr=nthr
)






