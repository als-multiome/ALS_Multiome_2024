
  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 

library(qs)
library(tidyverse)
library(DESeq2)



qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)




  ### 1.0 Load data ------------------------------------------------------------

DDS_8SVs_list <- qread(
  paste0(
    "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/SVA/1_15_SVs/", 
    "DDS_8SVs_list", 
    ".qrds"
  ), 
  nthr=nthr
)


DESeq_Results_8SVs_Index <- qread(
  paste0(
    "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase/SVA/1_15_SVs/", 
    "DESeq_Results_8SVs_Index", 
    ".qrds"
  ), 
  nthr=nthr
)


Peak_Features <- qread(
  paste0(
    "../Data/Annotations/Genomic_Regions/", 
    "peakAnno_df", 
    ".qrds"
  ), 
  nthr=nthr
)

ATAC_Peaks_Links <- qread(
  paste0(
    "../Data/Annotations/Genomic_Regions/", 
    "M0_ATAC_Peak_Links", 
    ".qrds"
  ), 
  nthr=nthr
)




  ### 2.0 Define AUX functions -------------------------------------------------

get_sign_DARs_df <- function(res, fdr=0.05){
  res <- data.frame(res)
  res <- res[!is.na(res$padj),]
  res <- res[res$padj<fdr,]
  return(res)
}

get_nDARs <- function(x, name){
  return(nrow(get_sign_DARs_df(results(x, alpha=0.05, lfcThreshold = 0, name = name))))
}




  ### 4.0 Generate Gene Universes for 'AllCase' comparisons --------------------

Universe_AllCase_ALS_List <- list()
Universe_AllCase_ALSFTD_List <- list() 
Universe_AllCase_ALS_ALSFTD_List <- list()



    ## 4.1 AllCells comparisons ------------------------------------------------

ind <- which(
  DESeq_Results_8SVs_Index$CellTypeLevel=="AllCells" & 
    DESeq_Results_8SVs_Index$Comparison=="All_Cases"
)


tmp_ALS <- as.data.frame(
  results(
   DDS_8SVs_list[[ind]], 
    name="Case_ALS_vs_HC", 
    alpha=0.05, 
    lfcThreshold = 0
  )
)

tmp_ALS <- tmp_ALS[tmp_ALS$baseMean >= min(tmp_ALS$baseMean[tmp_ALS$padj<0.05], na.rm = TRUE),]

tmp_ALS <- data.frame(
  ID2 = rownames(tmp_ALS)
)
tmp_ALS$ID <- Peak_Features$ID[match(tmp_ALS$ID2, Peak_Features$ID2)] 
Universe_AllCase_ALS_List[["AllCells"]] <- tmp_ALS 




tmp_ALSFTD <- as.data.frame(
  results(
    DDS_8SVs_list[[ind]], 
    name="Case_ALS_FTD_vs_HC", 
    alpha=0.05, 
    lfcThreshold = 0
  )
)

tmp_ALSFTD <- tmp_ALSFTD[tmp_ALSFTD$baseMean >= min(tmp_ALSFTD$baseMean[tmp_ALSFTD$padj<0.05], na.rm = TRUE),]

tmp_ALSFTD <- data.frame(
  ID2 = rownames(tmp_ALSFTD)
)
tmp_ALSFTD$ID <- Peak_Features$ID[match(tmp_ALSFTD$ID2, Peak_Features$ID2)] 
Universe_AllCase_ALSFTD_List[["AllCells"]] <- tmp_ALSFTD 


tmp_ALS_ALSFTD <- data.frame(
  ID2 = unique(
    c(
      tmp_ALS$ID2,
      tmp_ALSFTD$ID2
    )
  )
)
tmp_ALS_ALSFTD$ID <- Peak_Features$ID[match(tmp_ALS_ALSFTD$ID2, Peak_Features$ID2)] 
Universe_AllCase_ALS_ALSFTD_List[["AllCells"]] <- tmp_ALS_ALSFTD 


rm(tmp_ALS, tmp_ALSFTD, tmp_ALS_ALSFTD, ind)



    ## 4.2 WNN_L25 comparisons -------------------------------------------------

ind <- which(
 DESeq_Results_8SVs_Index$CellTypeLevel=="WNN_L25" & 
    DESeq_Results_8SVs_Index$Comparison=="All_Cases"
)

for(i in ind){
  
  tmp_ALS <- as.data.frame(
    results(
      DDS_8SVs_list[[i]], 
      name="Case_ALS_vs_HC", 
      alpha=0.05, 
      lfcThreshold = 0
    )
  )
  
  if(any(tmp_ALS$padj<0.05, na.rm = TRUE)){
    
    tmp_ALS <- tmp_ALS[tmp_ALS$baseMean >= min(tmp_ALS$baseMean[tmp_ALS$padj<0.05], na.rm = TRUE),]
  
    tmp_ALS <- data.frame(
      ID2 = rownames(tmp_ALS)
    )
    tmp_ALS$ID <- Peak_Features$ID[match(tmp_ALS$ID2, Peak_Features$ID2)] 
    
  } else {
    
    tmp_ALS <- NULL 
    
    
  }
  
  if(!is.null(tmp_ALS)) {
    Universe_AllCase_ALS_List[[DESeq_Results_8SVs_Index$CellType[i]]] <- tmp_ALS
  } else {
    Universe_AllCase_ALS_List[[DESeq_Results_8SVs_Index$CellType[i]]] <- "NA"
  }
  
  tmp_ALSFTD <- as.data.frame(
    results(
      DDS_8SVs_list[[i]], 
      name="Case_ALS_FTD_vs_HC", 
      alpha=0.05, 
      lfcThreshold = 0
    )
  )
  
  
  if(any(tmp_ALSFTD$padj<0.05, na.rm = TRUE)){
    
    tmp_ALSFTD <- tmp_ALSFTD[tmp_ALSFTD$baseMean >= min(tmp_ALSFTD$baseMean[tmp_ALSFTD$padj<0.05], na.rm = TRUE),]
    
    tmp_ALSFTD <- data.frame(
      ID2 = rownames(tmp_ALSFTD)
    )
    tmp_ALSFTD$ID <- Peak_Features$ID[match(tmp_ALSFTD$ID2, Peak_Features$ID2)] 
    
    
  } else { 
    
    tmp_ALSFTD <- NULL 
    
    
  }
  
  if(!is.null(tmp_ALSFTD)) {
    Universe_AllCase_ALSFTD_List[[DESeq_Results_8SVs_Index$CellType[i]]] <- tmp_ALSFTD   
  } else {
    Universe_AllCase_ALSFTD_List[[DESeq_Results_8SVs_Index$CellType[i]]] <- "NA" 
  }
  
  
  if(all(is.null(c(tmp_ALS, tmp_ALSFTD)))){
    
    tmp_ALS_ALSFTD <- "NA"
    
  } else {
    
    tmp_ALS_ALSFTD <- data.frame(
      ID2 = unique(
        c(
          tmp_ALS$ID2,
          tmp_ALSFTD$ID2
        )
      )
    ) 
    
    tmp_ALS_ALSFTD$ID <- Peak_Features$ID[match(tmp_ALS_ALSFTD$ID2, Peak_Features$ID2)] 
    
  }
  
  Universe_AllCase_ALS_ALSFTD_List[[DESeq_Results_8SVs_Index$CellType[i]]] <- tmp_ALS_ALSFTD 
  
  
  rm(tmp_ALS, tmp_ALSFTD, tmp_ALS_ALSFTD)
  
}
rm(ind, i)



    ## 4.3 Save data -----------------------------------------------------------

qsave(
  Universe_AllCase_ALS_List, 
  paste0(
    "../Data/Annotations/Genomic_Regions/", 
    "Universe_AllCase_ALS_List",
    ".qrds"
  ), 
  nthr=nthr
) 

qsave(
  Universe_AllCase_ALSFTD_List, 
  paste0(
    "../Data/Annotations/Genomic_Regions/", 
    "Universe_AllCase_ALSFTD_List",
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  Universe_AllCase_ALS_ALSFTD_List, 
  paste0(
    "../Data/Annotations/Genomic_Regions/", 
    "Universe_AllCase_ALS_ALSFTD_List",
    ".qrds"
  ), 
  nthr=nthr
)






