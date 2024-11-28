
  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 


library(qs) 
library(GenomicRanges)
library(DESeq2)
library(BiocParallel)



qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)




  ### 1.0 Load data ------------------------------------------------------------ 

DDS_8SVs_list <- qread(
  "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/SVA/1_15_SVs/DDS_8SVs_list.qrds", nthr=nthr
)

DESeq_Results_8SVs_ALS <- qread(
  "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/SVA/1_15_SVs/DESeq_Results_8SVs_ALS.qrds", nthr=nthr
)

DESeq_Results_8SVs_ALSFTD <- qread(
  "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/SVA/1_15_SVs/DESeq_Results_8SVs_ALSFTD.qrds", nthr=nthr
)

DESeq_Results_8SVs_Index <- qread(
  "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/SVA/1_15_SVs/DESeq_Results_8SVs_Index.qrds", nthr=nthr
)




  ### DARs WNN_L25 ALS, ALS-FTD L2FC Shrinkage ---------------------------------

register(MulticoreParam(44))

L2FC_Shrink_Results_8SVs_ALS <- list() 
L2FC_Shrink_Results_8SVs_ALSFTD <- list() 


for (i in 1:length(DDS_8SVs_list)){
  
  L2FC_Shrink_Results_8SVs_ALS[[i]] <- lfcShrink(
    DDS_8SVs_list[[i]], 
    coef = if(DESeq_Results_8SVs_Index$Comparison[i] == "All_Cases") "Case_ALS_vs_HC" else "Rand_R2_vs_R1", 
    res=DESeq_Results_8SVs_ALS[[i]], 
    type="apeglm", 
    parallel = TRUE
  ) 
  
  L2FC_Shrink_Results_8SVs_ALSFTD[[i]] <- lfcShrink(
    DDS_8SVs_list[[i]], 
    coef = if(DESeq_Results_8SVs_Index$Comparison[i] == "All_Cases") "Case_ALS_FTD_vs_HC" else "Rand_R3_vs_R1", 
    res=DESeq_Results_8SVs_ALSFTD[[i]], 
    type="apeglm", 
    parallel = TRUE
  ) 
  
}


qsave(
  L2FC_Shrink_Results_8SVs_ALS, 
  paste0(
    "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/SVA/1_15_SVs/", 
  "L2FC_Shrink_Results_8SVs_ALS", 
  ".qrds"
  ), 
  nthr=nthr
)

qsave(
  L2FC_Shrink_Results_8SVs_ALSFTD, 
  paste0(
    "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/SVA/1_15_SVs/", 
    "L2FC_Shrink_Results_8SVs_ALSFTD", 
    ".qrds"
  ), 
  nthr=nthr
)

