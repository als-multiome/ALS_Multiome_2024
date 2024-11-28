

  ### 0.0 Load libraries -------------------------------------------------------

source("~/ALS_Brain_Multiome.Rcfg")


library(qs) 
library(DESeq2)
library(BiocParallel)



qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)




  ### 1.0 Load data ------------------------------------------------------------

DDS_SVA_12sv_list <-qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/12_SVs/", 
    "DDS_SVA_12sv_list", 
    ".qrds"
  ), 
  nthr=nthr
)


DESeq_Results_12SVs_ALS <-qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/12_SVs/", 
    "DESeq_Results_12SVs_ALS", 
    ".qrds"
  ), 
  nthr=nthr
)

DESeq_Results_12SVs_ALSFTD <-qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/12_SVs/", 
    "DESeq_Results_12SVs_ALSFTD", 
    ".qrds"
  ), 
  nthr=nthr
)

DESeq_Results_12SVs_Index <-qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/12_SVs/", 
    "DESeq_Results_12SVs_Index", 
    ".qrds"
  ), 
  nthr=nthr
)




  ### 2.0 DGEs WNN_L25, ALS, ALSFTD L2FC Shrinkage -----------------------------

register(MulticoreParam(44))

L2FC_Shrink_Results_12SVs_ALS <- list() 
L2FC_Shrink_Results_12SVs_ALSFTD <- list() 


for (i in 1:length(DDS_SVA_12sv_list)){
  
  L2FC_Shrink_Results_12SVs_ALS[[i]] <- lfcShrink(
    DDS_SVA_12sv_list[[i]], 
    coef = if(DESeq_Results_12SVs_Index$Comparison[i] == "All_Cases") "Case_ALS_vs_HC" else "Rand_R2_vs_R1", 
    res=DESeq_Results_12SVs_ALS[[i]], 
    type="apeglm", 
    parallel = TRUE
  ) 
  
  L2FC_Shrink_Results_12SVs_ALSFTD[[i]] <- lfcShrink(
    DDS_SVA_12sv_list[[i]], 
    coef = if(DESeq_Results_12SVs_Index$Comparison[i] == "All_Cases") "Case_ALS_FTD_vs_HC" else "Rand_R3_vs_R1", 
    res=DESeq_Results_12SVs_ALSFTD[[i]], 
    type="apeglm", 
    parallel = TRUE
  ) 
  print(i)
}




  ### 3.0 Save Results ---------------------------------------------------------

qsave(
  L2FC_Shrink_Results_12SVs_ALS, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/12_SVs/", 
    "L2FC_Shrink_Results_12SVs_ALS", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  L2FC_Shrink_Results_12SVs_ALSFTD, 
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/12_SVs/", 
    "L2FC_Shrink_Results_12SVs_ALSFTD", 
    ".qrds"
  ), 
  nthr=nthr
)
