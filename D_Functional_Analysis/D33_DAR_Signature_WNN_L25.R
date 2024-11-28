
  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 

library(qs)  
library(tidyverse) 
library(DESeq2)
library(GenomicRanges)




qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)




  ### 1.0 Load data ------------------------------------------------------------

DESeq_Results_8SVs_ALS <- qread(
  paste0(
    "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/SVA/1_15_SVs/", 
    "DESeq_Results_8SVs_ALS", 
    ".qrds"
  ), 
  nthr=nthr
)

DESeq_Results_8SVs_ALSFTD <- qread(
  paste0(
    "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/SVA/1_15_SVs/", 
    "DESeq_Results_8SVs_ALSFTD", 
    ".qrds"
  ), 
  nthr=nthr
)

DESeq_Results_8SVs_Index <- qread(
  paste0(
    "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/SVA/1_15_SVs/", 
    "DESeq_Results_8SVs_Index", 
    ".qrds"
  ), 
  nthr=nthr
)


ATAC_Peak_Links <- qread(
  paste0(
    "../Data/Annotations/Genomic_Regions/", 
    "M0_ATAC_Peak_Links", 
    ".qrds"
  ), 
  nthr=nthr
)


Peaks <- qread(
  paste0(
    "../Data/Annotations/Genomic_Regions/", 
    "peakAnno_df", 
    ".qrds"
  ), 
  nthr=nthr
)




  ### 2.0 Define AUX functions -------------------------------------------------

get_sign_peaks_df <- function(res, fdr=0.05){
  res <- data.frame(res)
  res <- res[!is.na(res$padj),]
  res <- res[res$padj<fdr,]
  return(res)
}

get_sign_peaks_names <- function(res, fdr=0.05){
  res <- data.frame(res)
  res <- res[!is.na(res$padj),]
  res <- res[res$padj<fdr,]
  return(rownames(res))
}


get_same_N_rand_comparison_peaks_names <- function(res, N){
  res <- data.frame(res)
  res <- res[!is.na(res$pvalue),]
  res <- res[order(res$pvalue),]
  return(rownames(res)[1:N])
}




  ### 3.0 Get all DARs that are significant in WNN_L25 ------------------------

ind_WNN_L25_Case_Comparisons <- which(
  DESeq_Results_8SVs_Index$CellTypeLevel=="WNN_L25" & 
    DESeq_Results_8SVs_Index$Comparison=="All_Cases"
)

ind_WNN_L25_Rand_Comparisons <- which(
  DESeq_Results_8SVs_Index$CellTypeLevel=="WNN_L25" & 
    DESeq_Results_8SVs_Index$Comparison=="Rand"
)



    ## 3.1 ALS vs HC comparison ------------------------------------------------

All_DARs_WNN_L25_ALS <- unique(
  unlist(
    lapply(
      DESeq_Results_8SVs_ALS[ind_WNN_L25_Case_Comparisons], 
      FUN=get_sign_peaks_names
    )
  ) 
)



    ## 3.2 ALSFTD vs HC comparison ---------------------------------------------

All_DARs_WNN_L25_ALSFTD <- unique(
  unlist(
    lapply(
      DESeq_Results_8SVs_ALSFTD[ind_WNN_L25_Case_Comparisons], 
      FUN=get_sign_peaks_names
    )
  ) 
)



    ## 3.3 R2 vs R1 comparison -------------------------------------------------

Ns_DARs_ALS <- unlist(
  lapply(
    DESeq_Results_8SVs_ALS[ind_WNN_L25_Case_Comparisons], 
    FUN=function(x){
      return(
        length(
          get_sign_peaks_names(x)
        )
      )
    }
  )
)

i <- 0
All_DARs_WNN_L25_R2 <- unique(
  unlist(
    lapply(
      DESeq_Results_8SVs_ALS[ind_WNN_L25_Rand_Comparisons], 
      FUN=function(x){
        i <<- i + 1
        return(
          get_same_N_rand_comparison_peaks_names(
            res = x, 
            N=Ns_DARs_ALS[i]
          ) 
        )
      }
    )
  ) 
)
rm(Ns_DARs_ALS, i)



    ## 3.4 R3 vs R1 comparison -------------------------------------------------

Ns_DARs_ALSFTD <- unlist(
  lapply(
    DESeq_Results_8SVs_ALSFTD[ind_WNN_L25_Case_Comparisons], 
    FUN=function(x){
      return(
        length(
          get_sign_peaks_names(x)
        )
      )
    }
  )
)

i <- 0
All_DARs_WNN_L25_R3 <- unique(
  unlist(
    lapply(
      DESeq_Results_8SVs_ALSFTD[ind_WNN_L25_Rand_Comparisons], 
      FUN=function(x){
        i <<- i + 1
        return(
          get_same_N_rand_comparison_peaks_names(
            res = x, 
            N=Ns_DARs_ALSFTD[i]
          ) 
        )
      }
    )
  ) 
)
rm(Ns_DARs_ALSFTD, i) 



    ## 3.5 Sanity checks -------------------------------------------------------

table(All_DARs_WNN_L25_ALS %in% All_DARs_WNN_L25_ALSFTD)
table(All_DARs_WNN_L25_ALSFTD %in% All_DARs_WNN_L25_ALS)

table(All_DARs_WNN_L25_R2 %in% All_DARs_WNN_L25_R3)
table(All_DARs_WNN_L25_R3 %in% All_DARs_WNN_L25_R2)




  ### 4.0 Compile and save data ------------------------------------------------



    ## 4.1 ALS, ALSFTD Signature -----------------------------------------------

ALS_ALSFTD_WNN_L25_Signature <- data.frame(
  Peak = sort(
    unique(c(
      All_DARs_WNN_L25_ALS, 
      All_DARs_WNN_L25_ALSFTD
    ))
  )
)

table(ALS_ALSFTD_WNN_L25_Signature$Peak %in% ATAC_Peak_Links$peak)
table(unique(ATAC_Peak_Links$peak) %in% ALS_ALSFTD_WNN_L25_Signature$Peak)

table(ALS_ALSFTD_WNN_L25_Signature$Peak %in% Peaks$ID2)


ALS_ALSFTD_WNN_L25_Signature$ID2 <- ALS_ALSFTD_WNN_L25_Signature$Peak
ALS_ALSFTD_WNN_L25_Signature$ID <- Peaks$ID[match(ALS_ALSFTD_WNN_L25_Signature$Peak, Peaks$ID2)]


ALS_ALSFTD_WNN_L25_Signature$Signif_ALS <- ALS_ALSFTD_WNN_L25_Signature$Peak %in% All_DARs_WNN_L25_ALS
ALS_ALSFTD_WNN_L25_Signature$Signif_ALSFTD <- ALS_ALSFTD_WNN_L25_Signature$Peak %in% All_DARs_WNN_L25_ALSFTD
ALS_ALSFTD_WNN_L25_Signature$Signif_In <- "NA"
ALS_ALSFTD_WNN_L25_Signature$Signif_In[ALS_ALSFTD_WNN_L25_Signature$Signif_ALS & !ALS_ALSFTD_WNN_L25_Signature$Signif_ALSFTD] <- "ALS"
ALS_ALSFTD_WNN_L25_Signature$Signif_In[!ALS_ALSFTD_WNN_L25_Signature$Signif_ALS & ALS_ALSFTD_WNN_L25_Signature$Signif_ALSFTD] <- "ALSFTD"
ALS_ALSFTD_WNN_L25_Signature$Signif_In[ALS_ALSFTD_WNN_L25_Signature$Signif_ALS & ALS_ALSFTD_WNN_L25_Signature$Signif_ALSFTD] <- "Both"

table(ALS_ALSFTD_WNN_L25_Signature$Signif_ALS) 
table(ALS_ALSFTD_WNN_L25_Signature$Signif_ALSFTD)
table(ALS_ALSFTD_WNN_L25_Signature$Signif_In)

any(duplicated(ALS_ALSFTD_WNN_L25_Signature$Peak))



    ## 4.2 R3, R2 Signature ----------------------------------------------------

ALS_ALSFTD_WNN_L25_Rand_Signature <- data.frame(
  Peak = sort(
    unique(c(
      All_DARs_WNN_L25_R2, 
      All_DARs_WNN_L25_R3
    ))
  )
)

table(ALS_ALSFTD_WNN_L25_Rand_Signature$Peak %in% ATAC_Peak_Links$peak)
table(unique(ATAC_Peak_Links$peak) %in% ALS_ALSFTD_WNN_L25_Rand_Signature$Peak)

table(ALS_ALSFTD_WNN_L25_Rand_Signature$Peak %in% Peaks$ID2)


ALS_ALSFTD_WNN_L25_Rand_Signature$ID2 <- ALS_ALSFTD_WNN_L25_Rand_Signature$Peak
ALS_ALSFTD_WNN_L25_Rand_Signature$ID <- Peaks$ID[match(ALS_ALSFTD_WNN_L25_Rand_Signature$Peak, Peaks$ID2)]


ALS_ALSFTD_WNN_L25_Rand_Signature$Signif_R2 <- ALS_ALSFTD_WNN_L25_Rand_Signature$Peak %in% All_DARs_WNN_L25_R2
ALS_ALSFTD_WNN_L25_Rand_Signature$Signif_R3 <- ALS_ALSFTD_WNN_L25_Rand_Signature$Peak %in% All_DARs_WNN_L25_R3
ALS_ALSFTD_WNN_L25_Rand_Signature$Signif_In <- "NA"
ALS_ALSFTD_WNN_L25_Rand_Signature$Signif_In[ALS_ALSFTD_WNN_L25_Rand_Signature$Signif_R2 & !ALS_ALSFTD_WNN_L25_Rand_Signature$Signif_R3] <- "R2"
ALS_ALSFTD_WNN_L25_Rand_Signature$Signif_In[!ALS_ALSFTD_WNN_L25_Rand_Signature$Signif_R2 & ALS_ALSFTD_WNN_L25_Rand_Signature$Signif_R3] <- "R3"
ALS_ALSFTD_WNN_L25_Rand_Signature$Signif_In[ALS_ALSFTD_WNN_L25_Rand_Signature$Signif_R2 & ALS_ALSFTD_WNN_L25_Rand_Signature$Signif_R3] <- "Both"

table(ALS_ALSFTD_WNN_L25_Rand_Signature$Signif_R2) 
table(ALS_ALSFTD_WNN_L25_Rand_Signature$Signif_R3)
table(ALS_ALSFTD_WNN_L25_Rand_Signature$Signif_In)

any(duplicated(ALS_ALSFTD_WNN_L25_Rand_Signature$Peak))



    ## 4.3 Save data -----------------------------------------------------------

qsave(
  ALS_ALSFTD_WNN_L25_Signature, 
  paste0(
    "../Data/Annotations/Signatures/ATAC/", 
    "Signature_ATAC_ALS_ALSFTD_WNN_L25", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  ALS_ALSFTD_WNN_L25_Rand_Signature, 
  paste0(
    "../Data/Annotations/Signatures/ATAC/", 
    "Signature_ATAC_R2_R3_WNN_L25", 
    ".qrds"
  ), 
  nthr=nthr
)



