

  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 

library(qs)  
library(tidyverse) 
library(DESeq2)
library(GenomicRanges)




qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)




  ### 1.0 Load data ------------------------------------------------------------

DESeq_Results_12SVs_ALS <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/12_SVs/", 
    "DESeq_Results_12SVs_ALS", 
    ".qrds"
  ), 
  nthr=nthr
)

DESeq_Results_12SVs_ALSFTD <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/12_SVs/", 
    "DESeq_Results_12SVs_ALSFTD", 
    ".qrds"
  ), 
  nthr=nthr
)

DESeq_Results_12SVs_Index <- qread(
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
  nthr
)

GEX_Features_ENSEMBL_ENTREZ <- qread(
  paste0(
    "../Data/Annotations/Genes/", 
    "GEX_Features_ENSEMBL_ENTREZ", 
    ".qrds"
  ), 
  nthr
)




  ### 2.0 Define AUX functions -------------------------------------------------

get_sign_genes_df <- function(res, fdr=0.05){
  res <- data.frame(res)
  res <- res[!is.na(res$padj),]
  res <- res[res$padj<fdr,]
  return(res)
}

get_sign_genes_names <- function(res, fdr=0.05){
  res <- data.frame(res)
  res <- res[!is.na(res$padj),]
  res <- res[res$padj<fdr,]
  return(rownames(res))
}


get_same_N_rand_comparison_genes_names <- function(res, N){
  res <- data.frame(res)
  res <- res[!is.na(res$pvalue),]
  res <- res[order(res$pvalue),]
  return(rownames(res)[1:N])
}




  ### 3.0 Get all DEGs that are significant in AllCells ------------------------

ind_AllCells_Case_Comparisons <- which(
  DESeq_Results_12SVs_Index$CellTypeLevel=="AllCells" & 
    DESeq_Results_12SVs_Index$Comparison=="All_Cases"
)

ind_AllCells_Rand_Comparisons <- which(
  DESeq_Results_12SVs_Index$CellTypeLevel=="AllCells" & 
    DESeq_Results_12SVs_Index$Comparison=="Rand"
)



    ## 3.1 ALS vs HC comparison ------------------------------------------------

All_DEGs_AllCells_ALS <- unique(
  unlist(
    lapply(
      DESeq_Results_12SVs_ALS[ind_AllCells_Case_Comparisons], 
      FUN=get_sign_genes_names
    )
  ) 
)



    ## 3.2 ALSFTD vs HC comparison ---------------------------------------------

All_DEGs_AllCells_ALSFTD <- unique(
  unlist(
    lapply(
      DESeq_Results_12SVs_ALSFTD[ind_AllCells_Case_Comparisons], 
      FUN=get_sign_genes_names
    )
  ) 
)



    ## 3.3 R2 vs R1 comparison -------------------------------------------------

Ns_DEGs_ALS <- unlist(
  lapply(
    DESeq_Results_12SVs_ALS[ind_AllCells_Case_Comparisons], 
    FUN=function(x){
      return(
        length(
          get_sign_genes_names(x)
        )
      )
    }
  )
)

i <- 0
All_DEGs_AllCells_R2 <- unique(
  unlist(
    lapply(
      DESeq_Results_12SVs_ALS[ind_AllCells_Rand_Comparisons], 
      FUN=function(x){
        i <<- i + 1
        return(
          get_same_N_rand_comparison_genes_names(
            res = x, 
            N=Ns_DEGs_ALS[i]
          ) 
        )
      }
    )
  ) 
)
rm(Ns_DEGs_ALS, i)



    ## 3.4 R3 vs R1 comparison -------------------------------------------------

Ns_DEGs_ALSFTD <- unlist(
  lapply(
    DESeq_Results_12SVs_ALSFTD[ind_AllCells_Case_Comparisons], 
    FUN=function(x){
      return(
        length(
          get_sign_genes_names(x)
        )
      )
    }
  )
)

i <- 0
All_DEGs_AllCells_R3 <- unique(
  unlist(
    lapply(
      DESeq_Results_12SVs_ALSFTD[ind_AllCells_Rand_Comparisons], 
      FUN=function(x){
        i <<- i + 1
        return(
          get_same_N_rand_comparison_genes_names(
            res = x, 
            N=Ns_DEGs_ALSFTD[i]
          ) 
        )
      }
    )
  ) 
)
rm(Ns_DEGs_ALSFTD, i) 



    ## 3.5 Sanity checks -------------------------------------------------------

table(All_DEGs_AllCells_ALS %in% All_DEGs_AllCells_ALSFTD)
table(All_DEGs_AllCells_ALSFTD %in% All_DEGs_AllCells_ALS)

table(All_DEGs_AllCells_R2 %in% All_DEGs_AllCells_R3)
table(All_DEGs_AllCells_R3 %in% All_DEGs_AllCells_R2)




  ### 4.0 Compile and save data ------------------------------------------------



    ## 4.1 ALS, ALSFTD Signature -----------------------------------------------

ALS_ALSFTD_AllCells_Signature <- data.frame(
  Gene = sort(
    unique(c(
      All_DEGs_AllCells_ALS, 
      All_DEGs_AllCells_ALSFTD
    ))
  )
)

ALS_ALSFTD_AllCells_Signature$ID_10X <- ALS_ALSFTD_AllCells_Signature$Gene 
ALS_ALSFTD_AllCells_Signature$SYMBOL <- GEX_Features$SYMBOL[match(ALS_ALSFTD_AllCells_Signature$ID_10X, GEX_Features$ID_10X)]
ALS_ALSFTD_AllCells_Signature$GENETYPE <- GEX_Features$GENETYPE[match(ALS_ALSFTD_AllCells_Signature$ID_10X, GEX_Features$ID_10X)]


ALS_ALSFTD_AllCells_Signature$Signif_ALS <- ALS_ALSFTD_AllCells_Signature$Gene %in% All_DEGs_AllCells_ALS
ALS_ALSFTD_AllCells_Signature$Signif_ALSFTD <- ALS_ALSFTD_AllCells_Signature$Gene %in% All_DEGs_AllCells_ALSFTD
ALS_ALSFTD_AllCells_Signature$Signif_In <- "NA"
ALS_ALSFTD_AllCells_Signature$Signif_In[ALS_ALSFTD_AllCells_Signature$Signif_ALS & !ALS_ALSFTD_AllCells_Signature$Signif_ALSFTD] <- "ALS"
ALS_ALSFTD_AllCells_Signature$Signif_In[!ALS_ALSFTD_AllCells_Signature$Signif_ALS & ALS_ALSFTD_AllCells_Signature$Signif_ALSFTD] <- "ALSFTD"
ALS_ALSFTD_AllCells_Signature$Signif_In[ALS_ALSFTD_AllCells_Signature$Signif_ALS & ALS_ALSFTD_AllCells_Signature$Signif_ALSFTD] <- "Both"

table(ALS_ALSFTD_AllCells_Signature$Signif_ALS) 
table(ALS_ALSFTD_AllCells_Signature$Signif_ALSFTD)
table(ALS_ALSFTD_AllCells_Signature$Signif_In)

any(duplicated(ALS_ALSFTD_AllCells_Signature$Gene))



    ## 4.2 R3, R2 Signature ----------------------------------------------------

ALS_ALSFTD_AllCells_Rand_Signature <- data.frame(
  Gene = sort(
    unique(c(
      All_DEGs_AllCells_R2, 
      All_DEGs_AllCells_R3
    ))
  )
)

ALS_ALSFTD_AllCells_Rand_Signature$ID_10X <- ALS_ALSFTD_AllCells_Rand_Signature$Gene 
ALS_ALSFTD_AllCells_Rand_Signature$SYMBOL <- GEX_Features$SYMBOL[match(ALS_ALSFTD_AllCells_Rand_Signature$ID_10X, GEX_Features$ID_10X)]
ALS_ALSFTD_AllCells_Rand_Signature$GENETYPE <- GEX_Features$GENETYPE[match(ALS_ALSFTD_AllCells_Rand_Signature$ID_10X, GEX_Features$ID_10X)]


ALS_ALSFTD_AllCells_Rand_Signature$Signif_R2 <- ALS_ALSFTD_AllCells_Rand_Signature$Gene %in% All_DEGs_AllCells_R2
ALS_ALSFTD_AllCells_Rand_Signature$Signif_R3 <- ALS_ALSFTD_AllCells_Rand_Signature$Gene %in% All_DEGs_AllCells_R3
ALS_ALSFTD_AllCells_Rand_Signature$Signif_In <- "NA"
ALS_ALSFTD_AllCells_Rand_Signature$Signif_In[ALS_ALSFTD_AllCells_Rand_Signature$Signif_R2 & !ALS_ALSFTD_AllCells_Rand_Signature$Signif_R3] <- "R2"
ALS_ALSFTD_AllCells_Rand_Signature$Signif_In[!ALS_ALSFTD_AllCells_Rand_Signature$Signif_R2 & ALS_ALSFTD_AllCells_Rand_Signature$Signif_R3] <- "R3"
ALS_ALSFTD_AllCells_Rand_Signature$Signif_In[ALS_ALSFTD_AllCells_Rand_Signature$Signif_R2 & ALS_ALSFTD_AllCells_Rand_Signature$Signif_R3] <- "Both"

table(ALS_ALSFTD_AllCells_Rand_Signature$Signif_R2) 
table(ALS_ALSFTD_AllCells_Rand_Signature$Signif_R3)
table(ALS_ALSFTD_AllCells_Rand_Signature$Signif_In)

any(duplicated(ALS_ALSFTD_AllCells_Rand_Signature$Gene))



    ## 4.3 Save data -----------------------------------------------------------

qsave(
  ALS_ALSFTD_AllCells_Signature, 
  paste0(
    "../Data/Annotations/Signatures/RNA/", 
    "Signature_RNA_ALS_ALSFTD_AllCells", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  ALS_ALSFTD_AllCells_Rand_Signature, 
  paste0(
    "../Data/Annotations/Signatures/RNA/", 
    "Signature_RNA_R2_R3_AllCells", 
    ".qrds"
  ), 
  nthr=nthr
)
