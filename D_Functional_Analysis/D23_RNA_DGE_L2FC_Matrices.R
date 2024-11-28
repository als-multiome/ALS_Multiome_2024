
  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 


library(qs) 
library(DESeq2)
library(GenomicRanges)



qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)




  ### 1.0 Load Data ------------------------------------------------------------ 



    ## 1.1. Load Data ----------------------------------------------------------

RNA_L2FC_Shrink_Results_12SVs_ALS <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/12_SVs/", 
    "L2FC_Shrink_Results_12SVs_ALS", 
    ".qrds"
  ), 
  nthr=nthr
)

RNA_L2FC_Shrink_Results_12SVs_ALSFTD <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/12_SVs/", 
    "L2FC_Shrink_Results_12SVs_ALSFTD", 
    ".qrds"
  ), 
  nthr=nthr
)

RNA_L2FC_Shrink_Results_12SVs_Index <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/12_SVs/", 
    "DESeq_Results_12SVs_Index", 
    ".qrds"
  ), 
  nthr=nthr
) 



ALS_ALSFTD_WNN_L25_Signature <- qread(
  paste0(
    "../Data/Annotations/Signatures/RNA/", 
    "ALS_ALSFTD_WNN_L25_Signature", 
    ".qrds"
  ), 
  nthr=nthr
)
table(ALS_ALSFTD_WNN_L25_Signature$Signif_In)



    ## 1.2 Define AUX Functions ------------------------------------------------

Get_DEGs <- function(x, alpha=0.05, sign="All"){
  x <- x[!is.na(x$padj),] 
  return(
    if(sign=="All") rownames(x)[x$padj<alpha] else
      if(sign=="Up") rownames(x)[x$padj<alpha & x$log2FoldChange>0] else
        if(sign=="Down") rownames(x)[x$padj<alpha & x$log2FoldChange<0]
  )
} 



Overlap_Table <- function(x, y){
  
   m <- table(
      !unique(c(x,y)) %in% x, 
      !unique(c(x,y)) %in% y
    ) |> matrix(
      nrow = 2, 
      ncol = 2, 
      byrow = FALSE, 
      dimnames = list(
        "In X" = c("Yes", "No"), 
        "In Y" = c("Yes", "No") 
      ) 
    ) 
  
   names(dimnames(m)) <- c(
     paste0("In ", deparse(substitute(x))), 
     paste0("In ", deparse(substitute(y)))
   )
   
   return(m)
}




  ### 3.0 ALS, ALSFTD WNN_L25 L2FC matrices with L2FC Shrinkage ----------------

WNN_L25_Celltypes <- unique(
  RNA_L2FC_Shrink_Results_12SVs_Index$CellType[
    RNA_L2FC_Shrink_Results_12SVs_Index$CellTypeLevel=="WNN_L25"
  ]
)



    ## 3.1 ALS, ALS_FTD Signature ----------------------------------------------


      # ALS_ALSFTD Signature with AllCases comparisons -------------------------

Mtx_ALS_ALSFTD_Sign_AllCases <- data.frame(
  Gene = ALS_ALSFTD_WNN_L25_Signature$ID_10X
)

for (Celltype in WNN_L25_Celltypes){
  ind <- which(
    RNA_L2FC_Shrink_Results_12SVs_Index$CellTypeLevel == "WNN_L25" & 
      RNA_L2FC_Shrink_Results_12SVs_Index$CellType == Celltype & 
      RNA_L2FC_Shrink_Results_12SVs_Index$Comparison == "All_Cases"
  )
  
  Mtx_ALS_ALSFTD_Sign_AllCases[[paste0("ALS_", Celltype)]] <- RNA_L2FC_Shrink_Results_12SVs_ALS[[ind]]$log2FoldChange[
    match(
      Mtx_ALS_ALSFTD_Sign_AllCases$Gene, 
      rownames(RNA_L2FC_Shrink_Results_12SVs_ALS[[ind]])
    )
  ]
  
  Mtx_ALS_ALSFTD_Sign_AllCases[[paste0("ALSFTD_", Celltype)]] <- RNA_L2FC_Shrink_Results_12SVs_ALSFTD[[ind]]$log2FoldChange[
    match(
      Mtx_ALS_ALSFTD_Sign_AllCases$Gene, 
      rownames(RNA_L2FC_Shrink_Results_12SVs_ALSFTD[[ind]])
    )
  ]
}

rownames(Mtx_ALS_ALSFTD_Sign_AllCases) <- Mtx_ALS_ALSFTD_Sign_AllCases$Gene
Mtx_ALS_ALSFTD_Sign_AllCases <- Mtx_ALS_ALSFTD_Sign_AllCases[,-which(colnames(Mtx_ALS_ALSFTD_Sign_AllCases)=="Gene")]
Mtx_ALS_ALSFTD_Sign_AllCases <- as.matrix(Mtx_ALS_ALSFTD_Sign_AllCases)

rm(Celltype, ind)


      # ALS_ALSFTD Signature with Rand comparisons -----------------------------

Mtx_ALS_ALSFTD_Sign_Rand <- data.frame(
  Gene = ALS_ALSFTD_WNN_L25_Signature$ID_10X
)

for (Celltype in WNN_L25_Celltypes){
  ind <- which(
    RNA_L2FC_Shrink_Results_12SVs_Index$CellTypeLevel == "WNN_L25" & 
      RNA_L2FC_Shrink_Results_12SVs_Index$CellType == Celltype & 
      RNA_L2FC_Shrink_Results_12SVs_Index$Comparison == "Rand"
  )
  
  Mtx_ALS_ALSFTD_Sign_Rand[[paste0("R2_", Celltype)]] <- RNA_L2FC_Shrink_Results_12SVs_ALS[[ind]]$log2FoldChange[
    match(
      Mtx_ALS_ALSFTD_Sign_Rand$Gene, 
      rownames(RNA_L2FC_Shrink_Results_12SVs_ALS[[ind]])
    )
  ]
  
  Mtx_ALS_ALSFTD_Sign_Rand[[paste0("R3_", Celltype)]] <- RNA_L2FC_Shrink_Results_12SVs_ALSFTD[[ind]]$log2FoldChange[
    match(
      Mtx_ALS_ALSFTD_Sign_Rand$Gene, 
      rownames(RNA_L2FC_Shrink_Results_12SVs_ALSFTD[[ind]])
    )
  ]
}

rownames(Mtx_ALS_ALSFTD_Sign_Rand) <- Mtx_ALS_ALSFTD_Sign_Rand$Gene
Mtx_ALS_ALSFTD_Sign_Rand <- Mtx_ALS_ALSFTD_Sign_Rand[,-which(colnames(Mtx_ALS_ALSFTD_Sign_Rand)=="Gene")]
Mtx_ALS_ALSFTD_Sign_Rand <- as.matrix(Mtx_ALS_ALSFTD_Sign_Rand)
rm(Celltype, ind)



    ## 3.2 RandGenes Signature -------------------------------------------------

Rand_Genes <- list()

for(CellType in WNN_L25_Celltypes){
  
  ind <- which(
    RNA_L2FC_Shrink_Results_12SVs_Index$CellTypeLevel == "WNN_L25" & 
      RNA_L2FC_Shrink_Results_12SVs_Index$CellType == CellType & 
      RNA_L2FC_Shrink_Results_12SVs_Index$Comparison == "All_Cases"
  )
  
  nGenes_ALS = length(
    Get_DEGs(
      RNA_L2FC_Shrink_Results_12SVs_ALS[[ind]]
    )
  ) 
  
  nGenes_ALSFTD = length(
    Get_DEGs(
      RNA_L2FC_Shrink_Results_12SVs_ALSFTD[[ind]]
    )
  ) 
  
  ind <- which(
    RNA_L2FC_Shrink_Results_12SVs_Index$CellTypeLevel == "WNN_L25" & 
      RNA_L2FC_Shrink_Results_12SVs_Index$CellType == CellType & 
      RNA_L2FC_Shrink_Results_12SVs_Index$Comparison == "Rand"
  )
  
  res_R2 <- RNA_L2FC_Shrink_Results_12SVs_ALS[[ind]] 
  res_R2 <- res_R2[order(res_R2$pvalue),] 
  
  res_R3 <- RNA_L2FC_Shrink_Results_12SVs_ALSFTD[[ind]] 
  res_R3 <- res_R3[order(res_R3$pvalue),]
  
  Rand_Genes[[paste0("R2_",  CellType)]] <- rownames(res_R2)[1:nGenes_ALS] 
  Rand_Genes[[paste0("R3_",  CellType)]] <- rownames(res_R3)[1:nGenes_ALSFTD] 
  
  
}

R3_R2_WNN_L25_Signature <- unique(
  as.character(
    unlist(Rand_Genes)
  )
)

rm(Rand_Genes, CellType, ind, res_R2, res_R3, nGenes_ALS, nGenes_ALSFTD)

Overlap_Table(ALS_ALSFTD_WNN_L25_Signature$ID_10X, R3_R2_WNN_L25_Signature)


      #  RandGenes Signature with AllCases comparisons -------------------------

Mtx_RandGenes_Sign_AllCases <- data.frame(
  Gene = R3_R2_WNN_L25_Signature
)

for (Celltype in WNN_L25_Celltypes){
  ind <- which(
    RNA_L2FC_Shrink_Results_12SVs_Index$CellTypeLevel == "WNN_L25" & 
      RNA_L2FC_Shrink_Results_12SVs_Index$CellType == Celltype & 
      RNA_L2FC_Shrink_Results_12SVs_Index$Comparison == "All_Cases"
  )
  
  Mtx_RandGenes_Sign_AllCases[[paste0("ALS_", Celltype)]] <- RNA_L2FC_Shrink_Results_12SVs_ALS[[ind]]$log2FoldChange[
    match(
      Mtx_RandGenes_Sign_AllCases$Gene, 
      rownames(RNA_L2FC_Shrink_Results_12SVs_ALS[[ind]])
    )
  ]
  
  Mtx_RandGenes_Sign_AllCases[[paste0("ALSFTD_", Celltype)]] <- RNA_L2FC_Shrink_Results_12SVs_ALSFTD[[ind]]$log2FoldChange[
    match(
      Mtx_RandGenes_Sign_AllCases$Gene, 
      rownames(RNA_L2FC_Shrink_Results_12SVs_ALSFTD[[ind]])
    )
  ]
}
rownames(Mtx_RandGenes_Sign_AllCases) <- Mtx_RandGenes_Sign_AllCases$Gene
Mtx_RandGenes_Sign_AllCases <- Mtx_RandGenes_Sign_AllCases[,-which(colnames(Mtx_RandGenes_Sign_AllCases)=="Gene")]
Mtx_RandGenes_Sign_AllCases <- as.matrix(Mtx_RandGenes_Sign_AllCases)
rm(Celltype, ind)


      # RandGenes Signature with Rand comparisons ------------------------------

Mtx_RandGenes_Sign_Rand <- data.frame(
  Gene = R3_R2_WNN_L25_Signature
)

for (Celltype in WNN_L25_Celltypes){
  ind <- which(
    RNA_L2FC_Shrink_Results_12SVs_Index$CellTypeLevel == "WNN_L25" & 
      RNA_L2FC_Shrink_Results_12SVs_Index$CellType == Celltype & 
      RNA_L2FC_Shrink_Results_12SVs_Index$Comparison == "Rand"
  )
  
  Mtx_RandGenes_Sign_Rand[[paste0("R2_", Celltype)]] <- RNA_L2FC_Shrink_Results_12SVs_ALS[[ind]]$log2FoldChange[
    match(
      Mtx_RandGenes_Sign_Rand$Gene, 
      rownames(RNA_L2FC_Shrink_Results_12SVs_ALS[[ind]])
    )
  ]
  
  Mtx_RandGenes_Sign_Rand[[paste0("R3_", Celltype)]] <- RNA_L2FC_Shrink_Results_12SVs_ALSFTD[[ind]]$log2FoldChange[
    match(
      Mtx_RandGenes_Sign_Rand$Gene, 
      rownames(RNA_L2FC_Shrink_Results_12SVs_ALSFTD[[ind]])
    )
  ]
}
rownames(Mtx_RandGenes_Sign_Rand) <- Mtx_RandGenes_Sign_Rand$Gene
Mtx_RandGenes_Sign_Rand <- Mtx_RandGenes_Sign_Rand[,-which(colnames(Mtx_RandGenes_Sign_Rand)=="Gene")]
Mtx_RandGenes_Sign_Rand <- as.matrix(Mtx_RandGenes_Sign_Rand)
rm(Celltype, ind)




  ### 4.0 Save Data ------------------------------------------------------------

qsave(
  R3_R2_WNN_L25_Signature, 
  paste0(
    "../Data/Annotations/Signatures/RNA/", 
    "R3_R2_WNN_L25_Signature", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  Mtx_ALS_ALSFTD_Sign_AllCases, 
  paste0(
    "../Data/Annotations/Signatures/RNA/", 
    "Mtx_ALS_ALSFTD_WNN_L25_Signature_AllCases", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  Mtx_ALS_ALSFTD_Sign_Rand, 
  paste0(
    "../Data/Annotations/Signatures/RNA/", 
    "Mtx_ALS_ALSFTD_WNN_L25_Signature_Rand", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  Mtx_RandGenes_Sign_AllCases, 
  paste0(
    "../Data/Annotations/Signatures/RNA/", 
    "Mtx_RandGenes_WNN_L25_Signature_AllCases", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  Mtx_RandGenes_Sign_Rand, 
  paste0(
    "../Data/Annotations/Signatures/RNA/", 
    "Mtx_RandGenes_WNN_L25_Signature_Rand", 
    ".qrds"
  ), 
  nthr=nthr
)

