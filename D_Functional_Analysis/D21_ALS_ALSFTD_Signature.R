

  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 

library(qs)  
library(tidyverse) 
library(DESeq2)
library(GenomicRanges)
library(GenomeInfoDb)




qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)




  ### 1.0 Load data ------------------------------------------------------------

ALS_ALSFTD_WNN_L25_Signature <- qread(
  paste0(
    "../Data/Annotations/Signatures/RNA/", 
    "ALS_ALSFTD_WNN_L25_Signature", 
    ".qrds"
  ), 
  nthr=nthr
)

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


fisher.test(
  matrix(
    c(
      346, 5125, 861, 35367-346-5125-861
    ), 
    2, 2
    
  )
)

EpiSig_41Genes <- fread("../../Bulk_ATAC_ALS_PBMCs/EpiChrome_41Core_Genes.txt", data.table = FALSE, header=FALSE)

table(EpiSig_41Genes$V1 %in% EpiSig$Gene)
table(EpiSig_41Genes$V1 %in% ALS_ALSFTD_WNN_L25_Signature$Gene) 
table(ALS_ALSFTD_WNN_L25_Signature$Gene %in% EpiSig_41Genes$V1)
EpiSig_41Genes$V1[EpiSig_41Genes$V1 %in% ALS_ALSFTD_WNN_L25_Signature$Gene] 
fisher.test(
  matrix(
    c(
      34, 5437, 72, 35367-34-5437-72
    ), 
    2, 2
    
  )
)

ALS_Sign_510 <- ALS_ALSFTD_WNN_L25_Signature$SYMBOL[ALS_ALSFTD_WNN_L25_Signature$Signif_In=="Both"]
table(EpiSig_Genes %in% ALS_Sign_510)
table(ALS_Sign_510 %in% EpiSig_Genes)

fisher.test(
  matrix(
    c(
      31, 1176, 479, 35367-31-1176-479
    ), 
    2, 2
    
  )
)

table(EpiSig_41Genes$V1 %in% ALS_Sign_510)
table(ALS_Sign_510 %in% EpiSig_41Genes$V1)

EpiSig_41Genes$V1[EpiSig_41Genes$V1 %in% ALS_Sign_510]


fisher.test(
  matrix(
    c(
      4, 102, 506, 35367-4-102-506
    ), 
    2, 2
    
  )
)

  ### 3.0 Collect DGE stats for each comparison --------------------------------

ALS_Pseudobulk_Sign <- get_sign_genes_names(DESeq_Results_12SVs_ALS[[1]]) 
ALSFTD_Pseudobulk_Sign <- get_sign_genes_names(DESeq_Results_12SVs_ALSFTD[[1]])

ALS_ALSFTD_Pseudobulk_Sign <- unique(c(
  ALS_Pseudobulk_Sign, 
  ALSFTD_Pseudobulk_Sign
))

table(ALS_Pseudobulk_Sign %in% ALS_ALSFTD_WNN_L25_Signature$Gene)
table(ALS_ALSFTD_WNN_L25_Signature$Gene %in% ALS_Pseudobulk_Sign)

write.csv(ALS_ALSFTD_WNN_L25_Signature, "../ALS_Signa.csv", quote=FALSE)

table(ALS_Pseudobulk_Sign[!ALS_Pseudobulk_Sign %in% ALS_ALSFTD_WNN_L25_Signature$Gene] %in% markers)
table(ALS_Pseudobulk_Sign[ALS_Pseudobulk_Sign %in% ALS_ALSFTD_WNN_L25_Signature$Gene] %in% markers)



WNN_L25_HC_Markers <- qread("../Data/Annotations/CellType_Markers/WNN_L25_Markers_HC_Wilcox.qrds")
markers <- sort(unique(as.character(unlist(sapply(WNN_L25_HC_Markers, FUN=function(x){return(x$SYMBOL[x$p_val_adj<1e-10 & x$avg_log2FC>2])})))))



      # Collect WNN_L25 comparison indices -------------------------------------

ind_WNN_L25_Case_Comparisons <- which(
  DESeq_Results_12SVs_Index$CellTypeLevel=="WNN_L25" & 
    DESeq_Results_12SVs_Index$Comparison=="All_Cases"
)

ind_WNN_L25_Rand_Comparisons <- which(
  DESeq_Results_12SVs_Index$CellTypeLevel=="WNN_L25" & 
    DESeq_Results_12SVs_Index$Comparison=="Rand"
)



    ## 1.1 Collect DGE stats ---------------------------------------------------

for (ind in ind_WNN_L25_Case_Comparisons){
  
  ALS_ALSFTD_WNN_L25_Signature[[
    paste0(
      "LOG2FC_", 
      DESeq_Results_12SVs_Index$CellType[ind], 
      "_ALS"
    )
  ]] <- DESeq_Results_12SVs_ALS[[ind]]$log2FoldChange[match(
    ALS_ALSFTD_WNN_L25_Signature$Gene, rownames(DESeq_Results_12SVs_ALS[[ind]])
  )]  
  
  ALS_ALSFTD_WNN_L25_Signature[[
    paste0(
      "LOG2FC_", 
      DESeq_Results_12SVs_Index$CellType[ind], 
      "_ALSFTD"
    )
  ]] <- DESeq_Results_12SVs_ALSFTD[[ind]]$log2FoldChange[match(
    ALS_ALSFTD_WNN_L25_Signature$Gene, rownames(DESeq_Results_12SVs_ALSFTD[[ind]])
  )]  
  
  
  ALS_ALSFTD_WNN_L25_Signature[[
    paste0(
      "PVAL_", 
      DESeq_Results_12SVs_Index$CellType[ind], 
      "_ALS"
    )
  ]] <- DESeq_Results_12SVs_ALS[[ind]]$pvalue[match(
    ALS_ALSFTD_WNN_L25_Signature$Gene, rownames(DESeq_Results_12SVs_ALS[[ind]])
  )]  
  
  ALS_ALSFTD_WNN_L25_Signature[[
    paste0(
      "PVAL_", 
      DESeq_Results_12SVs_Index$CellType[ind], 
      "_ALSFTD"
    )
  ]] <- DESeq_Results_12SVs_ALSFTD[[ind]]$pvalue[match(
    ALS_ALSFTD_WNN_L25_Signature$Gene, rownames(DESeq_Results_12SVs_ALSFTD[[ind]])
  )] 
  
  
  ALS_ALSFTD_WNN_L25_Signature[[
    paste0(
      "QVAL_", 
      DESeq_Results_12SVs_Index$CellType[ind], 
      "_ALS"
    )
  ]] <- DESeq_Results_12SVs_ALS[[ind]]$padj[match(
    ALS_ALSFTD_WNN_L25_Signature$Gene, rownames(DESeq_Results_12SVs_ALS[[ind]])
  )]  
  
  ALS_ALSFTD_WNN_L25_Signature[[
    paste0(
      "QVAL_", 
      DESeq_Results_12SVs_Index$CellType[ind], 
      "_ALSFTD"
    )
  ]] <- DESeq_Results_12SVs_ALSFTD[[ind]]$padj[match(
    ALS_ALSFTD_WNN_L25_Signature$Gene, rownames(DESeq_Results_12SVs_ALSFTD[[ind]])
  )] 
  
}



    ## 1.2 Collect DGE stats for Random control --------------------------------

ALS_ALSFTD_WNN_L25_Signature_Random_Control <- ALS_ALSFTD_WNN_L25_Signature

for (ind in ind_WNN_L25_Rand_Comparisons){
  
  ALS_ALSFTD_WNN_L25_Signature_Random_Control[[
    paste0(
      "LOG2FC_", 
      DESeq_Results_12SVs_Index$CellType[ind], 
      "_ALS_Rand"
    )
  ]] <- DESeq_Results_12SVs_ALS[[ind]]$log2FoldChange[match(
    ALS_ALSFTD_WNN_L25_Signature_Random_Control$Gene, rownames(DESeq_Results_12SVs_ALS[[ind]])
  )]  
  
  ALS_ALSFTD_WNN_L25_Signature_Random_Control[[
    paste0(
      "LOG2FC_", 
      DESeq_Results_12SVs_Index$CellType[ind], 
      "_ALSFTD_Rand"
    )
  ]] <- DESeq_Results_12SVs_ALSFTD[[ind]]$log2FoldChange[match(
    ALS_ALSFTD_WNN_L25_Signature_Random_Control$Gene, rownames(DESeq_Results_12SVs_ALSFTD[[ind]])
  )]  
  
  
  ALS_ALSFTD_WNN_L25_Signature_Random_Control[[
    paste0(
      "PVAL_", 
      DESeq_Results_12SVs_Index$CellType[ind], 
      "_ALS_Rand"
    )
  ]] <- DESeq_Results_12SVs_ALS[[ind]]$pvalue[match(
    ALS_ALSFTD_WNN_L25_Signature_Random_Control$Gene, rownames(DESeq_Results_12SVs_ALS[[ind]])
  )]  
  
  ALS_ALSFTD_WNN_L25_Signature_Random_Control[[
    paste0(
      "PVAL_", 
      DESeq_Results_12SVs_Index$CellType[ind], 
      "_ALSFTD_Rand"
    )
  ]] <- DESeq_Results_12SVs_ALSFTD[[ind]]$pvalue[match(
    ALS_ALSFTD_WNN_L25_Signature_Random_Control$Gene, rownames(DESeq_Results_12SVs_ALSFTD[[ind]])
  )] 
  
  
  ALS_ALSFTD_WNN_L25_Signature_Random_Control[[
    paste0(
      "QVAL_", 
      DESeq_Results_12SVs_Index$CellType[ind], 
      "_ALS_Rand"
    )
  ]] <- DESeq_Results_12SVs_ALS[[ind]]$padj[match(
    ALS_ALSFTD_WNN_L25_Signature_Random_Control$Gene, rownames(DESeq_Results_12SVs_ALS[[ind]])
  )]  
  
  ALS_ALSFTD_WNN_L25_Signature_Random_Control[[
    paste0(
      "QVAL_", 
      DESeq_Results_12SVs_Index$CellType[ind], 
      "_ALSFTD_Rand"
    )
  ]] <- DESeq_Results_12SVs_ALSFTD[[ind]]$padj[match(
    ALS_ALSFTD_WNN_L25_Signature_Random_Control$Gene, rownames(DESeq_Results_12SVs_ALSFTD[[ind]])
  )] 
  
}



    ## 1.1 Collect DGE stats in long format for plotting -----------------------

# Correlation to Pseudobulk 




celltypes <- colnames(ALS_ALSFTD_WNN_L25_Signature)[grep("LOG2FC_", colnames(ALS_ALSFTD_WNN_L25_Signature))]
celltypes <- str_split(celltypes, "LOG2FC_", simplif=TRUE)[,2]
celltypes <- str_split(celltypes, "_ALS", simplify=TRUE)[,1]
celltypes <- unique(celltypes)

df <- NULL 
for (celltype in celltypes){
 
  if (is.null(df)){
    df = data.frame(
      Gene = ALS_ALSFTD_WNN_L25_Signature$Gene, 
      CellType = celltype, 
      LOG2FC_ALS = ALS_ALSFTD_WNN_L25_Signature[[paste0("LOG2FC_", celltype, "_ALS")]], 
      LOG2FC_ALSFTD = ALS_ALSFTD_WNN_L25_Signature[[paste0("LOG2FC_", celltype, "_ALSFTD")]], 
      SIGN_IN_BOTH = ALS_ALSFTD_WNN_L25_Signature[[paste0("QVAL_", celltype, "_ALS")]] < 0.05 & ALS_ALSFTD_WNN_L25_Signature[[paste0("QVAL_", celltype, "_ALSFTD")]] < 0.05, 
      SIGN_IN_ATLEASTONE = ALS_ALSFTD_WNN_L25_Signature[[paste0("QVAL_", celltype, "_ALS")]] < 0.05 | ALS_ALSFTD_WNN_L25_Signature[[paste0("QVAL_", celltype, "_ALSFTD")]] < 0.05
    )
  } else {
    df = rbind(
      df, 
      data.frame(
        Gene = ALS_ALSFTD_WNN_L25_Signature$Gene, 
        CellType = celltype, 
        LOG2FC_ALS = ALS_ALSFTD_WNN_L25_Signature[[paste0("LOG2FC_", celltype, "_ALS")]], 
        LOG2FC_ALSFTD = ALS_ALSFTD_WNN_L25_Signature[[paste0("LOG2FC_", celltype, "_ALSFTD")]], 
        SIGN_IN_BOTH = ALS_ALSFTD_WNN_L25_Signature[[paste0("QVAL_", celltype, "_ALS")]] < 0.05 & ALS_ALSFTD_WNN_L25_Signature[[paste0("QVAL_", celltype, "_ALSFTD")]] < 0.05, 
        SIGN_IN_ATLEASTONE = ALS_ALSFTD_WNN_L25_Signature[[paste0("QVAL_", celltype, "_ALS")]] < 0.05 | ALS_ALSFTD_WNN_L25_Signature[[paste0("QVAL_", celltype, "_ALSFTD")]] < 0.05
      )
    )
  }
}

df_rand <- NULL 
for (celltype in celltypes){
  
  if (is.null(df_rand)){
    df_rand = data.frame(
      Gene = ALS_ALSFTD_WNN_L25_Signature_Random_Control$Gene, 
      CellType = celltype, 
      LOG2FC_ALS_Rand = ALS_ALSFTD_WNN_L25_Signature_Random_Control[[paste0("LOG2FC_", celltype, "_ALS_Rand")]], 
      LOG2FC_ALSFTD_Rand = ALS_ALSFTD_WNN_L25_Signature_Random_Control[[paste0("LOG2FC_", celltype, "_ALSFTD_Rand")]]
    )
  } else {
    df_rand = rbind(
      df_rand, 
      data.frame(
        Gene = ALS_ALSFTD_WNN_L25_Signature_Random_Control$Gene, 
        CellType = celltype, 
        LOG2FC_ALS_Rand = ALS_ALSFTD_WNN_L25_Signature_Random_Control[[paste0("LOG2FC_", celltype, "_ALS_Rand")]], 
        LOG2FC_ALSFTD_Rand = ALS_ALSFTD_WNN_L25_Signature_Random_Control[[paste0("LOG2FC_", celltype, "_ALSFTD_Rand")]]
      )
    )
  }
} 

table(df$SIGN_IN_ATLEASTONE)
table(df$SIGN_IN_ATLEASTONE, df$SIGN_IN_BOTH)

df <- df[which(df$SIGN_IN_ATLEASTONE),]
plot(df$LOG2FC_ALS, df$LOG2FC_ALSFTD, xlim=c(-5,5), ylim=c(-5,5))
cor.test(df$LOG2FC_ALS, df$LOG2FC_ALSFTD, method="spearman")


df_rand <- df_rand[match(paste0(df$Gene, "_", df$CellType), paste0(df_rand$Gene, "_", df_rand$CellType)),]

plot(df_rand$LOG2FC_ALS_Rand, df_rand$LOG2FC_ALSFTD_Rand, xlim=c(-5,5), ylim=c(-5,5))
cor.test(df_rand$LOG2FC_ALS_Rand, df_rand$LOG2FC_ALSFTD_Rand, method="spearman")


df <- df[which(df$SIGN_IN_BOTH),]

plot(df$LOG2FC_ALS, df$LOG2FC_ALSFTD, xlim=c(-5,5), ylim=c(-5,5))
cor.test(df$LOG2FC_ALS, df$LOG2FC_ALSFTD, method="spearman")

df_rand <- df_rand[match(paste0(df$Gene, "_", df$CellType), paste0(df_rand$Gene, "_", df_rand$CellType)),]

plot(df_rand$LOG2FC_ALS_Rand, df_rand$LOG2FC_ALSFTD_Rand, xlim=c(-5,5), ylim=c(-5,5))
cor.test(df_rand$LOG2FC_ALS_Rand, df_rand$LOG2FC_ALSFTD_Rand, method="spearman")

ggplot(df) + 
  aes(LOG2FC_ALS, LOG2FC_ALSFTD, fill=CellType) + 
  geom_point(pch=21)

colnames(ALS_ALSFTD_WNN_L25_Signature)[grep("LOG2FC_", colnames(ALS_ALSFTD_WNN_L25_Signature))]

write.csv(sort(unique(df$Gene)), "../Teeest.csv", quote=FALSE)

data.frame(
  Gene = ALS_ALSFTD_WNN_L25_Signature$Gene, 
  CellType = 
)

str_split(
  str_split(
    colnames(ALS_ALSFTD_WNN_L25_Signature)[grep("LOG2FC_", colnames(ALS_ALSFTD_WNN_L25_Signature))], 
    "LOG2FC_", 
    simplify = TRUE
  )[,2], 
  "_ALS", 
  simplify = TRUE
)[,1]

  ### 4.0 Compile and save data ------------------------------------------------

ALS_ALSFTD_WNN_L25_Signature <- data.frame(
  Gene = sort(
    unique(c(
      All_DEGs_WNN_L25_ALS, 
      All_DEGs_WNN_L25_ALSFTD
    ))
  )
)

ALS_ALSFTD_WNN_L25_Signature$SYMBOL <- ALS_ALSFTD_WNN_L25_Signature$Gene 
ALS_ALSFTD_WNN_L25_Signature$ENSG <- features$ENSG[match(ALS_ALSFTD_WNN_L25_Signature$Gene, features$Symbol)] 
ALS_ALSFTD_WNN_L25_Signature$Biotype <- features$Gene_biotype[match(ALS_ALSFTD_WNN_L25_Signature$Gene, features$Symbol)] 
table(is.na(ALS_ALSFTD_WNN_L25_Signature$ENSG))

ALS_ALSFTD_WNN_L25_Signature$Signif_ALS <- ALS_ALSFTD_WNN_L25_Signature$Gene %in% All_DEGs_WNN_L25_ALS
ALS_ALSFTD_WNN_L25_Signature$Signif_ALSFTD <- ALS_ALSFTD_WNN_L25_Signature$Gene %in% All_DEGs_WNN_L25_ALSFTD
ALS_ALSFTD_WNN_L25_Signature$Signif_In <- "NA"
ALS_ALSFTD_WNN_L25_Signature$Signif_In[ALS_ALSFTD_WNN_L25_Signature$Signif_ALS & !ALS_ALSFTD_WNN_L25_Signature$Signif_ALSFTD] <- "ALS"
ALS_ALSFTD_WNN_L25_Signature$Signif_In[!ALS_ALSFTD_WNN_L25_Signature$Signif_ALS & ALS_ALSFTD_WNN_L25_Signature$Signif_ALSFTD] <- "ALSFTD"
ALS_ALSFTD_WNN_L25_Signature$Signif_In[ALS_ALSFTD_WNN_L25_Signature$Signif_ALS & ALS_ALSFTD_WNN_L25_Signature$Signif_ALSFTD] <- "Both"

table(ALS_ALSFTD_WNN_L25_Signature$Signif_ALS) 
table(ALS_ALSFTD_WNN_L25_Signature$Signif_ALSFTD)
table(ALS_ALSFTD_WNN_L25_Signature$Signif_In)

any(duplicated(ALS_ALSFTD_WNN_L25_Signature$Gene))

qsave(
  ALS_ALSFTD_WNN_L25_Signature, 
  paste0(
    "../Data/Annotations/Signatures/RNA/", 
    "ALS_ALSFTD_WNN_L25_Signature", 
    ".qrds"
  ), 
  nthr=nthr
)

