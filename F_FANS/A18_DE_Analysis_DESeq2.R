# A18_DE_Analysis_DESeq2.R




  ### 0.0 Load libraries ------------------------------------------------------- 

library(qs2)
library(tidyverse)
library(Seurat)
library(DESeq2)
library(BiocParallel)
library(apeglm)
  

  
  
  ### 1.0 Load data ------------------------------------------------------------

FANS <- qs_read(
    "../Data/FANS/SeuratObjects/FANS.Unsubsetted.qs2"
)

FANS.SubsetCells <- qs_read(
  "../Data/FANS/SeuratObjects/FANS.SubsetCells.qs2"
)

FANS.SubsetCounts <- qs_read(
  "../Data/FANS/SeuratObjects/FANS.SubsetCounts.qs2"
)

Samples <- qs_read(
  "../Data/FANS/Samples/Samples.qs2"
)

SampleClusters <- qs_read(
  "../Data/FANS/Samples/Sample_Clusters.qs2"
)




  ### 2.0 Define DESeq2 Workflow -----------------------------------------------

DESeq2_DE <- function(
    SeuratObject  
){
  
  SeuratObject$Pseudobulk_ID <- paste0(
    SeuratObject$Sample_donor, 
    "-", 
    SeuratObject$TDP43
  )
  AggrCounts <- AggregateExpression(
    SeuratObject, 
    assays = "RNA", 
    group.by = "Pseudobulk_ID", 
    slot = "counts"
  )
  AggrCounts <- AggrCounts$RNA 
  
  
  colData <- SeuratObject@meta.data %>% 
    group_by(Pseudobulk_ID) %>% 
    summarize(
      TDP43 = dplyr::first(TDP43), 
      Sex = dplyr::first(Sex), 
      PctMito = mean(percent.mt)
    )
  colData$Pseudobulk_ID <- str_replace_all(
    colData$Pseudobulk_ID, 
    "_", 
    "-"
  )
  colData <- as.data.frame(colData)
  rownames(colData) <- colData$Pseudobulk_ID
  
  all(colnames(AggrCounts) == rownames(colData)) 
  colData$TDP43 <- factor(colData$TDP43, levels = c("High", "Low"))
  colData$Sex <- factor(colData$Sex, levels = c("f", "m")) 
  colData$PctMito <- scale(colData$PctMito)
  
  dds <- DESeqDataSetFromMatrix(
    countData = AggrCounts, 
    colData = colData, 
    design = ~ Sex + PctMito + TDP43
    
  )
  DDS <- DESeq(dds)
  res <- results(DDS, alpha=0.05)  
  
  print(summary(res, alpha = 0.05))
  
  res_L2FC_Shrink <- lfcShrink(
    DDS, 
    coef = "TDP43_Low_vs_High", 
    res=res, 
    type="apeglm", 
    parallel = TRUE
  ) 
  
  res.df = data.frame(res)
  res.df$Gene = rownames(res.df)
  
  res_L2FC_Shrink.df = data.frame(res_L2FC_Shrink) 
  res_L2FC_Shrink.df$Gene = rownames(res_L2FC_Shrink.df)
  return(
    list(
      AggrCounts = AggrCounts, 
      colData = colData, 
      dds = dds, 
      DDS = DDS, 
      results = res, 
      results_L2FC_Shrink = res_L2FC_Shrink, 
      res = res.df, 
      res_L2FC_Shrink = res_L2FC_Shrink.df
    )
  )

}




  ### 3.0 DESeq2 DE Analysis ---------------------------------------------------



    ## 3.1 All Cells DE Analysis -----------------------------------------------

FANS_DE_DESeq2_Unsubsetted <- DESeq2_DE(
  FANS
)

FANS_DE_DESeq2_SubsetCells <- DESeq2_DE(
  FANS.SubsetCells
)

FANS_DE_DESeq2_SubsetCounts <- DESeq2_DE(
  FANS.SubsetCounts
) 


qs_save(
  FANS_DE_DESeq2_Unsubsetted, 
  "../Data/FANS/DE/FANS_DE_DESeq2_Unsubsetted.qs2"
)

qs_save(
  FANS_DE_DESeq2_SubsetCells, 
  "../Data/FANS/DE/FANS_DE_DESeq2_SubsetCells.qs2"
) 

qs_save(
  FANS_DE_DESeq2_SubsetCounts, 
  "../Data/FANS/DE/FANS_DE_DESeq2_SubsetCounts.qs2"
)



    ## 3.2 Exc_LINC00507 DE Analysis -------------------------------------------

LINC00507.Unsubsetted <- subset(FANS, subset = ID_WNN_L25_Predicted == "Exc_LINC00507")
LINC00507.SubsetCells <- subset(FANS.SubsetCells, subset = ID_WNN_L25_Predicted == "Exc_LINC00507")
LINC00507.SubsetCounts <- subset(FANS.SubsetCounts, subset = ID_WNN_L25_Predicted == "Exc_LINC00507")

LINC00507_DE_DESeq2_Unsubsetted <- DESeq2_DE(
  LINC00507.Unsubsetted
)

LINC00507_DE_DESeq2_SubsetCells <- DESeq2_DE(
  LINC00507.SubsetCells
)

LINC00507_DE_DESeq2_SubsetCounts <- DESeq2_DE(
  LINC00507.SubsetCounts
) 


qs_save(
  LINC00507_DE_DESeq2_Unsubsetted, 
  "../Data/FANS/DE/LINC00507_DE_DESeq2_Unsubsetted.qs2"
)

qs_save(
  LINC00507_DE_DESeq2_SubsetCells, 
  "../Data/FANS/DE/LINC00507_DE_DESeq2_SubsetCells.qs2"
) 

qs_save(
  LINC00507_DE_DESeq2_SubsetCounts, 
  "../Data/FANS/DE/LINC00507_DE_DESeq2_SubsetCounts.qs2"
)



    ## 3.3 Exc_THEMIS DE Analysis ----------------------------------------------

THEMIS.Unsubsetted <- subset(FANS, subset = ID_WNN_L25_Predicted == "Exc_THEMIS")
THEMIS.SubsetCells <- subset(FANS.SubsetCells, subset = ID_WNN_L25_Predicted == "Exc_THEMIS")
THEMIS.SubsetCounts <- subset(FANS.SubsetCounts, subset = ID_WNN_L25_Predicted == "Exc_THEMIS")

THEMIS_DE_DESeq2_Unsubsetted <- DESeq2_DE(
  THEMIS.Unsubsetted
)

THEMIS_DE_DESeq2_SubsetCells <- DESeq2_DE(
  THEMIS.SubsetCells
)

THEMIS_DE_DESeq2_SubsetCounts <- DESeq2_DE(
  THEMIS.SubsetCounts
) 


qs_save(
  THEMIS_DE_DESeq2_Unsubsetted, 
  "../Data/FANS/DE/THEMIS_DE_DESeq2_Unsubsetted.qs2"
)

qs_save(
  THEMIS_DE_DESeq2_SubsetCells, 
  "../Data/FANS/DE/THEMIS_DE_DESeq2_SubsetCells.qs2"
) 

qs_save(
  THEMIS_DE_DESeq2_SubsetCounts, 
  "../Data/FANS/DE/THEMIS_DE_DESeq2_SubsetCounts.qs2"
)



    ## 3.4 Exc_RORB DE Analysis ------------------------------------------------

RORB.Unsubsetted <- subset(FANS, subset = ID_WNN_L25_Predicted == "Exc_RORB")
RORB.SubsetCells <- subset(FANS.SubsetCells, subset = ID_WNN_L25_Predicted == "Exc_RORB")
RORB.SubsetCounts <- subset(FANS.SubsetCounts, subset = ID_WNN_L25_Predicted == "Exc_RORB")

RORB_DE_DESeq2_Unsubsetted <- DESeq2_DE(
  RORB.Unsubsetted
)

RORB_DE_DESeq2_SubsetCells <- DESeq2_DE(
  RORB.SubsetCells
)

RORB_DE_DESeq2_SubsetCounts <- DESeq2_DE(
  RORB.SubsetCounts
) 


qs_save(
  RORB_DE_DESeq2_Unsubsetted, 
  "../Data/FANS/DE/RORB_DE_DESeq2_Unsubsetted.qs2"
)

qs_save(
  RORB_DE_DESeq2_SubsetCells, 
  "../Data/FANS/DE/RORB_DE_DESeq2_SubsetCells.qs2"
) 

qs_save(
  RORB_DE_DESeq2_SubsetCounts, 
  "../Data/FANS/DE/RORB_DE_DESeq2_SubsetCounts.qs2"
)



    ## 3.5 Exc_FEZF2 DE Analysis -----------------------------------------------

FEZF2.Unsubsetted <- subset(FANS, subset = ID_WNN_L25_Predicted == "Exc_FEZF2")
FEZF2.SubsetCells <- subset(FANS.SubsetCells, subset = ID_WNN_L25_Predicted == "Exc_FEZF2")
FEZF2.SubsetCounts <- subset(FANS.SubsetCounts, subset = ID_WNN_L25_Predicted == "Exc_FEZF2")

FEZF2_DE_DESeq2_Unsubsetted <- DESeq2_DE(
  FEZF2.Unsubsetted
)

FEZF2_DE_DESeq2_SubsetCells <- DESeq2_DE(
  FEZF2.SubsetCells
)

FEZF2_DE_DESeq2_SubsetCounts <- DESeq2_DE(
  FEZF2.SubsetCounts
) 


qs_save(
  FEZF2_DE_DESeq2_Unsubsetted, 
  "../Data/FANS/DE/FEZF2_DE_DESeq2_Unsubsetted.qs2"
)

qs_save(
  FEZF2_DE_DESeq2_SubsetCells, 
  "../Data/FANS/DE/FEZF2_DE_DESeq2_SubsetCells.qs2"
) 

qs_save(
  FEZF2_DE_DESeq2_SubsetCounts, 
  "../Data/FANS/DE/FEZF2_DE_DESeq2_SubsetCounts.qs2"
)



    ## 3.6 Oligodendrocytes DE Analysis ----------------------------------------

Oligodendrocytes.Unsubsetted <- subset(FANS, subset = ID_WNN_L25_Predicted == "Oligodendrocytes")
Oligodendrocytes.SubsetCells <- subset(FANS.SubsetCells, subset = ID_WNN_L25_Predicted == "Oligodendrocytes")
Oligodendrocytes.SubsetCounts <- subset(FANS.SubsetCounts, subset = ID_WNN_L25_Predicted == "Oligodendrocytes")

Oligodendrocytes_DE_DESeq2_Unsubsetted <- DESeq2_DE(
  Oligodendrocytes.Unsubsetted
)

Oligodendrocytes_DE_DESeq2_SubsetCells <- DESeq2_DE(
  Oligodendrocytes.SubsetCells
)

Oligodendrocytes_DE_DESeq2_SubsetCounts <- DESeq2_DE(
  Oligodendrocytes.SubsetCounts
) 


qs_save(
  Oligodendrocytes_DE_DESeq2_Unsubsetted, 
  "../Data/FANS/DE/Oligodendrocytes_DE_DESeq2_Unsubsetted.qs2"
)

qs_save(
  Oligodendrocytes_DE_DESeq2_SubsetCells, 
  "../Data/FANS/DE/Oligodendrocytes_DE_DESeq2_SubsetCells.qs2"
) 

qs_save(
  Oligodendrocytes_DE_DESeq2_SubsetCounts, 
  "../Data/FANS/DE/Oligodendrocytes_DE_DESeq2_SubsetCounts.qs2"
)



    ## 3.7 Inh_SST DE Analysis -------------------------------------------------

SST.Unsubsetted <- subset(FANS, subset = ID_WNN_L25_Predicted == "Inh_SST")
SST.SubsetCells <- subset(FANS.SubsetCells, subset = ID_WNN_L25_Predicted == "Inh_SST")
SST.SubsetCounts <- subset(FANS.SubsetCounts, subset = ID_WNN_L25_Predicted == "Inh_SST")

SST_DE_DESeq2_Unsubsetted <- DESeq2_DE(
  SST.Unsubsetted
)

SST_DE_DESeq2_SubsetCells <- DESeq2_DE(
  SST.SubsetCells
)

SST_DE_DESeq2_SubsetCounts <- DESeq2_DE(
  SST.SubsetCounts
) 


qs_save(
  SST_DE_DESeq2_Unsubsetted, 
  "../Data/FANS/DE/SST_DE_DESeq2_Unsubsetted.qs2"
)

qs_save(
  SST_DE_DESeq2_SubsetCells, 
  "../Data/FANS/DE/SST_DE_DESeq2_SubsetCells.qs2"
) 

qs_save(
  SST_DE_DESeq2_SubsetCounts, 
  "../Data/FANS/DE/SST_DE_DESeq2_SubsetCounts.qs2"
)



      ## 3.8 Exc_FEZF2_NTNG1 DE Analysis ---------------------------------------

FEZF2_NTNG1.Unsubsetted <- subset(FANS, subset = ID_WNN_L4_Predicted == "Exc_FEZF2_NTNG1")

FEZF2_NTNG1_DE_DESeq2_Unsubsetted <- DESeq2_DE(
  FEZF2_NTNG1.Unsubsetted
) 

qs_save(
  FEZF2_NTNG1_DE_DESeq2_Unsubsetted, 
  "../Data/FANS/DE/DESeq2/WNN_L4/FEZF2_NTNG1_DE_DESeq2_Unsubsetted.qs2"
)

summary(FEZF2_NTNG1_DE_DESeq2_Unsubsetted$res_L2FC_Shrink$padj)
table(FEZF2_NTNG1_DE_DESeq2_Unsubsetted$res_L2FC_Shrink$padj < 0.3)
rownames(FEZF2_NTNG1_DE_DESeq2_Unsubsetted$res_L2FC_Shrink %>% filter(padj < 0.2))
