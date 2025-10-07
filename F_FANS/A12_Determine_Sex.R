# A9_Determine_Sex.R 

  ### 0.0 Load libraries ------------------------------------------------------- 

library(qs2)
library(tidyverse)
library(Seurat)
library(ggrepel)
library(data.table)

  
  
      
  ### 1.0 Load data ------------------------------------------------------------
  
Samples <- qs_read(
  "../Data/FANS/Samples/Samples.qs2"
)

FANS <- qs_read(
  "../Data/FANS/SeuratObjects/FANS_Merged_QCFiltered.qs2"
)

FANS_RNA_Features <- qs_read(
  "../Data/FANS/Annotations/FANS_RNA_Features.qs2"
)




  ### 2.0 Aggregate expression by sample cluster -------------------------------

FANS$CellId <- rownames(FANS@meta.data)
FANS$CB <- str_split(FANS$CellId,"_", simplify = TRUE)[,ncol(str_split(FANS$CellId,"_", simplify = TRUE))]

FANS$ID_Cluster <- paste0(FANS$ID, "_", FANS$Souporcell_Cluster)
FANS_Pseudobulk_Counts <- AggregateExpression(FANS, group.by = "ID_Cluster")
FANS_Pseudobulk_Counts <- FANS_Pseudobulk_Counts$RNA
colnames(FANS_Pseudobulk_Counts) 
rownames(FANS_Pseudobulk_Counts)

colSums(FANS_Pseudobulk_Counts)
FANS_Pseudobulk_Counts <- data.frame(
  apply(
    FANS_Pseudobulk_Counts, 
    MARGIN = 2, 
    FUN = function(x){
      return(
        x*1e6/sum(x, na.rm = TRUE)
      )
    }
  )
)




  ## 2.0 Calculate expression of ChrX and ChrY genes per Sample Cluster --------

ChrX_Genes <- c("XIST","TSIX")
ChrY_Genes <- FANS_RNA_Features$ID_10X[which(FANS_RNA_Features$CHR=="Y")]

ChrX_Genes_Counts <- FANS_Pseudobulk_Counts[rownames(FANS_Pseudobulk_Counts) %in% ChrX_Genes,]
ChrY_Genes_Counts <- FANS_Pseudobulk_Counts[rownames(FANS_Pseudobulk_Counts) %in% ChrY_Genes,]

ChrX_Genes_Counts_Total <- colSums(ChrX_Genes_Counts)
ChrY_Genes_Counts_Total <- colSums(ChrY_Genes_Counts)
all(names(ChrX_Genes_Counts_Total) == names(ChrY_Genes_Counts_Total))




  ### 3.0 Compile Sample Cluster data ------------------------------------------

df = data.frame(
  ID_Cluster= names(ChrX_Genes_Counts_Total), 
  Sample = str_split(
    names(ChrX_Genes_Counts_Total), 
    "\\.", 
    simplify=TRUE
  )[,1],
  Souporcell_Cluster = str_split(
    names(ChrX_Genes_Counts_Total), 
    "\\.", 
    simplify=TRUE
  )[,2],
  Sample_Cluster = str_replace_all(
    names(ChrX_Genes_Counts_Total), 
    "\\.", 
    "_"
  ), 
  ChrX_Genes_Counts_Total, 
  ChrY_Genes_Counts_Total
)

df$Sex = "NA" 
df$Sex[df$ChrX_Genes_Counts_Total>1000 & df$ChrY_Genes_Counts_Total<1000] <- "f"
df$Sex[df$ChrX_Genes_Counts_Total<1000 & df$ChrY_Genes_Counts_Total>1000] <- "m"
table(df$Sex, useNA = "always")


df %>% 
  ggplot() + 
  aes(ChrX_Genes_Counts_Total, ChrY_Genes_Counts_Total, fill=Sex, label=ID_Cluster) + 
  geom_point(pch=21, size=3) + 
  geom_text_repel() + 
  theme_classic()

FANS$Sex <- df$Sex[match(FANS$ID_Cluster, df$Sample_Cluster)]
table(FANS$Sex, FANS$TDP43)




  ### 3.0 Save data ------------------------------------------------------------

qs_save(
  FANS, 
  "../Data/FANS/SeuratObjects/FANS_Merged_QCFiltered_SexAdded.qs2"
)

qs_save(
  Samples, 
  "../Data/FANS/Samples/Samples.qs2"
)

qs_save(
  df, 
  "../Data/FANS/Samples/Sample_Clusters.qs2"
)


