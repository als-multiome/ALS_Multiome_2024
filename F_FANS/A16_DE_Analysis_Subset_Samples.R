# A16_DE_Analysis_Subset_Samples.R 



  ### 0.0 Load libraries ------------------------------------------------------- 

library(qs2)
library(tidyverse)
library(Seurat)
library(scuttle)
library(effsize)
library(effectsize)

  
  
  
  ### 1.0 Load data ------------------------------------------------------------

FANS <- qs_read(
    "../Data/FANS/SeuratObjects/FANS_Integrated_LabelsTransferred.qs2"
) 

DefaultAssay(FANS) <- "RNA"
FANS@assays$SCT<- NULL 

FANS
table(FANS$orig.ident)
table(FANS$Sample_donor)
table(FANS$Sample_donor, FANS$TDP43)

Samples <- qs_read(
  "../Data/FANS/Samples/Samples.qs2"
)

Sample_Clusters <- qs_read(
  "../Data/FANS/Samples/Sample_Clusters.qs2"
)

Idents(FANS) <- FANS$TDP43




  ### 2.0 Subset TDP43-Low and TDP43-High samples for similar seq.depth --------

FANS@meta.data %>% 
  ggplot() + 
  aes(TDP43, nCount_RNA, fill=TDP43) + 
  geom_boxplot(outlier.shape=NA) + 
  scale_y_continuous(trans="log10") + 
  theme_classic() 

wilcox.test(nCount_RNA ~ TDP43, data = FANS@meta.data)
nCount_TDP43_Anova <- aov(nCount_RNA ~TDP43, data = FANS@meta.data) 
summary(nCount_TDP43_Anova)
eta_squared(nCount_TDP43_Anova, ci = 0.95)
cohen.d(nCount_RNA ~TDP43, data = FANS@meta.data) 

summary(FANS$nCount_RNA[FANS$TDP43=="High"])
summary(FANS$nCount_RNA[FANS$TDP43=="Low"])



    ## 2.1 Subset by removing cells --------------------------------------------

ind <- which(FANS$TDP43 =="High" & FANS$nCount_RNA > 12000)
set.seed(142)
ind <- sample(ind, size = 2500, replace = FALSE)

FANS$Flag <- "Yes"
FANS$Flag[ind] <- "No"

FANS.SubsetCells <- subset(FANS, subset = Flag == "Yes")


FANS.SubsetCells@meta.data %>% 
  ggplot() + 
  aes(TDP43, nCount_RNA, fill=TDP43) + 
  geom_boxplot(outlier.shape=NA) + 
  scale_y_continuous(trans="log10") + 
  theme_classic() 

wilcox.test(nCount_RNA ~ TDP43, data = FANS.SubsetCells@meta.data)
nCount_TDP43_Anova <- aov(nCount_RNA ~TDP43, data =  FANS.SubsetCells@meta.data) 
summary(nCount_TDP43_Anova)
eta_squared(nCount_TDP43_Anova, ci = 0.95)
cohen.d(nCount_RNA ~TDP43, data =  FANS.SubsetCells@meta.data) 


summary(FANS.SubsetCells$nCount_RNA[FANS.SubsetCells$TDP43=="High"])
summary(FANS.SubsetCells$nCount_RNA[FANS.SubsetCells$TDP43=="Low"]) 



    ## 2.2 Subset by removing counts -------------------------------------------

FANS.High <- subset(FANS, subset = TDP43 == "High")
FANS.Low <- subset(FANS, subset = TDP43 == "Low")

counts = as.matrix(x = GetAssayData(object = FANS.High, assay = "RNA", slot = "counts"))
mat.downsampled <- downsampleMatrix(
  x = counts, 
  prop = 0.5, 
  bycol = TRUE, 
  sink = NULL
)

FANS.High <- CreateAssayObject(
  counts = mat.downsampled
)
FANS.High <- CreateSeuratObject(
  FANS.High
)


all(rownames(FANS.High@meta.data) %in% rownames(FANS@meta.data))
tmp <- FANS@meta.data[match(rownames(FANS.High@meta.data), rownames(FANS@meta.data)),] 
tmp <- tmp %>% 
  select(-c(nCount_RNA, nFeature_RNA)) 
FANS.High@meta.data <- FANS.High@meta.data %>% 
  select(-orig.ident)
all(rownames(FANS.High@meta.data)==rownames(tmp))
FANS.High@meta.data <- cbind(
  FANS.High@meta.data, 
  tmp
)
setequal(
  colnames(FANS.High@meta.data), 
  colnames(FANS.Low@meta.data)
)

FANS.High@meta.data <- FANS.High@meta.data[
  ,match(
    colnames(FANS.Low@meta.data), 
    colnames(FANS.High@meta.data)
  )
]

all(
  colnames(FANS.High@meta.data) ==
  colnames(FANS.Low@meta.data)
)

FANS.SubsetCounts <- merge(
  FANS.High, 
  FANS.Low
)

FANS.SubsetCounts@meta.data %>% 
  ggplot() + 
  aes(TDP43, nCount_RNA, fill=TDP43) + 
  geom_boxplot(outlier.shape=NA) + 
  scale_y_continuous(trans="log10") + 
  theme_classic() 

wilcox.test(nCount_RNA ~ TDP43, data = FANS.SubsetCounts@meta.data)
nCount_TDP43_Anova <- aov(nCount_RNA ~TDP43, data =  FANS.SubsetCounts@meta.data) 
summary(nCount_TDP43_Anova)
eta_squared(nCount_TDP43_Anova, ci = 0.95)
cohen.d(nCount_RNA ~TDP43, data =  FANS.SubsetCounts@meta.data) 



summary(FANS.SubsetCounts$nCount_RNA[FANS.SubsetCounts$TDP43=="High"])
summary(FANS.SubsetCounts$nCount_RNA[FANS.SubsetCounts$TDP43=="Low"]) 

all(
  table(
    FANS$Sample_donor, 
    FANS$TDP43
  ) == 
    
    table(
      FANS.SubsetCounts$Sample_donor, 
      FANS.SubsetCounts$TDP43
    )
)

FANS.SubsetCounts <- JoinLayers(FANS.SubsetCounts)




  ### 3.0 Save data ------------------------------------------------------------

qs_save(
  FANS, 
  "../Data/FANS/SeuratObjects/FANS.Unsubsetted.qs2"
)

qs_save(
  FANS.SubsetCells, 
  "../Data/FANS/SeuratObjects/FANS.SubsetCells.qs2"
)

qs_save(
  FANS.SubsetCounts, 
  "../Data/FANS/SeuratObjects/FANS.SubsetCounts.qs2"
)

