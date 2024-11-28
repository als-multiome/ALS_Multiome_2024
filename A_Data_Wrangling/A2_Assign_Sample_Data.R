  ### 0.0 Load libraries -------------------------------------------------------

source("~/ALS_Brain_Multiome.Rcfg") 

library(Seurat) 
library(Signac)
library(tidyverse)
library(qs)
library(data.table)





  ### 1.0 Load data ------------------------------------------------------------


M0 <- qread(
  paste0(
    "../Data/SeuratObjects/", 
    "M0", 
    ".qrds"
  ), 
  nthr = nthr
)


Sample_data <- fread(
  paste0(
    "../Data/Input/", 
    "Sample_data.txt"
  )
)

Sample_data[Sample_data==""] <- NA

M0@misc$Sample_data <- Sample_data 




  ### 2.0 Transfer IDs and Group variables from Sample data --------------------



    ## 2.1 Transfer Sample IDs -------------------------------------------------

all(unique(M0$WellCluster) %in% Sample_data$WellCluster)
setequal(unique(M0$WellCluster), Sample_data$WellCluster)
M0$ID <- Sample_data$ID[match(M0$WellCluster, Sample_data$WellCluster)]



    ## 2.2 Transfer sex labels -------------------------------------------------

M0$Sex <- Sample_data$Sex[match(M0$WellCluster, Sample_data$WellCluster)]
table(is.na(M0$sex))
table(is.na(M0$Sex))

table(M0$sex, M0$Sex)

table(M0$WellCluster[M0$Sex=="f" & M0$sex=="1"])
table(M0$WellCluster[M0$Sex=="m" & M0$sex=="0"])

M0@meta.data <- M0@meta.data[,-which(colnames(M0@meta.data)=="sex")]




    ## 2.3. Transfer Sample IDs ------------------------------------------------

M0$Sample <- Sample_data$sample[match(M0$WellCluster, Sample_data$WellCluster)]
table(is.na(M0$Sample))
# Incomplete 


    ## 2.4 Transfer SampleGroup labels -----------------------------------------

M0$Case <- Sample_data$SampleGroup[match(M0$WellCluster, Sample_data$WellCluster)]
table(is.na(M0$type))
table(is.na(M0$Case))

table(M0$type, M0$Case)

table(M0$WellCluster[M0$type=="ALS" & M0$Case == "ALSFTD"])
table(M0$Sample[M0$type=="ALS" & M0$Case == "ALSFTD"]) 

table(M0$WellCluster[M0$type=="ALS" & M0$Case=="Control"])
table(M0$Sample[M0$type=="ALS" & M0$Case=="Control"])

M0$Case <- str_replace_all(M0$Case, "ALSFTD", "ALS_FTD")
M0$Case <- str_replace_all(M0$Case, "Control", "HC")

table(M0$Case, M0$type)

M0@meta.data <- M0@meta.data[,-which(colnames(M0@meta.data)=="type")]



    ## 2.5 Transfer batch labels -----------------------------------------------

colnames(M0@meta.data)[colnames(M0@meta.data)=="Batch"] <- "batch"
M0$Batch <- Sample_data$Batch[match(M0$WellCluster, Sample_data$WellCluster)]
table(is.na(M0$batch))
table(is.na(M0$Batch))

table(M0$batch, M0$Batch)

table(M0$WellCluster[M0$Batch=="15"])

M0@meta.data <- M0@meta.data[,-which(colnames(M0@meta.data)=="batch")]



    ## 2.6 Transfer Cohort labels ----------------------------------------------

M0$Cohort <- Sample_data$Cohort[match(M0$WellCluster, Sample_data$WellCluster)]
table(is.na(M0$Cohort))



    ## 2.7 Transfer C9ORF72 labels ---------------------------------------------

M0$C9ORF72 <- Sample_data$C9ORF72[match(M0$WellCluster, Sample_data$WellCluster)]
table(is.na(M0$C9ORF72))

table(M0$Case, M0$C9ORF72)



    ## 2.8 Define ALS_FTD C9 Subgroups -----------------------------------------

M0$Case_Type <- M0$Case
M0$Case_Type[M0$C9ORF72==1] <- "C9_ALS_FTD"

table(is.na(M0$Case_Type))
table(M0$Case, M0$Case_Type)



    ## 2.9 Transfer Age labels -------------------------------------------------

M0$Age <- Sample_data$Age[match(M0$WellCluster, Sample_data$WellCluster)]
table(is.na(M0$Age))

summary(M0@misc$Sample_data$Age)

summary(M0@misc$Sample_data$Age[M0@misc$Sample_data$SampleGroup=="Control"])
summary(M0@misc$Sample_data$Age[M0@misc$Sample_data$SampleGroup=="ALS"]) 
summary(M0@misc$Sample_data$Age[(M0@misc$Sample_data$SampleGroup=="ALS_FTD" | M0@misc$Sample_data$SampleGroup=="ALSFTD") & M0@misc$Sample_data$C9ORF72==0])
summary(M0@misc$Sample_data$Age[(M0@misc$Sample_data$SampleGroup=="ALS_FTD" | M0@misc$Sample_data$SampleGroup=="ALSFTD") & M0@misc$Sample_data$C9ORF72==1])




  ### 3.0 Generate sample-wise cumulative cell metrics -------------------------



    ## 3.1 Add cumulative cell-type metrics ------------------------------------

M0$percent_mito <- M0@misc$GEX_QC$percent_mito
M0$atac_peak_region_fragments <- M0@misc$ATAC_QC$atac_peak_region_fragments 
M0$atac_fragments <- M0@misc$ATAC_QC$atac_fragments


M0@meta.data %>% 
  group_by(ID) %>% 
  summarise(
    nCells = length(gex_barcode), 
    nCount_RNA = sum(nCount_RNA), 
    mean_nCount_RNA = mean(nCount_RNA),  
    mean_nFeature_RNA = mean(nFeature_RNA), 
    nCount_ATAC = sum(nCount_ATAC), 
    mean_nCount_ATAC = mean(nCount_ATAC), 
    mean_nFeature_ATAC = sum(nFeature_ATAC), 
    pctMito_RNA = sum(percent_mito*nCount_RNA)/sum(nCount_RNA), 
    mean_pctMito_RNA = mean(percent_mito), 
    peakRegionFragments_ATAC = sum(atac_peak_region_fragments), 
    mean_peakRegionFragments_ATAC = mean(atac_peak_region_fragments), 
    FRIP_ATAC = sum(atac_peak_region_fragments)/sum(atac_fragments), 
    mean_FRIP_ATAC = mean(atac_peak_region_fragments/atac_fragments)
  ) %>% 
  data.frame() -> 
  Sample_metrics 



    ## 3.2 Add scaled sample-wise metrics --------------------------------------

Sample_metrics$pctMito_RNA_scaled <- scale(Sample_metrics$pctMito_RNA, center = TRUE)[,1]
Sample_metrics$peakRegionFragments_ATAC_scaled <- scale(Sample_metrics$peakRegionFragments_ATAC, center=TRUE)[,1] 

Sample_metrics$Case <- M0$Case[match(Sample_metrics$ID, M0$ID)]
Sample_metrics$Case_Type <- M0$Case_Type[match(Sample_metrics$ID, M0$ID)]




  ### 4.0 Add Sample_metrics to M0 object --------------------------------------

M0@misc$Sample_metrics <- Sample_metrics 




  ### 5.0 Export data ----------------------------------------------------------

qsave(
  Sample_metrics,
  paste0(
    "../Data/Input/", 
    "Sample_metrics", 
    ".qrds"
  ), 
  nthr = nthr
)

qsave(
  M0,
  paste0(
    "../Data/SeuratObjects/", 
    "M0", 
    ".qrds"
  ), 
  nthr = nthr
)

