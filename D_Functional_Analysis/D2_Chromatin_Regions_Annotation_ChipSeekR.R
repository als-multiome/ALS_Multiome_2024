

  ### 0.0 Load libraries -------------------------------------------------------

source("~/ALS_Brain_Multiome.Rcfg")

library(qs) 
library(Seurat)
library(Signac)
library(tidyverse)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86) 
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38) 
library(ggrepel)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(clusterProfiler)




  ### 1.0 Load data ------------------------------------------------------------

M0_ATAC <- qread(
  paste0(
    "../Data/SeuratObjects/", 
    "M0_ATAC", 
    ".qrds"
  ), 
  nthr=nthr
)


Chromosomes <- qread(
  paste0(
    "../Data/Annotations/Genomic_Regions/", 
    "Chromosomes", 
    ".qrds"
  ), 
  nthr=nthr
)

Peaks <- qread(
  paste0(
    "../Data/Annotations/Genomic_Regions/", 
    "Peaks", 
    ".qrds"
  ), 
  nthr=nthr
)





  ### 2.0 Annotate Peaks with ChIPseeker ---------------------------------------

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


peakAnno <- annotatePeak(
  peak = Peaks, 
  tssRegion=c(-3000, 3000),
  TxDb=txdb, 
  annoDb="org.Hs.eg.db"
)




  ### 3.0 Summarize annotations ------------------------------------------------

peakAnno_df <- as.data.frame(peakAnno)


peakAnno_df$Annotation_L2 <- peakAnno_df$annotation

peakAnno_df$Annotation_L2[grep("Exon", peakAnno_df$annotation)] <- "Exon"
peakAnno_df$Annotation_L2[grep("Intron", peakAnno_df$annotation)] <- "Intron"
peakAnno_df$Annotation_L2[grep("Promoter", peakAnno_df$annotation)] <- "Promoter"
table(peakAnno_df$Annotation_L2)


peakAnno_df$Annotation_L1 <- "NA"
peakAnno_df$Annotation_L1[peakAnno_df$Annotation_L2 %in% c(
    "3' UTR", 
    "5' UTR", 
    "Exon",  
    "Intron", 
    "Promoter" 
  )] <- "Genic"

peakAnno_df$Annotation_L1[peakAnno_df$Annotation_L2 %in% c(
  "Distal Intergenic", 
  "Downstream (<=300bp)"
)] <- "Intergenic" 


table(peakAnno_df$Annotation_L1)


peakAnno_df$Annotation_L2_Plain <- peakAnno_df$Annotation_L2

peakAnno_df$Annotation_L2_Plain[peakAnno_df$Annotation_L2=="3' UTR"] <- "Exon"
peakAnno_df$Annotation_L2_Plain[peakAnno_df$Annotation_L2=="5' UTR"] <- "Exon"
peakAnno_df$Annotation_L2_Plain[peakAnno_df$Annotation_L2=="Distal Intergenic"] <- "Intergenic"
peakAnno_df$Annotation_L2_Plain[peakAnno_df$Annotation_L2=="Downstream (<=300bp)"] <- "Downstream"




  ### 4.0 Save Data ------------------------------------------------------------

qsave(
  peakAnno, 
  paste0(
    "../Data/Annotations/Genomic_Regions/", 
    "peakAnno_csAnno", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  peakAnno_df, 
  paste0(
    "../Data/Annotations/Genomic_Regions/", 
    "peakAnno_df", 
    ".qrds"
  ), 
  nthr=nthr
)



