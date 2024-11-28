


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




  ### 1.0 Load data ------------------------------------------------------------

M0_ATAC <- qread(
  paste0(
    "../Data/SeuratObjects/", 
    "M0_ATAC", 
    ".qrds"
  ), 
  nthr=nthr
)




  ### 2.0 Collect Genomic Ranges -----------------------------------------------



    ## 2.1 Save as GRanges -----------------------------------------------------

Peaks <- M0_ATAC@assays$ATAC@ranges
message(
  paste0(
    "All ATAC features collected: ", 
    length(Peaks) == dim(M0_ATAC)[1]
  )
)

Peaks_unique <- paste0(
  Peaks@seqnames, 
  "_", 
  start(Peaks), 
  "_", 
  end(Peaks)
)

message(
  paste0(
    "All open chromatin regions have an unique ID: ",
    !(any(duplicated(Peaks_unique)) | any(is.na(Peaks_unique)))
  )
)

Peaks$ID <- Peaks_unique
Peaks$ID2 <- str_replace_all(
  Peaks$ID, 
  pattern = "_", 
  replacement = "-"
)


all(
  paste0(
    Peaks@seqnames, 
    "-", 
    start(Peaks), 
    "-", 
    end(Peaks)
  ) == 
  Peaks$ID2
)

Peaks$Chromosome <- seqnames(Peaks)



    ## 2.2 Calculate GC Bias ---------------------------------------------------

M0_ATAC <- RegionStats(
  M0_ATAC, 
  genome = BSgenome.Hsapiens.UCSC.hg38, 
  assay = "ATAC"
)




  ### 3.0 Chromome-wise Genomic Ranges Characteristics -------------------------



    ## 3.1 Collect Chromosome Info ---------------------------------------------

BSgenome.Hsapiens.NCBI.GRCh38 
BSgenome.Hsapiens.UCSC.hg38 

table(seqlevels(Peaks) %in% seqlevels(BSgenome.Hsapiens.UCSC.hg38))  

HS.seqinfo.Len <- c() 
HS.seqinfo.GC <- c() 

lapply(
  seqlevels(Peaks), 
  FUN=function(chr){
    HS.seqinfo.Len[chr] <<- BSgenome.Hsapiens.UCSC.hg38[[chr]] |> length() 
    tmp <- HS.seqinfo.GC[chr] <- BSgenome.Hsapiens.UCSC.hg38[[chr]] |> 
      alphabetFrequency() 
    HS.seqinfo.GC[chr] <<- sum(tmp[c("G", "C")]) / sum(tmp[c("A", "G", "T", "C")]) 
  }
)


barplot(HS.seqinfo.Len)
barplot(HS.seqinfo.GC*100, ylim=c(0,100))
barplot(HS.seqinfo.GC*100, ylim=c(30,50))

Chromosomes <- data.frame(
  Chromosome = names(HS.seqinfo.GC), 
  Length = seqlengths(BSgenome.Hsapiens.UCSC.hg38)[match(names(HS.seqinfo.GC), names(seqlengths(BSgenome.Hsapiens.UCSC.hg38)))], 
  Length2=HS.seqinfo.Len,
  GC_Content = HS.seqinfo.GC
)
all(Chromosomes$Length==Chromosomes$Length2)


rm(HS.seqinfo.GC, HS.seqinfo.Len)
Chromosomes <- Chromosomes[,-which(colnames(Chromosomes)=="Length2")]
all(rownames(Chromosomes)==Chromosomes$Chromosome)  



    ## 3.2  GC Content ---------------------------------------------------------

all(Peaks$ID2 == rownames(M0_ATAC@assays$ATAC@meta.features))
mcols(Peaks) <- cbind(mcols(Peaks), M0_ATAC@assays$ATAC@meta.features)
Peaks

all(width(Peaks) == Peaks$sequence.length)

PeakGC.tmp <- mcols(Peaks) |> 
  as.data.frame() %>% 
    group_by(Chromosome) %>% 
      summarize(PeakGC = mean(GC.percent))

all(Chromosomes$Chromosome==PeakGC.tmp$Chromosome)
Chromosomes$PeakGC <- PeakGC.tmp$PeakGC



  ## 3.3 Genomic (chromosomal) distribution ------------------------------------

table(Peaks$Chromosome)
all(rownames(Chromosomes)==names(table(Peaks$Chromosome)))
Chromosomes$N_Peaks <- table(Peaks$Chromosome)[match(rownames(Chromosomes), names(table(Peaks$Chromosome)))]
Chromosomes$Chromosome <- factor(Chromosomes$Chromosome, levels=Chromosomes$Chromosome)
Chromosomes$N_Peaks <- as.numeric(Chromosomes$N_Peaks)




  ### 4.0 ATAC Peak counts statistics ------------------------------------------


M0_TFIDF <- RunTFIDF(M0_ATAC)

nCounts_raw <- rowSums(M0_TFIDF@assays$ATAC$counts)
nCounts_data <- rowSums(M0_TFIDF@assays$ATAC$data)
cor.test(nCounts_raw, nCounts_data)
plot(nCounts_raw, nCounts_data, pch=".")
abline(0,1,col="red")
Peaks$TFIDFcounts <- nCounts_data


cor.test(nCounts_raw, Peaks$count) 
cor.test(nCounts_data, Peaks$TFIDFcounts) 

plot(density(Peaks$TFIDFcounts))
cor.test(Peaks$sequence.length, Peaks$TFIDFcounts)
cor.test(Peaks$sequence.length, Peaks$TFIDFcounts, method = "spearman")
plot(Peaks$sequence.length, Peaks$TFIDFcounts, pch=".")

summary(Peaks$TFIDFcounts)




  ### 5.0 Export Data ----------------------------------------------------------

qsave(
  Chromosomes, 
  paste0(
    "../Data/Annotations/Genomic_Regions/", 
    "Chromosomes", 
    ".qrds"
  ), 
  nthr=nthr
)


qsave(
  Chromosomes, 
  paste0(
    "../Data/Annotations/Genomic_Regions/", 
    "Chromosomes", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  Peaks, 
  paste0(
    "../Data/Annotations/Genomic_Regions/", 
    "Peaks", 
    ".qrds"
  ), 
  nthr=nthr
)


