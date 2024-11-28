


  ### 0.0 Load libraries -------------------------------------------------------

source("~/ALS_Brain_Multiome.Rcfg")

library(qs) 
library(Seurat)
library(Signac)
library(tidyverse)
library(GenomicRanges)
library(GenomeInfoDb) 
library(EnsDb.Hsapiens.v86) 
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38) 




  ### 1.0 Load data ------------------------------------------------------------

M0 <- qread(
  paste0(
    "../Data/SeuratObjects/", 
    "M0", 
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




  ### 2.0 Annotate Peaks with PeakLinks ----------------------------------------



    ### 2.1 Calculate region stats ---------------------------------------------

M0 <- RegionStats(
  M0, 
  genome = BSgenome.Hsapiens.UCSC.hg38, 
  assay = "ATAC"
)

chromosomes <- as.character(
  unique(
    seqnames(
      M0@assays$ATAC@ranges
    )
  )
)



    ### 2.2 Calculate Peak Links -----------------------------------------------

for (chr in chromosomes){ 

  message(
    paste0(
      "Linking peaks on ", 
      str_to_title(chr), 
      "... "
    )
  )
  message(Sys.time())
  
  genes <- unique(M0@assays$ATAC@annotation$gene_name[which(seqnames(M0@assays$ATAC@annotation) == chr)])
  genes <- genes[which(genes %in% rownames(M0@assays$RNA))] 
  M0 <- LinkPeaks(
    object=M0, 
    peak.assay = "ATAC", 
    expression.assay = "RNA", 
    genes.use=genes
  )
  message(
    paste0(
      "Saving links..."
    )
  ) 
  message(Sys.time())
  qsave(
    M0@assays$ATAC@links, 
    paste0(
      "../Data/Annotations/Genomic_Regions/", 
      "M0_ATAC_Links_", 
      chr,
      ".qrds"
    ), 
    nthr=nthr
  )
}



    ## 2.3 Merge results -------------------------------------------------------

Links_Files <- list.files(
  paste0(
    "../Data/Annotations/Genomic_Regions/"
  )
)

Links_Files <- grep("Links", Links_Files, value = TRUE) |> 
  str_sort(numeric = TRUE) 


Links <- qread(
  paste0(
    "../Data/Annotations/Genomic_Regions/", 
    Links_Files[1]
  )
)

for (i in 2:length(Links_Files)){
  tmp <- qread(
    paste0(
      "../Data/Annotations/Genomic_Regions/", 
      Links_Files[i]
    )
  ) 
  stopifnot(
    all(
      names(mcols(Links)) == names(mcols(tmp))
    )
  ) 
  Links <- c(Links, tmp) 
  rm(tmp)
}



    ## 2.4 Save results --------------------------------------------------------

qsave(
  Links, 
  paste0(
    "../Data/Annotations/Genomic_Regions/", 
    "M0_ATAC_Peak_Links", 
    ".qrds"
  ), 
  nthr=nthr
)

M0@assays$ATAC@links <- Links

qsave(
  M0, 
  paste0(
    "../Data/SeuratObjects/", 
    "M0", 
    ".qrds"
  ), 
  nthr=nthr
) 


file.remove(
  paste0(
    "../Data/Annotations/Genomic_Regions/", 
    Links_Files
  )
)

rm(i, chromosomes, Links_Files, chr, genes)

