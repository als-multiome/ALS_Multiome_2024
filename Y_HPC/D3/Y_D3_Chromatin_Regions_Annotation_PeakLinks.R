##### HEAD #####################################################################
###                                                                          ###
### *D3_Chromatin_Regions_Annotation_PeakLinks.R                             ###
### *v1.0                                                                    ###
### *VG17.06.24                                                              ###
### *Project: ALS_Brain_Multiome                                             ###
###                                                                          ###
### *Description:  Annotation of ATAC open chromatin regions for             ###
###   functional analysis                                                    ###  
### *Requires: M0 Seurat object                                              ###
###                                                                          ###    
### *Returns: Serialized chromatin region annotations                        ###
###                                                                          ###  
##### ENDHEAD ##################################################################




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



## 3.3 Merge results ---- 


Links_Files <- list.files(
  paste0(
    project_folder, 
    "Source/Data/GenomicRanges/"
  )
)

Links_Files <- grep("Links", Links_files, value = TRUE) |> 
  str_sort(numeric = TRUE) 


Links <- qread(
  paste0(
    project_folder, 
    "Source/Data/GenomicRanges/", 
    Links_Files[1]
  )
)

for (i in 2:length(Links_Files)){
  tmp <- qread(
    paste0(
      project_folder, 
      "Source/Data/GenomicRanges/", 
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



## 3.4 Save results ---- 

qsave(
  Links, 
  paste0(
    project_folder, 
    "Source/Data/GenomicRanges/", 
    "M0_ATAC_Links", 
    ".qrds"
  ), 
  nthr=nthr
)

M0@assays$ATAC@links <- Links
file.remove(
  paste0(
    project_folder, 
    "Source/Data/GenomicRanges/", 
    Links_Files
  )
)

rm(i, chromosomes, Links_Files)




### 4.0 EDA Peak Links ---- 

Links 
summary(Links$score)
table(Links$score>0)

all(Links$peak %in% Peaks$ID2)
table(Links$peak %in% Peaks$ID2)

ggplot(
  data.frame(
    Group=c("No link found", "Link(s) found"), 
    Pct = as.numeric(round(table(Peaks$ID2 %in% M0@assays$ATAC@links$peak) |> prop.table()*100, 1))
  )
) + 
  aes(x="ATAC Peaks", y=Pct, fill=Group) + 
  geom_bar(stat="identity", col="black") +  
  geom_text(
    aes(label=Group),
    position=position_stack(vjust=0.5),
    size=4,
    fontface="bold.italic",
    col="white",
    hjust=0.5
  ) + 
  scale_fill_manual(values=c("#00BB00", "#DD0000")) + 
  scale_y_continuous(expand=expand_scale(mult=c(0,0.05))) + 
  xlab("\nATAC Peaks") + 
  ylab("Percent %") + 
  theme_classic() + 
  theme(
    axis.title = element_text(size=12, color="black", face="bold.italic"),
    axis.text.x = element_blank(), 
    axis.text.y = element_text(color="black", face="bold.italic"),
    axis.line = element_line(size=1, color="black"), 
    legend.position = "Null"
  )

ggplot(
  data.frame(
    NPeaks=factor(
      names(sort(table(as.numeric(sort(table(M0@assays$ATAC@links$peak), decreasing = TRUE))), decreasing=TRUE)), 
      levels = names(sort(table(as.numeric(sort(table(M0@assays$ATAC@links$peak), decreasing = TRUE))), decreasing=TRUE))
    ), 
    N=as.numeric(prop.table(sort(table(as.numeric(sort(table(M0@assays$ATAC@links$peak), decreasing = TRUE))), decreasing=TRUE)))*100
  )
) + 
  aes(NPeaks, N, fill=NPeaks) + 
  geom_bar(stat="identity") + 
  scale_y_continuous(expand=expand_scale(mult=c(0,0.05))) + 
  xlab("\nNumber of links") + 
  ylab("Percent of peaks with links [%]\n") + 
  theme_classic() + 
  theme(
    axis.title = element_text(size=12, color="black", face="bold.italic"),
    axis.text = element_text(color="black", face="bold.italic"), 
    legend.position = "Null"
  )


tmp <- data.frame(round(prop.table(sort(table(as.numeric(sort(table(M0@assays$ATAC@links$peak), decreasing = TRUE))), decreasing=TRUE))*100, 1))
colnames(tmp) <- c("Number of links", "Percent of peaks")
tmp$`Number of peaks` <- sort(table(as.numeric(sort(table(M0@assays$ATAC@links$peak), decreasing = TRUE))), decreasing=TRUE)*100 
tmp <- tmp[,c(3,2,1)]
View(tmp)
rm(tmp)


ggplot(
  data.frame(
    Group=c("No link found", "Link(s) found"), 
    Pct = as.numeric(round(table(rownames(M0@assays$RNA) %in% M0@assays$ATAC@links$gene) |> prop.table()*100, 1))
  )
) + 
  aes(x="Genes", y=Pct, fill=Group) + 
  geom_bar(stat="identity", col="black") +  
  geom_text(
    aes(label=Group),
    position=position_stack(vjust=0.5),
    size=4,
    fontface="bold.italic",
    col="white",
    hjust=0.5
  ) + 
  scale_fill_manual(values=c("#00BB00", "#DD0000")) + 
  scale_y_continuous(expand=expand_scale(mult=c(0,0.05))) + 
  xlab("\nGenes") + 
  ylab("Percent %") + 
  theme_classic() + 
  theme(
    axis.title = element_text(size=12, color="black", face="bold.italic"),
    axis.text.x = element_blank(), 
    axis.text.y = element_text(color="black", face="bold.italic"),
    axis.line = element_line(size=1, color="black"), 
    legend.position = "Null"
  )



ggplot(
  data.frame(
    NPeaks=factor(
      names(sort(table(as.numeric(sort(table(M0@assays$ATAC@links$gene), decreasing = TRUE))), decreasing=TRUE)), 
      levels = names(sort(table(as.numeric(sort(table(M0@assays$ATAC@links$gene), decreasing = TRUE))), decreasing=TRUE))
    ), 
    N=as.numeric(prop.table(sort(table(as.numeric(sort(table(M0@assays$ATAC@links$gene), decreasing = TRUE))), decreasing=TRUE)))*100
  )
) + 
  aes(as.numeric(NPeaks), N, fill=NPeaks) + 
  geom_bar(stat="identity") + 
  scale_x_continuous(limits=c(0,70), expand=expand_scale(mult=c(0,0.05))) + 
  scale_y_continuous(expand=expand_scale(mult=c(0,0.05))) + 
  xlab("\nNumber of links") + 
  ylab("Percent of genes with links [%]\n") + 
  theme_classic() + 
  theme(
    axis.title = element_text(size=12, color="black", face="bold.italic"),
    axis.text = element_text(color="black", face="bold.italic"), 
    legend.position = "Null"
  )

tmp <- data.frame(round(prop.table(sort(table(as.numeric(sort(table(M0@assays$ATAC@links$gene), decreasing = TRUE))), decreasing=TRUE))*100, 1))
colnames(tmp) <- c("Number of links", "Percent of genes")
tmp$`Number of genes` <- sort(table(as.numeric(sort(table(M0@assays$ATAC@links$gene), decreasing = TRUE))), decreasing=TRUE)
tmp <- tmp[,c(3,2,1)]
View(tmp)
rm(tmp)



### 5.0 Save data ---- 

qsave(
  M0, 
  paste0(
    project_folder, 
    "Source/Data/SeuratObjects/", 
    "M0", 
    ".qrds"
  ), 
  nthr=nthr
)


