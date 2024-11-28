
  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 


library(qs) 
library(circlize)
library(ComplexHeatmap)
library(tidyverse)
library(dendextend)
library(RColorBrewer)
library(GenomicRanges)
library(eulerr)
library(clusterProfiler)
library(org.Hs.eg.db)



qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)




  ### 1.0 Load Data ------------------------------------------------------------ 



    ## 1.1 Load Data -----------------------------------------------------------

Universe_RNA_AllCase_ALS_List <- qread(
  paste0(
    "../Data/Annotations/Genes/", 
    "Universe_AllCase_ALS_List", 
    ".qrds"
  ), 
  nthr=nthr
)

Universe_RNA_AllCase_ALSFTD_List <- qread(
  paste0(
    "../Data/Annotations/Genes/", 
    "Universe_AllCase_ALSFTD_List", 
    ".qrds"
  ), 
  nthr=nthr
)


RNA_DESeq_Results_12SVs_ALS <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/12_SVs/", 
    "DESeq_Results_12SVs_ALS", 
    ".qrds"
  ), 
  nthr=nthr
)

RNA_DESeq_Results_12SVs_ALSFTD <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/12_SVs/", 
    "DESeq_Results_12SVs_ALSFTD", 
    ".qrds"
  ), 
  nthr=nthr
)


RNA_DESeq_Results_12SVs_Index <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/12_SVs/", 
    "DESeq_Results_12SVs_Index", 
    ".qrds"
  ), 
  nthr=nthr
)


Mtx_ALS_ALSFTD_WNN_L25_Signature_AllCases <- qread(
  paste0(
    "../Data/Annotations/Signatures/RNA/", 
    "Mtx_ALS_ALSFTD_WNN_L25_Signature_AllCases", 
    ".qrds"
  ), 
  nthr=nthr
)



Universe_ATAC_AllCase_ALS_List <- qread(
  paste0(
    "../Data/Annotations/Genomic_Regions/", 
    "Universe_AllCase_ALS_List", 
    ".qrds"
  ), 
  nthr=nthr
)

Universe_ATAC_AllCase_ALSFTD_List <- qread(
  paste0(
    "../Data/Annotations/Genomic_Regions/", 
    "Universe_AllCase_ALSFTD_List", 
    ".qrds"
  ), 
  nthr=nthr
)


ATAC_DESeq_Results_8SVs_Index <- qread(
  paste0(
    "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/SVA/1_15_SVs/", 
    "DESeq_Results_8SVs_Index", 
    ".qrds"
  ), 
  nthr=nthr
)

ATAC_DESeq_Results_8SVs_ALS <- qread(
  paste0(
    "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/SVA/1_15_SVs/", 
    "DESeq_Results_8SVs_ALS", 
    ".qrds"
  ), 
  nthr=nthr
)

ATAC_DESeq_Results_8SVs_ALSFTD <- qread(
  paste0(
    "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/SVA/1_15_SVs/", 
    "DESeq_Results_8SVs_ALSFTD", 
    ".qrds"
  ), 
  nthr=nthr
)

ATAC_Peak_Links <- qread(
  paste0(
    "../Data/Annotations/Genomic_Regions/", 
    "M0_ATAC_Peak_Links", 
    ".qrds"
  ), 
  nthr=nthr
)

FANS_DE_TDP43_MAST <- qread(
  paste0(
    "../Data/DE/FANS/", 
    "FANS_DE_TDP43_MAST", 
    ".qrds"
  ), 
  nthr=nthr
)


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


RNA_L2FC_Shrink_Results_12SVs_ALS <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/12_SVs/", 
    "L2FC_Shrink_Results_12SVs_ALS", 
    ".qrds"
  ), 
  nthr=nthr
)


ATAC_L2FC_Shrink_Results_8SVs_ALS <- qread(
  paste0(
    "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/SVA/1_15_SVs/", 
    "L2FC_Shrink_Results_8SVs_ALS", 
    ".qrds"
  ), 
  nthr=nthr
)



ATAC_L2FC_Shrink_Results_8SVs_ALSFTD <- qread(
  paste0(
    "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/SVA/1_15_SVs/", 
    "L2FC_Shrink_Results_8SVs_ALSFTD", 
    ".qrds"
  ), 
  nthr=nthr
)


GEX_Features <- qread(
  paste0(
    "../Data/Annotations/Genes/", 
    "GEX_Features", 
    ".qrds"
  ), 
  nthr=nthr
)

GEX_Features_ENSEMBL_ENTREZ <- qread(
  paste0(
    "../Data/Annotations/Genes/", 
    "GEX_Features_ENSEMBL_ENTREZ", 
    ".qrds"
  ), 
  nthr=nthr
)



ColDict_WNN_L25 <- setNames(
  object = readxl::read_xlsx(paste0(
    "../Data/Visualization/", 
    "ALS_Brain_Multiome_ColDicts", 
    ".xlsx"
  ), 
  sheet = "WNN_L25"
  )$Color, 
  nm = readxl::read_xlsx(paste0(
    "../Data/Visualization/", 
    "ALS_Brain_Multiome_ColDicts", 
    ".xlsx"
  ), 
  sheet = "WNN_L25"
  )$WNN_L25
)


ColDict_Case <- setNames(
  object = readxl::read_xlsx(paste0(
    "../Data/Visualization/", 
    "ALS_Brain_Multiome_ColDicts", 
    ".xlsx"
  ), 
  sheet = "Case"
  )$Color, 
  nm = readxl::read_xlsx(paste0(
    "../Data/Visualization/", 
    "ALS_Brain_Multiome_ColDicts", 
    ".xlsx"
  ), 
  sheet = "Case"
  )$Case
)



    ## 1.2 Define AUX functions ------------------------------------------------

get_DEGs <- function(x, alpha=0.05){
  x <- as.data.frame(x) 
  x <- x[!is.na(x$padj),]
  x <- x[x$padj < alpha,] 
  return(rownames(x))
  
} 

get_DARs <- function(x, alpha=0.05){
  x <- as.data.frame(x) 
  x <- x[!is.na(x$padj),]
  x <- x[x$padj < alpha,] 
  return(rownames(x))
  
} 

get_DARs_Genes <- function(x, Links=ATAC_Peak_Links, alpha=0.05){
  
  x <- as.data.frame(x) 
  x <- x[!is.na(x$padj),]
  x <- x[x$padj < alpha,] 
  
  x <- rownames(x) 
  Genes <- Links$gene[Links$peak %in% x]
  return(
    unique(
      Genes
    )
    )
}


plot_GSEA_as_Tree <- function(GSEA_Results, save=TRUE, name, width, height){
  
  geneLists <- lapply(GSEA_Results$core_enrichment, function(x) unlist(strsplit(as.character(x), "/")))
  names(geneLists) <- GSEA_Results$Description
  geneLists <- geneLists[!is.na(geneLists)]
  
  n <- length(geneLists)
  w <- matrix(NA, nrow = n, ncol = n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      u <- unlist(geneLists[i])
      v <- unlist(geneLists[j])
      w[i, j] <- length(intersect(u, v)) / length(unique(c(u, v)))
    }
  }
  
  
  Terms <- names(geneLists)
  
  rownames(w) <- Terms
  colnames(w) <- Terms
  
  
  hclust <- as.dist(1 - w) %>% hclust(method = "average")
  dend <- as.dendrogram(hclust)
  
  # Get the vertex IDs of the leaves (bottom level nodes)
  leaf_ids <- which(V(graph_dend)$leaf)
  
  # Assign cluster labels to the leaf nodes
  V(graph_dend)$cluster[leaf_ids] <- clusters  # Add the cluster variable to the leaf nodes
  
  library(tidygraph)
  tbl_graph(dend)
  a <- ggraph(dend, 'dendrogram', height = height) + 
    geom_edge_elbow() + 
    geom_node_point(aes(filter = leaf, color = label))  # Color only the leaves
  
  a$data$NES <- GSEA_Kegg_ByL2FC@result$NES[match(a$data$label, GSEA_Kegg_ByL2FC@result$Description)] 
  
  count_GSEA_CoreEnrichment_Genes <- function(x){
    if(is.na(x)) {return(NA)} else{
      return(sum(str_split(x, "/", simplify=TRUE)[1,] != ""))
    }
  }
  a$data$nGenes <- sapply(
    GSEA_Kegg_ByL2FC@result$core_enrichment[match(a$data$label, GSEA_Kegg_ByL2FC@result$Description)], 
    FUN = count_GSEA_CoreEnrichment_Genes
  )
  colFun=colorRamp2(
    c(
      -max(abs(a$data$NES), na.rm=T), 
      0, 
      max(abs(a$data$NES), na.rm=T)
    ), 
    c(
      "blue", 
      "#EEEEEE", 
      "red"
    )
  )
  
  library(ggrepel)
  ggraph(a$data, 'dendrogram', height = height) + 
    geom_edge_elbow(strength = 1, flipped = FALSE, width=1) + 
    geom_node_point(aes(filter = leaf, fill = NES, size = nGenes), pch=21, col="#000000") + 
    geom_node_text(aes(filter = leaf, label = label), nudge_y = 0.1, hjust = 0) + 
    scale_size(
      limits=c(
        min(a$data$nGenes, na.rm=TRUE)-0.25*(max(a$data$nGenes, na.rm=TRUE)-min(a$data$nGenes, na.rm=TRUE)), 
        max(a$data$nGenes, na.rm=TRUE)
      )
    ) + 
    scale_fill_gradientn(
      colours = c(
        "#0000FF", 
        "#05f7ff", 
        "white", 
        "#f5dc00",
        "#FF0000"
      ), 
      limits=c(
        -max(abs(a$data$NES), na.rm=TRUE), 
        max(abs(a$data$NES), na.rm=TRUE)
      )
      
    ) + coord_flip(xlim = c(0, 8), ylim=c(1, -3)) + scale_y_reverse() + 
    theme(
      panel.background = element_blank(), 
      
    ) 
  
  
  
  GSEA_Results_Direction <- GSEA_Results$NES<0
  
  mycolors <- sort(rainbow(20))[c(1, 20, 10, 11, 2, 19, 3, 12, 4, 13, 5, 14, 6, 15, 7, 16, 8, 17, 9, 18)]
  
  leafType <- as.factor(gsub(" .*", "", GSEA_Results_Direction[ix]))
  
  if (max(nchar(GSEA_Results[ix])) >= 1) { # if "Up regulated or Downregulated"; not "A", "B"
    # leafColors = c("green","red")  else  # mycolors # k-Means
    leafColors <- mycolors[1:2]
  } else { # convert c("B","D","E") to c(2, 4, 5)
    # leafType= as.factor( gsub(" .*","", enrichedTerms$Direction[ix] ) )
    leafType <- match(gsub(" .*", "", GSEA_Results$Direction[ix]), toupper(letters))
    
    
    
  }
  leafColors <- c("#FF0000", "#0000FF")
  
  leafSize <- as.numeric(abs(GSEA_Results$NES[ix])) # leaf size represent NES values
  leafSize <- .9 * (leafSize - min(leafSize)) / (max(leafSize) - min(leafSize) + 1e-50) + .1 # scale more aggressively
  
  
  print(
    dend %>%
      as.dendrogram(hang = -1) %>%
      set("leaves_pch", 19) %>% # type of marker
      set("leaves_cex", leafSize*3) %>% # Size
      set("leaves_col", leafColors[leafType]) %>% # up or down genes
      plot(horiz = TRUE, axes = FALSE)
  ) 
  
  if(save){
    pdf(
      name, 
      width = width, 
      height = height
    ) 
    dend %>%
      as.dendrogram(hang = -1) %>%
      set("leaves_pch", 19) %>% # type of marker
      set("leaves_cex", leafSize*3) %>% # Size
      set("leaves_col", leafColors[leafType]) %>% # up or down genes #leafColors
      plot(horiz = TRUE, axes = FALSE)
    dev.off()
  }
  
}

Return_GSEA_as_Tree <- function(GSEA_Results, name, save=TRUE, width, height){
  geneLists <- lapply(GSEA_Results$core_enrichment, function(x) unlist(strsplit(as.character(x), "/")))
  names(geneLists) <- GSEA_Results$Description
  geneLists <- geneLists[!is.na(geneLists)]
  
  n <- length(geneLists)
  w <- matrix(NA, nrow = n, ncol = n)
  
  for (i in 1:n) {
    for (j in i:n) {
      u <- unlist(geneLists[i])
      v <- unlist(geneLists[j])
      w[i, j] <- length(intersect(u, v)) / length(unique(c(u, v)))
    }
  }
  
  for (i in 1:n) {
    for (j in 1:(i - 1)) {
      w[i, j] <- w[j, i]
    }
  }
  
  if (0) {
    total_elements <- 30000
    n <- length(geneLists)
    w <- matrix(rep(0, n * n), nrow = n, ncol = n)
    
    for (i in 1:n) {
      for (j in (i + 1):n) {
        u <- unlist(geneLists[i])
        v <- unlist(geneLists[j])
        xx <- length(intersect(u, v))
        if (xx == 0) {
          next
        }
        mm <- length(u)
        nn <- total_elements - mm
        kk <- length(v)
        w[i, j] <- -sqrt(-phyper(xx - 1, mm, nn, kk, lower.tail = FALSE, log.p = TRUE))
      }
    }
    
    
    for (i in 1:n) {
      for (j in 1:(i - 1)) {
        w[i, j] <- w[j, i]
      }
    }
  }
  
  Terms <- paste(
    names(geneLists)
  )
  rownames(w) <- Terms
  colnames(w) <- Terms
  rightMargin=33
  par(mar = c(0, 0, 1, rightMargin)) 
  
  dend <- as.dist(1 - w) %>%
    hclust(method = "average")
  ix <- dend$order # permutated order of leaves
  
  GSEA_Results_Direction <- GSEA_Results$NES<0
  
  mycolors <- sort(rainbow(20))[c(1, 20, 10, 11, 2, 19, 3, 12, 4, 13, 5, 14, 6, 15, 7, 16, 8, 17, 9, 18)]
  
  leafType <- as.factor(gsub(" .*", "", GSEA_Results_Direction[ix]))
  
  if (max(nchar(GSEA_Results[ix])) >= 1) { # if "Up regulated or Downregulated"; not "A", "B"
    # leafColors = c("green","red")  else  # mycolors # k-Means
    leafColors <- mycolors[1:2]
  } else { # convert c("B","D","E") to c(2, 4, 5)
    # leafType= as.factor( gsub(" .*","", enrichedTerms$Direction[ix] ) )
    leafType <- match(gsub(" .*", "", GSEA_Results$Direction[ix]), toupper(letters))
    
    
    
  }
  leafColors <- c("#FF0000", "#0000FF")
  
  leafSize <- as.numeric(abs(GSEA_Results$NES[ix])) # leaf size represent NES values
  leafSize <- .9 * (leafSize - min(leafSize)) / (max(leafSize) - min(leafSize) + 1e-50) + .1 # scale more aggressively
  
  
  dend2 <- dend %>%
      as.dendrogram(hang = -1) %>%
      set("leaves_pch", 19) %>% # type of marker
      set("leaves_cex", leafSize*3) %>% # Size
      set("leaves_col", leafColors[leafType]) %>% # up or down genes
      
  return(dend2)
}




  ### 2.0 Define Sample clustering  ----------------------------------------------

mtx <- Mtx_ALS_ALSFTD_WNN_L25_Signature_AllCases


Case=rep("ALS", ncol(mtx))
Case[grep("ALSFTD", colnames(mtx))] <- "ALS_FTD"
CellTypes <- str_replace_all(str_replace_all(colnames(mtx), "ALS_", ""), "ALSFTD_", "")

cols <- unique(Case)
cols <- setNames(rep(NA, length(cols)), nm = cols)
cols <- ColDict_Case[match(names(cols), names(ColDict_Case))]

cols2 <- unique(CellTypes)
cols2 <- setNames(rep(NA, length(cols2)), nm = cols2)
cols2 <- ColDict_WNN_L25[match(names(cols2), names(ColDict_WNN_L25))]

cols3 <- c(cols, cols2)


row_ha  = rowAnnotation(
  Case=Case, 
  CellType=CellTypes, 
  col=list(Case=cols, CellType=cols2), 
  gp = gpar(col="black")  
)


heatmap <- Heatmap(
  matrix=t(mtx), 
  name = "Mtx", 
  clustering_distance_row = "pearson", 
  clustering_method_rows="complete", 
  clustering_distance_columns  = "manhattan", 
  clustering_method_columns = "average", 
  show_column_dend = FALSE, 
  show_column_names = FALSE, 
  left_annotation = row_ha, 
  col=colorRamp2(
    c(
      min(mtx, na.rm=T), 
      0, 
      max(mtx, na.rm=T)
    ), 
    c(
      "blue", 
      "#EEEEEE", 
      "red"
    )
  ), 
  border_gp = gpar(col = "black", lty = 1), 
  use_raster = FALSE
)

ht <- draw(heatmap) 


Sample_dend <- row_dend(ht)
Sample_hclust <- as.hclust(Sample_dend)

Genes_dend <- column_dend(ht)
Genes_hclust <- as.hclust(Genes_dend)

Sample_Ord <- setNames(Sample_hclust$order, nm = Sample_hclust$labels)
Genes_Ord <- setNames(Genes_hclust$order, nm = Genes_hclust$labels)

Sample_labels = Sample_hclust$labels[Sample_hclust$order] 
all(Sample_labels == names(Sample_Ord)[Sample_Ord]) 


rm(mtx, Case, CellTypes, cols, cols2, cols3, row_ha, heatmap, ht)




  ### 3.0 Plot Pct of transcriptomic changes -----------------------------------



    ## 3.1 RNA DEGs ALS --------------------------------------------------------

all(
  unique(
    RNA_DESeq_Results_12SVs_Index$CellType[
      RNA_DESeq_Results_12SVs_Index$CellTypeLevel=="WNN_L25"
    ]
  ) %in% 
    names(Universe_RNA_AllCase_ALS_List)
)

nDEGs_ALS <- setNames(
  RNA_DESeq_Results_12SVs_Index$All_q_0.05_ALS[
    RNA_DESeq_Results_12SVs_Index$CellTypeLevel == "WNN_L25" & 
      RNA_DESeq_Results_12SVs_Index$Comparison == "All_Cases"
  ], 
  nm = RNA_DESeq_Results_12SVs_Index$CellType[
    RNA_DESeq_Results_12SVs_Index$CellTypeLevel == "WNN_L25" & 
      RNA_DESeq_Results_12SVs_Index$Comparison == "All_Cases"
  ]
)

nGenes_Universe_ALS <- sapply(
    Universe_RNA_AllCase_ALS_List[2:length(Universe_RNA_AllCase_ALS_List)], 
    FUN = function(x){
      return(
        nrow(x)
      )
    }, USE.NAMES = TRUE
)


all(
  names(nDEGs_ALS) == names(nGenes_Universe_ALS)
)  

nGenes_Universe_ALS[sapply(nGenes_Universe_ALS, is.null)] <- 1

PctDEGs_ALS <- (nDEGs_ALS/unlist(nGenes_Universe_ALS))*100

rm(nGenes_Universe_ALS, nDEGs_ALS)



    ## 3.2 RNA DEGs ALSFTD -----------------------------------------------------

all(
  unique(
    RNA_DESeq_Results_12SVs_Index$CellType[
      RNA_DESeq_Results_12SVs_Index$CellTypeLevel=="WNN_L25"
    ]
  ) %in% 
    names(Universe_RNA_AllCase_ALSFTD_List)
)

nDEGs_ALSFTD <- setNames(
  RNA_DESeq_Results_12SVs_Index$All_q_0.05_ALSFTD[
    RNA_DESeq_Results_12SVs_Index$CellTypeLevel == "WNN_L25" & 
      RNA_DESeq_Results_12SVs_Index$Comparison == "All_Cases"
  ], 
  nm = RNA_DESeq_Results_12SVs_Index$CellType[
    RNA_DESeq_Results_12SVs_Index$CellTypeLevel == "WNN_L25" & 
      RNA_DESeq_Results_12SVs_Index$Comparison == "All_Cases"
  ]
)

nGenes_Universe_ALSFTD <- sapply(
  Universe_RNA_AllCase_ALSFTD_List[2:length(Universe_RNA_AllCase_ALSFTD_List)], 
  FUN = function(x){
    return(
      nrow(x)
    )
  }, USE.NAMES = TRUE
)


all(
  names(nDEGs_ALSFTD) == names(nGenes_Universe_ALSFTD)
)  

nGenes_Universe_ALSFTD[sapply(nGenes_Universe_ALSFTD, is.null)] <- 1

PctDEGs_ALSFTD <- (nDEGs_ALSFTD/unlist(nGenes_Universe_ALSFTD))*100

rm(nGenes_Universe_ALSFTD, nDEGs_ALSFTD)



    ## 3.3 Plot PctDEGs --------------------------------------------------------

all(names(PctDEGs_ALS) == names(PctDEGs_ALSFTD))
plot(PctDEGs_ALS, PctDEGs_ALSFTD) 

Pct_DEGs <- c(
  setNames(PctDEGs_ALSFTD, nm = paste0("ALSFTD_", names(PctDEGs_ALSFTD))),  
  setNames(PctDEGs_ALS, nm = paste0("ALS_", names(PctDEGs_ALS)))
)


mtx_Pct_DEGs <- matrix(
  Pct_DEGs, 
  nrow = length(Pct_DEGs), 
  ncol = 1, 
  byrow = FALSE
)
rownames(mtx_Pct_DEGs) <- names(Pct_DEGs) 
  

col_fun_Pct_DEGs = colorRamp2(c(0, quantile(mtx_Pct_DEGs)["25%"], median(mtx_Pct_DEGs), quantile(mtx_Pct_DEGs)["75%"], max(mtx_Pct_DEGs)), c("#EEEEEE", "#ffeda3", "#ffdc91", "#ff8000", "red")) 

circos.clear()
circos.heatmap(mtx_Pct_DEGs, col = col_fun_Pct_DEGs, cell.border = "#000000", cluster=FALSE)
circos.clear()

rm(PctDEGs_ALS, PctDEGs_ALSFTD)




  ### 4.0 Plot Pct of Chromatin changes ----------------------------------------



    ## 4.1 ATAC DARs ALS -------------------------------------------------------

all(
  unique(
    ATAC_DESeq_Results_8SVs_Index$CellType[
      ATAC_DESeq_Results_8SVs_Index$CellTypeLevel=="WNN_L25"
    ]
  ) %in% 
    names(Universe_ATAC_AllCase_ALS_List)
)

nDARs_ALS <- setNames(
  ATAC_DESeq_Results_8SVs_Index$All_q_0.05_ALS[
    ATAC_DESeq_Results_8SVs_Index$CellTypeLevel == "WNN_L25" & 
      ATAC_DESeq_Results_8SVs_Index$Comparison == "All_Cases"
  ], 
  nm = ATAC_DESeq_Results_8SVs_Index$CellType[
    ATAC_DESeq_Results_8SVs_Index$CellTypeLevel == "WNN_L25" & 
      ATAC_DESeq_Results_8SVs_Index$Comparison == "All_Cases"
  ]
)

nGenomicRegions_Universe_ALS <- sapply(
  Universe_ATAC_AllCase_ALS_List[2:length(Universe_ATAC_AllCase_ALS_List)], 
  FUN = function(x){
    return(
      nrow(x)
    )
  }, USE.NAMES = TRUE
)


all(
  names(nDARs_ALS) == names(nGenomicRegions_Universe_ALS)
)  

nGenomicRegions_Universe_ALS[sapply(nGenomicRegions_Universe_ALS, is.null)] <- 1

PctDARs_ALS <- (nDARs_ALS/unlist(nGenomicRegions_Universe_ALS))*100

rm(nGenomicRegions_Universe_ALS, nDARs_ALS)



    ## 4.2 ATAC DARs ALSFTD ----------------------------------------------------

all(
  unique(
    ATAC_DESeq_Results_8SVs_Index$CellType[
      ATAC_DESeq_Results_8SVs_Index$CellTypeLevel=="WNN_L25"
    ]
  ) %in% 
    names(Universe_ATAC_AllCase_ALSFTD_List)
)

nDARs_ALSFTD <- setNames(
  ATAC_DESeq_Results_8SVs_Index$All_q_0.05_ALSFTD[
    ATAC_DESeq_Results_8SVs_Index$CellTypeLevel == "WNN_L25" & 
      ATAC_DESeq_Results_8SVs_Index$Comparison == "All_Cases"
  ], 
  nm = ATAC_DESeq_Results_8SVs_Index$CellType[
    ATAC_DESeq_Results_8SVs_Index$CellTypeLevel == "WNN_L25" & 
      ATAC_DESeq_Results_8SVs_Index$Comparison == "All_Cases"
  ]
)

nGenomicRanges_Universe_ALSFTD <- sapply(
  Universe_ATAC_AllCase_ALSFTD_List[2:length(Universe_ATAC_AllCase_ALSFTD_List)], 
  FUN = function(x){
    return(
      nrow(x)
    )
  }, USE.NAMES = TRUE
)


all(
  names(nDARs_ALSFTD) == names(nGenomicRanges_Universe_ALSFTD)
)  

nGenomicRanges_Universe_ALSFTD[sapply(nGenomicRanges_Universe_ALSFTD, is.null)] <- 1

PctDARs_ALSFTD <- (nDARs_ALSFTD/unlist(nGenomicRanges_Universe_ALSFTD))*100

rm(nGenomicRanges_Universe_ALSFTD, nDARs_ALSFTD)



    ## 4.3 Plot PctDARs --------------------------------------------------------

all(names(PctDARs_ALS) == names(PctDARs_ALSFTD))
plot(PctDARs_ALS, PctDARs_ALSFTD) 

Pct_DARs <- c(
  setNames(PctDARs_ALSFTD, nm = paste0("ALSFTD_", names(PctDARs_ALSFTD))),  
  setNames(PctDARs_ALS, nm = paste0("ALS_", names(PctDARs_ALS)))
)


mtx_Pct_DARs <- matrix(
  Pct_DARs, 
  nrow = length(Pct_DARs), 
  ncol = 1, 
  byrow = FALSE
)
rownames(mtx_Pct_DARs) <- names(Pct_DARs) 


col_fun_Pct_DARs = colorRamp2(
  c(
    0, 
    min(mtx_Pct_DARs[mtx_Pct_DARs>0]), 
    quantile(mtx_Pct_DARs[mtx_Pct_DARs>0])["75%"], 
    median(mtx_Pct_DARs[mtx_Pct_DARs>0]), 
    quantile(mtx_Pct_DARs[mtx_Pct_DARs>0])["75%"], 
    max(mtx_Pct_DARs[mtx_Pct_DARs>0])
  ), 
  c(
    "#EEEEEE", 
    "#fff5cf", 
    "#ffeda3", 
    "#ffdc91", 
    "#ff8000", 
    "red"
  )
) 

circos.clear()
circos.heatmap(mtx_Pct_DARs, col = col_fun_Pct_DARs, cluster=FALSE, cell.border = "#000000")
circos.clear() 


mat = mtx_Pct_DARs
col_fun = col_fun_Pct_DARs
sectors = rep("All", length(mat))

mat_RNA = mtx_Pct_DEGs
col_fun_RNA = col_fun_Pct_DEGs
sectors = rep("All", length(mat_RNA))


circos.initialize("a", xlim = c(0, 24)) 
circos.track(ylim = c(0, 1), track.height=0.2, bg.border=NA, panel.fun = function(x, y) {
  
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    breaks = seq(xlim[1], xlim[2], by = 1)
    n_breaks = length(breaks)
    
  sector.index = CELL_META$sector.index
  m = mat_RNA
  
  col_mat = col_fun_RNA(m)
  nr = nrow(m)
  nc = ncol(m)
  for(i in 1:nc) {
    circos.rect(1:nr - 1, rep(nc - i, nr), 
                1:nr, rep(nc - i + 1, nr), 
                border = "#000000", col = col_mat[, i])
  }
})
circos.track(ylim = c(0, 1), track.height=0.2, bg.border=NA, panel.fun = function(x, y) {
  
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  breaks = seq(xlim[1], xlim[2], by = 1)
  n_breaks = length(breaks)
  
  sector.index = CELL_META$sector.index
  m = mat
  
  col_mat = col_fun(m)
  nr = nrow(m)
  nc = ncol(m)
  for(i in 1:nc) {
    circos.rect(1:nr - 1, rep(nc - i, nr), 
                1:nr, rep(nc - i + 1, nr), 
                border = "#000000", col = col_mat[, i])
  }
})

all(
  rownames(mtx_Pct_DEGs) == 
    rownames(mtx_Pct_DARs)
)

circos.clear()
circos.heatmap(mtx_Pct_DEGs, col = col_fun_Pct_DEGs, cluster=FALSE, cell.border = "#000000")
circos.heatmap(mtx_Pct_DARs, col = col_fun_Pct_DARs, cluster=FALSE, cell.border = "#000000")


circos.clear() 


  ### 5.0 Plot samples Case and Cell type --------------------------------------

n_labels=length(Sample_labels)

Circ_Sample_dend = as.dendrogram(Sample_hclust)

Circ_Sample_dend = color_branches(Circ_Sample_dend, k = 1, col = 1:1)
Circ_Sample_dend_height = attr(Circ_Sample_dend, "height")


circos.clear()
circos.par(
  cell.padding = c(0, 0, 0, 0), 
  track.margin = c(0.006, 0.006), 
  start.degree = 260, 
  gap.degree = 10
)
circos.initialize("a", xlim = c(0, n_labels)) 
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.3, 
             panel.fun = function(x, y) {
               for(i in seq_len(n_labels)) {
                 circos.text(i-0.5, 0, Sample_labels[i], adj = c(0, 0.5), 
                             facing = "clockwise", niceFacing = TRUE,
                             cex = 1)
               }
             })

circos.track(ylim = c(0, 1), track.height=0.05, bg.border=NA, panel.fun = function(x, y) {
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  breaks = seq(xlim[1], xlim[2], by = 1)
  n_breaks = length(breaks)
  circos.rect(breaks[-n_breaks], rep(ylim[1], n_breaks - 1),
              breaks[-1], rep(ylim[2], n_breaks - 1),
              col = ColDict_WNN_L25[match(str_replace_all(str_replace_all(Sample_labels, "ALS_", ""), "ALSFTD_", ""), names(ColDict_WNN_L25))], border = "#000000" )
}) 

circos.track(ylim = c(0, 1), track.height=0.05, bg.border=NA, panel.fun = function(x, y) {
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  breaks = seq(xlim[1], xlim[2], by = 1)
  n_breaks = length(breaks)
  circos.rect(breaks[-n_breaks], rep(ylim[1], n_breaks - 1),
              breaks[-1], rep(ylim[2], n_breaks - 1),
              col = ColDict_Case[match(str_replace_all(str_split(Sample_hclust$labels[Sample_hclust$order], "_", simplify=TRUE)[,1], "ALSFTD", "ALS_FTD"), names(ColDict_Case))], border = "#000000" )
})



circos.track(ylim = c(-0.1, Circ_Sample_dend_height), bg.border = NA, 
             track.height = 0.3, panel.fun = function(x, y) {
               circos.dendrogram(Circ_Sample_dend)
             }) 

circos.clear()
rm(Circ_Sample_dend, Circ_Sample_dend_height)




  ### 6.0 Plot ATAC-RNA Jaccard Index and overlap ------------------------------

ind = which(
  RNA_DESeq_Results_12SVs_Index$CellTypeLevel == "WNN_L25" & 
    RNA_DESeq_Results_12SVs_Index$Comparison == "All_Cases"
)

ATAC_JIs <- list() 
ATAC_SMCs <- list() 
ATAC_Fisher_Ests <- list() 
ATAC_Fisher_Pvals <- list() 
for (i in ind){
  
  int <- length(
    intersect(
      get_DEGs(RNA_DESeq_Results_12SVs_ALS[[i]]), 
      get_DARs_Genes(ATAC_DESeq_Results_8SVs_ALS[[i]])
    )
  ) 
  
  uni <- length(
    union(
      get_DEGs(RNA_DESeq_Results_12SVs_ALS[[i]]), 
      get_DARs_Genes(ATAC_DESeq_Results_8SVs_ALS[[i]])
    ) 
  ) 
  
  ATAC_JIs[[paste0(
    "ALS_", 
    RNA_DESeq_Results_12SVs_Index$CellType[i]
  )]] <- int/uni 
  
  
  tryCatch({
    
    ov1 = table(get_DEGs(RNA_DESeq_Results_12SVs_ALS[[i]]) %in% get_DARs_Genes(ATAC_DESeq_Results_8SVs_ALS[[i]]))
    ov2 = table(get_DARs_Genes(ATAC_DESeq_Results_8SVs_ALS[[i]]) %in% get_DEGs(RNA_DESeq_Results_12SVs_ALS[[i]]))
    
    val1 = as.numeric(ov1["TRUE"])
    val2 = as.numeric(ov1["FALSE"])
    val3 = as.numeric(ov2["FALSE"]) 
    
    
    val4 = length(
      unique(
        c(
          Universe_RNA_AllCase_ALS_List[[
            which(
              names(Universe_RNA_AllCase_ALS_List) == RNA_DESeq_Results_12SVs_Index$CellType[i]
            )
          ]]$ID_10X, 
          
          ATAC_Peak_Links$gene[
            ATAC_Peak_Links$peak %in% Universe_ATAC_AllCase_ALS_List[[
              which(
                names(Universe_ATAC_AllCase_ALS_List) == ATAC_DESeq_Results_8SVs_Index$CellType[i]
              )
            ]]$ID2
          ]
        )
      )
    ) - val1 - val2 -val3
    
    fisher.t <- fisher.test(
      matrix(c(val1, val2, val3, val4), 2, 2)
    )
    
    ATAC_SMCs[[
      paste0(
        "ALS_", 
        RNA_DESeq_Results_12SVs_Index$CellType[i]
      )
    ]] <- (val1 + val4)/sum(val1, val2, val3, val4)
    
    ATAC_Fisher_Ests[[
      paste0(
        "ALS_", 
        RNA_DESeq_Results_12SVs_Index$CellType[i]
      )
    ]] <- fisher.t$estimate
    
    ATAC_Fisher_Pvals[[
      paste0(
        "ALS_", 
        RNA_DESeq_Results_12SVs_Index$CellType[i]
      )
    ]] <- fisher.t$p.value 
    
  }, error = function(e){ 
    
    ATAC_SMCs[[
      paste0(
        "ALS_", 
        RNA_DESeq_Results_12SVs_Index$CellType[i]
      )
    ]] <<- 0
    
    ATAC_Fisher_Ests[[
      paste0(
        "ALS_", 
        RNA_DESeq_Results_12SVs_Index$CellType[i]
      )
    ]] <<- 1
    
    ATAC_Fisher_Pvals[[
      paste0(
        "ALS_", 
        RNA_DESeq_Results_12SVs_Index$CellType[i]
      )
    ]] <<- 1
    
  })
  
  
  
  suppressWarnings(
    rm(int, uni, ov1, ov2, val1, val2, val3, val4, fisher.t) 
  )
  
  
  
  int <- length(
    intersect(
      get_DEGs(RNA_DESeq_Results_12SVs_ALSFTD[[i]]), 
      get_DARs_Genes(ATAC_DESeq_Results_8SVs_ALSFTD[[i]])
    )
  ) 
  
  uni <- length(
    union(
      get_DEGs(RNA_DESeq_Results_12SVs_ALSFTD[[i]]), 
      get_DARs_Genes(ATAC_DESeq_Results_8SVs_ALSFTD[[i]])
    ) 
  ) 
  
  ATAC_JIs[[
    paste0(
      "ALSFTD_", 
      RNA_DESeq_Results_12SVs_Index$CellType[i]
    )
  ]] <- int/uni 
  
  
  tryCatch({
    
    ov1 = table(get_DEGs(RNA_DESeq_Results_12SVs_ALSFTD[[i]]) %in% get_DARs_Genes(ATAC_DESeq_Results_8SVs_ALSFTD[[i]]))
    ov2 = table(get_DARs_Genes(ATAC_DESeq_Results_8SVs_ALSFTD[[i]]) %in% get_DEGs(RNA_DESeq_Results_12SVs_ALSFTD[[i]]))
    
    val1 = as.numeric(ov1["TRUE"])
    val2 = as.numeric(ov1["FALSE"])
    val3 = as.numeric(ov2["FALSE"]) 
    
    
    val4 = length(
      unique(
        c(
          Universe_RNA_AllCase_ALSFTD_List[[
            which(
              names(Universe_RNA_AllCase_ALSFTD_List) == RNA_DESeq_Results_12SVs_Index$CellType[i]
            )
          ]]$ID_10X, 
          
          ATAC_Peak_Links$gene[
            ATAC_Peak_Links$peak %in% Universe_ATAC_AllCase_ALSFTD_List[[
              which(
                names(Universe_ATAC_AllCase_ALSFTD_List) == ATAC_DESeq_Results_8SVs_Index$CellType[i]
              )
            ]]$ID2
          ]
        )
      )
    ) - val1 - val2 -val3
    
    fisher.t <- fisher.test(
      matrix(c(val1, val2, val3, val4), 2, 2)
    )
    
    ATAC_SMCs[[
      paste0(
        "ALSFTD_", 
        RNA_DESeq_Results_12SVs_Index$CellType[i]
      )      
    ]] <- (val1 + val4)/sum(val1, val2, val3, val4)
    
    ATAC_Fisher_Ests[[
      paste0(
        "ALSFTD_", 
        RNA_DESeq_Results_12SVs_Index$CellType[i]
      )
    ]] <- fisher.t$estimate
    
    ATAC_Fisher_Pvals[[
      paste0(
        "ALSFTD_", 
        RNA_DESeq_Results_12SVs_Index$CellType[i]
      )
    ]] <- fisher.t$p.value
    
  }, error = function(e){
    
    ATAC_SMCs[[
      paste0(
        "ALSFTD_", 
        RNA_DESeq_Results_12SVs_Index$CellType[i]
      )      
    ]] <<- 0 
    
    ATAC_Fisher_Ests[[
      paste0(
        "ALSFTD_", 
        RNA_DESeq_Results_12SVs_Index$CellType[i]
      )
    ]] <<- 1
    
    ATAC_Fisher_Pvals[[
      paste0(
        "ALSFTD_", 
        RNA_DESeq_Results_12SVs_Index$CellType[i]
      )
    ]] <<- 1
    
  })
  
  
  suppressWarnings(rm(int, uni, ov1, ov2, val1, val2, val3, val4, fisher.t))
  
}

rm(i)
all(names(ATAC_JIs)==names(ATAC_SMCs))
all(names(ATAC_JIs)==names(ATAC_Fisher_Ests))
all(names(ATAC_Fisher_Ests)==names(ATAC_Fisher_Pvals))
ATAC_Fisher_Qvals <- p.adjust(ATAC_Fisher_Pvals, method = "BH")
sort(ATAC_Fisher_Qvals)

mtx <- unlist(ATAC_JIs)
mtx <- mtx[match(Sample_labels, names(mtx))]
all(names(mtx)==Sample_labels)

Mtx <- matrix(
  mtx, 
  nrow = length(mtx), 
  ncol = 1, 
  byrow = FALSE
)
rownames(Mtx) <- names(mtx) 
all(rownames(Mtx) == Sample_labels)


col_fun_ATAC_JIs = colorRamp2(
  c(
    0, 
    quantile(mtx, na.rm = TRUE)["25%"], 
    median(mtx, na.rm = TRUE), 
    quantile(mtx, na.rm = TRUE)["75%"], 
    1
  ), 
  c(
    "#EEEEEE", 
    "#aceafc", 
    "#4fbcdb", 
    "#0084ff", 
    "#0015ff"
  )
) 

circos.clear()
circos.par(
  cell.padding = c(0, 0, 0, 0), 
  track.margin = c(0.006, 0.006), 
  start.degree = 265, 
  gap.degree = 10
)
circos.initialize(
  "All", 
  xlim = c(0, n_labels)
)  

circos.track(
  ylim = c(0, 1), 
  track.height=0.05, 
  bg.border=NA, 
  panel.fun = function(x, y) {
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    breaks = seq(xlim[1], xlim[2], by = 1)
    n_breaks = length(breaks)
    circos.rect(
      breaks[-n_breaks], 
      rep(ylim[1], n_breaks - 1),
      breaks[-1], 
      rep(ylim[2], n_breaks - 1),
      col = col_fun_ATAC_JIs(Mtx)
    )
  }
) 

mtx <- unlist(ATAC_Fisher_Pvals)
mtx <- mtx[match(Sample_labels, names(mtx))]
mtx <- p.adjust(mtx, method = "BH")
all(names(mtx)==Sample_labels)

Mtx <- matrix(
  mtx, 
  nrow = length(mtx), 
  ncol = 1, 
  byrow = FALSE
)
rownames(Mtx) <- names(mtx) 
all(rownames(Mtx) == Sample_labels)

col_fun_ATAC_FisherPVals = colorRamp2(
  c(
    0, 
    0.05, 
    1
  ), 
  c(
    "#ff0400", 
    "#EEEEEE", 
    "#EEEEEE"
  )
) 

circos.track(
  ylim = c(0, 1), 
  track.height=0.05, 
  bg.border=NA, 
  panel.fun = function(x, y) {
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    breaks = seq(xlim[1], xlim[2], by = 1)
    n_breaks = length(breaks)
    circos.rect(
      breaks[-n_breaks], 
      rep(ylim[1], n_breaks - 1),
      breaks[-1], 
      rep(ylim[2], n_breaks - 1),
      col = col_fun_ATAC_FisherPVals(Mtx)
    )
  }
) 

unlist(ATAC_Fisher_Pvals)

rm(mtx, col_fun_Pct_DARs)




  ### 7.0 Plot FANS_TDP43-RNA Jaccard Index and overlap ------------------------

FANS_TDP43_Signature <- FANS_DE_TDP43_MAST$gene[FANS_DE_TDP43_MAST$fdr<0.05]

ind = which(
  RNA_DESeq_Results_12SVs_Index$CellTypeLevel == "WNN_L25" & 
    RNA_DESeq_Results_12SVs_Index$Comparison == "All_Cases"
)

FANS_TDP43_JIs <- list() 
FANS_TDP43_SMCs <- list() 
FANS_TDP43_Fisher_Ests <- list() 
FANS_TDP43_Fisher_Pvals <- list() 

for (i in ind){
  
  int <- length(
    intersect(
      get_DEGs(RNA_DESeq_Results_12SVs_ALS[[i]]), 
      FANS_TDP43_Signature
    )
  ) 
  
  uni <- length(
    union(
      get_DEGs(RNA_DESeq_Results_12SVs_ALS[[i]]), 
      FANS_TDP43_Signature
    ) 
  ) 
  
  FANS_TDP43_JIs[[paste0(
    "ALS_", 
    RNA_DESeq_Results_12SVs_Index$CellType[i]
  )]] <- int/uni 
  
  
  tryCatch({
    
    ov1 = table(get_DEGs(RNA_DESeq_Results_12SVs_ALS[[i]]) %in% FANS_TDP43_Signature)
    ov2 = table(FANS_TDP43_Signature %in% get_DEGs(RNA_DESeq_Results_12SVs_ALS[[i]]))
    
    val1 = as.numeric(ov1["TRUE"])
    val2 = as.numeric(ov1["FALSE"])
    val3 = as.numeric(ov2["FALSE"]) 
    
    
    val4 = length(
      unique(
        c(
          Universe_RNA_AllCase_ALS_List[[
            which(
              names(Universe_RNA_AllCase_ALS_List) == RNA_DESeq_Results_12SVs_Index$CellType[i]
            )
          ]]$ID_10X, 
          
          FANS_DE_TDP43_MAST$gene
        )
      )
    ) - val1 - val2 -val3
    
    fisher.t <- fisher.test(
      matrix(c(val1, val2, val3, val4), 2, 2)
    )
    
    FANS_TDP43_SMCs[[
      paste0(
        "ALS_", 
        RNA_DESeq_Results_12SVs_Index$CellType[i]
      )
    ]] <- (val1 + val4)/sum(val1, val2, val3, val4)
    
    FANS_TDP43_Fisher_Ests[[
      paste0(
        "ALS_", 
        RNA_DESeq_Results_12SVs_Index$CellType[i]
      )
    ]] <- fisher.t$estimate
    
    FANS_TDP43_Fisher_Pvals[[
      paste0(
        "ALS_", 
        RNA_DESeq_Results_12SVs_Index$CellType[i]
      )
    ]] <- fisher.t$p.value 
    
  }, error = function(e){ 
    
    FANS_TDP43_SMCs[[
      paste0(
        "ALS_", 
        RNA_DESeq_Results_12SVs_Index$CellType[i]
      )
    ]] <<- 0
    
    FANS_TDP43_Fisher_Ests[[
      paste0(
        "ALS_", 
        RNA_DESeq_Results_12SVs_Index$CellType[i]
      )
    ]] <<- 1
    
    FANS_TDP43_Fisher_Pvals[[
      paste0(
        "ALS_", 
        RNA_DESeq_Results_12SVs_Index$CellType[i]
      )
    ]] <<- 1
    
  })
  
  
  
  suppressWarnings(
    rm(int, uni, ov1, ov2, val1, val2, val3, val4, fisher.t) 
  )
  
  
  
  int <- length(
    intersect(
      get_DEGs(RNA_DESeq_Results_12SVs_ALSFTD[[i]]), 
      FANS_TDP43_Signature
    )
  ) 
  
  uni <- length(
    union(
      get_DEGs(RNA_DESeq_Results_12SVs_ALSFTD[[i]]), 
      FANS_TDP43_Signature
    ) 
  ) 
  
  FANS_TDP43_JIs[[
    paste0(
      "ALSFTD_", 
      RNA_DESeq_Results_12SVs_Index$CellType[i]
    )
  ]] <- int/uni 
  
  
  tryCatch({
    
    ov1 = table(get_DEGs(RNA_DESeq_Results_12SVs_ALSFTD[[i]]) %in% FANS_TDP43_Signature)
    ov2 = table(FANS_TDP43_Signature %in% get_DEGs(RNA_DESeq_Results_12SVs_ALSFTD[[i]]))
    
    val1 = as.numeric(ov1["TRUE"])
    val2 = as.numeric(ov1["FALSE"])
    val3 = as.numeric(ov2["FALSE"]) 
    
    
    val4 = length(
      unique(
        c(
          Universe_RNA_AllCase_ALSFTD_List[[
            which(
              names(Universe_RNA_AllCase_ALSFTD_List) == RNA_DESeq_Results_12SVs_Index$CellType[i]
            )
          ]]$ID_10X, 
          
          FANS_DE_TDP43_MAST$gene
        )
      )
    ) - val1 - val2 -val3
    
    fisher.t <- fisher.test(
      matrix(c(val1, val2, val3, val4), 2, 2)
    )
    
    FANS_TDP43_SMCs[[
      paste0(
        "ALSFTD_", 
        RNA_DESeq_Results_12SVs_Index$CellType[i]
      )      
    ]] <- (val1 + val4)/sum(val1, val2, val3, val4)
    
    FANS_TDP43_Fisher_Ests[[
      paste0(
        "ALSFTD_", 
        RNA_DESeq_Results_12SVs_Index$CellType[i]
      )
    ]] <- fisher.t$estimate
    
    FANS_TDP43_Fisher_Pvals[[
      paste0(
        "ALSFTD_", 
        RNA_DESeq_Results_12SVs_Index$CellType[i]
      )
    ]] <- fisher.t$p.value
    
  }, error = function(e){
    
    FANS_TDP43_SMCs[[
      paste0(
        "ALSFTD_", 
        RNA_DESeq_Results_12SVs_Index$CellType[i]
      )      
    ]] <<- 0 
    
    FANS_TDP43_Fisher_Ests[[
      paste0(
        "ALSFTD_", 
        RNA_DESeq_Results_12SVs_Index$CellType[i]
      )
    ]] <<- 1
    
    FANS_TDP43_Fisher_Pvals[[
      paste0(
        "ALSFTD_", 
        RNA_DESeq_Results_12SVs_Index$CellType[i]
      )
    ]] <<- 1
    
  })
  
  
  suppressWarnings(rm(int, uni, ov1, ov2, val1, val2, val3, val4, fisher.t))
  
}

rm(i)
all(names(FANS_TDP43_JIs)==names(FANS_TDP43_SMCs))
all(names(FANS_TDP43_JIs)==names(FANS_TDP43_Fisher_Ests))
all(names(FANS_TDP43_Fisher_Ests)==names(FANS_TDP43_Fisher_Pvals))
FANS_TDP43_Fisher_QVals <- p.adjust(FANS_TDP43_Fisher_Pvals, method = "BH")
sort(FANS_TDP43_Fisher_QVals)


mtx <- unlist(FANS_TDP43_JIs)
names(mtx) <- str_replace_all(names(mtx), ".odds ratio", "")
mtx <- mtx[match(Sample_labels, names(mtx))]
all(names(mtx)==Sample_labels)

Mtx <- matrix(
  mtx, 
  nrow = length(mtx), 
  ncol = 1, 
  byrow = FALSE
)
rownames(Mtx) <- names(mtx) 
all(rownames(Mtx) == Sample_labels)


col_fun_FANS_TDP43_JIs = colorRamp2(
  c(
    0, 
    quantile(mtx, na.rm = TRUE)["25%"], 
    median(mtx, na.rm = TRUE), 
    quantile(mtx, na.rm = TRUE)["75%"], 
    1
  ), 
  c(
    "#EEEEEE", 
    "#aceafc", 
    "#4fbcdb", 
    "#0084ff", 
    "#0015ff"
  )
) 

circos.clear()
circos.par(
  cell.padding = c(0, 0, 0, 0), 
  track.margin = c(0.006, 0.006), 
  start.degree = 265, 
  gap.degree = 10
)
circos.initialize(
  "All", 
  xlim = c(0, n_labels)
)  

circos.track(
  ylim = c(0, 1), 
  track.height=0.05, 
  bg.border=NA, 
  panel.fun = function(x, y) {
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    breaks = seq(xlim[1], xlim[2], by = 1)
    n_breaks = length(breaks)
    circos.rect(
      breaks[-n_breaks], 
      rep(ylim[1], n_breaks - 1),
      breaks[-1], 
      rep(ylim[2], n_breaks - 1),
      col = col_fun_FANS_TDP43_JIs(Mtx)
    )
  }
) 

mtx <- unlist(FANS_TDP43_Fisher_Pvals)
mtx <- mtx[match(Sample_labels, names(mtx))]
mtx <- p.adjust(mtx, method = "BH")
all(names(mtx)==Sample_labels)

Mtx <- matrix(
  mtx, 
  nrow = length(mtx), 
  ncol = 1, 
  byrow = FALSE
)
rownames(Mtx) <- names(mtx) 
all(rownames(Mtx) == Sample_labels)

col_fun_FANS_TDP43_FisherPVals = colorRamp2(
  c(
    0, 
    0.001, 
    1
  ), 
  c(
    "#ff0400", 
    "#EEEEEE", 
    "#EEEEEE"
  )
) 

circos.track(
  ylim = c(0, 1), 
  track.height=0.05, 
  bg.border=NA, 
  panel.fun = function(x, y) {
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    breaks = seq(xlim[1], xlim[2], by = 1)
    n_breaks = length(breaks)
    circos.rect(
      breaks[-n_breaks], 
      rep(ylim[1], n_breaks - 1),
      breaks[-1], 
      rep(ylim[2], n_breaks - 1),
      col = col_fun_FANS_TDP43_FisherPVals(Mtx)
    )
  }
) 

rm(mtx, col_fun_FANS_TDP43_JIs, col_fun_FANS_TDP43_FisherPVals)




  ### 8.0 Plot all layers together ---------------------------------------------



    ## 8.1 Initialize Circos plot ----------------------------------------------

circos.clear()
circos.par(
  cell.padding = c(0, 0, 0, 0), 
  track.margin = c(0.006, 0.006), 
  start.degree = 265, 
  gap.degree = 16
)
circos.initialize(
  "All", 
  xlim = c(0, n_labels)
)  

n_labels=length(Sample_labels)



    ## 8.2 Add labels ----------------------------------------------------------

Sample_labels.tmp <- Sample_labels %>% 
  str_replace_all("ALS_", "") %>% 
    str_replace_all("ALSFTD_", "") %>% 
      str_replace_all("TAFA1_VIP", "TAFA1/VIP") %>% 
        str_replace_all("LAMP5_PAX6", "LAMP5/PAX6") %>% 
          str_replace_all("_", " ")

circos.track(
  ylim = c(0, 1), 
  bg.border = NA, 
  track.height = 0.3, 
  panel.fun = function(x, y) {
    for(i in seq_len(n_labels)) {
      circos.text(
        i-0.5, 
        0, 
        Sample_labels.tmp[i], 
        adj = c(0, 0.5), 
        facing = "clockwise", 
        niceFacing = TRUE,
        cex = 1
      )
    }
  }
)

rm(Sample_labels.tmp)



    ## 8.3 Add FANS_TDP43 JIs and Fisher ---------------------------------------

mtx <- unlist(FANS_TDP43_JIs)
mtx <- mtx[match(Sample_labels, names(mtx))]
all(names(mtx)==Sample_labels)

Mtx <- matrix(
  mtx, 
  nrow = length(mtx), 
  ncol = 1, 
  byrow = FALSE
)
rownames(Mtx) <- names(mtx) 
all(rownames(Mtx) == Sample_labels)


col_fun_FANS_TDP43_JIs = colorRamp2(
  c(
    0, 
    0.1,
    0.2
  ), 
  c(
    "#EEEEEE", 
    "#2cb8b1", 
    "#106B67"
  )
) 
  
circos.track(
  ylim = c(0, 1), 
  track.height=0.025, 
  bg.border=NA, 
  panel.fun = function(x, y) {
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    breaks = seq(xlim[1], xlim[2], by = 1)
    n_breaks = length(breaks)
    circos.rect(
      breaks[-n_breaks], 
      rep(ylim[1], n_breaks - 1),
      breaks[-1], 
      rep(ylim[2], n_breaks - 1),
      col = col_fun_FANS_TDP43_JIs(Mtx)
    )
  }
) 

mtx <- unlist(FANS_TDP43_Fisher_Ests)
names(mtx) <- str_replace_all(
  names(mtx), 
  ".odds ratio", 
  ""
) 

mtx_Sign <- unlist(FANS_TDP43_Fisher_Pvals)
mtx_Sign <- p.adjust(mtx_Sign, method="BH")
all(names(mtx)==names(mtx_Sign))

mtx <- mtx[match(Sample_labels, names(mtx))]   
mtx_Sign <- mtx_Sign[match(Sample_labels, names(mtx_Sign))]   
all(names(mtx)==names(mtx_Sign))

all(names(mtx)==Sample_labels)

Mtx <- matrix(
  mtx, 
  nrow = length(mtx), 
  ncol = 1, 
  byrow = FALSE
)
rownames(Mtx) <- names(mtx) 
all(rownames(Mtx) == Sample_labels)

col_fun_FANS_TDP43_FisherEsts = colorRamp2(
  c(
    0, 
    1, 
    max(Mtx, na.rm=TRUE)
  ), 
  c(
    "#0000FF", 
    "#EEEEEE", 
    "#FF0000"
  )
) 

circos.track(
  ylim = c(0, 1), 
  track.height=0.04, 
  bg.border=NA, 
  panel.fun = function(x, y) {
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    breaks = seq(xlim[1], xlim[2], by = 1)
    n_breaks = length(breaks)
    circos.rect(
      breaks[-n_breaks], 
      rep(ylim[1], n_breaks - 1),
      breaks[-1], 
      rep(ylim[2], n_breaks - 1),
      col = col_fun_FANS_TDP43_FisherEsts(Mtx)
    )
  }
) 

rm(mtx, mtx_Sign, Mtx, col_fun_FANS_TDP43_JIs, col_fun_FANS_TDP43_FisherEsts)



    ## 8.4 Add ATAC-RNA JIs and Fisher -----------------------------------------

mtx <- unlist(ATAC_JIs)
mtx <- mtx[match(Sample_labels, names(mtx))]
all(names(mtx)==Sample_labels)

Mtx <- matrix(
  mtx, 
  nrow = length(mtx), 
  ncol = 1, 
  byrow = FALSE
)
rownames(Mtx) <- names(mtx) 
all(rownames(Mtx) == Sample_labels)

col_fun_ATAC_JIs = colorRamp2(
  c(
    0, 
    0.1,
    0.2
  ), 
  c(
    "#EEEEEE", 
    "#b00ceb", 
    "#7e249e"
  )
) 

circos.track(
  ylim = c(0, 1), 
  track.height=0.02, 
  bg.border=NA, 
  track.margin = c(0.006, 0.05), 
  panel.fun = function(x, y) {
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    breaks = seq(xlim[1], xlim[2], by = 1)
    n_breaks = length(breaks)
    circos.rect(
      breaks[-n_breaks], 
      rep(ylim[1], n_breaks - 1),
      breaks[-1], 
      rep(ylim[2], n_breaks - 1),
      col = col_fun_ATAC_JIs(Mtx)
    )
  }
) 


mtx <- unlist(ATAC_Fisher_Ests)
names(mtx) <- str_replace_all(
  names(mtx), 
  ".odds ratio", 
  ""
) 

mtx_Sign <- unlist(ATAC_Fisher_Pvals)
mtx_Sign <- p.adjust(mtx_Sign, method="BH")
all(names(mtx)==names(mtx_Sign))

mtx <- mtx[match(Sample_labels, names(mtx))]   
mtx_Sign <- mtx_Sign[match(Sample_labels, names(mtx_Sign))]   
all(names(mtx)==names(mtx_Sign))

all(names(mtx)==Sample_labels)

Mtx <- matrix(
  mtx, 
  nrow = length(mtx), 
  ncol = 1, 
  byrow = FALSE
)
rownames(Mtx) <- names(mtx) 
all(rownames(Mtx) == Sample_labels)

col_fun_ATAC_FisherEsts = colorRamp2(
  c(
    0, 
    1, 
    max(Mtx, na.rm=TRUE)
  ), 
  c(
    "#0000FF", 
    "#EEEEEE", 
    "#FF0000"
  )
) 


circos.track(
  ylim = c(0, 1), 
  track.height=0.04, 
  bg.border=NA, 
  panel.fun = function(x, y) {
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    breaks = seq(xlim[1], xlim[2], by = 1)
    n_breaks = length(breaks)
    circos.rect(
      breaks[-n_breaks], 
      rep(ylim[1], n_breaks - 1),
      breaks[-1], 
      rep(ylim[2], n_breaks - 1),
      col = col_fun_ATAC_FisherEsts(Mtx)
    )
  }
) 

rm(mtx, Mtx, col_fun_ATAC_JIs, col_fun_ATAC_FisherEsts)



    ## 8.6 Add Pct DARs track --------------------------------------------------

mtx_Pct_DARs <- matrix(
  Pct_DARs, 
  nrow = length(Pct_DARs), 
  ncol = 1, 
  byrow = FALSE
)
rownames(mtx_Pct_DARs) <- names(Pct_DARs) 
mtx_Pct_DARs <- as.matrix(mtx_Pct_DARs[match(Sample_labels, rownames(mtx_Pct_DARs)),])
all(rownames(mtx_Pct_DARs) == Sample_labels)


col_fun_Pct_DARs = colorRamp2(
  c(
    0, 
    quantile(mtx_Pct_DARs)["25%"], 
    median(mtx_Pct_DARs), 
    quantile(mtx_Pct_DARs)["75%"], 
    max(mtx_Pct_DARs)
  ), 
  c(
    "#EEEEEE", 
    "#c936ff", 
    "#c936ff", 
    "#ac17e3", 
    "#8f02c2"
  )
) 

circos.track(
  ylim = c(0, 1), 
  track.height=0.05, 
  bg.border=NA, 
  panel.fun = function(x, y) {
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    breaks = seq(xlim[1], xlim[2], by = 1)
    n_breaks = length(breaks)
    circos.rect(
      breaks[-n_breaks], 
      rep(ylim[1], n_breaks - 1),
      breaks[-1], 
      rep(ylim[2], n_breaks - 1),
      col = col_fun_Pct_DARs(mtx_Pct_DARs)
    )
  }
) 

rm(mtx_Pct_DARs, col_fun_Pct_DARs)



    ## 8.6 Add Pct DEGs track --------------------------------------------------

mtx_Pct_DEGs <- matrix(
  Pct_DEGs, 
  nrow = length(Pct_DEGs), 
  ncol = 1, 
  byrow = FALSE
)
rownames(mtx_Pct_DEGs) <- names(Pct_DEGs) 
mtx_Pct_DEGs <- as.matrix(mtx_Pct_DEGs[match(Sample_labels, rownames(mtx_Pct_DEGs)),])
all(rownames(mtx_Pct_DEGs) == Sample_labels)


col_fun_Pct_DEGs = colorRamp2(
  c(
    0, 
    quantile(mtx_Pct_DEGs)["25%"], 
    median(mtx_Pct_DEGs), 
    quantile(mtx_Pct_DEGs)["75%"], 
    max(mtx_Pct_DEGs)
  ), 
  c(
    "#EEEEEE", 
    "#ffeda3", 
    "#ffdc91", 
    "#ff8000", 
    "red"
  )
) 

original_breaks <- c(0, 0.3129688, 1.3396221, 2.7930195, 14.0533877)  # Discrete values

# Define the colors corresponding to each discrete break
# One less color than breaks because the intervals are between breaks
colors <- c("#EEEEEE", 
            "#ffeda3", 
            "#ffdc91", 
            "#ff8000", 
            "red")

hist(Pct_DEGs)
plot(density(Pct_DEGs))
# Create a color mapping function using discrete breaks
color_fun <- colorRamp2(original_breaks, colors)

# Create a discrete legend for the intervals
lgd <- Legend(col_fun = color_fun, 
              title = "Discrete Logarithmic Scale", 
              at = original_breaks,  # Log-transformed breaks
              labels = original_breaks)  # Show original values as labels

# Draw the legend
draw(lgd)


# Draw the legend
draw(lgd)


circos.track(
  ylim = c(0, 1), 
  track.height=0.05, 
  bg.border=NA, 
  track.margin=c(0.006, 0.05), 
  panel.fun = function(x, y) {
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    breaks = seq(xlim[1], xlim[2], by = 1)
    n_breaks = length(breaks)
    circos.rect(
      breaks[-n_breaks], 
      rep(ylim[1], n_breaks - 1),
      breaks[-1], 
      rep(ylim[2], n_breaks - 1),
      col = col_fun_Pct_DEGs(mtx_Pct_DEGs)
    )
  }
) 

rm(mtx_Pct_DEGs, col_fun_Pct_DEGs)






    ## 8.7 Add Sample label track ----------------------------------------------

circos.track(
  ylim = c(0, 1), 
  track.height=0.05, 
  bg.border=NA, 
  panel.fun = function(x, y) {
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    breaks = seq(xlim[1], xlim[2], by = 1)
    n_breaks = length(breaks)
    circos.rect(
      breaks[-n_breaks], rep(ylim[1], n_breaks - 1),
      breaks[-1], rep(ylim[2], n_breaks - 1),
      col = ColDict_Case[
        match(
          str_replace_all(
            str_split(Sample_hclust$labels[Sample_hclust$order], "_", simplify=TRUE)[,1], 
            "ALSFTD", 
            "ALS_FTD"
          ), 
          names(ColDict_Case)
        )
      ], 
      border = "#000000" 
    )
  }
)



    ## 8.8 Add Cell type label track -------------------------------------------

circos.track(
  ylim = c(0, 1), 
  track.height=0.03, 
  bg.border=NA, 
  panel.fun = function(x, y){
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    breaks = seq(xlim[1], xlim[2], by = 1)
    n_breaks = length(breaks)
    circos.rect(
      breaks[-n_breaks], rep(ylim[1], n_breaks - 1),
      breaks[-1], rep(ylim[2], n_breaks - 1),
      col = ColDict_WNN_L25[
        match(
          str_replace_all(
            str_replace_all(
              Sample_labels, 
              "ALS_", 
              ""), 
            "ALSFTD_", ""
          ), 
          names(ColDict_WNN_L25)
        )
      ], 
      border = "#000000" 
    )
  }
)



    ## 8.9 Add dendrogram ------------------------------------------------------

Circ_Sample_dend = as.dendrogram(Sample_hclust)
Circ_Sample_dend_height = attr(Circ_Sample_dend, "height")
Circ_Sample_dend <- set(Circ_Sample_dend, "branches_lwd", 2) 

circos.track(
  ylim = c(-0.1, Circ_Sample_dend_height), 
  bg.border = NA, 
  track.height = 0.3, 
  panel.fun = function(x, y) {
    circos.dendrogram(
      Circ_Sample_dend
    )
  }
) 

circos.clear() 
rm(n_labels, Circ_Sample_dend, Circ_Sample_dend_height)



    ## 8.10 Save plot with and without labels ----------------------------------

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig7/", 
    "Transcriptomic_Changes_CircosPlot_woLabels", 
    ".pdf"
  ), 
  width = 8.83, 
  height = 8.83,
  units = "in"
)
  
ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig7/", 
    "Transcriptomic_Changes_CircosPlot_withLabels", 
    ".pdf"
  ), 
  width = 8.83, 
  height = 8.83,
  units = "in"
)




  ### 9.0 Plot legends --------------------------------------------------------



    ## 9.1 Legend Pct DEGs -----------------------------------------------------

ggplot() + 
  
  geom_rect(
    
    dat = data.frame(
      x1 = c(1,2,3,4,5), 
      x2 = c(2,3,4,5,6), 
      y1 = c(1,1,1,1,1), 
      y2 = c(2,2,2,2,2),
      Group=LETTERS[1:5]
    ), 
    aes(
      xmin=x1, 
      xmax=x2, 
      ymin=y1, 
      ymax=y2, 
      fill=Group
    ), 
    
    size = 1.6
  ) + 
  
  
  geom_rect(
    
    dat = data.frame(
      x1 = c(1), 
      x2 = c(6), 
      y1 = c(1), 
      y2 = c(2)
    ), 
    aes(
      xmin=x1, 
      xmax=x2, 
      ymin=y1, 
      ymax=y2
    ), 
    
    col="#000000", 
    fill="#00000000", 
    size = 1.6
    
  ) + 
  
  scale_x_continuous(
    breaks= c(1.5, 2.5, 3.5, 4.5, 5.5), 
    labels = c("0", "0.3", "1.3", "2.8", "14.1" )
  ) + 
  scale_y_continuous(limits=c(1,10), expand=c(0,0,0,0)) + 
  scale_fill_manual(
    values = c(
      "#EEEEEE", 
      "#ffeda3", 
      "#ffdc91", 
      "#ff8000", 
      "red", 
      "#FFFFFF00"
    )
  ) + 
  theme(
    legend.position = "Null", 
    panel.grid = element_blank(), 
    panel.background = element_blank(), 
    axis.text.y = element_blank(), 
    axis.ticks.x = element_line(color="#888888", size=2), 
    axis.ticks.y = element_blank(), 
    axis.ticks.length=unit(-0.4, "cm"), 
    axis.text.x = element_text(size = 14, face="bold.italic", color="#000000"), 
    axis.line = element_blank()
  )


ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig7/", 
    "Legend_DEGs_Pct", 
    ".pdf"
  ), 
  width = 4.5, 
  height = 4.5,
  units = "in"
)



    ## 9.2 Legend Pct DARs -----------------------------------------------------

ggplot() + 
  
  geom_rect(
    
    dat = data.frame(
      x1 = c(1,2,3), 
      x2 = c(2,3,4), 
      y1 = c(1,1,1), 
      y2 = c(2,2,2),
      Group=LETTERS[1:3]
    ), 
    aes(
      xmin=x1, 
      xmax=x2, 
      ymin=y1, 
      ymax=y2, 
      fill=Group
    ), 
    
    size = 1.6
  ) + 
  
  
  geom_rect(
    
    dat = data.frame(
      x1 = c(1), 
      x2 = c(4), 
      y1 = c(1), 
      y2 = c(2)
    ), 
    aes(
      xmin=x1, 
      xmax=x2, 
      ymin=y1, 
      ymax=y2
    ), 
    
    col="#000000", 
    fill="#00000000", 
    size = 1.6
    
  ) + 
  
  scale_x_continuous(
    breaks= c(1.5, 2.5, 3.5), 
    labels = c("0", "0.2", "9.8" )
  ) + 
  scale_y_continuous(limits=c(1,10), expand=c(0,0,0,0)) + 
  scale_fill_manual(
    values = c(
      c(
        "#EEEEEE", 
        "#ac17e3", 
        "#8f02c2"
      )
    )
  ) + 
  theme(
    legend.position = "Null", 
    panel.grid = element_blank(), 
    panel.background = element_blank(), 
    axis.text.y = element_blank(), 
    axis.ticks.x = element_line(color="#888888", size=2), 
    axis.ticks.y = element_blank(), 
    axis.ticks.length=unit(-0.4, "cm"), 
    axis.text.x = element_text(size = 14, face="bold.italic", color="#000000"), 
    axis.line = element_blank()
  )


ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig7/", 
    "Legend_DARs_Pct", 
    ".pdf"
  ), 
  width = 4.5, 
  height = 4.5,
  units = "in"
)





    ## 9.3 Legend ATAC JIs -----------------------------------------------------

ggplot() + 
  
  geom_rect(
    
    dat = data.frame(
      x1 = c(1,2,3), 
      x2 = c(2,3,4), 
      y1 = c(1,1,1), 
      y2 = c(2,2,2),
      Group=LETTERS[1:3]
    ), 
    aes(
      xmin=x1, 
      xmax=x2, 
      ymin=y1, 
      ymax=y2, 
      fill=Group
    ), 
    
    size = 1.6
  ) + 
  
  
  geom_rect(
    
    dat = data.frame(
      x1 = c(1), 
      x2 = c(4), 
      y1 = c(1), 
      y2 = c(2)
    ), 
    aes(
      xmin=x1, 
      xmax=x2, 
      ymin=y1, 
      ymax=y2
    ), 
    
    col="#000000", 
    fill="#00000000", 
    size = 1.6
    
  ) + 
  
  scale_x_continuous(
    breaks= c(1.5, 2.5, 3.5), 
    labels = c("0", "0.1", "0.2" )
  ) + 
  scale_y_continuous(limits=c(1,10), expand=c(0,0,0,0)) + 
  scale_fill_manual(
    values = c(
      c(
        "#EEEEEE", 
        "#b00ceb", 
        "#7e249e"
      )
    )
  ) + 
  theme(
    legend.position = "Null", 
    panel.grid = element_blank(), 
    panel.background = element_blank(), 
    axis.text.y = element_blank(), 
    axis.ticks.x = element_line(color="#888888", size=2), 
    axis.ticks.y = element_blank(), 
    axis.ticks.length=unit(-0.4, "cm"), 
    axis.text.x = element_text(size = 14, face="bold.italic", color="#000000"), 
    axis.line = element_blank()
  )


ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig7/", 
    "Legend_ATAC_JIs", 
    ".pdf"
  ), 
  width = 4.5, 
  height = 4.5,
  units = "in"
)





    ## 9.4 Legend ATAC Fisher --------------------------------------------------

ggplot() + 
  
  geom_rect(
    
    dat = data.frame(
      x1 = c(1,2,3), 
      x2 = c(2,3,4), 
      y1 = c(1,1,1), 
      y2 = c(2,2,2),
      Group=LETTERS[1:3]
    ), 
    aes(
      xmin=x1, 
      xmax=x2, 
      ymin=y1, 
      ymax=y2, 
      fill=Group
    ), 
    
    size = 1.6
  ) + 
  
  
  geom_rect(
    
    dat = data.frame(
      x1 = c(1), 
      x2 = c(4), 
      y1 = c(1), 
      y2 = c(2)
    ), 
    aes(
      xmin=x1, 
      xmax=x2, 
      ymin=y1, 
      ymax=y2
    ), 
    
    col="#000000", 
    fill="#00000000", 
    size = 1.6
    
  ) + 
  
  scale_x_continuous(
    breaks= c(1.5, 2.5, 3.5), 
    labels = c("0", "1", "30" )
  ) + 
  scale_y_continuous(limits=c(1,10), expand=c(0,0,0,0)) + 
  scale_fill_manual(
    values = c(
      c(
        "#0000FF", 
        "#EEEEEE", 
        "#FF0000"
      )
    )
  ) + 
  theme(
    legend.position = "Null", 
    panel.grid = element_blank(), 
    panel.background = element_blank(), 
    axis.text.y = element_blank(), 
    axis.ticks.x = element_line(color="#888888", size=2), 
    axis.ticks.y = element_blank(), 
    axis.ticks.length=unit(-0.4, "cm"), 
    axis.text.x = element_text(size = 14, face="bold.italic", color="#000000"), 
    axis.line = element_blank()
  )


ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig7/", 
    "Legend_ATAC_Fisher", 
    ".pdf"
  ), 
  width = 4.5, 
  height = 4.5,
  units = "in"
)




    ## 9.5 Legend FANS_TDP43 JIs -----------------------------------------------

ggplot() + 
  
  geom_rect(
    
    dat = data.frame(
      x1 = c(1,2,3), 
      x2 = c(2,3,4), 
      y1 = c(1,1,1), 
      y2 = c(2,2,2),
      Group=LETTERS[1:3]
    ), 
    aes(
      xmin=x1, 
      xmax=x2, 
      ymin=y1, 
      ymax=y2, 
      fill=Group
    ), 
    
    size = 1.6
  ) + 
  
  
  geom_rect(
    
    dat = data.frame(
      x1 = c(1), 
      x2 = c(4), 
      y1 = c(1), 
      y2 = c(2)
    ), 
    aes(
      xmin=x1, 
      xmax=x2, 
      ymin=y1, 
      ymax=y2
    ), 
    
    col="#000000", 
    fill="#00000000", 
    size = 1.6
    
  ) + 
  
  scale_x_continuous(
    breaks= c(1.5, 2.5, 3.5), 
    labels = c("0", "0.1", "0.2" )
  ) + 
  scale_y_continuous(limits=c(1,10), expand=c(0,0,0,0)) + 
  scale_fill_manual(
    values = c(
      c(
        "#EEEEEE", 
        "#2cb8b1", 
        "#106B67"
      )
    )
  ) + 
  theme(
    legend.position = "Null", 
    panel.grid = element_blank(), 
    panel.background = element_blank(), 
    axis.text.y = element_blank(), 
    axis.ticks.x = element_line(color="#888888", size=2), 
    axis.ticks.y = element_blank(), 
    axis.ticks.length=unit(-0.4, "cm"), 
    axis.text.x = element_text(size = 14, face="bold.italic", color="#000000"), 
    axis.line = element_blank()
  )


ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig7/", 
    "Legend_FANS_TDP43_JIs", 
    ".pdf"
  ), 
  width = 4.5, 
  height = 4.5,
  units = "in"
)



    ## 9.4 Legend FANS_TDP43 Fisher --------------------------------------------

ggplot() + 
  
  geom_rect(
    
    dat = data.frame(
      x1 = c(1,2,3), 
      x2 = c(2,3,4), 
      y1 = c(1,1,1), 
      y2 = c(2,2,2),
      Group=LETTERS[1:3]
    ), 
    aes(
      xmin=x1, 
      xmax=x2, 
      ymin=y1, 
      ymax=y2, 
      fill=Group
    ), 
    
    size = 1.6
  ) + 
  
  
  geom_rect(
    
    dat = data.frame(
      x1 = c(1), 
      x2 = c(4), 
      y1 = c(1), 
      y2 = c(2)
    ), 
    aes(
      xmin=x1, 
      xmax=x2, 
      ymin=y1, 
      ymax=y2
    ), 
    
    col="#000000", 
    fill="#00000000", 
    size = 1.6
    
  ) + 
  
  scale_x_continuous(
    breaks= c(1.5, 2.5, 3.5), 
    labels = c("0", "1", "4" )
  ) + 
  scale_y_continuous(limits=c(1,10), expand=c(0,0,0,0)) + 
  scale_fill_manual(
    values = c(
      c(
        "#0000FF", 
        "#EEEEEE", 
        "#FF0000"
      )
    )
  ) + 
  theme(
    legend.position = "Null", 
    panel.grid = element_blank(), 
    panel.background = element_blank(), 
    axis.text.y = element_blank(), 
    axis.ticks.x = element_line(color="#888888", size=2), 
    axis.ticks.y = element_blank(), 
    axis.ticks.length=unit(-0.4, "cm"), 
    axis.text.x = element_text(size = 14, face="bold.italic", color="#000000"), 
    axis.line = element_blank()
  )


ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig7/", 
    "Legend_FANS_TDP43_Fisher", 
    ".pdf"
  ), 
  width = 4.5, 
  height = 4.5,
  units = "in"
)




  ### 10.0 Plot overlaps of single cell types -----------------------------------



    ## 10.1 ALSFTD Oligodendrocytes ---------------------------------------------

DEGs <- get_DEGs(
  RNA_DESeq_Results_12SVs_ALSFTD[[
    which(
      RNA_DESeq_Results_12SVs_Index$CellTypeLevel == "WNN_L25" & 
        RNA_DESeq_Results_12SVs_Index$CellType == "Oligodendrocytes" & 
        RNA_DESeq_Results_12SVs_Index$Comparison == "All_Cases"
    )
  ]]
)

DARs_Genes <- get_DARs_Genes(
  ATAC_DESeq_Results_8SVs_ALSFTD[[
    which(
      ATAC_DESeq_Results_8SVs_Index$CellTypeLevel == "WNN_L25" & 
        ATAC_DESeq_Results_8SVs_Index$CellType == "Oligodendrocytes" & 
        ATAC_DESeq_Results_8SVs_Index$Comparison == "All_Cases"
    )
  ]], 
  Links = ATAC_Peak_Links
)

fit <- euler(
  combinations = c(
    "RNA" = length(DEGs), 
    "ATAC" = length(DARs_Genes), 
    "TDP43" = length(FANS_TDP43_Signature),
    "RNA&TDP43" = length(intersect(DEGs, FANS_TDP43_Signature)),
    "ATAC&TDP43" = length(intersect(DARs_Genes, FANS_TDP43_Signature)),
    "RNA&ATAC" = length(intersect(DEGs, DARs_Genes)), 
    "ATAC&RNA&TDP43" = length(intersect(intersect(DEGs, DARs_Genes), FANS_TDP43_Signature))
  )
)
plot(fit)

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig7/", 
    "ALSFTD_OLG_Overlap_Venn_diagram", 
    ".pdf"
  ), 
  width = 4.5, 
  height = 4.5,
  units = "in"
)



rm(DEGs, DARs_Genes, fit)



    ## 10.2 ALS Oligodendrocytes ------------------------------------------------

DEGs <- get_DEGs(
  RNA_DESeq_Results_12SVs_ALS[[
    which(
      RNA_DESeq_Results_12SVs_Index$CellTypeLevel == "WNN_L25" & 
        RNA_DESeq_Results_12SVs_Index$CellType == "Oligodendrocytes" & 
        RNA_DESeq_Results_12SVs_Index$Comparison == "All_Cases"
    )
  ]]
)

DARs_Genes <- get_DARs_Genes(
  ATAC_DESeq_Results_8SVs_ALS[[
    which(
      ATAC_DESeq_Results_8SVs_Index$CellTypeLevel == "WNN_L25" & 
        ATAC_DESeq_Results_8SVs_Index$CellType == "Oligodendrocytes" & 
        ATAC_DESeq_Results_8SVs_Index$Comparison == "All_Cases"
    )
  ]], 
  Links = ATAC_Peak_Links
)

fit <- euler(
  combinations = c(
    "RNA" = length(DEGs), 
    "ATAC" = length(DARs_Genes), 
    "TDP43" = length(FANS_TDP43_Signature),
    "RNA&TDP43" = length(intersect(DEGs, FANS_TDP43_Signature)),
    "ATAC&TDP43" = length(intersect(DARs_Genes, FANS_TDP43_Signature)),
    "RNA&ATAC" = length(intersect(DEGs, DARs_Genes)), 
    "ATAC&RNA&TDP43" = length(intersect(intersect(DEGs, DARs_Genes), FANS_TDP43_Signature))
  )
)
plot(fit)

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig7/", 
    "ALS_OLG_Overlap_Venn_diagram", 
    ".pdf"
  ), 
  width = 4.5, 
  height = 4.5,
  units = "in"
)

rm(DEGs, DARs_Genes, fit)



    ## 10.3 ALSFTD Exc RORB -----------------------------------------------------

DEGs <- get_DEGs(
  RNA_DESeq_Results_12SVs_ALSFTD[[
    which(
      RNA_DESeq_Results_12SVs_Index$CellTypeLevel == "WNN_L25" & 
        RNA_DESeq_Results_12SVs_Index$CellType == "Exc_RORB" & 
        RNA_DESeq_Results_12SVs_Index$Comparison == "All_Cases"
    )
  ]]
)

DARs_Genes <- get_DARs_Genes(
  ATAC_DESeq_Results_8SVs_ALSFTD[[
    which(
      ATAC_DESeq_Results_8SVs_Index$CellTypeLevel == "WNN_L25" & 
        ATAC_DESeq_Results_8SVs_Index$CellType == "Exc_RORB" & 
        ATAC_DESeq_Results_8SVs_Index$Comparison == "All_Cases"
    )
  ]], 
  Links = ATAC_Peak_Links
)

fit <- euler(
  combinations = c(
    "RNA" = length(DEGs), 
    "ATAC" = length(DARs_Genes), 
    "TDP43" = length(FANS_TDP43_Signature),
    "RNA&TDP43" = length(intersect(DEGs, FANS_TDP43_Signature)),
    "ATAC&TDP43" = length(intersect(DARs_Genes, FANS_TDP43_Signature)),
    "RNA&ATAC" = length(intersect(DEGs, DARs_Genes)), 
    "ATAC&RNA&TDP43" = length(intersect(intersect(DEGs, DARs_Genes), FANS_TDP43_Signature))
  )
)
plot(fit)

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig7/", 
    "ALSFTD_ExcRORB_Overlap_Venn_diagram", 
    ".pdf"
  ), 
  width = 4.5, 
  height = 4.5,
  units = "in"
)

rm(DEGs, DARs_Genes, fit)



    ## 10.4 ALSFTD Exc LINC00507 ------------------------------------------------

DEGs <- get_DEGs(
  RNA_DESeq_Results_12SVs_ALSFTD[[
    which(
      RNA_DESeq_Results_12SVs_Index$CellTypeLevel == "WNN_L25" & 
        RNA_DESeq_Results_12SVs_Index$CellType == "Exc_LINC00507" & 
        RNA_DESeq_Results_12SVs_Index$Comparison == "All_Cases"
    )
  ]]
)

DARs_Genes <- get_DARs_Genes(
  ATAC_DESeq_Results_8SVs_ALSFTD[[
    which(
      ATAC_DESeq_Results_8SVs_Index$CellTypeLevel == "WNN_L25" & 
        ATAC_DESeq_Results_8SVs_Index$CellType == "Exc_LINC00507" & 
        ATAC_DESeq_Results_8SVs_Index$Comparison == "All_Cases"
    )
  ]], 
  Links = ATAC_Peak_Links
)

fit <- euler(
  combinations = c(
    "RNA" = length(DEGs), 
    "ATAC" = length(DARs_Genes), 
    "TDP43" = length(FANS_TDP43_Signature),
    "RNA&TDP43" = length(intersect(DEGs, FANS_TDP43_Signature)),
    "ATAC&TDP43" = length(intersect(DARs_Genes, FANS_TDP43_Signature)),
    "RNA&ATAC" = length(intersect(DEGs, DARs_Genes)), 
    "ATAC&RNA&TDP43" = length(intersect(intersect(DEGs, DARs_Genes), FANS_TDP43_Signature))
  )
)
plot(fit)

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig7/", 
    "ALSFTD_ExcLINC00507_Overlap_Venn_diagram", 
    ".pdf"
  ), 
  width = 4.5, 
  height = 4.5,
  units = "in"
)

rm(DEGs, DARs_Genes, fit)



    ## 10.5 ALS Exc RORB --------------------------------------------------------

DEGs <- get_DEGs(
  RNA_DESeq_Results_12SVs_ALS[[
    which(
      RNA_DESeq_Results_12SVs_Index$CellTypeLevel == "WNN_L25" & 
        RNA_DESeq_Results_12SVs_Index$CellType == "Exc_RORB" & 
        RNA_DESeq_Results_12SVs_Index$Comparison == "All_Cases"
    )
  ]]
)

DARs_Genes <- get_DARs_Genes(
  ATAC_DESeq_Results_8SVs_ALS[[
    which(
      ATAC_DESeq_Results_8SVs_Index$CellTypeLevel == "WNN_L25" & 
        ATAC_DESeq_Results_8SVs_Index$CellType == "Exc_RORB" & 
        ATAC_DESeq_Results_8SVs_Index$Comparison == "All_Cases"
    )
  ]], 
  Links = ATAC_Peak_Links
)

fit <- euler(
  combinations = c(
    "RNA" = length(DEGs), 
    "ATAC" = length(DARs_Genes), 
    "TDP43" = length(FANS_TDP43_Signature),
    "RNA&TDP43" = length(intersect(DEGs, FANS_TDP43_Signature)),
    "ATAC&TDP43" = length(intersect(DARs_Genes, FANS_TDP43_Signature)),
    "RNA&ATAC" = length(intersect(DEGs, DARs_Genes)), 
    "ATAC&RNA&TDP43" = length(intersect(intersect(DEGs, DARs_Genes), FANS_TDP43_Signature))
  )
)
plot(fit)

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig7/", 
    "ALS_ExcRORB_Overlap_Venn_diagram", 
    ".pdf"
  ), 
  width = 4.5, 
  height = 4.5,
  units = "in"
)

rm(DEGs, DARs_Genes, fit)



    ## 10.7 ALS Exc LINC00507 ---------------------------------------------------

DEGs <- get_DEGs(
  RNA_DESeq_Results_12SVs_ALS[[
    which(
      RNA_DESeq_Results_12SVs_Index$CellTypeLevel == "WNN_L25" & 
        RNA_DESeq_Results_12SVs_Index$CellType == "Exc_LINC00507" & 
        RNA_DESeq_Results_12SVs_Index$Comparison == "All_Cases"
    )
  ]]
)

DARs_Genes <- get_DARs_Genes(
  ATAC_DESeq_Results_8SVs_ALS[[
    which(
      ATAC_DESeq_Results_8SVs_Index$CellTypeLevel == "WNN_L25" & 
        ATAC_DESeq_Results_8SVs_Index$CellType == "Exc_LINC00507" & 
        ATAC_DESeq_Results_8SVs_Index$Comparison == "All_Cases"
    )
  ]], 
  Links = ATAC_Peak_Links
)

fit <- euler(
  combinations = c(
    "RNA" = length(DEGs), 
    "ATAC" = length(DARs_Genes), 
    "TDP43" = length(FANS_TDP43_Signature),
    "RNA&TDP43" = length(intersect(DEGs, FANS_TDP43_Signature)),
    "ATAC&TDP43" = length(intersect(DARs_Genes, FANS_TDP43_Signature)),
    "RNA&ATAC" = length(intersect(DEGs, DARs_Genes)), 
    "ATAC&RNA&TDP43" = length(intersect(intersect(DEGs, DARs_Genes), FANS_TDP43_Signature))
  )
)
plot(fit)

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig7/", 
    "ALS_ExcLINC00507_Overlap_Venn_diagram", 
    ".pdf"
  ), 
  width = 4.5, 
  height = 4.5,
  units = "in"
)

rm(DEGs, DARs_Genes, fit) 



    ## 10.8 ALS Exc FEZF2 -------------------------------------------------------

DEGs <- get_DEGs(
  RNA_DESeq_Results_12SVs_ALS[[
    which(
      RNA_DESeq_Results_12SVs_Index$CellTypeLevel == "WNN_L25" & 
        RNA_DESeq_Results_12SVs_Index$CellType == "Exc_FEZF2" & 
        RNA_DESeq_Results_12SVs_Index$Comparison == "All_Cases"
    )
  ]]
)

DARs_Genes <- get_DARs_Genes(
  ATAC_DESeq_Results_8SVs_ALS[[
    which(
      ATAC_DESeq_Results_8SVs_Index$CellTypeLevel == "WNN_L25" & 
        ATAC_DESeq_Results_8SVs_Index$CellType == "Exc_FEZF2" & 
        ATAC_DESeq_Results_8SVs_Index$Comparison == "All_Cases"
    )
  ]], 
  Links = ATAC_Peak_Links
)

fit <- euler(
  combinations = c(
    "RNA" = length(DEGs), 
    "ATAC" = length(DARs_Genes), 
    "TDP43" = length(FANS_TDP43_Signature),
    "RNA&TDP43" = length(intersect(DEGs, FANS_TDP43_Signature)),
    "ATAC&TDP43" = length(intersect(DARs_Genes, FANS_TDP43_Signature)),
    "RNA&ATAC" = length(intersect(DEGs, DARs_Genes)), 
    "ATAC&RNA&TDP43" = length(intersect(intersect(DEGs, DARs_Genes), FANS_TDP43_Signature))
  )
)
plot(fit)

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig7/", 
    "ALS_ExcFEZF2_Overlap_Venn_diagram", 
    ".pdf"
  ), 
  width = 4.5, 
  height = 4.5,
  units = "in"
)

rm(DEGs, DARs_Genes, fit)
  ###
  ### 11.0 Plot Correlation of TDP-43 and RNA ----------------------------------



    ## 11.1 ALSFTD Exc LINC00507 -----------------------------------------------

ind <- which(
  RNA_DESeq_Results_12SVs_Index$CellTypeLevel == "WNN_L25" & 
    RNA_DESeq_Results_12SVs_Index$CellType == "Exc_LINC00507" & 
    RNA_DESeq_Results_12SVs_Index$Comparison == "All_Cases"
) 

DEGs <- get_DEGs(
  RNA_DESeq_Results_12SVs_ALSFTD[[ind]]
)

res <- RNA_L2FC_Shrink_Results_12SVs_ALSFTD[[ind]]

df <- data.frame(
  Gene = DEGs[DEGs %in% FANS_TDP43_Signature]
)

df$L2FC_RNA <- res$log2FoldChange[match(df$Gene, rownames(res))]

df$L2FC_TDP43 <- FANS_DE_TDP43_MAST$model_log2FC[match(df$Gene, FANS_DE_TDP43_MAST$gene)]
plot(df$L2FC_TDP43, df$L2FC_RNA)
cor.test(df$L2FC_TDP43, df$L2FC_RNA, method="spearman")

write.csv(df$Gene, "../ALSLINC_TDPRNA_All.csv", quote=FALSE)
write.csv(df$Gene[df$L2FC_RNA<0 & df$L2FC_TDP43<0], "../ALSLINC_TDPRNA_BothDecreased.csv", quote=FALSE)
write.csv(df$Gene[df$L2FC_RNA>0 & df$L2FC_TDP43<0], "../ALSLINC_TDPRNA_DecreasedTDP_Increased_RNA.csv", quote=FALSE)


rm(ind, DEGs, res, df)




  ### 12.0 Plot Correlation of ATAC and RNA ------------------------------------



    ## 12.1 ALSFTD Exc LINC00507 -----------------------------------------------

ind_ATAC <- which(
  ATAC_DESeq_Results_8SVs_Index$CellTypeLevel == "WNN_L25" & 
    ATAC_DESeq_Results_8SVs_Index$CellType == "Exc_LINC00507" & 
    ATAC_DESeq_Results_8SVs_Index$Comparison == "All_Cases"
) 

ind_RNA <- which(
  RNA_DESeq_Results_12SVs_Index$CellTypeLevel == "WNN_L25" & 
    RNA_DESeq_Results_12SVs_Index$CellType == "Exc_LINC00507" & 
    RNA_DESeq_Results_12SVs_Index$Comparison == "All_Cases"
) 

DEGs <- get_DEGs(
  RNA_DESeq_Results_12SVs_ALSFTD[[ind_RNA]]
) 

DARs <- get_DARs(
  ATAC_DESeq_Results_8SVs_ALSFTD[[ind_ATAC]]
)
ind_Link_Gene_Pairs <- which(ATAC_Peak_Links$gene %in% DEGs)

Peak_Gene_Pairs <- data.frame(
  peak = ATAC_Peak_Links$peak[ind_Link_Gene_Pairs], 
  gene = ATAC_Peak_Links$gene[ind_Link_Gene_Pairs]
)
Peak_Gene_Pairs <- Peak_Gene_Pairs[
  Peak_Gene_Pairs$peak %in% DARs,
]

res_RNA <- RNA_L2FC_Shrink_Results_12SVs_ALSFTD[[ind_RNA]]
res_ATAC <- ATAC_L2FC_Shrink_Results_8SVs_ALSFTD[[ind_ATAC]]

Peak_Gene_Pairs$ATAC_L2FC <- res_ATAC$log2FoldChange[match(Peak_Gene_Pairs$peak, rownames(res_ATAC))]
Peak_Gene_Pairs$RNA_L2FC <- res_RNA$log2FoldChange[match(Peak_Gene_Pairs$gene, rownames(res_RNA))]

plot(Peak_Gene_Pairs$ATAC_L2FC, Peak_Gene_Pairs$RNA_L2FC)
cor.test(Peak_Gene_Pairs$ATAC_L2FC, Peak_Gene_Pairs$RNA_L2FC)

write.csv(unique(Peak_Gene_Pairs$gene), "ALSFTDLINC_ATACRNA_Genes.csv", quote=FALSE)

rm(ind_ATAC, ind_RNA, ind_Link_Gene_Pairs, DGEs, DARs, res_RNA, res_ATAC)




  ### 14. GSEA DEG, DAR, TDP43 Genes -------------------------------------------



    ## 14.1 FANS TDP-43 Pathology signature ------------------------------------


#Universe <- data.frame(
#  ID_10X = FANS_DE_TDP43_MAST$gene
#)

Universe <- data.frame(
  ID_10X = FANS_DE_Universe # FANS_DE_Universe list is generated further below 
)
Universe$SYMBOL <- GEX_Features$SYMBOL[match(Universe$ID_10X, GEX_Features$ID_10X)]
Universe$Accession <- GEX_Features$Accession[match(Universe$ID_10X, GEX_Features$ID_10X)]

bitr <- bitr(
  geneID = Universe$SYMBOL, 
  fromType = "SYMBOL", 
  toType = "ENTREZID", 
  OrgDb = "org.Hs.eg.db"
)

Universe$ENTREZ <- bitr$ENTREZID[match(Universe$SYMBOL, bitr$SYMBOL)]
universe = Universe$ENTREZ
universe = universe[!is.na(universe)]

FANS_TDP43_Signature <- FANS_DE_TDP43_MAST$gene[FANS_DE_TDP43_MAST$fdr<0.05 ] 

Genes <- data.frame(
  Gene = FANS_TDP43_Signature
) 
Genes$ID_10X <- Genes$Gene
Genes$SYMBOL <- GEX_Features$SYMBOL[match(Genes$ID_10X, GEX_Features$ID_10X)]
all(Genes$ID_10X %in% Universe$ID_10X)
table(Genes$ID_10X %in% Universe$ID_10X)

Genes$ENTREZID <- Universe$ENTREZ[match(Genes$ID_10X, Universe$ID_10X)]
Genes$Accession <- GEX_Features$Accession[match(Genes$ID_10X, GEX_Features$ID_10X)]

universe = Universe$ENTREZ
universe = universe[!is.na(universe)]


      # Enrichment KEGG 

table(Genes$Accession)

genelist <- Genes$ENTREZID
genelist <- genelist[!is.na(genelist)]


Enrich_Kegg <- enrichKEGG(
  gene = genelist, 
  organism = "hsa", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05, 
  universe = universe
)
View(as.data.frame(Enrich_Kegg))
EnrichGO <- enrichGO(genelist, universe = universe, org.Hs.eg.db)
View(as.data.frame(EnrichGO))
rm(genelist)


# GSEA KEGG 

genelist <- FANS_DE_TDP43_MAST$model_log2FC[
  match(Genes$ID_10X, FANS_DE_TDP43_MAST$gene)
]
names(genelist) <- Genes$ENTREZID

genelist <- genelist[!is.na(names(genelist))]
genelist <- sort(genelist, decreasing = TRUE)

GSEA_Kegg_ByL2FC <- gseKEGG(  
  geneList=genelist, 
  organism = "hsa", 
  eps=0, pvalueCutoff = 0.5
)

View(as.data.frame(GSEA_Kegg_ByL2FC))

rm(Universe, universe, Genes, genelist, Enrich_Kegg, Enrich_GO, GSEA_Kegg_ByL2FC)


FANS <- qread("../Data/SeuratObjects/F2.qrds")
FANS
rownames(FANS)
Assays(FANS)
library(Seurat)
FANS
summary(c(FANS_DE_TDP43_MAST$avg_CPM_disease_1, FANS_DE_TDP43_MAST$avg_CPM_disease_2)) 
FANS_CPM <- NormalizeData(FANS, normalization.method = "RC", scale.factor = 1e6) 

table(FANS$TDP43) 
FANS_DE_TDP43_MAST$BaseMean <- (
  FANS_DE_TDP43_MAST$avg_CPM_disease_2*sum(FANS$TDP43=="TDP43_High") + 
  FANS_DE_TDP43_MAST$avg_CPM_disease_1*sum(FANS$TDP43=="TDP43_Low")
  )/length(FANS$TDP43)
  
summary(
  c(
    FANS_DE_TDP43_MAST$BaseMean[FANS_DE_TDP43_MAST$fdr<0.05 & abs(FANS_DE_TDP43_MAST$model_log2FC)>0.5]
  )
) 

table(
  rowMeans(FANS_CPM) < 
  min(FANS_DE_TDP43_MAST$BaseMean[FANS_DE_TDP43_MAST$fdr<0.05 & abs(FANS_DE_TDP43_MAST$model_log2FC)>0.5])    
)

FANS_DE_Universe <- rownames(FANS_CPM)[
  which(
    rowMeans(FANS_CPM) >= 
      min(FANS_DE_TDP43_MAST$BaseMean[FANS_DE_TDP43_MAST$fdr<0.05 & abs(FANS_DE_TDP43_MAST$model_log2FC)>0.5])    
  )
]


    ## 14.2 ALSFTD Oligodendrocytes Genes --------------------------------------

DEGs = get_DEGs(
  RNA_DESeq_Results_12SVs_ALS[[
    which(
      RNA_DESeq_Results_12SVs_Index$CellTypeLevel == "WNN_L25" & 
        RNA_DESeq_Results_12SVs_Index$CellType == "Oligodendrocytes" & 
          RNA_DESeq_Results_12SVs_Index$Comparison == "All_Cases"
    )
  ]]
)

DAR_Genes <- get_DARs_Genes(
  x=ATAC_DESeq_Results_8SVs_ALS[[
    which(
      ATAC_DESeq_Results_8SVs_Index$CellTypeLevel == "WNN_L25" & 
        ATAC_DESeq_Results_8SVs_Index$CellType == "Oligodendrocytes" & 
        ATAC_DESeq_Results_8SVs_Index$Comparison == "All_Cases"
    )
  ]]
)
table(DEGs %in% FANS_TDP43_Signature)
table(DEGs %in% DAR_Genes)

rm(DEGs, DAR_Genes)



    ## 14.3 ALS ExcRORB Genes --------------------------------------------------


      # Table of gene overlaps -------------------------------------------------

DEGs = get_DEGs(
  RNA_DESeq_Results_12SVs_ALS[[
    which(
      RNA_DESeq_Results_12SVs_Index$CellTypeLevel == "WNN_L25" & 
        RNA_DESeq_Results_12SVs_Index$CellType == "Exc_RORB" & 
        RNA_DESeq_Results_12SVs_Index$Comparison == "All_Cases"
    )
  ]]
)

DAR_Genes <- get_DARs_Genes(
  x=ATAC_DESeq_Results_8SVs_ALS[[
    which(
      ATAC_DESeq_Results_8SVs_Index$CellTypeLevel == "WNN_L25" & 
        ATAC_DESeq_Results_8SVs_Index$CellType == "Exc_RORB" & 
        ATAC_DESeq_Results_8SVs_Index$Comparison == "All_Cases"
    )
  ]]
)
table(DEGs %in% FANS_TDP43_Signature)
table(DEGs %in% DAR_Genes)

rm(DEGs, DAR_Genes)


      # ALS Exc_RORB DEGs TDP-43 Signature overlap -----------------------------

Universe <- Universe_RNA_AllCase_ALS_List[[
  "Exc_RORB"
]]

Universe$Accession <- GEX_Features$Accession[match(Universe$ID_10X, GEX_Features$ID_10X)]

bitr <- bitr(
  geneID = Universe$SYMBOL, 
  fromType = "SYMBOL", 
  toType = "ENTREZID", 
  OrgDb = "org.Hs.eg.db"
)

Universe$ENTREZ <- bitr$ENTREZID[match(Universe$SYMBOL, bitr$SYMBOL)]
universe = Universe$ENTREZ
universe = universe[!is.na(universe)]

Genes = get_DEGs(
  RNA_DESeq_Results_12SVs_ALS[[
    which(
      RNA_DESeq_Results_12SVs_Index$CellTypeLevel == "WNN_L25" & 
        RNA_DESeq_Results_12SVs_Index$CellType == "Exc_RORB" & 
        RNA_DESeq_Results_12SVs_Index$Comparison == "All_Cases"
    )
  ]]
)
Genes = Genes[Genes %in% FANS_TDP43_Signature]

Genes <- data.frame(
  Gene = Genes
) 
Genes$ID_10X <- Genes$Gene
Genes$SYMBOL <- GEX_Features$SYMBOL[match(Genes$ID_10X, GEX_Features$ID_10X)]
all(Genes$ID_10X %in% Universe$ID_10X)
table(Genes$ID_10X %in% Universe$ID_10X)

Genes$ENTREZID <- Universe$ENTREZ[match(Genes$ID_10X, Universe$ID_10X)]
Genes$Accession <- GEX_Features$Accession[match(Genes$ID_10X, GEX_Features$ID_10X)]

universe = Universe$ENTREZ
universe = universe[!is.na(universe)]


      # Enrichment KEGG 

table(Genes$Accession)

genelist <- Genes$ENTREZID
genelist <- genelist[!is.na(genelist)]


Enrich_Kegg <- enrichKEGG(
  gene = genelist, 
  organism = "hsa", 
  pvalueCutoff = 1, 
  qvalueCutoff = 1, 
  universe = universe
)
View(as.data.frame(Enrich_Kegg))


      # Enrichment GO 

EnrichGO <- enrichGO(genelist, universe = universe, org.Hs.eg.db)
View(as.data.frame(EnrichGO))
rm(genelist)



      # GSEA KEGG 

res <- RNA_L2FC_Shrink_Results_12SVs_ALS[[
  which(
    RNA_DESeq_Results_12SVs_Index$CellTypeLevel == "WNN_L25" & 
      RNA_DESeq_Results_12SVs_Index$CellType == "Exc_RORB" & 
      RNA_DESeq_Results_12SVs_Index$Comparison == "All_Cases"
  ) 
]] %>% 
  as.data.frame()  


genelist <- res$log2FoldChange[
  match(Genes$ID_10X, rownames(res))
]
names(genelist) <- Genes$ENTREZID

genelist <- genelist[!is.na(names(genelist))]
genelist <- sort(genelist, decreasing = TRUE)

GSEA_Kegg_ByL2FC <- gseKEGG(  
  geneList=genelist, 
  organism = "hsa", 
  eps=0, pvalueCutoff = 1, minGSSize = 1
)

View(as.data.frame(GSEA_Kegg_ByL2FC))

a <- Return_GSEA_as_Tree(
  GSEA_Kegg_ByL2FC
)

library(ggraph)
ggraph(a, 'dendrogram', height = height) + 
  geom_edge_elbow() + 
  geom_node_point(aes(col=))


plot_GSEA_as_Tree(
  GSEA_Kegg_ByL2FC, 
  name=paste0(
    "../Data/Visualization/Figures/Fig7/", 
    "GSEA_ALS_ExcRORB_DEGs_FANS_Tree", 
    ".pdf"
  ), 
  save = TRUE, width=11.07, height=4.23
)

## 
Genes = get_DEGs(
  RNA_DESeq_Results_12SVs_ALSFTD[[
    which(
      RNA_DESeq_Results_12SVs_Index$CellTypeLevel == "WNN_L25" & 
        RNA_DESeq_Results_12SVs_Index$CellType == "Exc_RORB" & 
        RNA_DESeq_Results_12SVs_Index$Comparison == "All_Cases"
    )
  ]]
)

Genes = Genes[!Genes %in% FANS_TDP43_Signature]

Genes <- data.frame(
  Gene = Genes
) 
Genes$ID_10X <- Genes$Gene
Genes$SYMBOL <- GEX_Features$SYMBOL[match(Genes$ID_10X, GEX_Features$ID_10X)]
all(Genes$ID_10X %in% Universe$ID_10X)
table(Genes$ID_10X %in% Universe$ID_10X)

Genes$ENTREZID <- Universe$ENTREZ[match(Genes$ID_10X, Universe$ID_10X)]
Genes$Accession <- GEX_Features$Accession[match(Genes$ID_10X, GEX_Features$ID_10X)]

universe = Universe$ENTREZ
universe = universe[!is.na(universe)]


# Enrichment KEGG 

table(Genes$Accession)

genelist <- Genes$ENTREZID
genelist <- genelist[!is.na(genelist)]


Enrich_Kegg <- enrichKEGG(
  gene = genelist, 
  organism = "hsa", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05, 
  universe = universe
)
View(as.data.frame(Enrich_Kegg))
EnrichGO <- enrichGO(genelist, universe = universe, org.Hs.eg.db)
View(as.data.frame(EnrichGO))
rm(genelist)


# GSEA KEGG 

res <- RNA_L2FC_Shrink_Results_12SVs_ALSFTD[[
  which(
    RNA_DESeq_Results_12SVs_Index$CellTypeLevel == "WNN_L25" & 
      RNA_DESeq_Results_12SVs_Index$CellType == "Exc_RORB" & 
      RNA_DESeq_Results_12SVs_Index$Comparison == "All_Cases"
  ) 
]] %>% 
  as.data.frame()  


genelist <- res$log2FoldChange[
  match(Genes$ID_10X, rownames(res))
]
names(genelist) <- Genes$ENTREZID

genelist <- genelist[!is.na(names(genelist))]
genelist <- sort(genelist, decreasing = TRUE)

GSEA_Kegg_ByL2FC <- gseKEGG(  
  geneList=genelist, 
  organism = "hsa", 
  eps=0, pvalueCutoff = 0.5
)

View(as.data.frame(GSEA_Kegg_ByL2FC))
rm(Universe, universe, Genes, genelist, Enrich_Kegg, Enrich_GO, GSEA_Kegg_ByL2FC)



    ## 14.3 ALSFTD ExcLINC00507 Genes -----------------------------------------------


      # Table of gene overlaps -------------------------------------------------

DEGs = get_DEGs(
  RNA_DESeq_Results_12SVs_ALS[[
    which(
      RNA_DESeq_Results_12SVs_Index$CellTypeLevel == "WNN_L25" & 
        RNA_DESeq_Results_12SVs_Index$CellType == "Exc_LINC00507" & 
        RNA_DESeq_Results_12SVs_Index$Comparison == "All_Cases"
    )
  ]]
)

DAR_Genes <- get_DARs_Genes(
  x=ATAC_DESeq_Results_8SVs_ALS[[
    which(
      ATAC_DESeq_Results_8SVs_Index$CellTypeLevel == "WNN_L25" & 
        ATAC_DESeq_Results_8SVs_Index$CellType == "Exc_LINC00507" & 
        ATAC_DESeq_Results_8SVs_Index$Comparison == "All_Cases"
    )
  ]]
)
table(DEGs %in% FANS_TDP43_Signature)
table(DEGs %in% DAR_Genes)


      # ALS Exc_RORB DEGs TDP-43 Signature overlap -----------------------------

Universe <- Universe_RNA_AllCase_ALS_List[[
  "Exc_LINC00507"
]]

Universe$Accession <- GEX_Features$Accession[match(Universe$ID_10X, GEX_Features$ID_10X)]

bitr <- bitr(
  geneID = Universe$SYMBOL, 
  fromType = "SYMBOL", 
  toType = "ENTREZID", 
  OrgDb = "org.Hs.eg.db"
)

Universe$ENTREZ <- bitr$ENTREZID[match(Universe$SYMBOL, bitr$SYMBOL)]
universe = Universe$ENTREZ
universe = universe[!is.na(universe)]

Genes = get_DEGs(
  RNA_DESeq_Results_12SVs_ALS[[
    which(
      RNA_DESeq_Results_12SVs_Index$CellTypeLevel == "WNN_L25" & 
        RNA_DESeq_Results_12SVs_Index$CellType == "Exc_LINC00507" & 
        RNA_DESeq_Results_12SVs_Index$Comparison == "All_Cases"
    )
  ]]
)
Genes = Genes[Genes %in% FANS_TDP43_Signature]

Genes <- data.frame(
  Gene = Genes
) 
Genes$ID_10X <- Genes$Gene
Genes$SYMBOL <- GEX_Features$SYMBOL[match(Genes$ID_10X, GEX_Features$ID_10X)]
all(Genes$ID_10X %in% Universe$ID_10X)
table(Genes$ID_10X %in% Universe$ID_10X)

Genes$ENTREZID <- Universe$ENTREZ[match(Genes$ID_10X, Universe$ID_10X)]
Genes$Accession <- GEX_Features$Accession[match(Genes$ID_10X, GEX_Features$ID_10X)]

universe = Universe$ENTREZ
universe = universe[!is.na(universe)]


# Enrichment KEGG 

table(Genes$Accession)

genelist <- Genes$ENTREZID
genelist <- genelist[!is.na(genelist)]


Enrich_Kegg <- enrichKEGG(
  gene = genelist, 
  organism = "hsa", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05, 
  universe = universe
)

View(as.data.frame(Enrich_Kegg))


# Enrichment GO 

EnrichGO <- enrichGO(genelist, universe = universe, org.Hs.eg.db)
View(as.data.frame(EnrichGO))
rm(genelist)



# GSEA KEGG 

res <- RNA_L2FC_Shrink_Results_12SVs_ALS[[
  which(
    RNA_DESeq_Results_12SVs_Index$CellTypeLevel == "WNN_L25" & 
      RNA_DESeq_Results_12SVs_Index$CellType == "Exc_LINC00507" & 
      RNA_DESeq_Results_12SVs_Index$Comparison == "All_Cases"
  ) 
]] %>% 
  as.data.frame()  



genelist <- res$log2FoldChange[
  match(Genes$ID_10X, rownames(res))
]
names(genelist) <- Genes$ENTREZID

genelist <- genelist[!is.na(names(genelist))]
genelist <- sort(genelist, decreasing = TRUE)

GSEA_Kegg_ByL2FC <- gseKEGG(  
  geneList=genelist, 
  minGSSize = 5, 
  maxGSSize = Inf,
  organism = "hsa", 
  eps=0, pvalueCutoff = 0.05
)

View(as.data.frame(GSEA_Kegg_ByL2FC)) 

GSEA_GO_ByL2FC <- gseGO(  
  OrgDb = org.Hs.eg.db, 
  geneList=genelist, 
  eps=0, pvalueCutoff = 0.05
)

View(as.data.frame(GSEA_GO_ByL2FC))


a <- Return_GSEA_as_Tree(
  GSEA_Kegg_ByL2FC
)

library(ggraph)
ggraph(a, 'dendrogram', height = height) + 
  geom_edge_elbow() + 
  geom_node_point(aes(col=))


plot_GSEA_as_Tree(
  GSEA_Kegg_ByL2FC, 
  name=paste0(
    "../Data/Visualization/Figures/Fig7/", 
    "GSEA_ALS_ExcRORB_DEGs_FANS_Tree", 
    ".pdf"
  ), 
  save = TRUE, width=11.07, height=4.23
)

## 
Genes = get_DEGs(
  RNA_DESeq_Results_12SVs_ALS[[
    which(
      RNA_DESeq_Results_12SVs_Index$CellTypeLevel == "WNN_L25" & 
        RNA_DESeq_Results_12SVs_Index$CellType == "Exc_FEZF2" & 
        RNA_DESeq_Results_12SVs_Index$Comparison == "All_Cases"
    )
  ]]
)

Genes = Genes[!Genes %in% FANS_TDP43_Signature]

Genes <- data.frame(
  Gene = Genes
) 
Genes$ID_10X <- Genes$Gene
Genes$SYMBOL <- GEX_Features$SYMBOL[match(Genes$ID_10X, GEX_Features$ID_10X)]
all(Genes$ID_10X %in% Universe$ID_10X)
table(Genes$ID_10X %in% Universe$ID_10X)

Genes$ENTREZID <- Universe$ENTREZ[match(Genes$ID_10X, Universe$ID_10X)]
Genes$Accession <- GEX_Features$Accession[match(Genes$ID_10X, GEX_Features$ID_10X)]

universe = Universe$ENTREZ
universe = universe[!is.na(universe)]


# Enrichment KEGG 

table(Genes$Accession)

genelist <- Genes$ENTREZID
genelist <- genelist[!is.na(genelist)]


Enrich_Kegg <- enrichKEGG(
  gene = genelist, 
  organism = "hsa", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05, 
  universe = universe
)
View(as.data.frame(Enrich_Kegg))
EnrichGO <- enrichGO(genelist, universe = universe, org.Hs.eg.db)
View(as.data.frame(EnrichGO))
rm(genelist)


# GSEA KEGG 

res <- RNA_L2FC_Shrink_Results_12SVs_ALS[[
  which(
    RNA_DESeq_Results_12SVs_Index$CellTypeLevel == "WNN_L25" & 
      RNA_DESeq_Results_12SVs_Index$CellType == "Exc_FEZF2" & 
      RNA_DESeq_Results_12SVs_Index$Comparison == "All_Cases"
  ) 
]] %>% 
  as.data.frame()  


genelist <- res$log2FoldChange[
  match(Genes$ID_10X, rownames(res))
]
names(genelist) <- Genes$ENTREZID

genelist <- genelist[!is.na(names(genelist))]
genelist <- sort(genelist, decreasing = TRUE)

GSEA_Kegg_ByL2FC <- gseKEGG(  
  geneList=genelist, 
  organism = "hsa", 
  eps=0, pvalueCutoff = 0.5
)

View(as.data.frame(GSEA_Kegg_ByL2FC))
rm(Universe, universe, Genes, genelist, Enrich_Kegg, Enrich_GO, GSEA_Kegg_ByL2FC)

DE_MAST_ALSvsHC_WNN_L25_Exc_LINC00507 <- qread("../Data/DE/WNN/Single_Cell/MAST/ALS/WNN_L25/DE_MAST_ALSvsHC_WNN_L25_Exc_LINC00507.qrds")
DE_MAST_ALSvsHC_WNN_L25_Exc_FEZF2 <- qread("../Data/DE/WNN/Single_Cell/MAST/ALS/WNN_L25/DE_MAST_ALSvsHC_WNN_L25_Exc_FEZF2.qrds")

table(DE_MAST_ALSvsHC_WNN_L25_Exc_LINC00507$fdr < 0.05)
table(DE_MAST_ALSvsHC_WNN_L25_Exc_FEZF2$fdr < 0.05)

DE_MAST_ALSvsHC_WNN_L25_Exc_LINC00507_DEGs <- DE_MAST_ALSvsHC_WNN_L25_Exc_LINC00507$gene[
  DE_MAST_ALSvsHC_WNN_L25_Exc_LINC00507$fdr<0.05
]
DE_MAST_ALSvsHC_WNN_L25_Exc_FEZF2_DEGs <- DE_MAST_ALSvsHC_WNN_L25_Exc_FEZF2$gene[
  DE_MAST_ALSvsHC_WNN_L25_Exc_FEZF2$fdr<0.05
]

table(DE_MAST_ALSvsHC_WNN_L25_Exc_LINC00507_DEGs %in% DEGs)
table(DEGs %in% DE_MAST_ALSvsHC_WNN_L25_Exc_LINC00507_DEGs) 

table(DE_MAST_ALSvsHC_WNN_L25_Exc_FEZF2_DEGs %in% DEGs)
table(DEGs %in% DE_MAST_ALSvsHC_WNN_L25_Exc_FEZF2_DEGs)

BothMethods <- data.frame(
  #Gene = DEGs[DEGs %in% DE_MAST_ALSvsHC_WNN_L25_Exc_LINC00507_DEGs]
  #Gene =DE_MAST_ALSvsHC_WNN_L25_Exc_LINC00507_DEGs
  Gene = DEGs 
  #Gene = unique(c(DEGs, DE_MAST_ALSvsHC_WNN_L25_Exc_LINC00507_DEGs))
)
BothMethods$L2FC_DESeq <- res$log2FoldChange[match(BothMethods$Gene, rownames(res))]  
BothMethods$L2FC_MAST <- DE_MAST_ALSvsHC_WNN_L25_Exc_LINC00507$model_log2FC[match(BothMethods$Gene, DE_MAST_ALSvsHC_WNN_L25_Exc_LINC00507$gene)]  
plot(BothMethods$L2FC_DESeq, BothMethods$L2FC_MAST)

BothMethods <- BothMethods[BothMethods$L2FC_DESeq < -0.1 & BothMethods$L2FC_MAST < -0.1,]
BothMethods$ID_10X <- BothMethods$Gene
BothMethods$SYMBOL <- GEX_Features$SYMBOL[match(BothMethods$ID_10X, GEX_Features$ID_10X)]
all(BothMethods$ID_10X %in% Universe$ID_10X)
table(BothMethods$ID_10X %in% Universe$ID_10X)

BothMethods$ENTREZID <- Universe$ENTREZ[match(BothMethods$ID_10X, Universe$ID_10X)]
BothMethods$Accession <- GEX_Features$Accession[match(BothMethods$ID_10X, GEX_Features$ID_10X)]

genelist <- res$log2FoldChange[
  match(BothMethods$Gene, rownames(res))
]
names(genelist) <- BothMethods$ENTREZID

genelist <- DE_MAST_ALSvsHC_WNN_L25_Exc_LINC00507$model_log2FC[
  match(BothMethods$Gene, DE_MAST_ALSvsHC_WNN_L25_Exc_LINC00507$gene)
]
names(genelist) <- BothMethods$ENTREZID

genelist <- genelist[!is.na(names(genelist))]
genelist <- sort(genelist, decreasing = TRUE)

GSEA_Kegg_ByL2FC <- gseKEGG(  
  geneList=genelist, 
  organism = "hsa", 
  eps=0, pvalueCutoff = 0.05
)

View(as.data.frame(GSEA_Kegg_ByL2FC))

ALS_Signature <- qread("../Data/Annotations/Signatures/RNA/ALS_ALSFTD_WNN_L25_Signature.qrds")

table(ALS_Signature$ID_10X %in% FANS_TDP43_Signature)  
table(FANS_TDP43_Signature %in% ALS_Signature$ID_10X)  

table(ALS_Signature$Gene %in% ATAC_Peak_Links$gene)
table(unique(ATAC_Peak_Links$gene) %in% ALS_Signature$Gene)


