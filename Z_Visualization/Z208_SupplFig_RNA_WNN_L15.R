

  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 

library(qs)
library(tidyverse)
library(DESeq2)



qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)


M0_RNA <- qread("../Data/SeuratObjects/M0_RNA.qrds", nthr)
Meta <- M0_RNA@meta.data %>% 
  group_by(ID) %>% 
  summarise(Case = dplyr::first(Case), CaseType = dplyr::first(Case_Type), Sex = dplyr::first(Sex))

table(Meta$CaseType, Meta$Sex)




  ### 1.0 Load data ------------------------------------------------------------



SVAs_15SVs_Index <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/1_15_SVs/SVAs_15SVs_Index.qrds"
  )
)

DE_C9 <- qread("../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9_NonC9/SVA/12_SVs/DDS_SVA_12sv_list.qrds", nthr=nthr)
C9_DDS1 <- DE_C9[[1]]
C9_DDS2 <- DE_C9[[2]]

resultsNames(C9_DDS1)
resultsNames(C9_DDS2)

design(C9_DDS1) <- ~SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + Case_Type
C9_DDS1 <- DESeq(C9_DDS1)

resultsNames(C9_DDS1)
resultsNames(C9_DDS2)

C9_Res1 <- results(C9_DDS1, name="Case_Type_ALS_FTD_vs_C9_ALS_FTD", alpha=0.05, lfcThreshold = 0) 
C9_Res2 <- results(C9_DDS2, name="Rand_R2_vs_R1", alpha=0.05, lfcThreshold = 0) 

C9_Res1 %>% 
  summary(alpha=0.05) 

C9_Res2 %>% 
  summary(alpha=0.05) 


SVAs_15SVs_Index %>% 
  filter(CellTypeLevel=="WNN_L15") %>% 
  pull(CellType) %>% 
  unique() 

SVA_DDS_List_Case_Glia <- list() 
SVA_DDS_List_Rand_Glia <- list()  

SVA_DDS_List_Case_Exc_Neurons <- list() 
SVA_DDS_List_Rand_Exc_Neurons <- list()  

SVA_DDS_List_Case_Inh_Neurons <- list() 
SVA_DDS_List_Rand_Inh_Neurons <- list()  

NoCov <- qread(
  paste0(
  "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/",
  "DDS_DESEqed_list", 
  ".qrds"
  ), 
  nthr=nthr
) 

NoCov_Index <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/", 
    "DESeq_Results_Index", 
    ".qrds" 
  ), 
  nthr=nthr
)



which(
  SVAs_15SVs_Index$CellTypeLevel=="WNN_L15" & 
)


SVA_DDS_List_Case <- list() 
SVA_DDS_List_Rand <- list() 

for (i in 1:15){
  tmp <- qread(
    paste0(
      "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/1_15_SVs/DDS_", 
      i, 
      "SVs_list_DESeqed", 
      ".qrds"
    ),
    nthr=nthr
  )
  SVA_DDS_List_Case[[paste0(i, "SVs")]] <- tmp[[1]]
  SVA_DDS_List_Rand[[paste0(i, "SVs")]] <- tmp[[2]] 
  rm(tmp)
}

SVA_DDS_List_Case[["NoCov"]] <- NoCov[[which(
  NoCov_Index$CellTypeLevel=="AllCells" & 
    NoCov_Index$CellType=="AllCells" & 
    NoCov_Index$Comparison=="All_Cases")
]]

SVA_DDS_List_Rand[["NoCov"]] <- NoCov[[which(
  NoCov_Index$CellTypeLevel=="AllCells" & 
    NoCov_Index$CellType=="AllCells" & 
    NoCov_Index$Comparison=="Rand")
]]

rm(NoCov, NoCov_Index)

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
  
GEX_Features <- qread(
  paste0(
    "../Data/Annotations/Genes/", 
    "GEX_Features", 
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

get_nDEGs <- function(x, name){
  return(nrow(get_sign_genes_df(results(x, alpha=0.05, lfcThreshold = 0, name = name))))
}

plot_GSEA_as_Tree <- function(GSEA_Results, name, save=TRUE, width, height){
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



  ### 3.0 Plot SVA DEG numbers -------------------------------------------------

nDEGs <- data.frame(
  ALSvsHC = c(
    sapply(SVA_DDS_List_Case, FUN=get_nDEGs, name="Case_ALS_vs_HC"), 
    sapply(SVA_DDS_List_Rand, FUN=get_nDEGs, name="Rand_R2_vs_R1")
  ), 
  Group = rep(c("ALSvsHC", "Random control"), each = length(SVA_DDS_List_Case))
)
   
  
nDEGs$SVs <- names(SVA_DDS_List_Case)
nDEGs$SVs <- factor(nDEGs$SVs, levels=unique(nDEGs$SVs))
  
ggplot(nDEGs) + 
  geom_line(aes(group=Group, col=Group)) + 
  aes(SVs, ALSvsHC, fill=Group) + 
  geom_point(size=4, pch=21, col="#000000") + 
  scale_x_discrete(labels=c(0:15)) + 
  scale_color_manual(values=c("ALSvsHC"=ColDict_Case[["ALS"]], "Random control"="#AAAAAA")) + 
  scale_fill_manual(values=c("ALSvsHC"=ColDict_Case[["ALS"]], "Random control"="#AAAAAA")) + 
  theme_classic() + 
  theme(
    axis.title.x = element_blank(), 
    axis.line = element_line(linewidth = 0.8, color="#000000"), 
    axis.ticks = element_line(linewidth = 0.8, color="#000000"), 
    axis.ticks.length=unit(4, "points"), 
    axis.text.x = element_text(size=10, color="#000000", angle=45, hjust=1, face="bold"), 
    axis.text.y = element_text(size=10, color="#000000", face="bold.italic")
  )

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/SupplFig_RNA_Complete_Pseudobulk/", 
    "DEGs_SVA_1_15SVs_ALS", 
    ".pdf"
  ),
  width = 150, 
  height = 110, 
  units="mm"
)  




nDEGs <- data.frame(
  ALSFTDvsHC = c(
    sapply(SVA_DDS_List_Case, FUN=get_nDEGs, name="Case_ALS_FTD_vs_HC"), 
    sapply(SVA_DDS_List_Rand, FUN=get_nDEGs, name="Rand_R3_vs_R1")
  ), 
  Group = rep(c("ALSFTDvsHC", "Random control"), each = length(SVA_DDS_List_Case))
)


nDEGs$SVs <- names(SVA_DDS_List_Case)
nDEGs$SVs <- factor(nDEGs$SVs, levels=unique(nDEGs$SVs))

ggplot(nDEGs) + 
  geom_line(aes(group=Group, col=Group)) + 
  aes(SVs, ALSFTDvsHC, fill=Group) + 
  geom_point(size=4, pch=21, col="#000000") + 
  scale_x_discrete(labels=c(0:15)) + 
  scale_color_manual(values=c("ALSFTDvsHC"=ColDict_Case[["ALS_FTD"]], "Random control"="#AAAAAA")) + 
  scale_fill_manual(values=c("ALSFTDvsHC"=ColDict_Case[["ALS_FTD"]], "Random control"="#AAAAAA")) + 
  theme_classic() + 
  theme(
    axis.title.x = element_blank(), 
    axis.line = element_line(linewidth = 0.8, color="#000000"), 
    axis.ticks = element_line(linewidth = 0.8, color="#000000"), 
    axis.ticks.length=unit(4, "points"), 
    axis.text.x = element_text(size=10, color="#000000", angle=45, hjust=1, face="bold"), 
    axis.text.y = element_text(size=10, color="#000000", face="bold.italic")
  )

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/SupplFig_RNA_Complete_Pseudobulk/", 
    "DEGs_SVA_1_15SVs_ALSFTD", 
    ".pdf"
  ),
  width = 150, 
  height = 110, 
  units="mm"
) 

DDS_12SVs_Case <- SVA_DDS_List_Case[["12SVs"]]
DDS_12SVs_Rand <- SVA_DDS_List_Rand[["12SVs"]]

rm(nDEGs, SVA_DDS_List_Case, SVA_DDS_List_Rand, SVAs_15SVs_Index, i)




  ### 4.0 Volcano Plots --------------------------------------------------------



    ## 4.1 ALS DESeq -----------------------------------------------------------

ALS_DESeq <- as.data.frame(
  results(DDS_12SVs_Case, name="Case_ALS_vs_HC", alpha=0.05, lfcThreshold = 0)
)
ALS_DESeq <- ALS_DESeq[!is.na(ALS_DESeq$padj),]

ALS_DESeq$Type = "NA"
ALS_DESeq$Type[ALS_DESeq$log2FoldChange<0 & ALS_DESeq$padj<0.05] <- "Down"
ALS_DESeq$Type[ALS_DESeq$log2FoldChange>0 & ALS_DESeq$padj<0.05] <- "Up"
table(ALS_DESeq$Type)

ggplot() + 
  aes(log2FoldChange, -log10(padj)) + 
  geom_hline(yintercept = c(-log10(0.05)), col="#AAAAAA", linewidth=1) + 
  ggrastr::geom_point_rast(aes(log2FoldChange, -log10(padj)), dat=ALS_DESeq[ALS_DESeq$Type=="NA",],alpha=0.2, col="#000000", fill="#000000", pch=21, size=3) + 
  ggrastr::geom_point_rast(aes(log2FoldChange, -log10(padj)), dat=ALS_DESeq[ALS_DESeq$Type=="Down",], pch=21, col="#000000", fill="#4444FF",alpha=0.4, size=3) + 
  ggrastr::geom_point_rast(aes(log2FoldChange,-log10(padj)), dat=ALS_DESeq[ALS_DESeq$Type=="Up",], pch=21, col="#000000", fill="#FF0000",alpha=0.4, size=3) +
  theme_classic() + 
  theme(
    axis.title.x = element_blank(), 
    axis.line = element_line(linewidth = 0.8, color="#000000"), 
    axis.ticks = element_line(linewidth = 0.8, color="#000000"), 
    axis.ticks.length=unit(4, "points"), 
    axis.text.x = element_text(size=10, color="#000000", face="bold"), 
    axis.text.y = element_text(size=10, color="#000000", face="bold.italic")
  )

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/SupplFig_RNA_Complete_Pseudobulk/", 
    "ALS_DESeq_VolcanoPlot", 
    ".pdf"
  ),
  width = 140, 
  height = 180, 
  units="mm"
) 


get_nDEGs(DDS_12SVs_Case, name="Case_ALS_vs_HC")
get_nDEGs(DDS_12SVs_Rand, name="Rand_R2_vs_R1")



    ## 4.2 ALSFTD DESeq --------------------------------------------------------

ALSFTD_DESeq <- as.data.frame(
  results(DDS_12SVs_Case, name="Case_ALS_FTD_vs_HC", alpha=0.05, lfcThreshold = 0)
)
ALSFTD_DESeq <- ALSFTD_DESeq[!is.na(ALSFTD_DESeq$padj),]

ALSFTD_DESeq$Type = "NA"
ALSFTD_DESeq$Type[ALSFTD_DESeq$log2FoldChange<0 & ALSFTD_DESeq$padj<0.05] <- "Down"
ALSFTD_DESeq$Type[ALSFTD_DESeq$log2FoldChange>0 & ALSFTD_DESeq$padj<0.05] <- "Up"
table(ALSFTD_DESeq$Type)

ggplot() + 
  aes(log2FoldChange, -log10(padj)) + 
  geom_hline(yintercept = c(-log10(0.05)), col="#AAAAAA", linewidth=1) + 
  ggrastr::geom_point_rast(aes(log2FoldChange, -log10(padj)), dat=ALSFTD_DESeq[ALSFTD_DESeq$Type=="NA",],alpha=0.2, col="#000000", fill="#000000", pch=21, size=3) + 
  ggrastr::geom_point_rast(aes(log2FoldChange, -log10(padj)), dat=ALSFTD_DESeq[ALSFTD_DESeq$Type=="Down",], pch=21, col="#000000", fill="#4444FF",alpha=0.4, size=3) + 
  ggrastr::geom_point_rast(aes(log2FoldChange,-log10(padj)), dat=ALSFTD_DESeq[ALSFTD_DESeq$Type=="Up",], pch=21, col="#000000", fill="#FF0000",alpha=0.4, size=3) +
  scale_x_continuous(limits=c(-3, 3)) + 
  theme_classic() + 
  theme(
    axis.title.x = element_blank(), 
    axis.line = element_line(linewidth = 0.8, color="#000000"), 
    axis.ticks = element_line(linewidth = 0.8, color="#000000"), 
    axis.ticks.length=unit(4, "points"), 
    axis.text.x = element_text(size=10, color="#000000", face="bold"), 
    axis.text.y = element_text(size=10, color="#000000", face="bold.italic")
  )

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/SupplFig_RNA_Complete_Pseudobulk/", 
    "ALSFTD_DESeq_VolcanoPlot", 
    ".pdf"
  ),
  width = 140, 
  height = 180, 
  units="mm"
) 


get_nDEGs(DDS_12SVs_Case, name="Case_ALS_FTD_vs_HC")
get_nDEGs(DDS_12SVs_Rand, name="Rand_R3_vs_R1") 




  ### 5.0 Overlap of ALS with ALS-FTD ------------------------------------------

ALS_Genes <- rownames(
  get_sign_genes_df(
    results(
      DDS_12SVs_Case, 
      name = "Case_ALS_vs_HC", 
      alpha=0.05, 
      lfcThreshold = 0
    ), 
    fdr=0.05
  )
)

ALSFTD_Genes <- rownames(
  get_sign_genes_df(
    results(
      DDS_12SVs_Case, 
      name = "Case_ALS_FTD_vs_HC", 
      alpha=0.05, 
      lfcThreshold = 0
    ), 
    fdr=0.05
  )
)

table(ALS_Genes %in% ALSFTD_Genes)
table(ALSFTD_Genes %in% ALS_Genes)

tmp <- as.data.frame(results(DDS_12SVs_Case, name="Case_ALS_vs_HC", alpha=0.05, lfcThreshold = 0))   
tmp <- tmp[!is.na(tmp$padj),]
ALS_All_Genes <- rownames(tmp)
table(is.na(tmp$padj))

tmp <- as.data.frame(results(DDS_12SVs_Case, name="Case_ALS_FTD_vs_HC", alpha=0.05, lfcThreshold = 0))  
tmp <- tmp[!is.na(tmp$padj),]
ALSFTD_All_Genes <- rownames(tmp)
table(is.na(tmp$padj))

rm(tmp)
AllGenes <- unique(c(ALS_All_Genes, ALSFTD_All_Genes))
length(AllGenes)


min(ALSFTD_DESeq$baseMean[ALSFTD_DESeq$padj<0.05])

tmp <- as.data.frame(results(DDS_12SVs_Case, name="Case_ALS_vs_HC", alpha=0.05, lfcThreshold = 0))   
tmp <- tmp[tmp$baseMean>=min(ALS_DESeq$baseMean[ALS_DESeq$padj<0.05]),]
ALS_All_EquallyExpressed_Genes <- rownames(tmp)
table(is.na(tmp$padj))

tmp <- as.data.frame(results(DDS_12SVs_Case, name="Case_ALS_FTD_vs_HC", alpha=0.05, lfcThreshold = 0))   
tmp <- tmp[tmp$baseMean>=min(ALSFTD_DESeq$baseMean[ALSFTD_DESeq$padj<0.05]),]
ALSFTD_All_EquallyExpressed_Genes <- rownames(tmp)
table(is.na(tmp$padj))

AllEquallyExpressedGenes <- unique(c(ALS_All_EquallyExpressed_Genes, ALSFTD_All_EquallyExpressed_Genes))

table(
  AllGenes %in% 
  AllEquallyExpressedGenes
) 

table(
  AllEquallyExpressedGenes %in% 
    AllGenes
)


fisher.test(
  matrix(
    c(
      180, 
      487, 
      3647, 
      25067-180-487-3647
    ), 
    2, 2
  )
)
rm(AllGenes, ALS_All_Genes, ALSFTD_All_Genes)


  

  ### 6.0 ALS, ALS-FTD Correlation plot ----------------------------------------

All_Genes <- unique(c(ALS_Genes, ALSFTD_Genes))

All_Genes <- data.frame(
  Gene = All_Genes, 
  ALS_L2FC = ALS_DESeq$log2FoldChange[match(All_Genes, rownames(ALS_DESeq))], 
  ALSFTD_L2FC = ALSFTD_DESeq$log2FoldChange[match(All_Genes, rownames(ALSFTD_DESeq))], 
  ALS_PVAL = ALS_DESeq$pvalue[match(All_Genes, rownames(ALS_DESeq))], 
  ALSFTD_PVAL = ALSFTD_DESeq$pvalue[match(All_Genes, rownames(ALSFTD_DESeq))], 
  ALS_QVAL = ALS_DESeq$padj[match(All_Genes, rownames(ALS_DESeq))], 
  ALSFTD_QVAL = ALSFTD_DESeq$padj[match(All_Genes, rownames(ALSFTD_DESeq))], 
  ALS_STAT = ALS_DESeq$stat[match(All_Genes, rownames(ALS_DESeq))], 
  ALSFTD_STAT = ALSFTD_DESeq$stat[match(All_Genes, rownames(ALSFTD_DESeq))]
) 
All_Genes$SignIn <- "None"
All_Genes$SignIn[All_Genes$ALS_QVAL<0.05] <- "ALS"
All_Genes$SignIn[All_Genes$ALSFTD_QVAL<0.05] <- "ALSFTD"
All_Genes$SignIn[All_Genes$ALS_QVAL<0.05 & All_Genes$ALSFTD_QVAL<0.05] <- "BOTH"
table(All_Genes$SignIn)

ALS_Genes[ALS_Genes %in% ALSFTD_Genes]

ggplot() + 
  geom_hline(yintercept = 0, col="#AAAAAA") + 
  geom_vline(xintercept = 0, col="#AAAAAA") + 
  geom_abline(slope = 1, intercept = 0, col="#FF000055") + 
  ggrastr::geom_point_rast(aes(ALS_L2FC, ALSFTD_L2FC), data=All_Genes[All_Genes$SignIn=="ALS",], size=1, pch=21, col="#000000", fill="#891D14", alpha=0.5) + 
  ggrastr::geom_point_rast(aes(ALS_L2FC, ALSFTD_L2FC), data=All_Genes[All_Genes$SignIn=="ALSFTD",], size=1, pch=21, col="#000000", fill="#1E709B", alpha=0.5) + 
  ggrastr::geom_point_rast(aes(ALS_L2FC, ALSFTD_L2FC), data=All_Genes[All_Genes$SignIn=="BOTH",], size=2, pch=21, col="#000000", fill="#FFF98C") + 
  scale_x_continuous(limits=c(-3,3)) + 
  scale_y_continuous(limits=c(-3,3)) + 
  xlab("\nALS vs. HC [ fold change, log2 ]") +  
  ylab("ALS-FTD vs. HC [ fold change, log2 ]\n") + 
  theme_classic() + 
  theme(
    axis.title = element_text(size=12, color="#000000", face="bold.italic"), 
    axis.line = element_line(linewidth = 0.8, color="#000000"), 
    axis.ticks = element_line(linewidth = 0.8, color="#000000"), 
    axis.ticks.length=unit(4, "points"), 
    axis.text.x = element_text(size=10, color="#000000", face="bold.italic"), 
    axis.text.y = element_text(size=10, color="#000000", face="bold.italic")
  )

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/SupplFig_RNA_Complete_Pseudobulk/", 
    "Correlation_Plot_ALS_ALSFTD_DEGs", 
    ".pdf"
  ),
  width = 100, 
  height = 100, 
  units="mm"
) 

cor.test(
  All_Genes$ALS_L2FC, 
  All_Genes$ALSFTD_L2FC
)

cor.test(
  All_Genes$ALS_L2FC[All_Genes$SignIn=="BOTH"], 
  All_Genes$ALSFTD_L2FC[All_Genes$SignIn=="BOTH"]
)

table(
  All_Genes$ALS_L2FC[All_Genes$SignIn=="BOTH"]<0, 
  All_Genes$ALSFTD_L2FC[All_Genes$SignIn=="BOTH"]<0 
)
173/180




  ### 7.0 Plot GENETYPE Annotation ---------------------------------------------

Genetypes <- data.frame(
  GENETYPE=unique(GEX_Features$GENETYPE)
)


All_Features <- data.frame(
  ID_10X = rownames(DDS_12SVs_Case)
) 
All_Features <- inner_join(
  All_Features, 
  GEX_Features, 
  by=c("ID_10X"="ID_10X")
)

Genetypes <- full_join(
  full_join(
    Genetypes, 
    as.data.frame(prop.table(table(All_Features$GENETYPE))), 
    by=c("GENETYPE"="Var1")
  ) %>% rename(All_Features_Pct=Freq), 
  as.data.frame(table(All_Features$GENETYPE)), 
  by=c("GENETYPE"="Var1")
) %>% rename(All_Features_N=Freq)
rm(All_Features)



All_Expressed <- data.frame(
  ID_10X = AllEquallyExpressedGenes
) 
All_Expressed <- inner_join(
  All_Expressed, 
  GEX_Features, 
  by=c("ID_10X"="ID_10X")
)

Genetypes <- full_join(
  full_join(
    Genetypes, 
    as.data.frame(prop.table(table(All_Expressed$GENETYPE))), 
    by=c("GENETYPE"="Var1")
  ) %>% rename(Expressed_Features_Pct=Freq), 
  as.data.frame(table(All_Expressed$GENETYPE)), 
  by=c("GENETYPE"="Var1")
) %>% rename(Expressed_Features_N=Freq)
rm(All_Expressed)




All_DEGs <- data.frame(
  ID_10X = All_Genes$Gene
) 
All_DEGs <- inner_join(
  All_DEGs, 
  GEX_Features, 
  by=c("ID_10X"="ID_10X")
)

Genetypes <- full_join(
  full_join(
    Genetypes, 
    as.data.frame(prop.table(table(All_DEGs$GENETYPE))), 
    by=c("GENETYPE"="Var1"), keep = FALSE
  ) %>% rename(All_DEGs_Features_Pct=Freq), 
  as.data.frame(table(All_DEGs$GENETYPE)), 
  by=c("GENETYPE"="Var1")
) %>% rename(All_DEGs_Features_N=Freq)
rm(All_DEGs)



ALS_DEGs <- data.frame(
  ID_10X = All_Genes$Gene[All_Genes$SignIn %in% c("ALS", "BOTH")]
) 
ALS_DEGs <- inner_join(
  ALS_DEGs, 
  GEX_Features, 
  by=c("ID_10X"="ID_10X")
)

Genetypes <- full_join(
  full_join(
    Genetypes, 
    as.data.frame(prop.table(table(ALS_DEGs$GENETYPE))), 
    by=c("GENETYPE"="Var1")
  ) %>% rename(ALS_DEGs_Features_Pct=Freq), 
  as.data.frame(table(ALS_DEGs$GENETYPE)), 
  by=c("GENETYPE"="Var1")
) %>% rename(ALS_DEGs_Features_N=Freq)
rm(ALS_DEGs)



ALSFTD_DEGs <- data.frame(
  ID_10X = All_Genes$Gene[All_Genes$SignIn %in% c("ALSFTD", "BOTH")]
) 
ALSFTD_DEGs <- inner_join(
  ALSFTD_DEGs, 
  GEX_Features, 
  by=c("ID_10X"="ID_10X")
)

Genetypes <- full_join(
  full_join(
    Genetypes, 
    as.data.frame(prop.table(table(ALSFTD_DEGs$GENETYPE))), 
    by=c("GENETYPE"="Var1")
  ) %>% rename(ALSFTD_DEGs_Features_Pct=Freq), 
  as.data.frame(table(ALSFTD_DEGs$GENETYPE)), 
  by=c("GENETYPE"="Var1")
) %>% rename(ALSFTD_DEGs_Features_N=Freq)
rm(ALSFTD_DEGs)



Common_DEGs <- data.frame(
  ID_10X = All_Genes$Gene[All_Genes$SignIn %in% c("BOTH")]
) 
Common_DEGs <- inner_join(
  Common_DEGs, 
  GEX_Features, 
  by=c("ID_10X"="ID_10X")
)

Genetypes <- full_join(
  full_join(
    Genetypes, 
    as.data.frame(prop.table(table(Common_DEGs$GENETYPE))), 
    by=c("GENETYPE"="Var1")
  ) %>% rename(Common_DEGs_Features_Pct=Freq), 
  as.data.frame(table(Common_DEGs$GENETYPE)), 
  by=c("GENETYPE"="Var1")
) %>% rename(Common_DEGs_Features_N=Freq)
rm(Common_DEGs)

Genetypes[is.na(Genetypes)] <- 0

Genetypes2 <- inner_join(
  
  Genetypes %>% 
    dplyr::select(GENETYPE, ends_with("Features_N")) %>% 
    pivot_longer(!GENETYPE, names_transform = function(x){return(str_replace_all(x, "_Features_N", ""))}, names_to = "Group", values_to = "N"), 
  
  Genetypes %>% 
    dplyr::select(GENETYPE, ends_with("Features_Pct")) %>% 
    pivot_longer(!GENETYPE, names_transform = function(x){return(str_replace_all(x, "_Features_Pct", ""))}, names_to = "Group", values_to = "Pct"), 

  by=c("Group"="Group", "GENETYPE"="GENETYPE")
)


Genetypes %>% 
  dplyr::select(GENETYPE, ends_with("Features_N")) %>% 
  pivot_longer(!GENETYPE, names_transform = function(x){return(str_replace_all(x, "_Features_N", ""))}, names_to = "Group", values_to = "N")

ggplot(Genetypes2 %>% filter(Group %in% c("All", "Expressed", "All_DEGs"))) + 
  aes(factor(Group, levels=unique(Genetypes2$Group)[c(1,2,3)]), N, fill=as.character(GENETYPE)) + 
  geom_bar(stat="identity", position="fill", col="#000000", lwd=0.5) + 
  scale_y_continuous(expand=c(0,0,0,0)) + 
  scale_fill_manual(values=c(c("#0072B2", "#56B4E9", "#009E73",
                               "#F0E442", "#E69F00", "#D55E00", "#CC79A7"))) + 
  theme_classic() + 
  theme(
    axis.title.x = element_blank(), 
    axis.line = element_line(linewidth = 0.8, color="#000000"), 
    axis.ticks = element_line(linewidth = 0.8, color="#000000"), 
    axis.ticks.length=unit(4, "points"), 
    axis.text.x = element_text(size=10, color="#000000", angle=45, hjust=1, face="bold"), 
    axis.text.y = element_text(size=10, color="#000000", face="bold.italic")
  )

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/SupplFig_RNA_Complete_Pseudobulk/", 
    "Genetypes_Barplot_All_DEGs", 
    ".pdf"
  ),
  width = 120, 
  height = 120, 
  units="mm"
)  


rm(Genetypes, Genetypes2)




  ### 8.0 ALS, ALS-FTD DEGs GSEA -----------------------------------------------

Universe <- data.frame(
  ID_10X = AllEquallyExpressedGenes
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

All_Genes$ID_10X <- All_Genes$Gene
All_Genes$SYMBOL <- GEX_Features$SYMBOL[match(All_Genes$ID_10X, GEX_Features$ID_10X)]
all(All_Genes$ID_10X %in% Universe$ID_10X)
All_Genes$ENTREZID <- Universe$ENTREZ[match(All_Genes$ID_10X, Universe$ID_10X)]
All_Genes$Accession <- GEX_Features$Accession[match(All_Genes$ID_10X, GEX_Features$ID_10X)]

universe = Universe$ENTREZ
universe = universe[!is.na(universe)]



    ## 8.1 ALS DEGs ------------------------------------------------------------


      # Enrichment KEGG 

table(All_Genes$SignIn, All_Genes$Accession)

genelist <- All_Genes$ENTREZID[
  All_Genes$SignIn %in% c("ALS", "BOTH") 
  ]
genelist <- genelist[!is.na(genelist)]


Enrich_Kegg_ALS_DEGs <- enrichKEGG(
  gene =genelist, 
  organism = "hsa", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05, 
  universe = universe
)

rm(genelist)


      # GSEA KEGG 

genelist <- All_Genes$ALS_STAT[
  All_Genes$SignIn %in% c("ALS", "BOTH") & 
    !is.na(All_Genes$ENTREZID)
]
names(genelist) <- All_Genes$ENTREZID[
  All_Genes$SignIn %in% c("ALS", "BOTH") & 
    !is.na(All_Genes$ENTREZID)
]
genelist <- genelist[!is.na(names(genelist))]
genelist <- sort(genelist, decreasing = TRUE)

GSEA_Kegg_ALS_DEGs_BySTAT <- gseKEGG(  
  geneList=genelist, 
  organism = "hsa", 
  eps=0
)

rm(genelist)




    ## 8.2 ALSFTD DEGs ---------------------------------------------------------


      # Enrichment KEGG 

table(All_Genes$SignIn, All_Genes$Accession)

genelist <- All_Genes$ENTREZID[
  All_Genes$SignIn %in% c("ALSFTD", "BOTH") 
]
genelist <- genelist[!is.na(genelist)]


Enrich_Kegg_ALSFTD_DEGs <- enrichKEGG(
  gene =genelist, 
  organism = "hsa", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05, 
  universe = universe
)

rm(genelist)


      # GSEA KEGG 

genelist <- All_Genes$ALSFTD_STAT[
  All_Genes$SignIn %in% c("ALSFTD", "BOTH") & 
    !is.na(All_Genes$ENTREZID)
]
names(genelist) <- All_Genes$ENTREZID[
  All_Genes$SignIn %in% c("ALSFTD", "BOTH") & 
    !is.na(All_Genes$ENTREZID)
]
genelist <- genelist[!is.na(names(genelist))]
genelist <- sort(genelist, decreasing = TRUE)

GSEA_Kegg_ALSFTD_DEGs_BySTAT <- gseKEGG(  
  geneList=genelist, 
  organism = "hsa", 
  eps=0
)

rm(genelist)



    ## 8.3 ALS + ALSFTD DEGs --------------------------------------------------- 


      # Enrichment KEGG 

table(All_Genes$SignIn, All_Genes$Accession)

genelist <- All_Genes$ENTREZID
genelist <- genelist[!is.na(genelist)]


Enrich_Kegg_All_DEGs <- enrichKEGG(
  gene =genelist, 
  organism = "hsa", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05, 
  universe = universe
)

rm(genelist)


      # GSEA KEGG 

genelist <- rowMeans(cbind(
    All_Genes$ALS_STAT, 
    All_Genes$ALSFTD_STAT
  )
)
names(genelist) <- All_Genes$ENTREZID

genelist <- genelist[!is.na(names(genelist))]
genelist <- sort(genelist, decreasing = TRUE)

GSEA_Kegg_All_DEGs_BySTAT <- gseKEGG(  
  geneList=genelist, 
  organism = "hsa", 
  eps=0
)

rm(genelist)



    ## 8.4 ALS & ALSFTD DEGs --------------------------------------------------- 


      # Enrichment KEGG 

table(All_Genes$SignIn, All_Genes$Accession)

genelist <- All_Genes$ENTREZID[
  All_Genes$SignIn %in% c("BOTH") 
]
genelist <- genelist[!is.na(genelist)]


Enrich_Kegg_Common_DEGs <- enrichKEGG(
  gene =genelist, 
  organism = "hsa", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05, 
  universe = universe
)

rm(genelist)


      # GSEA KEGG 

genelist <- rowMeans(
  cbind(
    All_Genes$ALS_STAT, 
    All_Genes$ALSFTD_STAT
  )
)
names(genelist) <- All_Genes$ENTREZID
genelist <- genelist[
  All_Genes$SignIn %in% c("BOTH") 
]

genelist <- genelist[!is.na(names(genelist))]
genelist <- sort(genelist, decreasing = TRUE)

GSEA_Kegg_Common_DEGs_BySTAT <- gseKEGG(  
  geneList=genelist, 
  organism = "hsa", 
  eps=0
)

rm(genelist)



    ## 8.5 GSEA plots ----------------------------------------------------------

plot_GSEA_as_Tree(
  GSEA_Kegg_All_DEGs_BySTAT, 
   name=paste0(
      "../Data/Visualization/Figures/SupplFig_RNA_Complete_Pseudobulk/", 
      "GSEA_All_DEGs_Tree", 
      ".pdf"
      ), 
  save = TRUE, width=3, height=6
  )
