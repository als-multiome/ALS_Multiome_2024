

### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 


library(qs) 
library(tidyverse) 
library(data.table)
library(ggrastr)
library(GenomicRanges)
library(ChIPseeker)
library(DESeq2)
library(BiocParallel)
library(sva)
library(UpSetR)  



qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)




  ### 1.0 Load data ------------------------------------------------------------ 

M0_ATAC <- qread(
  paste0(
    "../Data/SeuratObjects/", 
    "M0_ATAC", 
    ".qrds"
  ), 
  nthr=nthr
)

M0_RNA <- qread(
  paste0(
    "../Data/SeuratObjects/", 
    "M0_RNA", 
    ".qrds"
  ), 
  nthr=nthr
)


ATAC_Links <- qread(
  paste0(
    "../Data/Annotations/Genomic_Regions/", 
    "M0_ATAC_Peak_Links", 
    ".qrds"
  ), 
  nthr=nthr
)

Peaks <- qread(
  paste0(
    "../Data/Annotations/Genomic_Regions/", 
    "peakAnno_df", 
    ".qrds"
  ), 
  nthr=nthr
)



DDS_8SVs_list <- qread(
  "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/SVA/1_15_SVs/DDS_8SVs_list.qrds", nthr=nthr
)

DESeq_Results_8SVs_ALS <- qread(
  "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/SVA/1_15_SVs/DESeq_Results_8SVs_ALS.qrds", nthr=nthr
)

DESeq_Results_8SVs_ALSFTD <- qread(
  "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/SVA/1_15_SVs/DESeq_Results_8SVs_ALSFTD.qrds", nthr=nthr
)

DESeq_Results_8SVs_Index <- qread(
  "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/SVA/1_15_SVs/DESeq_Results_8SVs_Index.qrds", nthr=nthr
)



L2FC_Shrink_Results_8SVs_ALS <- qread(
  paste0(
    "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/SVA/1_15_SVs/", 
    "L2FC_Shrink_Results_8SVs_ALS", 
    ".qrds"
  ), 
  nthr=nthr
)


L2FC_Shrink_Results_8SVs_ALSFTD <- qread(
  paste0(
    "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/SVA/1_15_SVs/", 
    "L2FC_Shrink_Results_8SVs_ALSFTD", 
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

RNA_L2FC_Shrink_Results_12SVs_Index <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/12_SVs/", 
    "DESeq_Results_12SVs_Index", 
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

ColDict_Modality <- setNames(
  object = readxl::read_xlsx(paste0(
    "../Data/Visualization/", 
    "ALS_Brain_Multiome_ColDicts", 
    ".xlsx"
  ), 
  sheet = "Modality"
  )$Color, 
  nm = readxl::read_xlsx(paste0(
    "../Data/Visualization/", 
    "ALS_Brain_Multiome_ColDicts", 
    ".xlsx"
  ), 
  sheet = "Modality"
  )$Modality
)



Universe_AllCase_ALS_List <- qread(
  paste0(
    "../Data/Annotations/Genes/", 
    "Universe_AllCase_ALS_List", 
    ".qrds"
  ), 
  nthr=nthr
)

Universe_AllCase_ALSFTD_List <- qread(
  paste0(
    "../Data/Annotations/Genes/", 
    "Universe_AllCase_ALSFTD_List", 
    ".qrds"
  ), 
  nthr=nthr
)



    ## 1.1 Define AUX functions ------------------------------------------------

Get_DARs <- function(x, alpha=0.05, sign="All"){
  x <- x[!is.na(x$padj),] 
  return(
    if(sign=="All") rownames(x)[x$padj<alpha] else
    if(sign=="Up") rownames(x)[x$padj<alpha & x$log2FoldChange>0] else
    if(sign=="Down") rownames(x)[x$padj<alpha & x$log2FoldChange<0]
  )
}




  ### 2.0 Barplot peaks and genes with link ------------------------------------

data.frame(
  Cat=factor(c("Peaks", "Peaks", "Genes", "Genes"), levels=c("Peaks", "Genes")), 
  Group=factor(c("Link", "NoLink", "GeneLink", "NoLink"), levels=c("NoLink", "Link", "GeneLink")), 
  N=c(
    length(unique(ATAC_Links$peak)), 
    nrow(M0_ATAC) - length(unique(ATAC_Links$peak)), 
    length(unique(ATAC_Links$gene)), 
    nrow(M0_RNA) - length(unique(ATAC_Links$gene))
  )
) %>% 
  ggplot() + 
  aes(x=Cat, y=N, fill=Group) + 
  geom_col(col="#000000", position="fill", width=0.75) + 
  scale_y_continuous(expand=c(0,0,0.05,0)) + 
  scale_fill_manual(values=c("NoLink"="#E8E6DF55", "Link"="#9E7BB0", "GeneLink"="#E3C491")) + 
  theme_classic() 
  
ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/Fig_ATAC/", 
    "Peaks_Genes_Links_Barplot", 
    ".pdf"
  ),
  width = 100, 
  height = 100, 
  units="mm"
)  

table(rownames(M0_ATAC) %in% ATAC_Links$peak)
table(rownames(M0_RNA) %in% ATAC_Links$gene)




  ### 3.0 Barplot Peaks annotation ---------------------------------------------

data.frame(
  table(
    Peaks$Annotation_L2_Plain, 
    Peaks$ID2 %in% ATAC_Links$peak
  )
) %>% 
  ggplot() + 
  aes(
    x=plyr::revalue(factor(Var2), c("TRUE"="Link Peaks", "FALSE"="All Peaks")), 
    y=Freq, 
    fill=factor(Var1, levels=c("Intergenic", "Promoter", "Exon", "Intron", "Downstream"))
  ) + 
  geom_col(position="fill", col="#000000") + 
  scale_y_continuous(expand=c(0,0,0.05,0)) + 
  scale_fill_manual(values=c("#FFA94D", "#0d0887", "#7e03a8", "#cc4778", "#f89540")) + 
  theme_classic() + 
  theme(
    legend.title = element_blank(), 
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size=14, color="#000000", face="bold.italic"), 
    axis.text = element_text(color="#000000", face="bold.italic", size=12), 
    legend.text = element_text(size=12, color="#000000", face="bold.italic") 
    
  )
  
ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/Fig_ATAC/", 
    "LinksPeaks_Annotations_Barplot", 
    ".pdf"
  ),
  width = 100, 
  height = 100, 
  units="mm"
)  




  ### 4.0 Barplots # DARs ------------------------------------------------------



    ## 4.1 DAR #N summary ------------------------------------------------------

ind <- which(
  DESeq_Results_8SVs_Index$CellTypeLevel=="WNN_L25" & 
    DESeq_Results_8SVs_Index$CellType %in% c(
      "Exc_RORB", 
      "Exc_LINC00507", 
      "Oligodendrocytes"
    ) & 
    DESeq_Results_8SVs_Index$Comparison=="All_Cases"
)
names(ind) <- DESeq_Results_8SVs_Index$CellType[ind]


Results_WNN_L25 <- data.frame(
  CellType=names(ind), 
  Up_q_0.05_ALS = DESeq_Results_8SVs_Index$Up_q_0.05_ALS[ind], 
  Down_q_0.05_ALS = DESeq_Results_8SVs_Index$Down_q_0.05_ALS[ind], 
  Up_q_0.05_ALSFTD = DESeq_Results_8SVs_Index$Up_q_0.05_ALSFTD[ind], 
  Down_q_0.05_ALSFTD = DESeq_Results_8SVs_Index$Down_q_0.05_ALSFTD[ind]
)



    ## 4.2 ALS 8 SVs Barplots --------------------------------------------------

ggplot(
  Results_WNN_L25
) + 
  aes(CellType, Up_q_0.05_ALS, fill=CellType) + 
  geom_col(col="#000000", fill="#FF5744", width=0.8) + 
  scale_y_continuous(expand=c(0,0,0.2,0), limits=c(0,40)) + 
  theme_classic() + 
  theme(
    axis.text.x = element_text(angle=43, hjust=1.0), 
    panel.grid.major.y = element_line(color="#DDDDDD")
  ) 

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig5/", 
    "DA_ATAC_LinkPeaks_WNN_L25_nDEG_Up_ALS",
    ".pdf"
  ), 
  width = 1.92, 
  height = 2.24,  
  units="in"
)


ggplot(
  Results_WNN_L25
) + 
  aes(CellType, Down_q_0.05_ALS, fill=CellType) + 
  geom_col(fill="#3354E2", col="#000000", width=0.8) + 
  scale_y_reverse(expand=c(0.1,0,0,0), limits=c(40,0)) + 
  theme_classic() + 
  theme(
    axis.text.x = element_text(angle=43, hjust=1.0), 
    panel.grid.major.y = element_line(color="#DDDDDD")
  ) 

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig5/", 
    "DA_ATAC_LinkPeaks_WNN_L25_nDEG_Down_ALS",
    ".pdf"
  ), 
  width = 1.92, 
  height = 2.24,  
  units="in"
)



    ## 4.3 ALSFTD 8 SVs Barplots -----------------------------------------------

ggplot(
  Results_WNN_L25
) + 
  aes(CellType, Up_q_0.05_ALSFTD, fill=CellType) + 
  geom_col(col="#000000", fill="#FF5744", width=0.8) + 
  scale_y_continuous(expand=c(0,0,0.2,0), limits=c(0,2000)) + 
  theme_classic() + 
  theme(
    axis.text.x = element_text(angle=43, hjust=1.0), 
    panel.grid.major.y = element_line(color="#DDDDDD")
  ) 

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig5/", 
    "DA_ATAC_LinkPeaks_WNN_L25_nDEG_Up_ALSFTD",
    ".pdf"
  ), 
  width = 1.92, 
  height = 2.24,  
  units="in"
)


ggplot(
  Results_WNN_L25
) + 
  aes(CellType, Down_q_0.05_ALSFTD, fill=CellType) + 
  geom_col(fill="#3354E2", col="#000000", width=0.8) + 
  scale_y_reverse(expand=c(0.1,0,0,0), limits=c(2000,0)) + 
  theme_classic() + 
  theme(
    axis.text.x = element_text(angle=43, hjust=1.0), 
    panel.grid.major.y = element_line(color="#DDDDDD")
  ) 

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig5/", 
    "DA_ATAC_LinkPeaks_WNN_L25_nDEG_Down_ALSFTD",
    ".pdf"
  ), 
  width = 1.92, 
  height = 2.24,  
  units="in"
)




  ### 5.0 Upset Plots ----------------------------------------------------------



    ## 5.1 ALS DARs WNN_L25 ----------------------------------------------------

DARs_ALS_Exc_LINC00507 <- Get_DARs(DESeq_Results_8SVs_ALS[[ind["Exc_LINC00507"]]])
DARs_ALS_Exc_RORB <- Get_DARs(DESeq_Results_8SVs_ALS[[ind["Exc_RORB"]]])
DARs_ALS_OLG <- Get_DARs(DESeq_Results_8SVs_ALS[[ind["Oligodendrocytes"]]])

ALS_DARs <- list(
  "Exc LINC00507" = DARs_ALS_Exc_LINC00507, 
  "Exc RORB" = DARs_ALS_Exc_RORB, 
  "Oligodendrocytes" = DARs_ALS_OLG
)


upset(
  fromList(ALS_DARs), 
  order.by = c("degree"), 
  sets = c("Exc LINC00507", "Exc RORB", "Oligodendrocytes"), 
  keep.order = TRUE, 
  group.by = "degree", 
  empty.intersections = TRUE, 
  point.size = 4, 
  matrix.dot.alpha = 0.8, 
  text.scale = 2, 
  line.size = 1, 
  queries = list(
    list(
      query = intersects, 
      params = list(
        "Exc RORB"
      ), 
      color = ColDict_WNN_L25["Exc_RORB"], 
      active = T
    ), 
    list(
      query = intersects, 
      params = list(
        "Exc LINC00507"
      ), 
      color = ColDict_WNN_L25["Exc_LINC00507"], 
      active = T
    ), 
    list(
      query = intersects, 
      params = list(
        c("Oligodendrocytes")
      ), 
      color = ColDict_WNN_L25["Oligodendrocytes"], 
      active = T
    ) 
  )
)

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/Fig_ATAC/", 
    "DARS_WNN_L25_UpsetPlot_ALS", 
    ".pdf"
  ),
  width = 207, 
  height = 153, 
  units="mm"
)


rm(DARs_ALS_Exc_LINC00507, DARs_ALS_Exc_RORB, DARs_ALS_OLG, ALS_DARs) 



    ## 5.2 ALSFTD DARs WNN_L25 -------------------------------------------------

DARs_ALSFTD_Exc_LINC00507 <- Get_DARs(DESeq_Results_8SVs_ALSFTD[[ind["Exc_LINC00507"]]])
DARs_ALSFTD_Exc_RORB <- Get_DARs(DESeq_Results_8SVs_ALSFTD[[ind["Exc_RORB"]]])
DARs_ALSFTD_OLG <- Get_DARs(DESeq_Results_8SVs_ALSFTD[[ind["Oligodendrocytes"]]])

ALSFTD_DARs <- list(
  "Exc LINC00507" = DARs_ALSFTD_Exc_LINC00507, 
  "Exc RORB" = DARs_ALSFTD_Exc_RORB, 
  "OLG" = DARs_ALSFTD_OLG
)


upset(
  fromList(ALSFTD_DARs), 
  order.by = c("degree"), 
  sets = c("Exc LINC00507", "Exc RORB", "OLG"), 
  keep.order = TRUE, 
  group.by = "degree", 
  empty.intersections = TRUE, 
  point.size = 4, 
  matrix.dot.alpha = 0.8, 
  text.scale = 2, 
  line.size = 1, 
  queries = list(
    list(
      query = intersects, 
      params = list(
        "Exc RORB"
      ), 
      color = ColDict_WNN_L25["Exc_RORB"], 
      active = T
    ), 
    list(
      query = intersects, 
      params = list(
        "Exc LINC00507"
      ), 
      color = ColDict_WNN_L25["Exc_LINC00507"], 
      active = T
    ), 
    list(
      query = intersects, 
      params = list(
        c("OLG")
      ), 
      color = ColDict_WNN_L25["Oligodendrocytes"], 
      active = T
    ) 
  )
)

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/Fig_ATAC/", 
    "DARS_WNN_L25_UpsetPlot_ALSFTD", 
    ".pdf"
  ),
  width = 207, 
  height = 153, 
  units="mm"
)

rm(DARs_ALSFTD_Exc_LINC00507, DARs_ALSFTD_Exc_RORB, DARs_ALSFTD_OLG, ALSFTD_DARs) 




  ### 6.0 DAR WNN_L25 CellType correlation plots -------------------------------



    ## 6.1 Define plotting function --------------------------------------------

plot_cor_WNN_L25_CellTypes <- function(
    ResultsList, 
    ResultsListIndex, 
    CellTypeLevel = "WNN_L25", 
    CellType1, 
    CellType2, 
    Comparison="Case", 
    alpha=0.05, 
    ColDict=NULL
    ){
  
  if(Comparison=="Case"){
    
    ind1 <- which(
      ResultsListIndex$CellTypeLevel == CellTypeLevel & 
        ResultsListIndex$CellType == CellType1 & 
          ResultsListIndex$Comparison=="All_Cases"
    ) 
    
    ind2 <- which(
      ResultsListIndex$CellTypeLevel == CellTypeLevel & 
        ResultsListIndex$CellType == CellType2 & 
        ResultsListIndex$Comparison=="All_Cases"
    )
    
    DARs_CellType1 <- Get_DARs(ResultsList[[ind1]], alpha = alpha)
    DARs_CellType2 <- Get_DARs(ResultsList[[ind2]], alpha = alpha)
    
    df <- data.frame(
      Peak = unique(
        c(
          DARs_CellType1, 
          DARs_CellType2
        )
      )
    ) 
    
    df$L2FC_CellType1 <- ResultsList[[ind1]]$log2FoldChange[match(
      df$Peak, 
      rownames(ResultsList[[ind1]])
    )]
    
    df$L2FC_CellType2 <- ResultsList[[ind2]]$log2FoldChange[match(
      df$Peak, 
      rownames(ResultsList[[ind2]])
    )] 
    
    df$PAdj_CellType1 <- ResultsList[[ind1]]$padj[match(
      df$Peak, 
      rownames(ResultsList[[ind1]])
    )]
    
    df$PAdj_CellType2 <- ResultsList[[ind2]]$padj[match(
      df$Peak, 
      rownames(ResultsList[[ind2]])
    )]
    
  } else {
    
    if(Comparison=="Rand"){
      
      ind1 <- which(
        ResultsListIndex$CellTypeLevel == CellTypeLevel & 
          ResultsListIndex$CellType == CellType1 & 
          ResultsListIndex$Comparison=="All_Cases"
      )  
      
      ind2 <- which(
        ResultsListIndex$CellTypeLevel == CellTypeLevel & 
          ResultsListIndex$CellType == CellType2 & 
          ResultsListIndex$Comparison=="All_Cases"
      )
      
      DARs_CellType1 <- Get_DARs(ResultsList[[ind1]], alpha = alpha)
      DARs_CellType2 <- Get_DARs(ResultsList[[ind2]], alpha = alpha)
      
      ind1 <- which(
        ResultsListIndex$CellTypeLevel == CellTypeLevel & 
          ResultsListIndex$CellType == CellType1 & 
          ResultsListIndex$Comparison=="Rand"
      )  
      
      ind2 <- which(
        ResultsListIndex$CellTypeLevel == CellTypeLevel & 
          ResultsListIndex$CellType == CellType2 & 
          ResultsListIndex$Comparison=="Rand"
      )
      
      res1 <- ResultsList[[ind1]] 
      res1 <- res1[order(res1$pvalue),] 
      DARs_CellType1 <- rownames(res1)[1:length(DARs_CellType1)]
      
      res2 <- ResultsList[[ind2]] 
      res2 <- res2[order(res2$pvalue),] 
      DARs_CellType2 <- rownames(res2)[1:length(DARs_CellType2)] 
      
      
      df <- data.frame(
        Peak = unique(
          c(
            DARs_CellType1, 
            DARs_CellType2
          )
        )
      ) 
      
      df$L2FC_CellType1 <- ResultsList[[ind1]]$log2FoldChange[match(
        df$Peak, 
        rownames(ResultsList[[ind1]])
      )]
      
      df$L2FC_CellType2 <- ResultsList[[ind2]]$log2FoldChange[match(
        df$Peak, 
        rownames(ResultsList[[ind2]])
      )] 
      
      df$PAdj_CellType1 <- ResultsList[[ind1]]$padj[match(
        df$Peak, 
        rownames(ResultsList[[ind1]])
      )]
      
      df$PAdj_CellType2 <- ResultsList[[ind2]]$padj[match(
        df$Peak, 
        rownames(ResultsList[[ind2]])
      )]
      
    } 
    
  }
  
  
  df$Sign <- "None"
  df$Sign[df$PAdj_CellType1 < alpha] <- CellType1
  df$Sign[df$PAdj_CellType2 < alpha] <- CellType2 
  df$Sign[df$PAdj_CellType1 < alpha & df$PAdj_CellType2 < alpha] <- "Both"
  df$Sign <- factor(df$Sign)
  lim=max(abs(c(df$L2FC_CellType1, df$L2FC_CellType2)))
 
  cor <- cor.test(df$L2FC_CellType1, df$L2FC_CellType2) 
  
 p1 <- 
    ggplot() + 
    geom_hline(yintercept = 0, col = "#CCCCCC", size = 1) + 
    geom_vline(xintercept = 0, col = "#CCCCCC", size = 1) + 
    (if(cor$p.value<0.05) geom_smooth(aes(L2FC_CellType1, L2FC_CellType2), data = df, method = "lm", col="#238BFA", alpha=0.6, fill= "#CCCCCC44")) +  
    geom_point(aes(x = L2FC_CellType1, y = L2FC_CellType2, fill=Sign), data=df[df$Sign=="None",], size = 4, pch = 21, alpha = 0.6) + 
    geom_point(aes(x = L2FC_CellType1, y = L2FC_CellType2, fill=Sign), data=df[df$Sign==CellType1,], size = 4, pch = 21, alpha = 0.6) + 
    geom_point(aes(x = L2FC_CellType1, y = L2FC_CellType2, fill=Sign), data=df[df$Sign==CellType2,], size = 4, pch = 21, alpha = 0.6) + 
    geom_point(aes(x = L2FC_CellType1, y = L2FC_CellType2, fill=Sign), data=df[df$Sign=="Both",], size = 5, pch = 21, alpha = 1) + 
    scale_x_continuous(limits=c(-lim, lim)) +  
    scale_y_continuous(limits=c(-lim, lim)) +  
    (if(!is.null(ColDict)) {scale_fill_manual(values=ColDict)}) + 
    xlab(paste0("\n", "L2FC ", CellType1)) + 
    ylab(paste0("L2FC ", CellType2, "\n")) + 
    theme_classic() + 
    theme(
      panel.background = element_rect(color="#000000", size=1.2), 
      axis.line = element_blank(), 
      axis.title = element_text(face="bold.italic", color="#000000"), 
      axis.text = element_text(face="bold.italic", color="#000000"), 
      legend.box = element_blank(),  
      legend.box.background = element_blank(), 
      legend.background = element_blank(), 
      legend.key = element_blank(), 
      legend.text = element_text(face="bold.italic", color="#000000"), 
      legend.position = "Null"
    )
   
 message(paste0("#N Peaks: ", nrow(df)))
 message(paste0("Estimate: ", round(cor$estimate, 2)))
 message(paste0("PVal: ", signif(cor$p.value, 2)))
 message(paste0("95CI: ", round(cor$conf.int, 2)[1], "-", round(cor$conf.int, 2)[2])) 

 p1 + annotate(
   geom = "text", 
   x = -0.99*lim, 
   y = 0.95*lim, 
   label = paste0(
   paste0(
     "r = ", 
     round(cor$estimate, 2), 
     " (", 
     round(cor$conf.int, 2)[1], 
     " - ", 
     round(cor$conf.int, 2)[2], 
     ")"
   ), 
   "\n", 
   (if(cor$p.value >= 0.05){paste0("n.s.")}),
   (if(cor$p.value < 0.05 & cor$p.value >= 0.01){paste0("*p<0.05")}), 
   (if(cor$p.value < 0.01 & cor$p.value >= 0.001){paste0("**p<0.01")}), 
   (if(cor$p.value < 0.001 & cor$p.value >= 0.0001){paste0("***p<0.01")}), 
   (if(cor$p.value < 0.0001){paste0("****p<", signif(cor$p.value, 1))})
  ), 
  hjust=0, 
  vjust=1, 
  fontface=4,
  color="#000000"
 )
}

ColDict_Cor <- c(ColDict_WNN_L25, "Both"="#FBFCCC", "None"="#DDDDDD")



    ## 6.2 DAR WNN_L25 ALS L2FC Shrinkage --------------------------------------

plot_cor_WNN_L25_CellTypes(
  L2FC_Shrink_Results_8SVs_ALS, 
  DESeq_Results_8SVs_Index, 
  CellType1 = "Exc_LINC00507", 
  CellType2 = "Exc_RORB", 
  Comparison="Case", 
  alpha=0.05, 
  ColDict = ColDict_Cor
)

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/Fig_ATAC/", 
    "LinkPeaks_Correlation_ALS_ExcLINC00507_ExcRORB", 
    ".pdf"
  ),
  width = 100, 
  height = 100, 
  units="mm"
)  


plot_cor_WNN_L25_CellTypes(
  L2FC_Shrink_Results_8SVs_ALS, 
  DESeq_Results_8SVs_Index, 
  CellType1 = "Exc_LINC00507", 
  CellType2 = "Oligodendrocytes", 
  Comparison="Case", 
  alpha=0.05, 
  ColDict = ColDict_Cor
)

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/Fig_ATAC/", 
    "LinkPeaks_Correlation_ALS_ExcLINC00507_OLG", 
    ".pdf"
  ),
  width = 100, 
  height = 100, 
  units="mm"
)  


plot_cor_WNN_L25_CellTypes(
  L2FC_Shrink_Results_8SVs_ALS, 
  DESeq_Results_8SVs_Index, 
  CellType1 = "Exc_RORB", 
  CellType2 = "Oligodendrocytes", 
  Comparison="Case", 
  alpha=0.05, 
  ColDict = ColDict_Cor
)

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/Fig_ATAC/", 
    "LinkPeaks_Correlation_ALS_ExcRORB_OLG", 
    ".pdf"
  ),
  width = 100, 
  height = 100, 
  units="mm"
)  



    ## 6.3 DAR WNN_L25 ALS_FTD L2FC Shrinkage ----------------------------------

plot_cor_WNN_L25_CellTypes(
  L2FC_Shrink_Results_8SVs_ALSFTD, 
  DESeq_Results_8SVs_Index, 
  CellType1 = "Exc_LINC00507", 
  CellType2 = "Exc_RORB", 
  Comparison="Case", 
  alpha=0.05, 
  ColDict = ColDict_Cor
)

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/Fig_ATAC/", 
    "LinkPeaks_Correlation_ALSFTD_ExcLINC00507_ExcRORB", 
    ".pdf"
  ),
  width = 100, 
  height = 100, 
  units="mm"
)  


plot_cor_WNN_L25_CellTypes(
  L2FC_Shrink_Results_8SVs_ALSFTD, 
  DESeq_Results_8SVs_Index, 
  CellType1 = "Exc_LINC00507", 
  CellType2 = "Oligodendrocytes", 
  Comparison="Case", 
  alpha=0.05, 
  ColDict = ColDict_Cor
)

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/Fig_ATAC/", 
    "LinkPeaks_Correlation_ALSFTD_ExcLINC00507_OLG", 
    ".pdf"
  ),
  width = 100, 
  height = 100, 
  units="mm"
)  


plot_cor_WNN_L25_CellTypes(
  L2FC_Shrink_Results_8SVs_ALSFTD, 
  DESeq_Results_8SVs_Index, 
  CellType1 = "Exc_RORB", 
  CellType2 = "Oligodendrocytes", 
  Comparison="Case", 
  alpha=0.05, 
  ColDict = ColDict_Cor
)

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/Fig_ATAC/", 
    "LinkPeaks_Correlation_ALSFTD_ExcRORB_OLG", 
    ".pdf"
  ),
  width = 100, 
  height = 100, 
  units="mm"
)  



  ### 7.0 DAR WNN_L25 Disease correlation plots --------------------------------



    ## 7.1 Define plotting function --------------------------------------------

plot_cor_Case <- function(
    ResultsList1, 
    ResultsList2,
    ResultsListIndex, 
    CellTypeLevel = "WNN_L25", 
    CellType, 
    Comparison="Case", 
    Case1, 
    Case2,
    alpha=0.05, 
    ColDict=NULL
){
  
  if(Comparison=="Case"){
    
    ind <- which(
      ResultsListIndex$CellTypeLevel == CellTypeLevel & 
        ResultsListIndex$CellType == CellType & 
        ResultsListIndex$Comparison=="All_Cases"
    ) 
    
    
    DARs_Case1 <- Get_DARs(ResultsList1[[ind]], alpha = alpha)
    DARs_Case2 <- Get_DARs(ResultsList2[[ind]], alpha = alpha)
    
    df <- data.frame(
      Peak = unique(
        c(
          DARs_Case1, 
          DARs_Case2
        )
      )
    ) 
    
    df$L2FC_Case1 <- ResultsList1[[ind]]$log2FoldChange[match(
      df$Peak, 
      rownames(ResultsList1[[ind]])
    )]
    
    df$L2FC_Case2 <- ResultsList2[[ind]]$log2FoldChange[match(
      df$Peak, 
      rownames(ResultsList2[[ind]])
    )] 
    
    df$PAdj_Case1 <- ResultsList1[[ind]]$padj[match(
      df$Peak, 
      rownames(ResultsList1[[ind]])
    )]
    
    df$PAdj_Case2 <- ResultsList2[[ind]]$padj[match(
      df$Peak, 
      rownames(ResultsList2[[ind]])
    )]
    
  } else {
    
    if(Comparison=="Rand"){
      
      ind <- which(
        ResultsListIndex$CellTypeLevel == CellTypeLevel & 
          ResultsListIndex$CellType == CellType & 
          ResultsListIndex$Comparison=="All_Cases"
      )
      
      DARs_Case1 <- Get_DARs(ResultsList1[[ind]], alpha = alpha)
      DARs_Case2 <- Get_DARs(ResultsList2[[ind]], alpha = alpha)
      

      ind <- which(
        ResultsListIndex$CellTypeLevel == CellTypeLevel & 
          ResultsListIndex$CellType == CellType & 
          ResultsListIndex$Comparison=="Rand"
      )
      
      res1 <- ResultsList1[[ind]] 
      res1 <- res1[order(res1$pvalue),] 
      DARs_Case1 <- rownames(res1)[1:length(DARs_Case1)]
      
      res2 <- ResultsList2[[ind]] 
      res2 <- res2[order(res2$pvalue),] 
      DARs_Case2 <- rownames(res2)[1:length(DARs_Case2)] 
      
      
      df <- data.frame(
        Peak = unique(
          c(
            DARs_Case1, 
            DARs_Case2
          )
        )
      ) 
      
      df$L2FC_Case1 <- ResultsList1[[ind]]$log2FoldChange[match(
        df$Peak, 
        rownames(ResultsList1[[ind]])
      )]
      
      df$L2FC_Case2 <- ResultsList2[[ind]]$log2FoldChange[match(
        df$Peak, 
        rownames(ResultsList2[[ind]])
      )] 
      
      df$PAdj_Case1 <- ResultsList1[[ind]]$padj[match(
        df$Peak, 
        rownames(ResultsList1[[ind]])
      )]
      
      df$PAdj_Case2 <- ResultsList2[[ind]]$padj[match(
        df$Peak, 
        rownames(ResultsList2[[ind]])
      )]
      
    } 
    
  }
  
  
  df$Sign <- "None"
  df$Sign[df$PAdj_Case1 < alpha] <- Case1
  df$Sign[df$PAdj_Case2 < alpha] <- Case2 
  df$Sign[df$PAdj_Case1 < alpha & df$PAdj_Case2 < alpha] <- "Both"
  df$Sign <- factor(df$Sign)
  lim=max(abs(c(df$L2FC_Case1, df$L2FC_Case2)))
  
  cor <- cor.test(df$L2FC_Case1, df$L2FC_Case2) 
  
  p1 <- 
    ggplot() + 
    geom_hline(yintercept = 0, col = "#CCCCCC", size = 1) + 
    geom_vline(xintercept = 0, col = "#CCCCCC", size = 1) + 
    (if(cor$p.value<0.05) geom_smooth(aes(L2FC_Case1, L2FC_Case2), data = df, method = "lm", col="#238BFA", alpha=0.6, fill= "#CCCCCC44")) +  
    geom_point(aes(x = L2FC_Case1, y = L2FC_Case2, fill=Sign), data=df[df$Sign=="None",], size = 4, pch = 21, alpha = 0.6) + 
    geom_point(aes(x = L2FC_Case1, y = L2FC_Case2, fill=Sign), data=df[df$Sign==Case1,], size = 4, pch = 21, alpha = 0.6) + 
    geom_point(aes(x = L2FC_Case1, y = L2FC_Case2, fill=Sign), data=df[df$Sign==Case2,], size = 4, pch = 21, alpha = 0.6) + 
    geom_point(aes(x = L2FC_Case1, y = L2FC_Case2, fill=Sign), data=df[df$Sign=="Both",], size = 5, pch = 21, alpha = 1) + 
    scale_x_continuous(limits=c(-lim, lim)) +  
    scale_y_continuous(limits=c(-lim, lim)) +  
    (if(!is.null(ColDict)) {scale_fill_manual(values=ColDict)}) + 
    xlab(paste0("\n", "L2FC ", Case1)) + 
    ylab(paste0("L2FC ", Case2, "\n")) + 
    theme_classic() + 
    theme(
      panel.background = element_rect(color="#000000", size=1.2), 
      axis.line = element_blank(), 
      axis.title = element_text(face="bold.italic", color="#000000"), 
      axis.text = element_text(face="bold.italic", color="#000000"), 
      legend.box = element_blank(),  
      legend.box.background = element_blank(), 
      legend.background = element_blank(), 
      legend.key = element_blank(), 
      legend.text = element_text(face="bold.italic", color="#000000"), 
      legend.position = "Null"
    )
  
  message(paste0("#N Peaks: ", nrow(df)))
  message(paste0("Estimate: ", round(cor$estimate, 2)))
  message(paste0("PVal: ", signif(cor$p.value, 2)))
  message(paste0("95CI: ", round(cor$conf.int, 2)[1], "-", round(cor$conf.int, 2)[2])) 
  
  p1 + annotate(
    geom = "text", 
    x = -0.99*lim, 
    y = 0.95*lim, 
    label = paste0(
      paste0(
        "r = ", 
        round(cor$estimate, 2), 
        " (", 
        round(cor$conf.int, 2)[1], 
        " - ", 
        round(cor$conf.int, 2)[2], 
        ")"
      ), 
      "\n", 
      (if(cor$p.value >= 0.05){paste0("n.s.")}),
      (if(cor$p.value < 0.05 & cor$p.value >= 0.01){paste0("*p<0.05")}), 
      (if(cor$p.value < 0.01 & cor$p.value >= 0.001){paste0("**p<0.01")}), 
      (if(cor$p.value < 0.001 & cor$p.value >= 0.0001){paste0("***p<0.01")}), 
      (if(cor$p.value < 0.0001){paste0("****p<", signif(cor$p.value, 1))})
    ), 
    hjust=0, 
    vjust=1, 
    fontface=4,
    color="#000000"
  )
}

ColDict_Cor <- c(ColDict_Case, "Both"="#FBFCCC", "None"="#DDDDDD")



    ## 7.2 DAR WNN_L25 ALS ALSFTD L2FC Shrinkage -------------------------------

plot_cor_Case(
  ResultsList1 = L2FC_Shrink_Results_8SVs_ALS, 
  ResultsList2 = L2FC_Shrink_Results_8SVs_ALSFTD, 
  ResultsListIndex = DESeq_Results_8SVs_Index, 
  Case1 = "ALS", 
  Case2 = "ALS_FTD", 
  CellType="Exc_LINC00507", 
  Comparison="Case", 
  alpha=0.05, 
  ColDict = ColDict_Cor
)

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/Fig_ATAC/", 
    "LinkPeaks_Correlation_ExcLINC00507_ALS_ALSFTD", 
    ".pdf"
  ),
  width = 100, 
  height = 100, 
  units="mm"
)  


plot_cor_Case(
  ResultsList1 = L2FC_Shrink_Results_8SVs_ALS, 
  ResultsList2 = L2FC_Shrink_Results_8SVs_ALSFTD, 
  ResultsListIndex = DESeq_Results_8SVs_Index, 
  Case1 = "ALS", 
  Case2 = "ALS_FTD", 
  CellType="Exc_RORB", 
  Comparison="Case", 
  alpha=0.05, 
  ColDict = ColDict_Cor
)

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/Fig_ATAC/", 
    "LinkPeaks_Correlation_ExcRORB_ALS_ALSFTD", 
    ".pdf"
  ),
  width = 100, 
  height = 100, 
  units="mm"
)  


plot_cor_Case(
  ResultsList1 = L2FC_Shrink_Results_8SVs_ALS, 
  ResultsList2 = L2FC_Shrink_Results_8SVs_ALSFTD, 
  ResultsListIndex = DESeq_Results_8SVs_Index, 
  Case1 = "ALS", 
  Case2 = "ALS_FTD", 
  CellType="Oligodendrocytes", 
  Comparison="Case", 
  alpha=0.05, 
  ColDict = ColDict_Cor
)

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/Fig_ATAC/", 
    "LinkPeaks_Correlation_OLG_ALS_ALSFTD", 
    ".pdf"
  ),
  width = 100, 
  height = 100, 
  units="mm"
)   




  ### 8.0 DARs Peaks annotation barplots ---------------------------------------



  ## 8.1 Define plotting function --------------------------------------------

plot_DAR_Annotation <- function(
    ResultsList1, 
    ResultsList2,
    ResultsListIndex, 
    CellTypeLevel = "WNN_L25", 
    CellType, 
    Comparison="Case", 
    Case1, 
    Case2, 
    PeaksAnnotation, 
    alpha=0.05, 
    ColDict=NULL
){
  
  if(Comparison=="Case"){
    
    ind <- which(
      ResultsListIndex$CellTypeLevel == CellTypeLevel & 
        ResultsListIndex$CellType == CellType & 
        ResultsListIndex$Comparison=="All_Cases"
    ) 
    
    
    DARs_Case1 <- Get_DARs(ResultsList1[[ind]], alpha = alpha)
    DARs_Case2 <- Get_DARs(ResultsList2[[ind]], alpha = alpha)
    
    df <- data.frame(
      Peak = unique(
        c(
          DARs_Case1, 
          DARs_Case2
        )
      )
    ) 
    
    df$Annotation_L2_Plain <- Peaks$Annotation_L2_Plain[match(
      df$Peak, 
      Peaks$ID2
    )]
    
    
  } else {
    
    if(Comparison=="Rand"){
      
      ind <- which(
        ResultsListIndex$CellTypeLevel == CellTypeLevel & 
          ResultsListIndex$CellType == CellType & 
          ResultsListIndex$Comparison=="All_Cases"
      )
      
      DARs_Case1 <- Get_DARs(ResultsList1[[ind]], alpha = alpha)
      DARs_Case2 <- Get_DARs(ResultsList2[[ind]], alpha = alpha)
      
      
      ind <- which(
        ResultsListIndex$CellTypeLevel == CellTypeLevel & 
          ResultsListIndex$CellType == CellType & 
          ResultsListIndex$Comparison=="Rand"
      )
      
      res1 <- ResultsList1[[ind]] 
      res1 <- res1[order(res1$pvalue),] 
      DARs_Case1 <- rownames(res1)[1:length(DARs_Case1)]
      
      res2 <- ResultsList2[[ind]] 
      res2 <- res2[order(res2$pvalue),] 
      DARs_Case2 <- rownames(res2)[1:length(DARs_Case2)] 
      
      
      df <- data.frame(
        Peak = unique(
          c(
            DARs_Case1, 
            DARs_Case2
          )
        )
      )
      
      df$Annotation_L2_Plain <- Peaks$Annotation_L2_Plain[match(
        df$Peak, 
        Peaks$ID2
      )]
 
    } 
    
  }
  
  p1 <- data.frame(
    table(
      df$Annotation_L2_Plain
    )
  ) %>% 
  ggplot() + 
    aes(
      x="All", 
      y=Freq, 
      fill=factor(Var1, levels=c("Intergenic", "Promoter", "Exon", "Intron", "Downstream"))
    ) + 
    geom_col(position="fill", col="#000000") + 
    scale_y_continuous(expand=c(0,0,0.05,0)) + 
    scale_fill_manual(values=c("#FFA94D", "#0d0887", "#7e03a8", "#cc4778", "#f89540")) + 
    theme_classic() + 
    theme(
      legend.title = element_blank(), 
      axis.title.x = element_blank(), 
      axis.title.y = element_text(size=14, color="#000000", face="bold.italic"), 
      axis.text = element_text(color="#000000", face="bold.italic", size=12), 
      legend.text = element_text(size=12, color="#000000", face="bold.italic") 
      
    )
  
  a <- round(prop.table(table(df$Annotation_L2_Plain))*100, 2)
  
  for (i in 1:length(a)){
    
    message(paste0(names(a)[i], ": ", a[i], " %"))
  
  }
  
  p1
  
}



    ## 8.2 Plot ALS/ALS-FTD Peaks Annotation -----------------------------------

plot_DAR_Annotation(
  ResultsList1 = L2FC_Shrink_Results_8SVs_ALS, 
  ResultsList2 = L2FC_Shrink_Results_8SVs_ALSFTD, 
  ResultsListIndex = DESeq_Results_8SVs_Index, 
  CellType="Exc_LINC00507", 
  Comparison="Case", 
  PeaksAnnotation = Peaks, 
  alpha=0.05, 
  ColDict = ColDict_Cor
)

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/Fig_ATAC/", 
    "LinkPeaks_DARs_ALS_ALSFTD_Annotations_ExcLINC00507_Barplot", 
    ".pdf"
  ),
  width = 100, 
  height = 100, 
  units="mm"
)  


plot_DAR_Annotation(
  ResultsList1 = L2FC_Shrink_Results_8SVs_ALS, 
  ResultsList2 = L2FC_Shrink_Results_8SVs_ALSFTD, 
  ResultsListIndex = DESeq_Results_8SVs_Index, 
  CellType="Exc_RORB", 
  Comparison="Case", 
  PeaksAnnotation = Peaks, 
  alpha=0.05, 
  ColDict = ColDict_Cor
)

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/Fig_ATAC/", 
    "LinkPeaks_DARs_ALS_ALSFTD_Annotations_ExcRORB_Barplot", 
    ".pdf"
  ),
  width = 100, 
  height = 100, 
  units="mm"
)  


plot_DAR_Annotation(
  ResultsList1 = L2FC_Shrink_Results_8SVs_ALS, 
  ResultsList2 = L2FC_Shrink_Results_8SVs_ALSFTD, 
  ResultsListIndex = DESeq_Results_8SVs_Index, 
  CellType="Oligodendrocytes", 
  Comparison="Case", 
  PeaksAnnotation = Peaks, 
  alpha=0.05, 
  ColDict = ColDict_Cor
)


ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/Fig_ATAC/", 
    "LinkPeaks_DARs_ALS_ALSFTD_Annotations_Oligodendrocytes_Barplot", 
    ".pdf"
  ),
  width = 100, 
  height = 100, 
  units="mm"
)  

rm(i, ind)




  ### 9.0 Plot ATAC/RNA Overlap ------------------------------------------------



    ## 9.1 Define plotting function --------------------------------------------

plot_ATAC_RNA_Correlation <- function(
    
    ATAC_RNA_Links, 
    ResultsList_ATAC,
    ResultsList_ATAC_Index, 
    ResultsList_RNA,
    ResultsList_RNA_Index, 
    CellTypeLevel = "WNN_L25", 
    CellType, 
    Comparison="All_Cases", 
    alpha=0.05, 
    ColDict=NULL, 
    Gene_Universe, 
    pt.size=3, 
    returnData=FALSE, 
    text.size=12
){
  
  if(Comparison=="All_Cases"){
    
    ind_ATAC <- which(
      ResultsList_ATAC_Index$CellTypeLevel == "WNN_L25" & 
        ResultsList_ATAC_Index$CellType == CellType & 
        ResultsList_ATAC_Index$Comparison == Comparison
    ) 
    
    atac_res <-  as.data.frame(ResultsList_ATAC[[ind_ATAC]])
      
    
    ind_RNA <- which(
      ResultsList_RNA_Index$CellTypeLevel == "WNN_L25" & 
        ResultsList_RNA_Index$CellType == CellType & 
        ResultsList_RNA_Index$Comparison == Comparison
    )
    
    
    rna_res <-  as.data.frame(
      ResultsList_RNA[[ind_RNA]]
    )
      
    
    DARs <- rownames(atac_res)[which(atac_res$padj<alpha)] 
    DEGs <- rownames(rna_res)[which(rna_res$padj<alpha)]  

  } else {
    
    if(Comparison=="Rand"){
      
      ind_ATAC <- which(
        ResultsList_ATAC_Index$CellTypeLevel == "WNN_L25" & 
          ResultsList_ATAC_Index$CellType == CellType & 
          ResultsList_ATAC_Index$Comparison == "All_Cases"
      ) 
      
      atac_res <-  as.data.frame(ResultsList_ATAC[[ind_ATAC]])
      
      
      ind_RNA <- which(
        ResultsList_RNA_Index$CellTypeLevel == "WNN_L25" & 
          ResultsList_RNA_Index$CellType == CellType & 
          ResultsList_RNA_Index$Comparison == "All_Cases"
      )
      
      
      rna_res <-  as.data.frame(
        ResultsList_RNA[[ind_RNA]]
      )
      
      
      DARs <- rownames(atac_res)[which(atac_res$padj<alpha)] 
      DEGs <- rownames(rna_res)[which(rna_res$padj<alpha)]  
      
      ind_ATAC <- which(
        ResultsList_ATAC_Index$CellTypeLevel == "WNN_L25" & 
          ResultsList_ATAC_Index$CellType == CellType & 
          ResultsList_ATAC_Index$Comparison == Comparison
      ) 
      
      ind_RNA <- which(
        ResultsList_RNA_Index$CellTypeLevel == "WNN_L25" & 
          ResultsList_RNA_Index$CellType == CellType & 
          ResultsList_RNA_Index$Comparison == Comparison
      )
      
      atac_res <- ResultsList_ATAC[[ind_ATAC]] 
      atac_res <- atac_res[order(atac_res$pvalue),] 
      DARs <- rownames(atac_res)[1:length(DARs)]
      
      rna_res <- ResultsList_RNA[[ind_RNA]] 
      rna_res <- rna_res[order(rna_res$pvalue),] 
      DEGs <- rownames(rna_res)[1:length(DEGs)]
      
    }
  }
      

  df <- data.frame(
    peak = ATAC_RNA_Links$peak[ATAC_RNA_Links$peak %in% DARs], 
    gene = ATAC_RNA_Links$gene[ATAC_RNA_Links$peak %in% DARs]
  )

  df$Peak_L2FC <- atac_res$log2FoldChange[match(df$peak, rownames(atac_res))]
  df$Gene_L2FC <- rna_res$log2FoldChange[match(df$gene, rownames(rna_res))]
  
  coef1 = table(df$gene %in% DEGs)["TRUE"]
  coef2 = table(df$gene %in% DEGs)["FALSE"]
  coef3 = table(DEGs %in% df$gene)["FALSE"]
  coef4 = length(unique(c(unique(ATAC_RNA_Links$gene), Gene_Universe[[CellType]]$ID_10X)))
  
  fisher = fisher.test(
    matrix(
      c(
        coef1, 
        coef2, 
        coef3, 
        coef4
      ), 
      2, 
      2
    )
  )

  df$DEG <- FALSE
  df$DEG[df$gene %in% DEGs] <- TRUE
  
  
  lim=max(abs(c(df$Peak_L2FC, df$Gene_L2FC)))
  
  cor <- cor.test(df$Peak_L2FC[df$DEG], df$Gene_L2FC[df$DEG]) 

  p1 <- 
    ggplot() + 
    geom_hline(yintercept = 0, col = "#CCCCCC", size = 1) + 
    geom_vline(xintercept = 0, col = "#CCCCCC", size = 1) + 
    (if(cor$p.value<0.05 & coef1>10) geom_smooth(aes(Peak_L2FC, Gene_L2FC), data = df[df$DEG,], method = "lm", col="#F0C000", alpha=0.6, fill= "#CCCCCC44")) +  
    geom_point(aes(x = Peak_L2FC, y = Gene_L2FC, fill=DEG), data=df[!df$DEG,], size = pt.size, pch = 21, alpha = 0.6) + 
    geom_point(aes(x = Peak_L2FC, y = Gene_L2FC, fill=DEG), data=df[df$DEG,], size = pt.size, pch = 21, alpha = 0.6) + 
    scale_x_continuous(limits=c(-lim, lim)) +  
    scale_y_continuous(limits=c(-lim, lim)) +  
    (if(!is.null(ColDict)) {scale_fill_manual(values=ColDict)}) + 
    xlab(paste0("\n", "L2FC ATAC")) + 
    ylab(paste0("L2FC RNA", "\n")) + 
    theme_classic() + 
    theme(
      panel.background = element_rect(color="#000000", size=1.2), 
      axis.line = element_blank(), 
      axis.title = element_text(face="bold.italic", color="#000000", size=text.size), 
      axis.text = element_text(face="bold.italic", color="#000000", size=text.size), 
      legend.box = element_blank(),  
      legend.box.background = element_blank(), 
      legend.background = element_blank(), 
      legend.key = element_blank(), 
      legend.text = element_text(face="bold.italic", color="#000000", size=text.size), 
      legend.position = "Null", 
      axis.ticks =  element_line(color="#000000", linewidth = 1.2), 
      axis.ticks.length = unit(6, "points")
    )
  
  tryCatch({
    message(paste0("#N Peak-gene links: ", nrow(df)))
    message(paste0("#N Unique peaks: ", length(unique(df$peak))))
    message(paste0("#N Unique genes: ", length(unique(df$gene))))
    message(paste0("\n"))
    
    message(paste0("#PEARSON Estimate: ", round(cor$estimate, 2)))
    message(paste0("#PEARSON PVal: ", signif(cor$p.value, 2)))
    message(paste0("#PEARSON 95CI: ", round(cor$conf.int, 2)[1], "-", round(cor$conf.int, 2)[2])) 
    
    message(paste0("\n"))
    message(paste0("# FISHER OddsRatio: ", round(fisher$estimate,2)))
    message(paste0("# FISHER PVal: ", signif(fisher$p.value,2)))
    message(paste0("# FISHER 95CI: ", round(fisher$conf.int, 2)[1], "-", round(fisher$conf.int, 2)[2]))
  }, error=function(e){})
  
  
  
  p2 <- p1 
  
  tryCatch({
    p2 <- p2 + annotate(
      geom = "text", 
      x = -0.99*lim, 
      y = 0.95*lim, 
      label = paste0(
        paste0(
          "r = ", 
          round(cor$estimate, 2), 
          "\n(", 
          round(cor$conf.int, 2)[1], 
          " - ", 
          round(cor$conf.int, 2)[2], 
          ")"
        ), 
        "\n", 
        (if(cor$p.value >= 0.05){paste0("n.s.")}),
        (if(cor$p.value < 0.05 & cor$p.value >= 0.01){paste0("*p<0.05")}), 
        (if(cor$p.value < 0.01 & cor$p.value >= 0.001){paste0("**p<0.01")}), 
        (if(cor$p.value < 0.001 & cor$p.value >= 0.0001){paste0("***p<0.001")}), 
        (if(cor$p.value < 0.0001){paste0("****p<0.001")})
      ), 
      hjust=0, 
      vjust=1, 
      fontface=4,
      color="#F0C000", 
      size=6
    )}, error=function(e){}
  )

  plot(p2)
  if(returnData){
    return(df)
  }
}



    ## 9.2 ATAC RNA Correlation ALS --------------------------------------------

plot_ATAC_RNA_Correlation(
  ATAC_RNA_Links = ATAC_Links, 
  ResultsList_ATAC = L2FC_Shrink_Results_8SVs_ALS,
  ResultsList_ATAC_Index = DESeq_Results_8SVs_Index, 
  ResultsList_RNA = RNA_L2FC_Shrink_Results_12SVs_ALS,
  ResultsList_RNA_Index = RNA_L2FC_Shrink_Results_12SVs_Index, 
  CellTypeLevel = "WNN_L25", 
  CellType = "Exc_LINC00507", 
  Comparison="All_Cases", 
  alpha=0.05, 
  ColDict=c("FALSE"=as.character(ColDict_Modality["ATAC"]), "TRUE"=as.character(ColDict_Modality["RNA"])), 
  Gene_Universe = Universe_AllCase_ALS_List, 
  pt.size=4, 
  returnData = FALSE, 
  text.size=14
)

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/Fig_ATAC/", 
    "ATAC_RNA_Correlation_ALS_ExcLINC00507", 
    ".pdf"
  ),
  width = 100, 
  height = 100, 
  units="mm"
)  



plot_ATAC_RNA_Correlation(
  ATAC_RNA_Links = ATAC_Links, 
  ResultsList_ATAC = L2FC_Shrink_Results_8SVs_ALS,
  ResultsList_ATAC_Index = DESeq_Results_8SVs_Index, 
  ResultsList_RNA = RNA_L2FC_Shrink_Results_12SVs_ALS,
  ResultsList_RNA_Index = RNA_L2FC_Shrink_Results_12SVs_Index, 
  CellTypeLevel = "WNN_L25", 
  CellType = "Exc_RORB", 
  Comparison="All_Cases", 
  alpha=0.05, 
  ColDict=c("FALSE"=as.character(ColDict_Modality["ATAC"]), "TRUE"=as.character(ColDict_Modality["RNA"])), 
  Gene_Universe = Universe_AllCase_ALS_List, 
  pt.size=4, 
  returnData = FALSE, 
  text.size=14
)

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/Fig_ATAC/", 
    "ATAC_RNA_Correlation_ALS_ExcRORB", 
    ".pdf"
  ),
  width = 100, 
  height = 100, 
  units="mm"
)  



plot_ATAC_RNA_Correlation(
  ATAC_RNA_Links = ATAC_Links, 
  ResultsList_ATAC = L2FC_Shrink_Results_8SVs_ALS,
  ResultsList_ATAC_Index = DESeq_Results_8SVs_Index, 
  ResultsList_RNA = RNA_L2FC_Shrink_Results_12SVs_ALS,
  ResultsList_RNA_Index = RNA_L2FC_Shrink_Results_12SVs_Index, 
  CellTypeLevel = "WNN_L25", 
  CellType = "Oligodendrocytes", 
  Comparison="All_Cases", 
  alpha=0.05, 
  ColDict=c("FALSE"=as.character(ColDict_Modality["ATAC"]), "TRUE"=as.character(ColDict_Modality["RNA"])), 
  Gene_Universe = Universe_AllCase_ALS_List, 
  pt.size=3, 
  returnData = FALSE, 
  text.size=14
)

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/Fig_ATAC/", 
    "ATAC_RNA_Correlation_ALS_Oligodendrocytes", 
    ".pdf"
  ),
  width = 100, 
  height = 100, 
  units="mm"
)  



    ## 9.3 ATAC RNA Correlation ALSFTD -----------------------------------------

plot_ATAC_RNA_Correlation(
  ATAC_RNA_Links = ATAC_Links, 
  ResultsList_ATAC = L2FC_Shrink_Results_8SVs_ALSFTD,
  ResultsList_ATAC_Index = DESeq_Results_8SVs_Index, 
  ResultsList_RNA = RNA_L2FC_Shrink_Results_12SVs_ALSFTD,
  ResultsList_RNA_Index = RNA_L2FC_Shrink_Results_12SVs_Index, 
  CellTypeLevel = "WNN_L25", 
  CellType = "Exc_LINC00507", 
  Comparison="All_Cases", 
  alpha=0.05, 
  ColDict=c("FALSE"=as.character(ColDict_Modality["ATAC"]), "TRUE"=as.character(ColDict_Modality["RNA"])), 
  Gene_Universe = Universe_AllCase_ALS_List, 
  pt.size=3, 
  returnData = FALSE, 
  text.size=14
)

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/Fig_ATAC/", 
    "ATAC_RNA_Correlation_ALSFTD_ExcLINC00507", 
    ".pdf"
  ),
  width = 100, 
  height = 100, 
  units="mm"
)  



plot_ATAC_RNA_Correlation(
  ATAC_RNA_Links = ATAC_Links, 
  ResultsList_ATAC = L2FC_Shrink_Results_8SVs_ALSFTD,
  ResultsList_ATAC_Index = DESeq_Results_8SVs_Index, 
  ResultsList_RNA = RNA_L2FC_Shrink_Results_12SVs_ALSFTD,
  ResultsList_RNA_Index = RNA_L2FC_Shrink_Results_12SVs_Index, 
  CellTypeLevel = "WNN_L25", 
  CellType = "Exc_RORB", 
  Comparison="All_Cases", 
  alpha=0.05, 
  ColDict=c("FALSE"=as.character(ColDict_Modality["ATAC"]), "TRUE"=as.character(ColDict_Modality["RNA"])), 
  Gene_Universe = Universe_AllCase_ALS_List, 
  pt.size=3, 
  returnData = FALSE, 
  text.size=14
)

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/Fig_ATAC/", 
    "ATAC_RNA_Correlation_ALSFTD_ExcRORB", 
    ".pdf"
  ),
  width = 100, 
  height = 100, 
  units="mm"
)  



plot_ATAC_RNA_Correlation(
  ATAC_RNA_Links = ATAC_Links, 
  ResultsList_ATAC = L2FC_Shrink_Results_8SVs_ALSFTD,
  ResultsList_ATAC_Index = DESeq_Results_8SVs_Index, 
  ResultsList_RNA = RNA_L2FC_Shrink_Results_12SVs_ALSFTD,
  ResultsList_RNA_Index = RNA_L2FC_Shrink_Results_12SVs_Index, 
  CellTypeLevel = "WNN_L25", 
  CellType = "Oligodendrocytes", 
  Comparison="All_Cases", 
  alpha=0.05, 
  ColDict=c("FALSE"=as.character(ColDict_Modality["ATAC"]), "TRUE"=as.character(ColDict_Modality["RNA"])), 
  Gene_Universe = Universe_AllCase_ALS_List, 
  pt.size=3, 
  returnData = FALSE, 
  text.size=14
)

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/Fig_ATAC/", 
    "ATAC_RNA_Correlation_ALSFTD_Oligodendrocytes", 
    ".pdf"
  ),
  width = 100, 
  height = 100, 
  units="mm"
)  





