
  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 

library(qs) 
library(tidyverse) 
library(data.table)
library(ggrastr)
library(APAlyzer)

  

qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)




  ### 1.0 Load data ------------------------------------------------------------ 

F2 <- qread(
  paste0(
    "../Data/SeuratObjects/", 
    "F2", 
    ".qrds"
  ), 
  nthr=nthr
)

FANS_DE_TDP43_Wilcox <- qread(
  paste0(
    "../Data/DE/FANS/", 
    "FANS_DE_TDP43_Wilcox", 
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

APA_Results <- qread(
  paste0(
    "../Data/APA/DaPars1_TDP43Pos_vs_TDP43Neg/", 
    "APA_Results", 
    ".qrds"
  ), 
  nthr=nthr
)

APA_Results <- APA_Results[!is.na(APA_Results$PDUI_Group_diff),]
APA_Results$Sign <- FALSE
APA_Results$Sign[APA_Results$adjusted.P_val<0.0001] <- TRUE
APA_Results$Threshold <- FALSE
APA_Results$Threshold[abs(APA_Results$PDUI_Group_diff)>=0.2] <- TRUE

table(APA_Results$Threshold, APA_Results$Sign)

APA_Results$SignVolcano <- APA_Results$Threshold & APA_Results$Sign

APA_Results$adjusted.P_val[APA_Results$adjusted.P_val==0] <- min(APA_Results$adjusted.P_val[!APA_Results$adjusted.P_val==0])


APAlyzer_Results <- qread(
  paste0(
    "../Data/APA/APAlyzer_TFP43Neg_vs_TDP43Pos/", 
    "Res_3UTR_TDPNeg_vs_TDPPos", 
    ".qrds"
  ), 
  nthr=nthr
)

ColDict_WNN_L25 <- setNames(
  object = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "WNN_L25")$Color, 
  nm = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "WNN_L25")$WNN_L25
)


ALS_ALSFTD_AllCells_Signature <- qread(
  paste0(
    "../Data/Annotations/Signatures/RNA/", 
    "ALS_ALSFTD_AllCells_Signature", 
    ".qrds"
  ), 
  nthr=nthr
)

ALS_ALSFTD_WNN_L25_Signature <- qread(
  paste0(
    "../Data/Annotations/Signatures/RNA/", 
    "ALS_ALSFTD_WNN_L25_Signature", 
    ".qrds"
  ), 
  nthr=nthr
)


DESeq_Results_12SVs_Disease_AllCells_list <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/Disease/SVA/12_SVs/", 
    "DESeq_Results_12SVs_Disease_AllCells_list", 
    ".qrds"
  ), 
  nthr=nthr
)

DESeq_Results_12SVs_Disease_AllCells_Index <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/Disease/SVA/12_SVs/", 
    "DESeq_Results_12SVs_Disease_AllCells_Index", 
    ".qrds"
  ), 
  nthr=nthr
)


DESeq_Results_12SVs_Disease_WNN_L25_list <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/Disease/SVA/12_SVs/", 
    "DESeq_Results_12SVs_Disease_WNN_L25_list", 
    ".qrds"
  ), 
  nthr=nthr
)

DESeq_Results_12SVs_Disease_WNN_L25_Index <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/Disease/SVA/12_SVs/", 
    "DESeq_Results_12SVs_Disease_WNN_L25_Index", 
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


DE_ALS_AllCells_MAST <- qread("../Data/DE/WNN/Single_Cell/MAST/ALS/AllCells/DE_MAST_ALSvsHC_AllCells_AllCells.qrds")
DE_ALSFTD_AllCells_MAST <- qread("../Data/DE/WNN/Single_Cell/MAST/ALSFTD/AllCells/DE_MAST_ALSFTDvsHC_AllCells_AllCells.qrds")
table(DE_ALS_AllCells_MAST$fdr<0.05)
table(DE_ALSFTD_AllCells_MAST$fdr<0.05)

table(DE_ALS_AllCells_MAST$gene[DE_ALS_AllCells_MAST$fdr<0.05] %in% ALS_ALSFTD_AllCells_Signature$ID_10X)
table(DE_ALS_AllCells_MAST$gene[DE_ALS_AllCells_MAST$fdr<0.05] %in% ALS_ALSFTD_WNN_L25_Signature$ID_10X)

table(DE_ALSFTD_AllCells_MAST$gene[DE_ALSFTD_AllCells_MAST$fdr<0.05] %in% ALS_ALSFTD_AllCells_Signature$ID_10X)
table(DE_ALSFTD_AllCells_MAST$gene[DE_ALSFTD_AllCells_MAST$fdr<0.05] %in% ALS_ALSFTD_WNN_L25_Signature$ID_10X)

table(FANS_DE_TDP43_MAST$gene[FANS_DE_TDP43_MAST$fdr<0.05] %in% ALS_ALSFTD_AllCells_Signature$ID_10X)
table(FANS_DE_TDP43_MAST$gene[FANS_DE_TDP43_MAST$fdr<0.05] %in% ALS_ALSFTD_WNN_L25_Signature$ID_10X) 
#!Check Richtung ALS TDP DE 
table(FANS_DE_TDP43_MAST$gene[FANS_DE_TDP43_MAST$fdr<0.05] %in% DE_ALS_AllCells_MAST$gene[DE_ALS_AllCells_MAST$fdr<0.05])
table(FANS_DE_TDP43_MAST$gene[FANS_DE_TDP43_MAST$fdr<0.05] %in% DE_ALSFTD_AllCells_MAST$gene[DE_ALS_AllCells_MAST$fdr<0.05])

MAST_LINC_ALS <- qread("../Data/DE/WNN/Single_Cell/MAST/ALS/WNN_L25/DE_MAST_ALSvsHC_WNN_L25_Exc_LINC00507.qrds")
MAST_LINC_ALSFTD <- qread("../Data/DE/WNN/Single_Cell/MAST/ALSFTD/WNN_L25/")
table(MAST_LINC_ALS$fdr<0.05)

table(MAST_LINC_ALS$gene[MAST_LINC_ALS$fdr<0.05] %in% FANS_DE_TDP43_MAST$gene[FANS_DE_TDP43_MAST$fdr<0.05])
Common <- MAST_LINC_ALS$gene[MAST_LINC_ALS$gene[MAST_LINC_ALS$fdr<0.05] %in% FANS_DE_TDP43_MAST$gene[FANS_DE_TDP43_MAST$fdr<0.05]]
Common <- data.frame(
  Gene=Common, 
  DE_L2FC = MAST_LINC_ALS$avg_log2FC[match(Common, MAST_LINC_ALS$gene)], 
  FANS_L2FC = FANS_DE_TDP43_MAST$avg_log2FC[match(Common, FANS_DE_TDP43_MAST$gene)]
)
plot(Common$DE_L2FC, Common$FANS_L2FC)
abline(h=0, col="red")
abline(v=0, col="red")
cor.test(Common$DE_L2FC, Common$FANS_L2FC)

genelist <- FANS_DE_TDP43_MAST$model_log2FC[abs(FANS_DE_TDP43_MAST$model_log2FC) > log2(1.5) & FANS_DE_TDP43_MAST$fdr <0.1]
names(genelist) <- FANS_DE_TDP43_MAST$gene[abs(FANS_DE_TDP43_MAST$model_log2FC) > log2(1.5) & FANS_DE_TDP43_MAST$fdr <0.1]


names(genelist) <- GEX_Features_ENSEMBL_ENTREZ$ENTREZID[match(names(genelist), GEX_Features_ENSEMBL_ENTREZ$ID_10X)]
genelist <- genelist[!is.na(names(genelist))]
genelist <- sort(genelist, decreasing = TRUE)
gseKegg <- gseKEGG(genelist, organism = "hsa", pvalueCutoff = 0.5)
enrKegg <- enrichKEGG(names(genelist), organism = "hsa")
enrKegg <- enrKegg@result
enrKegg <- enrKegg[enrKegg$qvalue<0.05,]


  ### 2.0 M0 Label Transfer Cell-Type proportions ------------------------------

table(F2$Predicted_ID_M0_WNN_L25)

ggplot(F2@meta.data) + 
  aes(y=TDP43 , fill=Predicted_ID_M0_WNN_L25) + 
  geom_bar(position="fill", col= "#000000") + 
  scale_x_continuous(expand=c(0,0,0.1,0.1)) + 
  scale_fill_manual(values=ColDict_WNN_L25) + 
  theme_classic()


ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/Fig_APA/", 
    "F2_M0_CellType_Proportions", 
    ".pdf"
  ),
  width = 420, 
  height = 200, 
  units="mm"
)  




  ### 3.0 FANS DE Wilcox Volcano Plot ------------------------------------------

FANS_DE_TDP43_Wilcox$Sign <- FALSE
FANS_DE_TDP43_Wilcox$Sign[abs(FANS_DE_TDP43_Wilcox$avg_log2FC) > log2(1.5) & FANS_DE_TDP43_Wilcox$p_val_adj <0.0001] <- TRUE
FANS_DE_TDP43_Wilcox$TDP43_Reg <- "DN"
FANS_DE_TDP43_Wilcox$TDP43_Reg[FANS_DE_TDP43_Wilcox$avg_log2FC>0] <- "UP"


ggplot() + 
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), col="#CCCCCC") + 
  geom_hline(yintercept = -log10(0.0001), col="#CCCCCC") + 
  ggrastr::geom_point_rast(aes(avg_log2FC, -log10(p_val_adj)), dat=FANS_DE_TDP43_Wilcox[!FANS_DE_TDP43_Wilcox$Sign,], alpha=0.2, col="#000000", fill="#000000", pch=21, size=3) + 
  ggrastr::geom_point_rast(aes(avg_log2FC, -log10(p_val_adj)), dat=FANS_DE_TDP43_Wilcox[FANS_DE_TDP43_Wilcox$Sign & !FANS_DE_TDP43_Wilcox$TDP43_Reg=="DN",], pch=21, col="#000000", fill="#4444FF",alpha=0.4, size=3) + 
  ggrastr::geom_point_rast(aes(avg_log2FC, -log10(p_val_adj)), dat=FANS_DE_TDP43_Wilcox[FANS_DE_TDP43_Wilcox$Sign & !FANS_DE_TDP43_Wilcox$TDP43_Reg=="UP",], pch=21, col="#000000", fill="#FF0000",alpha=0.4, size=3) + 
  scale_x_continuous(limits=c(-12,12)) + 
  scale_y_continuous(expand=c(0,4,0.1,0)) + 
  theme_classic()

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/Fig_APA/", 
    "FANS_DE_TDP43_Wilcox_Volcano", 
    ".pdf"
  ),
  width = 160, 
  height = 200, 
  units="mm"
)  




  ### 4.0 FANS DE MAST Volcano Plot ------------------------------------------

FANS_DE_TDP43_MAST$Sign <- FALSE
FANS_DE_TDP43_MAST$Sign[abs(FANS_DE_TDP43_MAST$model_log2FC) > log2(1.5) & FANS_DE_TDP43_MAST$fdr <0.1] <- TRUE
FANS_DE_TDP43_MAST$TDP43_Reg <- "DN"
FANS_DE_TDP43_MAST$TDP43_Reg[FANS_DE_TDP43_MAST$model_log2FC>0] <- "UP"


ggplot() + 
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), col="#CCCCCC") + 
  geom_hline(yintercept = -log10(0.1), col="#CCCCCC") + 
  ggrastr::geom_point_rast(aes(model_log2FC, -log10(fdr)), dat=FANS_DE_TDP43_MAST[!FANS_DE_TDP43_MAST$Sign,], alpha=0.2, col="#000000", fill="#000000", pch=21, size=3) + 
  ggrastr::geom_point_rast(aes(model_log2FC, -log10(fdr)), dat=FANS_DE_TDP43_MAST[FANS_DE_TDP43_MAST$Sign & !FANS_DE_TDP43_MAST$TDP43_Reg=="DN",], pch=21, col="#000000", fill="#4444FF",alpha=0.4, size=3) + 
  ggrastr::geom_point_rast(aes(model_log2FC, -log10(fdr)), dat=FANS_DE_TDP43_MAST[FANS_DE_TDP43_MAST$Sign & !FANS_DE_TDP43_MAST$TDP43_Reg=="UP",], pch=21, col="#000000", fill="#FF0000",alpha=0.4, size=3) + 
  scale_x_continuous(limits=c(-5,5)) +  
  scale_y_continuous(limits=c(0,4)) + 
  theme_classic()

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/Fig_APA/", 
    "FANS_DE_TDP43_MAST_Volcano", 
    ".pdf"
  ),
  width = 160, 
  height = 200, 
  units="mm"
)  




  ### 5.0 FANS APA TDP43PosTDP43Neg DaPars Volcano Plot ------------------------

ggplot() + 
  geom_vline(xintercept = c(-0.2, 0.2), col="#CCCCCC") + 
  ggrastr::geom_point_rast(aes(-PDUI_Group_diff, -log10(adjusted.P_val)), dat=APA_Results[!APA_Results$SignVolcano,],alpha=0.2, col="#000000", fill="#000000", pch=21, size=3) + 
  ggrastr::geom_point_rast(aes(-PDUI_Group_diff, -log10(adjusted.P_val)), dat=APA_Results[APA_Results$SignVolcano & APA_Results$PDUI_Group_diff>0,], pch=21, col="#000000", fill="#4444FF",alpha=0.4, size=3) + 
  ggrastr::geom_point_rast(aes(-PDUI_Group_diff, -log10(adjusted.P_val)), dat=APA_Results[APA_Results$SignVolcano & APA_Results$PDUI_Group_diff<0,], pch=21, col="#000000", fill="#FF0000",alpha=0.4, size=3) + 
  scale_x_continuous(limits=c(-1,1)) + 
  scale_y_continuous(expand=c(0,0,0.1,0)) + 
  theme_classic()

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/Fig_APA/", 
    "FANS_APA_DaPars_TDP43Neg_TDP43Pos_Volcano", 
    ".pdf"
  ),
  width = 160, 
  height = 200, 
  units="mm"
)  




  ### 6.0 FANS APA TDP43PosTDP43Neg APAlyzer Volcano Plot   --------------------

APAlyzer_Results$Sign <- FALSE
APAlyzer_Results$Sign[abs(APAlyzer_Results$RED)>1 & APAlyzer_Results$p_adj<0.0001] <- TRUE
table(APAlyzer_Results$Sign)

APAlyzer_Results$DoubleSign <- APAlyzer_Results$APAreg %in% c("DN", "UP") & APAlyzer_Results$Sign
table(APAlyzer_Results$DoubleSign)

ggplot() + 
  geom_vline(xintercept = c(-1, 1), col="#CCCCCC") + 
  ggrastr::geom_point_rast(aes(RED, -log10(p_adj)), dat=APAlyzer_Results[!APAlyzer_Results$DoubleSign,],alpha=0.2, col="#000000", fill="#000000", pch=21, size=3) + 
  ggrastr::geom_point_rast(aes(RED, -log10(p_adj)), dat=APAlyzer_Results[APAlyzer_Results$DoubleSign & APAlyzer_Results$APAreg == "DN",], pch=21, col="#000000", fill="#4444FF",alpha=0.4, size=3) + 
  ggrastr::geom_point_rast(aes(RED, -log10(p_adj)), dat=APAlyzer_Results[APAlyzer_Results$DoubleSign & APAlyzer_Results$APAreg == "UP",], pch=21, col="#000000", fill="#FF0000",alpha=0.4, size=3) + 
  scale_x_continuous(limits=c(-10,10)) + 
  scale_y_continuous(expand=c(0,4,0.1,0)) + 
  theme_classic()

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/Fig_APA/", 
    "FANS_APA_APAlyzer_TDP43Neg_TDP43Pos_Volcano", 
    ".pdf"
  ),
  width = 160, 
  height = 200, 
  units="mm"
)  




  ### 7.0 FANS APA TDP43PosTDP43Neg Barplot ------------------------------------

APA_Sign_Sorted <- APA_Results[APA_Results$SignVolcano & abs(APA_Results$PDUI_Group_diff)>0.2,]

common <- APA_Sign_Sorted$Gene[APA_Sign_Sorted$Symbol %in% APAlyzer_Results$gene_symbol[APAlyzer_Results$DoubleSign]]
common <- data.frame(
  Transcript=common
)
common$Gene = APA_Sign_Sorted$Symbol[match(common$Transcript,APA_Sign_Sorted$Gene)] 
common$DaPars=APA_Sign_Sorted$PDUI_Group_diff[match(common$Gene, APA_Sign_Sorted$Symbol)]
common$APAlyzer=APAlyzer_Results$APAreg[match(common$Gene, APAlyzer_Results$gene_symbol)]

common <- common[((common$DaPars<0 & common$APAlyzer=="UP") | (common$DaPars>0 & common$APAlyzer=="DN")),]
common$PDUI <- APA_Sign_Sorted$PDUI_Group_diff[match(common$Transcript, APA_Sign_Sorted$Gene)]
common$DaPars_PAdj <- APA_Sign_Sorted$adjusted.P_val[match(common$Transcript, APA_Sign_Sorted$Gene)]
common$RED <- APAlyzer_Results$RED[match(common$Gene, APAlyzer_Results$gene_symbol)]
common$APAlyzer_PAdj <- APAlyzer_Results$p_adj[match(common$Gene, APAlyzer_Results$gene_symbol)]
plot(common$PDUI, common$RED)
plot(common$PDUI, common$RED)

common <- common[order(common$PDUI, decreasing=FALSE),]
common$Transcript <- factor(common$Transcript, levels=c(common$Transcript))
common$Direction <- common$PDUI < 0
common$Direction[common$Direction] <- "Longer"
common$Direction[common$Direction==FALSE] <- "Shorter"

ggplot(common) + 
  aes(y=Transcript, x=-PDUI, fill=Direction) + 
  geom_bar(stat="identity") + 
  scale_x_continuous(limits=c(-1,1)) + 
  scale_y_discrete(limits=rev) + 
  scale_fill_manual(values=c("Shorter"="#0000FF", "Longer"="#FF0000")) + 
  theme_classic() + 
  theme(
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(), 
    legend.position = "Null"
  ) 

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/Fig_APA/", 
    "FANS_APA_Common_TDP43Neg_TDP43Pos_Barplot", 
    ".pdf"
  ),
  width = 160, 
  height = 200, 
  units="mm"
)  

table(common$PDUI<0)




  ### 8.0 Overlap FANS DE & FANS APA -------------------------------------------

FANS_DE_MAST_Sign <- FANS_DE_TDP43_MAST$gene[abs(FANS_DE_TDP43_MAST$model_log2FC) > log2(1.5) & FANS_DE_TDP43_MAST$fdr <0.1]
FANS_APA_APAlyzer_Sign <- APAlyzer_Results$gene_symbol[abs(APAlyzer_Results$RED)>1 & APAlyzer_Results$p_adj<0.0001]

table(FANS_DE_MAST_Sign %in% FANS_APA_APAlyzer_Sign)
table(FANS_APA_APAlyzer_Sign %in% FANS_DE_MAST_Sign)

table(FANS_DE_MAST_Sign %in% ALS_ALSFTD_AllCells_Signature$SYMBOL)
table(ALS_ALSFTD_AllCells_Signature$SYMBOL %in% FANS_DE_MAST_Sign)

table(FANS_DE_MAST_Sign %in% ALS_ALSFTD_WNN_L25_Signature$SYMBOL)
table(ALS_ALSFTD_WNN_L25_Signature$SYMBOL %in% FANS_DE_MAST_Sign)

fisher.test(
  matrix(
    c(889, 1024, 4582, 24000-1024-889-4582), 
    2, 2 
  )
)

APA_WNNL25_FANSDE_Common <- intersect(
  FANS_APA_APAlyzer_Sign[FANS_APA_APAlyzer_Sign %in% FANS_DE_MAST_Sign], 
  ALS_ALSFTD_WNN_L25_Signature$SYMBOL
)

APA_WNNL25_FANSDE_Common <- data.frame(
  SYMBOL = APA_WNNL25_FANSDE_Common
)
APA_WNNL25_FANSDE_Common$FANS_DE_L2FC <- FANS_DE_TDP43_MAST$model_log2FC[match(APA_WNNL25_FANSDE_Common$SYMBOL, FANS_DE_TDP43_MAST$gene)]
tmp <- as.data.frame(DESeq_Results_12SVs_Disease_AllCells_list[[1]])
tmp$ID_10X <- rownames(tmp)
tmp$SYMBOL <- GEX_Features$SYMBOL[match(tmp$ID_10X, GEX_Features$ID_10X)]
APA_WNNL25_FANSDE_Common$WNNL5_DE_L2FC <- tmp$log2FoldChange[match(APA_WNNL25_FANSDE_Common$SYMBOL, tmp$SYMBOL)]

plot(APA_WNNL25_FANSDE_Common$FANS_DE_L2FC, APA_WNNL25_FANSDE_Common$WNNL5_DE_L2FC)

  ### 9.0 Overlap FANS DE & ALS_ALSFTD Signatures ------------------------------

Overlap_FANS_DE_AllCells_DEGs <- ALS_ALSFTD_AllCells_Signature$SYMBOL[ALS_ALSFTD_AllCells_Signature$SYMBOL %in% FANS_DE_MAST_Sign]

ind <- which(
  DESeq_Results_12SVs_Disease_AllCells_Index$CellTypeLevel=="AllCells" & 
    DESeq_Results_12SVs_Disease_AllCells_Index$Comparison=="Disease"
)

tmp <- as.data.frame(DESeq_Results_12SVs_Disease_AllCells_list[[ind]])
tmp$ID_10X <- rownames(tmp)
tmp$SYMBOL <- GEX_Features$SYMBOL[match(tmp$ID_10X, GEX_Features$ID_10X)]

Overlap_FANS_DE_AllCells_DEGs <- data.frame(
  SYMBOL = Overlap_FANS_DE_AllCells_DEGs, 
  FANS_DE_L2FC = FANS_DE_TDP43_MAST$model_log2FC[match(Overlap_FANS_DE_AllCells_DEGs, FANS_DE_TDP43_MAST$gene)], 
  DE_AllCells_ALS_L2FC = tmp$log2FoldChange[match(Overlap_FANS_DE_AllCells_DEGs, tmp$SYMBOL)]
)

plot(Overlap_FANS_DE_AllCells_DEGs$FANS_DE_L2FC, Overlap_FANS_DE_AllCells_DEGs$DE_AllCells_ALS_L2FC)
cor.test(Overlap_FANS_DE_AllCells_DEGs$FANS_DE_L2FC, Overlap_FANS_DE_AllCells_DEGs$DE_AllCells_ALS_L2FC)

ind <- which(
  DESeq_Results_12SVs_Disease_WNN_L25_Index$CellTypeLevel=="WNN_L25" & 
    DESeq_Results_12SVs_Disease_AllCells_Index$Comparison=="Disease"
)


Overlap_FANS_DE_WNN_L25_DEGs_Genes <- ALS_ALSFTD_WNN_L25_Signature$SYMBOL[ALS_ALSFTD_WNN_L25_Signature$SYMBOL %in% FANS_DE_MAST_Sign]
for(i in ind){
  
  tmp <- as.data.frame(DESeq_Results_12SVs_Disease_WNN_L25_list[[i]])
  tmp$ID_10X <- rownames(tmp)
  tmp$SYMBOL <- GEX_Features$SYMBOL[match(tmp$ID_10X, GEX_Features$ID_10X)]
  
  Overlap_FANS_DE_WNN_L25_DEGs <- data.frame(
    SYMBOL = Overlap_FANS_DE_WNN_L25_DEGs_Genes, 
    FANS_DE_L2FC = FANS_DE_TDP43_MAST$model_log2FC[match(Overlap_FANS_DE_WNN_L25_DEGs_Genes, FANS_DE_TDP43_MAST$gene)], 
    DE_WNN_L25_ALS_L2FC = tmp$log2FoldChange[match(Overlap_FANS_DE_WNN_L25_DEGs_Genes, tmp$SYMBOL)]
  )
  
  print(DESeq_Results_12SVs_Disease_WNN_L25_Index$CellType[i])
  plot(Overlap_FANS_DE_WNN_L25_DEGs$FANS_DE_L2FC, Overlap_FANS_DE_WNN_L25_DEGs$DE_WNN_L25_ALS_L2FC)
  print(cor.test(Overlap_FANS_DE_WNN_L25_DEGs$FANS_DE_L2FC, Overlap_FANS_DE_WNN_L25_DEGs$DE_WNN_L25_ALS_L2FC) )
  readline("Next? ")
}
library(clusterProfiler) 
library(org.Hs.eg.db)

table(FANS_DE_WNN_L25_DEGs %in% FANS_APA_APAlyzer_Sign)

geneList <- GEX_Features_ENSEMBL_ENTREZ$ENTREZID[match(FANS_DE_WNN_L25_DEGs, GEX_Features_ENSEMBL_ENTREZ$SYMBOL)]
go_Overlap_FANS_DE_WNN_L25_DE <- enrichGO(
  geneList, org.Hs.eg.db
)
go_Overlap_FANS_DE_WNN_L25_DE <- go_Overlap_FANS_DE_WNN_L25_DE@result
go_Overlap_FANS_DE_WNN_L25_DE <- go_Overlap_FANS_DE_WNN_L25_DE[go_Overlap_FANS_DE_WNN_L25_DE$qvalue<0.05,]

kegg_Overlap_FANS_DE_WNN_L25_DE <- enrichKEGG(
  geneList, organism="hsa"
)
kegg_Overlap_FANS_DE_WNN_L25_DE <- kegg_Overlap_FANS_DE_WNN_L25_DE@result
kegg_Overlap_FANS_DE_WNN_L25_DE <- kegg_Overlap_FANS_DE_WNN_L25_DE[kegg_Overlap_FANS_DE_WNN_L25_DE$qvalue<0.05,]

# Add Universe! 
# Plot ALS/ALSFTD LOG2FCs! 
# Plot in which cell Types sign & maybe spatial 


