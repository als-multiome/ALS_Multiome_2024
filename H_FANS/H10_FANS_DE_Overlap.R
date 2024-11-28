
  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 
 
library(qs)  

library(Seurat)
library(SeuratObject) 
library(MAST)
library(zinbwave)
library(patchwork)
library(tidyverse)




qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)




  ### 1.0 Load data ------------------------------------------------------------

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

table(FANS_DE_TDP43_MAST$gene %in% FANS_DE_TDP43_Wilcox$Gene)
table(FANS_DE_TDP43_Wilcox$Gene %in% FANS_DE_TDP43_MAST$gene) 

table(FANS_DE_TDP43_MAST$gene[FANS_DE_TDP43_MAST$fdr<0.05] %in% FANS_DE_TDP43_Wilcox$Gene[FANS_DE_TDP43_Wilcox$p_val_adj<0.05])
table(FANS_DE_TDP43_Wilcox$Gene[FANS_DE_TDP43_Wilcox$p_val_adj<0.05] %in% FANS_DE_TDP43_MAST$gene[FANS_DE_TDP43_MAST$fdr<0.05])

WNN_L25_HC_Markers <- qread(
  paste0(
    "../Data/Annotations/CellType_Markers/WNN_L25_Markers_HC_Wilcox.qrds"
  )
)

Markers <- unique(unlist(sapply(WNN_L25_HC_Markers, FUN=function(x){return(x$SYMBOL[x$p_val_adj<1e-10 & x$avg_log2FC>5])})))
table(FANS_DE_TDP43_MAST$gene[FANS_DE_TDP43_MAST$fdr<0.05] %in% Markers)
table(FANS_DE_TDP43_Wilcox$Gene[FANS_DE_TDP43_Wilcox$p_val_adj<0.05] %in% Markers)

MAST <- FANS_DE_TDP43_MAST$gene[FANS_DE_TDP43_MAST$fdr<0.05]
table(MAST %in% ALS_Signature$ID_10X)
table(ALS_Signature$ID_10X %in% MAST)

All_ALS_Signature <- All_ALS_Signature$ID_10X

table(MAST %in% All_ALS_Signature)
table(All_ALS_Signature %in% MAST)
fisher.test(
  matrix(
    c(
      924, 989, 4547, 24000-924-989-4547
    ), 2, 2 
  )
)

fisher.test(
  matrix(
    c(
      924, 989, 4547, 14000-924-989-4547
    ), 2, 2 
  )
)

fisher.test(
  matrix(
    c(
      499, 1414, 3815, 24000-499-1414-3815
    ), 2, 2 
  )
)

fisher.test(
  matrix(
    c(
      499, 1414, 3815, 14000-499-1414-3815
    ), 2, 2 
  )
)


WILCOX <- FANS_DE_TDP43_Wilcox$Gene[FANS_DE_TDP43_Wilcox$p_val_adj<0.05] 
Common <- MAST[MAST %in% WILCOX]
MAST <- MAST[!MAST %in% Common]
WILCOX <- WILCOX[!WILCOX %in% Common]

table(Markers %in% MAST)
table(Markers %in% WILCOX)
table(Markers %in% Common)

table(MAST %in% Markers)
table(WILCOX %in% Markers)
table(Common %in% Markers)

All <- unique(unlist(sapply(WNN_L25_HC_Markers, FUN = function(x){return(x$SYMBOL)})))
table(FANS_DE_TDP43_MAST$gene %in% All)
table(FANS_DE_TDP43_Wilcox$Gene %in% All)

length(unique(c(FANS_DE_TDP43_MAST$gene, FANS_DE_TDP43_Wilcox$Gene, All)))

fisher.test(
  matrix(
    c(
      0,
      985, 
      307, 
      24320-985-307
      ), 
    2, 2
  )
)

fisher.test(
  matrix(
    c(
      20,
      965, 
      5771, 
      24320-20-965-5771
    ), 
    2, 2
  )
)

fisher.test(
  matrix(
    c(
      10,
      1931, 
      1596, 
      24320-10-1931-1596
    ), 
    2, 2
  )
)

All_ALS_Signature <- qread("../Data/Annotations/Signatures/RNA/ALS_ALSFTD_AllCells_Signature.qrds")
ALS_Signature <- qread("../Data/Annotations/Signatures/RNA/ALS_ALSFTD_WNN_L25_Signature.qrds")


table()

library(clusterProfiler)
TDP43_MAST_Signature <- FANS_DE_TDP43_MAST$gene[FANS_DE_TDP43_MAST$fdr<0.05]
TDP43_MAST_Signature

geneList <- FANS_DE_TDP43_MAST$model_log2FC[FANS_DE_TDP43_MAST$fdr<0.05]
names(geneList) <- unique(GEX_Features_ENSEMBL_ENTREZ$ENTREZID[match(TDP43_MAST_Signature, GEX_Features_ENSEMBL_ENTREZ$ID_10X)])
geneList <- geneList[!is.na(geneList)]
geneList <- sort(geneList, decreasing=TRUE)

gsea_kegg <- gseKEGG(
  geneList, organism = 'hsa', eps = 0
)

ALS_Signature$In_TDP43_Signature <- ALS_Signature$ID_10X %in% TDP43_MAST_Signature
table(
  ALS_Signature$Signif_In, 
  ALS_Signature$In_TDP43_Signature
) %>% prop.table(margin=2)

enrich_kegg <- enrichKEGG(names(geneList), organism = 'hsa')

ALS_df <- get_sign_genes_df(DESeq_Results_12SVs_ALS[[1]])
ALSFTD_df <- get_sign_genes_df(DESeq_Results_12SVs_ALSFTD[[1]])

plot(
  ALS_df$log2FoldChange, 
  FANS_DE_TDP43_MAST$avg_log2FC[match(rownames(ALS_df), FANS_DE_TDP43_MAST$gene)]
) 

common <- rownames(ALS_df)[rownames(ALS_df) %in% FANS_DE_TDP43_MAST$gene[FANS_DE_TDP43_MAST$fdr<0.05]]
common <- data.frame(
  Gene=common, 
  ALS_L2FC = ALS_df$log2FoldChange[match(common, rownames(ALS_df))], 
  TDP43_L2FC = FANS_DE_TDP43_MAST$avg_log2FC[match(common, FANS_DE_TDP43_MAST$gene)]
)

common2 <- rownames(ALSFTD_df)[rownames(ALSFTD_df) %in% FANS_DE_TDP43_MAST$gene[FANS_DE_TDP43_MAST$fdr<0.05]]
common2 <- data.frame(
  Gene=common2, 
  ALSFTD_L2FC = ALSFTD_df$log2FoldChange[match(common2, rownames(ALSFTD_df))], 
  TDP43_L2FC = FANS_DE_TDP43_MAST$avg_log2FC[match(common2, FANS_DE_TDP43_MAST$gene)]
)



plot(common$ALS_L2FC, common$TDP43_L2FC)

genelist = common$Gene[common$ALS_L2FC <0 & common$TDP43_L2FC<0]
genelist2 = common$Gene[common$ALS_L2FC >0 & common$TDP43_L2FC<0]

go <- enrichKEGG(GEX_Features_ENSEMBL_ENTREZ$ENTREZID[match(genelist, GEX_Features_ENSEMBL_ENTREZ$ID_10X)], organism = 'hsa', qvalueCutoff = 0.01)
go2 <- enrichKEGG(GEX_Features_ENSEMBL_ENTREZ$ENTREZID[match(genelist2, GEX_Features_ENSEMBL_ENTREZ$ID_10X)], organism = 'hsa', qvalueCutoff = 0.01)

genelist = common$ALS_L2FC[common$ALS_L2FC <0 & common$TDP43_L2FC<0]
names(genelist) =GEX_Features_ENSEMBL_ENTREZ$ENTREZID[match(common$Gene[common$ALS_L2FC <0 & common$TDP43_L2FC<0], GEX_Features_ENSEMBL_ENTREZ$ID_10X)]
genelist <- genelist[!is.na(names(genelist))]
genelist <- sort(genelist, decreasing=TRUE)


genelist2 = common$ALS_L2FC[common$ALS_L2FC >0 & common$TDP43_L2FC<0]
names(genelist2) =GEX_Features_ENSEMBL_ENTREZ$ENTREZID[match(common$Gene[common$ALS_L2FC >0 & common$TDP43_L2FC<0], GEX_Features_ENSEMBL_ENTREZ$ID_10X)]
genelist2 <- genelist2[!is.na(names(genelist2))]
genelist2 <- sort(genelist2, decreasing=TRUE)

gse <- gseKEGG(genelist, 'hsa')
gse2 <- gseKEGG(genelist2, 'hsa', pvalueCutoff = NULL, scoreType="pos")


genelist3 = common2$ALSFTD_L2FC[common2$ALSFTD_L2FC <0 & common2$TDP43_L2FC<0]
names(genelist3) =GEX_Features_ENSEMBL_ENTREZ$ENTREZID[match(common2$Gene[common2$ALSFTD_L2FC <0 & common2$TDP43_L2FC<0], GEX_Features_ENSEMBL_ENTREZ$ID_10X)]
genelist3 <- genelis3t[!is.na(names(genelist3))]
genelist3 <- sort(genelist3, decreasing=TRUE)


genelist4 = common2$ALSFTD_L2FC[common2$ALSFTD_L2FC >0 & common2$TDP43_L2FC<0]
names(genelist4) =GEX_Features_ENSEMBL_ENTREZ$ENTREZID[match(common2$Gene[common2$ALSFTD_L2FC >0 & common2$TDP43_L2FC<0], GEX_Features_ENSEMBL_ENTREZ$ID_10X)]
genelist4 <- genelist4[!is.na(names(genelist4))]
genelist4 <- sort(genelist4, decreasing=TRUE)

gse3 <- gseKEGG(genelist3, 'hsa', pvalueCutoff = 0.5)
gse4 <- gseKEGG(genelist4, 'hsa', pvalueCutoff = NULL, scoreType="pos")
plot(
  ALSFTD_df$log2FoldChange, 
  FANS_DE_TDP43_MAST$avg_log2FC[match(rownames(ALSFTD_df), FANS_DE_TDP43_MAST$gene)]
)

#common2 <- rownames(ALSFTD_df)[rownames(ALSFTD_df) %in% FANS_DE_TDP43_MAST$gene[FANS_DE_TDP43_MAST$fdr<0.05]]
#common2 <- data.frame(
#  Gene=common2, 
#  ALSFTD_L2FC = ALSFTD_df$log2FoldChange[match(common2, rownames(ALSFTD_df))], 
#  TDP43_L2FC = FANS_DE_TDP43_MAST$avg_log2FC[match(common2, FANS_DE_TDP43_MAST$gene)]
#)


plot(common$ALS_L2FC, common$TDP43_L2FC)
table(APA_Genes$Symbol %in% common$Gene)
common$APA <- common$Gene %in% APA_Genes$Symbol
ggplot(common) + 
  aes(ALS_L2FC, TDP43_L2FC, fill=APA) + 
  geom_point(pch=21, size=2)


  ### 1.0 AAAAAA ---------------------------------------------------------------

TDP43_Signature_NoModel <- qread(
  paste0(
    "../Data/DE/FANS/", 
    "FANS_DE_TDP43_Wilcox", 
    ".qrds"
  ), 
  nthr=nthr
)
#Here instead of TDP43_Signature_NoModel, use Joao's Pseudobulk results (can't do because no samples?)
TDP43_Signature_NoModel <- TDP43_Signature_NoModel[TDP43_Signature_NoModel$p_val_adj<0.05 & TDP43_Signature_NoModel$avg_log2FC>0.5,]
markers <- sort(unique(as.character(unlist(sapply(WNN_L25_HC_Markers, FUN=function(x){return(x$SYMBOL[x$p_val_adj<1e-5 & x$avg_log2FC>2])})))))

prop.table(table(TDP43_Signature_NoModel$Gene %in% markers))
prop.table(table(TDP43_Signature$gene %in% markers))

table(ALS_Signature$Gene %in% FANS_DE_TDP43_MAST$gene)
table(FANS_DE_TDP43_MAST$gene %in% ALS_Signature$Gene)

TDP43_Signature <- FANS_DE_TDP43_MAST[FANS_DE_TDP43_MAST$fdr<0.1 & abs(FANS_DE_TDP43_MAST$avg_log2FC) > 0.5,]

table(ALS_Signature$Gene %in% TDP43_Signature$gene)
table(TDP43_Signature$gene %in% ALS_Signature$Gene)

fisher.test(
  matrix(
    c(
      924, 4547, 989, 35367-924-4547-989
    ), 
    2, 2
    
  )
)

genelist <- ATAC_TDP43_RNA

## Here extract the common genes and calculate with MAST the fold change with 
## the same formula as TDP43 MAST with all cells to get one fold change per gene 


ALS_Signature.bkp <- ALS_Signature
ALS_Signature <- ALS_Signature.bkp 


ALS_Signature <- ALS_Signature[ALS_Signature$QVAL_Exc_RORB_ALSFTD<0.05,]
TDP43_Signature <- FANS_DE_TDP43_MAST[FANS_DE_TDP43_MAST$fdr<0.05 & abs(FANS_DE_TDP43_MAST$avg_log2FC) > 0.5,] 


common <- ALS_Signature$Gene[ALS_Signature$Gene %in% TDP43_Signature$gene]
common <- data.frame(
  Gene = common, 
  LOG2FC_ALS_Signature = ALS_Signature$LOG2FC_Exc_RORB_ALSFTD[match(common, ALS_ALSFTD_WNN_L25_Signature$Gene)], 
  LOG2FC_TDP43_Signature = TDP43_Signature$model_log2FC[match(common, TDP43_Signature$gene)]
)

plot(common$LOG2FC_TDP43_Signature, common$LOG2FC_ALS_Signature)
cor.test(common$LOG2FC_TDP43_Signature, common$LOG2FC_ALS_Signature)

a <- common$Gene[common$LOG2FC_ALS_Signature<0 & common$LOG2FC_TDP43_Signature<0]
a <- a[!is.na(a)]
a <- unique(a)
write.csv(a, "../a.csv", quote = FALSE)
table(F2$Predicted_ID_M0_WNN_L25, F2$TDP43)

df2 <- df[df$SIGN_IN_BOTH,]
table(TDP43_Signature$gene %in% df$Gene)
df$Gene[df$Gene %in% TDP43_Signature$gene] %in% a

F2_LINC00507 <- subset(F2, subset=predicted.id=="Exc_LINC00507")

F2_OLG <- subset(F2, subset=predicted.id=="Oligodendrocytes")

RORB <- data.table::fread("../RORB_ALSFTD_Common.csv", data.table=FALSE)
table(RORB$Gene %in% ALS_Sign$Gene)
RORB$Gene[!RORB$Gene %in% ALS_Sign$Gene]
table(RORB$Gene %in% Res_MAST_Manual$gene)
table(RORB$Gene %in% Res_MAST_Manual$gene[Res_MAST_Manual$fdr<0.1])
ALS_Sign$InTDP43Sign <- FALSE
ALS_Sign$InTDP43Sign[ALS_Sign$Gene %in% Res_MAST_Manual$gene[Res_MAST_Manual$fdr<0.1]] <- TRUE
table(ALS_Sign$InTDP43Sign)

ALS_Sign$InATAC <- FALSE
ALS_Sign$InATAC[ALS_Sign$Gene %in% RORB$Gene] <- TRUE
table(ALS_Sign$InATAC)
ALS_Sign$InTDP43Sign2 <- "TDP43 No"
ALS_Sign$InTDP43Sign2[ALS_Sign$InTDP43Sign] <- "TDP43 Yes"

table(ALS_Sign$InTDP43Sign2, ALS_Sign$InATAC)
fisher.test(matrix(c(4312,235,833,91), 2, 2))
TDP_Linc_Markers <- Seurat::FindMarkers(F2_LINC00507, ident.1 = "TDP43_Low", ident.2 = "TDP43_High", test.use="MAST", latent.vars = c("Sex"))
TDP_OLG_Markers <- Seurat::FindMarkers(F2_OLG, ident.1 = "TDP43_Low", ident.2 = "TDP43_High", test.use="MAST", latent.vars = c("Sex"))

TDP_Linc_Markers$Gene <- rownames(TDP_Linc_Markers )
View(TDP_Linc_Markers)

TDP_OLG_Markers$Gene <- rownames(TDP_OLG_Markers) 
View(TDP_OLG_Markers)

table(
  TDP_Linc_Markers$Gene[TDP_Linc_Markers$p_val_adj<0.001] %in% 
    TDP_OLG_Markers$Gene[TDP_OLG_Markers$p_val_adj<0.05]
)

table(Markers_TDP43_Sex_CellType$Gene[Markers_TDP43_Sex_CellType$p_val_adj<0.0001] %in% 
        Res_MAST_Manual$gene[Res_MAST_Manual$fdr<0.1])

table(Res_MAST_Manual$gene[Res_MAST_Manual$fdr<0.1] %in% Markers_TDP43_Sex_CellType$Gene[Markers_TDP43_Sex_CellType$p_val_adj<0.0001]
        )

table(Res_MAST_Manual$gene %in% Markers_TDP43_Sex_CellType$Gene)

library(xlsx)
library(readxl)
tmp <- readxl::excel_sheets("../Data/Input/TDP43_Cao_and_Scotter/TableS2.xlsx")[7]
Cao_TDPNeg_Down <- readxl::read_xlsx("../Data/Input/TDP43_Cao_and_Scotter/TableS2.xlsx", tmp)
summary(Cao_TDPNeg_Down$padj)
colnames(Cao_TDPNeg_Down)[1] <- "Gene"
table(Cao_TDPNeg_Down$Gene %in% Markers_TDP43_Sex_CellType$Gene)
table(Cao_TDPNeg_Down$Gene %in% Markers_TDP43$Gene)
table(Cao_TDPNeg_Down$Gene %in% Res_MAST_Manual$gene[Res_MAST_Manual$fdr<0.1])

tmp <- readxl::excel_sheets("../Data/Input/TDP43_Cao_and_Scotter/TableS3.xlsx")[7]
Cao_TDPNeg_Up <- readxl::read_xlsx("../Data/Input/TDP43_Cao_and_Scotter/TableS3.xlsx", tmp)
colnames(Cao_TDPNeg_Up)[1] <- "Gene" 
summary(Cao_TDPNeg_Down$padj)
table(Cao_TDPNeg_Up$Gene %in% Markers_TDP43_Sex_CellType$Gene)
table(Cao_TDPNeg_Up$Gene %in% Res_MAST_Manual$gene[Res_MAST_Manual$fdr<0.1])
table(Res_MAST_Manual$fdr<0.1)

BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
tx_by_gene <- transcriptsBy(txdb, by="gene")
gene_lens <- max(width(tx_by_gene))
gene_lens <- data.frame(
  ENSG = paste0("ENSG", names(gene_lens)), 
  lgn = as.numeric(gene_lens)
)

library(EnsDb.Hsapiens.v86)
genes <- genes(EnsDb.Hsapiens.v86) 
width(genes)
genes
features <- data.frame(
  ENSG = rownames(LIBD_Spatial_Seurat)
)

grep("\\.", features$ENS)

features$Gene_Name <- genes$gene_name[match(features$ENSG, genes$gene_id)]
features$Symbol <- genes$symbol[match(features$ENSG, genes$gene_id)]

all(features$Gene_Name==features$Symbol, na.rm=TRUE)
features$Gene_biotype = genes$gene_biotype[match(features$ENSG, genes$gene_id)]

features$Key <- features$ENSG
features$Key[!is.na(features$Symbol)] <- features$Symbol[!is.na(features$Symbol)]
length(grep("ENSG", features$Key, value=TRUE))==sum(is.na(features$Symbol))

grep("_", features$Key, value = TRUE)
features$Key[duplicated(features$Key)] |> sort()

Res_MAST_Manual$GeneLength <- width(genes)[match(Res_MAST_Manual$gene, genes$symbol)]
table(is.na(Res_MAST_Manual$GeneLength))
table(is.na(Res_MAST_Manual$GeneLength), Res_MAST_Manual$fdr<0.1)
Res_MAST_Manual$Sign <- "NotSign"
Res_MAST_Manual$Sign[Res_MAST_Manual$fdr < 0.1] <- "Sign"
table(Res_MAST_Manual$Sign)
boxplot(Res_MAST_Manual$GeneLength ~ Res_MAST_Manual$Sign, log="y")
wilcox.test(Res_MAST_Manual$GeneLength ~ Res_MAST_Manual$Sign, log="y")

TDP43_iCLIP <- readxl::read_xlsx("../Data/Input/TDP43_iCLIP_Tollervey/SupplTable3.xlsx", col_names = FALSE)
colnames(TDP43_iCLIP)[1] <- "Gene"
table(TDP43_iCLIP$Gene %in% Res_MAST_Manual$gene[Res_MAST_Manual$Sign])
table(TDP43_iCLIP$Gene %in% Markers_TDP43_Sex_CellType$Gene[Markers_TDP43_Sex_CellType$p_val_adj<0.0001])
table(TDP43_iCLIP$Gene %in% Markers_TDP43_Sex_CellType$Gene)
table(Markers_TDP43_Sex_CellType$Gene[Markers_TDP43_Sex_CellType$p_val_adj<0.0001] %in% TDP43_iCLIP$Gene)

AllM0Genes <- Features(M0_RNA)
AllM0Genes <- data.frame(
  Gene=AllM0Genes
)
AllM0Genes$GeneLength <- width(genes)[match(AllM0Genes$Gene, genes$symbol)]
AllM0Genes$TDP_MAST_Sign <- AllM0Genes$Gene %in% Res_MAST_Manual$gene[Res_MAST_Manual$fdr<0.1]
boxplot(AllM0Genes$GeneLength ~ AllM0Genes$TDP_MAST_Sign, log="y")

AllFANSGenes <- data.frame(
  Gene = Features(F2)
)
AllFANSGenes$GeneLength <- width(genes)[match(AllFANSGenes$Gene, genes$symbol)]
AllFANSGenes$TDP_MAST <- AllFANSGenes$Gene %in% Res_MAST_Manual$gene
AllFANSGenes$TDP_MAST_Sign <- AllFANSGenes$Gene %in% Res_MAST_Manual$gene[Res_MAST_Manual$fdr<0.1]
AllFANSGenes$Rand <- AllFANSGenes$Gene %in% (AllFANSGenes$Gene[sample(1:nrow(AllFANSGenes), size=table(AllFANSGenes$TDP_MAST_Sign)["TRUE"]
, replace = FALSE)])
boxplot(AllFANSGenes$GeneLength ~ AllFANSGenes$TDP_MAST, log="y")
boxplot(AllFANSGenes$GeneLength ~ AllFANSGenes$TDP_MAST_Sign, log="y")
boxplot(AllFANSGenes$GeneLength ~ AllFANSGenes$Rand, log="y")

AllFANSGenes$TDP_Markers <- AllFANSGenes$Gene %in% Markers_TDP43_Sex_CellType$Gene
AllFANSGenes$TDP_Markers_Sign <- AllFANSGenes$Gene %in% Markers_TDP43_Sex_CellType$Gene[Markers_TDP43_Sex_CellType$p_val_adj<0.0001]
boxplot(AllFANSGenes$GeneLength ~ AllFANSGenes$TDP_Markers, log="y")
boxplot(AllFANSGenes$GeneLength ~ AllFANSGenes$TDP_Markers_Sign, log="y")

Res_Des_DEG <- qread("../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/12_SVs/DESeq_Results_12SVs_ALSFTD.qrds", nthr=nthr)
res <- Res_Des_DEG[[1]]
AllFANSGenes$DEG_DES_PSdo_All <- AllFANSGenes$Gene %in% rownames(res)
AllFANSGenes$DEG_DES_PSdo_All_Measures <- AllFANSGenes$Gene %in% rownames(res)[!is.na(res$padj)]
res <- res[!is.na(res$padj),]
AllFANSGenes$DEG_DES_PSdo_Sign <- AllFANSGenes$Gene %in% rownames(res)[res$padj < 0.05]
boxplot(AllFANSGenes$GeneLength ~ AllFANSGenes$DEG_DES_PSdo_All, log="y")
table(AllFANSGenes$DEG_DES_PSdo_All)
boxplot(AllFANSGenes$GeneLength ~ AllFANSGenes$DEG_DES_PSdo_All_Measures, log="y")
boxplot(AllFANSGenes$GeneLength ~ AllFANSGenes$DEG_DES_PSdo_Sign, log="y")

res$Sign <- res$padj < 0.05
boxplot(res$baseMean ~ res$Sign, log="y")
res$Gene <- rownames(res)
res$GeneLength <- width(genes)[match(res$Gene, genes$symbol)]
plot(res$baseMean ~ res$GeneLength, log="xy")

cor.test(res$baseMean, res$GeneLength, method="spearman")

res_Bulk_Julia_PBMCs <- data.table::fread("../../Bulk_ATAC_ALS_PBMCs/RNA_ALS_vs_HC_none_0p01.txt")
cor.test(res_Bulk_Julia_PBMCs$AveExpr, res_Bulk_Julia_PBMCs$adj.P.Val)
res_Bulk_Julia_PBMCs$Sign <- res_Bulk_Julia_PBMCs$adj.P.Val<0.05
res_Bulk_Julia_PBMCs$GeneLength <- width(genes)[match(res_Bulk_Julia_PBMCs$ENSEMBL, genes$gene_id)]
cor.test(res_Bulk_Julia_PBMCs$GeneLength, res_Bulk_Julia_PBMCs$AveExpr)
cor.test(res_Bulk_Julia_PBMCs$GeneLength, res_Bulk_Julia_PBMCs$P.Value)
boxplot(res_Bulk_Julia_PBMCs$GeneLength ~ res_Bulk_Julia_PBMCs$Sign, log="y")

res_Bulk_Julia_PBMCs2 <- readRDS("../../Bulk_ATAC_ALS_PBMCs/res_ATAC.rds")
res_Bulk_Julia_PBMCs2 <- as.data.frame(res_Bulk_Julia_PBMCs2)
cor.test(res_Bulk_Julia_PBMCs$baseMean, res_Bulk_Julia_PBMCs$padj)
res_Bulk_Julia_PBMCs2$ENSEMBL <- res_Bulk_Julia_PBMCs$ENSEMBL[match(rowames(res_Bulk_Julia_PBMCs2), res_Bulk_Julia_PBMCs$)]
res_Bulk_Julia_PBMCs$Sign <- res_Bulk_Julia_PBMCs$padj <0.05
res
res_Bulk_Julia_PBMCs$GeneLength <- width(genes)[match(res_Bulk_Julia_PBMCs$, genes$symbol)]
cor.test(res_Bulk_Julia_PBMCs$GeneLength, res_Bulk_Julia_PBMCs$baseMean)
cor.test(res_Bulk_Julia_PBMCs$GeneLength, res_Bulk_Julia_PBMCs$P.Value)
boxplot(res_Bulk_Julia_PBMCs$GeneLength ~ res_Bulk_Julia_PBMCs$Sign, log="y")


####
M0_RNA <- AddModuleScore(
  M0_RNA, 
  features=list(Res_MAST_Manual$gene[Res_MAST_Manual$fdr<0.1 & Res_MAST_Manual$model_log2FC>0]), 
  name="TDP43LowScore_Up"
)

VlnPlot(M0_RNA, feature=c("TDP43LowScore_Up1"), group.by="Case", pt.size=0)
write.csv(Res_MAST_Manual[,c(1:8)], "../Res_TDP43_MAST_WieInLi.csv", quote=FALSE)
ggplot(M0_RNA@meta.data$TDP)
TDP_Linc_Markers$Gene[TDP_Linc_Markers$p_val_adj<0.001][TDP_Linc_Markers$Gene[TDP_Linc_Markers$p_val_adj<0.001] %in% 
  TDP_OLG_Markers$Gene[TDP_OLG_Markers$p_val_adj<0.05]]

View(Markers_TDP43_Sex_CellType) 
summary(Markers_TDP43_Sex_CellType$avg_log2FC)
table(Markers_TDP43_Sex_CellType$p_val_adj<0.001) 

M0_RNA <- AddModuleScore(M0_RNA, features=list(
  Markers_TDP43_Sex_CellType$Gene[
    Markers_TDP43_Sex_CellType$p_val_adj<0.001 & 
    Markers_TDP43_Sex_CellType$avg_log2FC> 0.5
    ]), 
                         name = "Score_TDP43_Markeres_Sex_CellType2_")

M0_RNA <- AddModuleScore(M0_RNA, features=list(
  TDP_Linc_Markers$Gene[
    TDP_Linc_Markers$p_val_adj<0.001 & 
      TDP_Linc_Markers$avg_log2FC> 0.5
  ]), 
  name = "TDP_Linc_Markers_") 

M0_RNA <- AddModuleScore(M0_RNA, features=list(
  TDP_Linc_Markers$Gene[
    TDP_Linc_Markers$p_val_adj<0.001 & 
      TDP_Linc_Markers$avg_log2FC < 0
  ]), 
  name = "TDP_Linc_Markers_Down_")

VlnPlot(M0_RNA, features=c("Score_TDP43_Markeres_Sex_CellType_1"), group.by="WNN_L25", pt.size=0)
VlnPlot(M0_RNA, features=c("Score_TDP43_Markeres_Sex_CellType_1"), group.by="Case", pt.size=0)
VlnPlot(M0_RNA, features=c("Score_TDP43_Markeres_Sex_CellType2_1"), group.by="Case", pt.size=0)
M0_
VlnPlot(M0_RNA, features=c("TDP_Linc_Markers_1"), group.by="Case", pt.size=0)
VlnPlot(M0_RNA, features=c("TDP_Linc_Markers_1"), group.by="WNN_L25", pt.size=0)

library(gghalves)
df = M0_RNA@meta.data
ggplot(df[df$Case %in% c("ALS_FTD", "HC"),]) +
  geom_half_violin(
    aes(x = WNN_L25, y = TDP_Linc_Markers_1, split = Case, fill = Case),
    position = "identity"
  )

ggplot(df[df$Case %in% c("ALS_FTD", "HC"),]) +
  geom_half_violin(
    aes(x = WNN_L25, y = TDP_Linc_Markers_Down_1, split = Case, fill = Case),
    position = "identity"
  )


############ ############ M0_RNA.anchors


FC3 <- NormalizeData(FC3)
FC3 <- ScaleData(FC3)
DefaultAssay(M0) <- "RNA"
M0 <- ScaleData(M0)
VariableFeatures(M0)
M0 <- RunPCA(M0, features = VariableFeatures(M0))

M0.anchors <- FindTransferAnchors(reference = M0, query = FC3, dims = 1:30,
                                  reference.reduction = "pca")

predictions <- TransferData(anchorset =M0.anchors, refdata = M0$RNA_L4, dims = 1:30)
FC3 <- AddMetaData(FC3, metadata = predictions)

M0.anchors_FC3 = M0.anchors 
rm(M0.anchors, predictions)

FC2 <- Read10X(
  "/home/veselin/NAS/Bioinformatic_Data/Cellranger_Output/Human/scRNA/scRNA_FACS_Cortex_Sarah_210623/FC2_cellranger/outs/filtered_feature_bc_matrix/"
)
FC2 <- CreateSeuratObject(FC2)
FC2 <- NormalizeData(FC2)
FC2 <- ScaleData(FC2)

M0.anchors <- FindTransferAnchors(reference = M0, query = FC2, dims = 1:30,
                                  reference.reduction = "pca")

predictions <- TransferData(anchorset =M0.anchors, refdata = M0$RNA_L4, dims = 1:30)
FC2 <- AddMetaData(FC2, metadata = predictions) 

M0.anchors_FC2 = M0.anchors 
rm(M0.anchors, predictions)


FC1 <- Read10X(
  "/home/veselin/NAS/Bioinformatic_Data/Cellranger_Output/Human/scRNA/scRNA_FACS_Cortex_Sarah_210623/Re_Run/FC1/outs/filtered_feature_bc_matrix/"
)
FC1 <- CreateSeuratObject(FC1)
FC1 <- NormalizeData(FC1)
FC1 <- ScaleData(FC1)

M0.anchors <- FindTransferAnchors(reference = M0, query = FC1, dims = 1:30,
                                  reference.reduction = "pca")

predictions <- TransferData(anchorset = M0.anchors, refdata = M0$RNA_L4, dims = 1:30, k.weight = 30)
FC1 <- AddMetaData(FC1, metadata = predictions) 

M0.anchors_FC1 = M0.anchors 
rm(M0.anchors, predictions) 

sort(round(prop.table(table(FC1$predicted.id))*100, 1), decreasing=TRUE)
sort(round(prop.table(table(FC2$predicted.id))*100, 1), decreasing=TRUE)
sort(round(prop.table(table(FC3$predicted.id))*100, 1), decreasing=TRUE)

predictions <- TransferData(anchorset = M0.anchors_FC1, refdata = M0$RNA_L5, dims = 1:30, k.weight = 30)
FC1_L5 <- AddMetaData(FC1, metadata = predictions) 

predictions <- TransferData(anchorset = M0.anchors_FC2, refdata = M0$RNA_L5, dims = 1:30, k.weight = 30)
FC2_L5 <- AddMetaData(FC2, metadata = predictions) 

predictions <- TransferData(anchorset = M0.anchors_FC3, refdata = M0$RNA_L5, dims = 1:30, k.weight = 30)
FC3_L5 <- AddMetaData(FC3, metadata = predictions) 

sort(round(prop.table(table(FC1_L5$predicted.id))*100, 1), decreasing=TRUE)
sort(round(prop.table(table(FC2_L5$predicted.id))*100, 1), decreasing=TRUE)
sort(round(prop.table(table(FC3_L5$predicted.id))*100, 1), decreasing=TRUE)

APA_Genes <- data.table::fread("../Data/APA/DaPars1_TDP43Pos_vs_TDP43Neg/APA_Sign_Results.csv", data.table=FALSE)

LINC_TDP_Sign <-  TDP_Linc_Markers[TDP_Linc_Markers$p_val_adj<0.001,]
LINC_TDP_Sign <- LINC_TDP_Sign$Gene[LINC_TDP_Sign$avg_log2FC>0]
write.csv(LINC_TDP_Sign, "../LINC_TDP_Sign.csv", quote=FALSE)
table(F2$predicted.id, F2$TDP43)

APA_EA <- LINC_TDP_Sign[LINC_TDP_Sign %in% APA_Genes$Symbol]
grep("NRXN3", APA_EA)

M0_RNA <- AddModuleScore(M0_RNA, features=list(APA_EA), ctrl=500, name="APA_EA_Score_")
VlnPlot(M0_RNA, features=c("APA_EA_Score_1"), group.by="WNN_L25", pt.size=0) 

df <- M0_RNA@meta.data
ggplot(df[df$Case %in% c("ALS", "HC") & df$WNN_L15=="Exc_Neurons",]) +
  geom_half_violin(
    aes(x = WNN_L4, y = APA_EA_Score_1, split = Case, fill = Case),
    position = "identity"
  ) + theme(
    axis.text.x = element_text(angle=45, hjust=1)
  )



