

  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 

library(qs)
library(tidyverse)
library(data.table)
library("spatialLIBD")
library(ExperimentHub)
library(Seurat)
library(EnsDb.Hsapiens.v86) 
library(gghalves)


qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)




  ### 1.0 Load data ------------------------------------------------------------

spe <- fetch_data(type = "spe", eh = ExperimentHub())

spe

M0_RNA <- qread(
  paste0(
    "../Data/SeuratObjects/", 
    "M0_RNA", 
    ".qrds"
  ), 
  nthr=nthr
)

libd_layer_colors2 <- libd_layer_colors
libd_layer_colors2["WM"] <- "#A4A4A4"
libd_layer_colors3 <- libd_layer_colors2
libd_layer_colors3["WM"] <- "#515151"

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




  ### 2.0 Generate a SeuratObject for the spatial data -------------------------

vis_counts <- spe@assays@data$counts
all(colnames(vis_counts)==rownames(spe@colData))
colnames(vis_counts) <- spe@colData$key

vis_logcounts <- spe@assays@data$logcounts
all(colnames(vis_logcounts)==rownames(spe@colData))
colnames(vis_logcounts) <- spe@colData$key

metadata <- spe@colData
metadata$Original_Barcode <- rownames(metadata)
rownames(metadata) <- metadata$key


Spatial_Counts <- CreateAssayObject(counts=vis_counts)
Spatial_LogCounts <- CreateAssayObject(counts=vis_logcounts)
all(colnames(Spatial_Counts)==colnames(Spatial_LogCounts))
all(rownames(Spatial_Counts)==rownames(Spatial_LogCounts))
all(colnames(Spatial_Counts)==rownames(metadata))


Visium <- CreateSeuratObject(Spatial_Counts)
all(rownames(Visium@meta.data)==rownames(metadata))

Visium@meta.data <- cbind(Visium@meta.data, metadata)
Visium@assays$SpatialCounts <- Visium@assays$RNA
Visium@assays$SpatialLogCounts <- Spatial_LogCounts
DefaultAssay(Visium) <- "SpatialLogCounts"
Visium@assays$RNA <- NULL 

rm(vis_counts, vis_logcounts, Spatial_Counts, Spatial_LogCounts, metadata) 




  ### 3.0 Generate a dictionary for ENSG IDs and HUGO Gene Symbols -------------

genes <- genes(EnsDb.Hsapiens.v86) 
features <- data.frame(
  ENSG = rownames(vis_counts)
)

grep("\\.", features$ENS)

features$Gene_Name <- genes$gene_name[match(features$ENSG, genes$gene_id)]
features$Symbol <- genes$symbol[match(features$ENSG, genes$gene_id)]

all(features$Gene_Name==features$Symbol, na.rm=TRUE)
features$Gene_biotype = genes$gene_biotype[match(features$ENSG, genes$gene_id)]
table(features$Gene_biotype)

features$Key <- features$ENSG
features$Key[!is.na(features$Symbol)] <- features$Symbol[!is.na(features$Symbol)]
length(grep("ENSG", features$Key, value=TRUE))==sum(is.na(features$Symbol))



grep("_", features$Key, value = TRUE)
features$Key[duplicated(features$Key)] |> sort()

rm(genes)




  ### 4.0 Generate CellType-specific signatures in the M0 dataset --------------

WNN_L25_Markers <- list()  
Idents(M0_RNA) <- M0_RNA$WNN_L25

for (CellType in unique(M0_RNA$WNN_L25)){
  tryCatch({
    WNN_L25_Markers[[CellType]] <- FindMarkers(
      M0_RNA, 
      ident.1 = CellType, 
      test.use="wilcox"
    ) 
  }, error=function(e){
    WNN_L25_Markers[[CellType]] <- "NA"
  }
  )
}
rm(CellType)
qsave(
  WNN_L25_Markers, 
  paste0(
    "../Data/Annotations/CellType_Markers/", 
    "WNN_L25_Markers_Wilcox", 
    ".qrds"
    ), 
  nthr=nthr
)

lapply(
  WNN_L25_Markers, 
  FUN=function(x){
    sum(x$p_val_adj < 1e-10 & x$avg_log2FC > 2.32)
  }
)




  ### 5.0 Plot GeneType Signatures ---------------------------------------------

CellTypeMarkerFeatures_ENSG <- list()
CellTypeMarkerFeatures_SYMB <- list()
DefaultAssay(Visium) <- "SpatialLogCounts"
sample="151673" 
df <- data.frame(
  WNN_L15=M0_RNA$WNN_L15,
  WNN_L25=M0_RNA$WNN_L25
)


for (CellType in names(WNN_L25_Markers)){
  fs <- WNN_L25_Markers[[CellType]] 
  fs <- rownames(fs)[fs$p_val_adj < 1e-10 & fs$avg_log2FC > 2.32]
  fs <- features$ENSG[match(fs, features$Symbol)]
  fs <- fs[!is.na(fs)]
  fs_s <- features$Symbol[match(fs, features$ENSG)] 
  
  CellTypeMarkerFeatures_ENSG[[CellType]] <- fs 
  CellTypeMarkerFeatures_SYMB[[CellType]] <- fs_s  
  
  Visium <- AddModuleScore(Visium, features=list(fs), name=paste0(CellType, "_ModuleScore_50fts_"), ctrl=50 , assay = "SpatialLogCounts")
  M0_RNA <- AddModuleScore(M0_RNA, features=list(fs_s), name=paste0(CellType, "_ModuleScore_50fts_"), ctrl=50 , assay = "RNA")
  
  spe@colData[[paste0(CellType, "_ModuleScore_50fts_1")]] <- base::scale(Visium[[paste0(CellType, "_ModuleScore_50fts_1")]])[,1]
 
  b <- vis_gene(
    spe = spe[,!is.na(spe$layer_guess_reordered)],
    sampleid = sample,
    geneid = paste0(CellType, "_ModuleScore_50fts_1"), 
    spatial = TRUE, 
    point_size = 1.2
  )
  
  b <- b + NoLegend() + theme(title=element_blank())                                                                               
  ggsave(plot = b, filename = paste0("../", CellType, "VisiumPlot.pdf"), dpi=300, width = 1200, height = 1200, units = "px") 
  ggsave(plot = b, filename = paste0("../", CellType, "VisiumPlot.tiff"), dpi=300, width = 1200, height = 1200, units = "px") 
  
  if(CellType %in% c("Astrocytes", "Oligodendrocytes", "OPC", "Microglia")){
    df$tmp <- df$WNN_L25
    df$tmp[df$tmp==CellType] <- str_replace_all(CellType, "_", " ") 
    df$tmp[df$WNN_L15=="Glia" & df$tmp!=str_replace_all(CellType, "_", " ") ] <- "Other Glia"
    df$tmp[df$WNN_L15=="Inh_Neurons"] <- "Inh Neurons"
    df$tmp[df$WNN_L15=="Exc_Neurons"] <- "Exc Neurons"
    df$tmp <- factor(df$tmp, levels=c(str_replace_all(CellType, "_", " ") , "Other Glia", "Exc Neurons", "Inh Neurons"))
  }
  
  if(CellType %in% c("Exc_RORB", "Exc_THEMIS", "Exc_FEZF2", "Exc_LINC00507")){
    df$tmp <- df$WNN_L25
    df$tmp[df$tmp==CellType] <- str_replace_all(CellType, "_", " ") 
    df$tmp[df$WNN_L15=="Exc_Neurons" & df$tmp!=str_replace_all(CellType, "_", " ")] <- "Other Exc Neurons"
    df$tmp[df$WNN_L15=="Inh_Neurons"] <- "Inh Neurons"
    df$tmp[df$WNN_L15=="Glia"] <- "Glia" 
    df$tmp <- factor(df$tmp, levels=c(str_replace_all(CellType, "_", " "), "Other Exc Neurons", "Inh Neurons", "Glia"))
  }

 
  
  if(CellType %in% c("Inh_TAFA1_VIP", "Inh_PVALB", "Inh_SST", "Inh_LAMP5_PAX6")){
    df$tmp <- df$WNN_L25
    df$tmp[df$tmp==CellType] <- str_replace_all(CellType, "_", " ")
    df$tmp[df$WNN_L15=="Inh_Neurons" & df$tmp!=str_replace_all(CellType, "_", " ")] <- "Other Inh Neurons"
    df$tmp[df$WNN_L15=="Exc_Neurons"] <- "Exc Neurons"
    df$tmp[df$WNN_L15=="Glia"] <- "Glia"
    df$tmp <- factor(df$tmp, levels=c(str_replace_all(CellType, "_", " "), "Other Inh Neurons", "Exc Neurons", "Glia"))
   }
  
  
  df$Score <- M0_RNA[[paste0(CellType, "_ModuleScore_50fts_1")]][,1] 
  fill_cols <- c(
    setNames(ColDict_WNN_L25, nm=str_replace_all(names(ColDict_WNN_L25), "_", " ")), 
    "Exc Neurons"="#AA0000", 
    "Other Exc Neurons"="#AA0000", 
    "Inh Neurons"="#0000AA", 
    "Other Inh Neurons"="#0000AA", 
    "Glia"="#00AA00", 
    "Other Glia"="#00AA00"
  )
  
  a <- ggplot(df) + 
    aes(tmp, scale(Score), fill=tmp) + 
    geom_hline(yintercept = 0, size=1, col="#00000066", linetype=2) + 
    geom_boxplot(fatten=3, col="#000000", outlier.shape = NA, ) + 
    scale_fill_manual(values = fill_cols) + 
    ylab("Gene signature score\n") + 
    theme_classic()   + 
    theme(
      axis.title.x = element_blank(), 
      axis.title.y = element_text(color="#000000", face="bold.italic"), 
      axis.text.x = element_text(color="#000000", angle=45, hjust = 1, face="bold"), 
      axis.text.y = element_text(color="#000000", face="bold.italic"), 
      legend.position = "Null"
    )
  
  ggsave(plot = a, filename = paste0("../", CellType, "BoxPlot.pdf"), dpi=300, width = 1200, height = 1200, units = "px") 
  ggsave(plot = a, filename = paste0("../", CellType, "BoxPlot.tiff"), dpi=300, width = 800, height = 1200, units = "px") 
  
  df$tmp <- "NA" 
  df$Score <- NA 
  rm (fs, fs_s, a, b)
  
  df2 <- data.frame(
    Sample = Visium$sample_id, 
    Subject = Visium$subject, 
    Position = Visium$position, 
    Layer = Visium$layer_guess_reordered, 
    Score = Visium[[paste0(CellType, "_ModuleScore_50fts_1")]][,1]
  )
  
  df2 <- df2[!is.na(df2$Layer),]
  
  c <- df2 %>% 
    group_by(Sample, Layer, Subject, Position) %>% 
    summarize(Score=mean(Score), Subject=dplyr::first(Subject), Position=dplyr::first(Position)) %>% 
    ggplot() + 
    geom_hline(yintercept = 0, size=1, col="#00000066", linetype=2) + 
    geom_boxplot(aes(Layer, Score, fill=Layer), fatten=3, col="#000000", outlier.shape=NA, show.legend = FALSE) + 
    geom_jitter(aes(Layer, Score, shape=Subject), fill="#00000044", size=3, alpha=0.6, width=0.15, col="#000000FF") + 
    scale_shape_manual(values=c(21,22,24)) +  
    scale_fill_manual(values=libd_layer_colors2) + 
    ylab("Gene signature score\n") + 
    theme_classic() + 
    theme(
      axis.text.x = element_text(color="#000000", angle=45, hjust=1, face="bold"), 
      axis.title.x = element_blank(), 
      axis.text.y = element_text(color="#000000", face = "bold.italic"), 
      axis.title.y = element_text(color="#000000", face = "bold.italic"), 
      legend.position =  "Null"
    )
  
  ggsave(plot = c, filename = paste0("../", CellType, "_Score_Layers_BoxPlot.pdf"), dpi=300, width = 1200, height = 1200, units = "px") 
  ggsave(plot = c, filename = paste0("../", CellType, "_Score_Layers_BoxPlot.tiff"), dpi=300, width = 1200, height = 1200, units = "px") 
  
  rm(df2, c)  
}










###### 
Markers_LINC00507.bkp <- Markers_LINC00507
Markers_LINC00507$Gene <- rownames(Markers_LINC00507)

Markers_LINC00507 <- Markers_LINC00507[Markers_LINC00507$p_val_adj == 0,]
Markers_LINC00507 <- Markers_LINC00507[Markers_LINC00507$avg_log2FC > 2.5,]
Markers_LINC00507$Gene <- rownames(Markers_LINC00507)

ENS_LINC00507_Markers <- features$ENSG[match(Markers_LINC00507$Gene, features$Symbol)]
table(is.na(ENS_LINC00507_Markers))
Visium <- AddModuleScore(Visium, features=list(ENS_LINC00507_Markers), name = "LINC00507_Score_Counts")
Visium <- AddModuleScore(Visium, features=list(ENS_LINC00507_Markers), assay="SpatialLogCounts", name = "LINC00507_Score_LogCounts")

VlnPlot(Visium, features = "LINC00507_Score_Counts1", group.by = "layer_guess_reordered_short", pt.size=0) + 
VlnPlot(Visium, features = "LINC00507_Score_LogCounts1", group.by = "layer_guess_reordered_short", pt.size=0)

Idents(M0_RNA) <- M0_RNA$WNN_L25
Markers_LINC00507 <- FindMarkers(M0_RNA, ident.1 = "Exc_LINC00507")
Markers_LINC00507.bkp <- Markers_LINC00507
Markers_LINC00507$Gene <- rownames(Markers_LINC00507)

Markers_LINC00507 <- Markers_LINC00507[Markers_LINC00507$p_val_adj == 0,]
Markers_LINC00507 <- Markers_LINC00507[Markers_LINC00507$avg_log2FC > 2.5,]
Markers_LINC00507$Gene <- rownames(Markers_LINC00507)

ENS_LINC00507_Markers <- features$ENSG[match(Markers_LINC00507$Gene, features$Symbol)]
table(is.na(ENS_LINC00507_Markers))
Visium <- AddModuleScore(Visium, features=list(ENS_LINC00507_Markers), name = "LINC00507_Score_Counts")
VlnPlot(Visium, features = "LINC00507_Score_Counts1", group.by = "layer_guess_reordered_short", pt.size=0)

ALS_GWAS_1 <- fread("../../../EFO_0001357_associations_export.tsv", data.table=FALSE)
ALS_GWAS_1 <- unique(unlist(str_split(ALS_GWAS_1$mappedGenes, ",")))

ALS_GWAS_1_ENSG <- features$ENSG[match(ALS_GWAS_1, features$Symbol)]
ALS_GWAS_1_ENSG <- ALS_GWAS_1_ENSG[!is.na(ALS_GWAS_1_ENSG)]


Visium <- AddModuleScore(Visium, features = list(ALS_GWAS_1_ENSG), assay="SpatialLogCounts", name="ALS_GWAS_1_LogCounts_Score")
VlnPlot(Visium, features=c("ALS_GWAS_1_Score1"), group.by="layer_guess_reordered_short", pt.size=0) + 
VlnPlot(Visium, features=c("ALS_GWAS_1_LogCounts_Score1"), group.by="layer_guess_reordered_short", pt.size=0)


ALS_GWAS_2 <- fread("../../../MONDO_0004976_associations_export.tsv", data.table=FALSE)
ALS_GWAS_2 <- unique(unlist(str_split(ALS_GWAS_2$mappedGenes, ",")))

ALS_GWAS_2_ENSG <- features$ENSG[match(ALS_GWAS_2, features$Symbol)]
ALS_GWAS_2_ENSG <- ALS_GWAS_2_ENSG[!is.na(ALS_GWAS_2_ENSG)]

Visium <- AddModuleScore(Visium, features = list(ALS_GWAS_2_ENSG), name="ALS_GWAS_2_Score")
Visium <- AddModuleScore(Visium, features = list(ALS_GWAS_2_ENSG), assay="SpatialLogCounts", name="ALS_GWAS_2_LogCounts_Score")

VlnPlot(Visium, features=c("ALS_GWAS_2_Score1"), group.by="layer_guess_reordered_short", pt.size=0) + 
VlnPlot(Visium, features=c("ALS_GWAS_2_LogCounts_Score1"), group.by="layer_guess_reordered_short", pt.size=0) 

DefaultAssay(M0_RNA)

Markers_Oligodendrocytes <- FindMarkers(M0_RNA, ident.1 = "Oligodendrocytes")
Markers_Oligodendrocytes.bkp <- Markers_Oligodendrocytes
Markers_Oligodendrocytes$Gene <- rownames(Markers_Oligodendrocytes)

Markers_Oligodendrocytes <- Markers_Oligodendrocytes[Markers_Oligodendrocytes$p_val_adj == 0,]
Markers_Oligodendrocytes <- Markers_Oligodendrocytes[Markers_Oligodendrocytes$avg_log2FC > 2.5,]
Markers_Oligodendrocytes$Gene <- rownames(Markers_Oligodendrocytes)

Oligodendrocytes_Markers <- features$ENSG[match(Markers_Oligodendrocytes$Gene, features$Symbol)]
table(is.na(Oligodendrocytes_Markers))


Visium <- AddModuleScore(Visium, features=list(Oligodendrocytes_Markers), name = "Oligodendrocytes_Score_Counts")
VlnPlot(Visium, features = "Oligodendrocytes_Score_Counts1", group.by = "layer_guess_reordered_short", pt.size=0)
DefaultAssay(Visium)

DefaultAssay(Visium) <- "SpatialLogCounts"
Visium <- AddModuleScore(Visium, features=list(Oligodendrocytes_Markers), assay = "SpatialLogCounts", name = "Oligodendrocytes_Score_LogCounts")

VlnPlot(Visium, features = "Oligodendrocytes_Score_LogCounts1", group.by = "layer_guess_reordered_short", pt.size=0)
DefaultAssay(Visium) <- "SpatialCounts"

Idents(M0_RNA) <- M0_RNA$WNN_L4
Markers_Exc_LINC00507_FREM3 <- FindMarkers(M0_RNA, ident.1 = "Exc_LINC00507_FREM3")
Markers_Exc_LINC00507_FREM3.bkp <- Markers_Exc_LINC00507_FREM3
Markers_Exc_LINC00507_FREM3$Gene <- rownames(Markers_Exc_LINC00507_FREM3)

Markers_Exc_LINC00507_FREM3 <- Markers_Exc_LINC00507_FREM3[Markers_Exc_LINC00507_FREM3$p_val_adj == 0,]
Markers_Exc_LINC00507_FREM3 <- Markers_Exc_LINC00507_FREM3[Markers_Exc_LINC00507_FREM3$avg_log2FC > 2.5,]
Markers_Exc_LINC00507_FREM3$Gene <- rownames(Markers_Exc_LINC00507_FREM3)

Exc_LINC00507_FREM3_Markers <- features$ENSG[match(Markers_Exc_LINC00507_FREM3$Gene, features$Symbol)]
table(is.na(Exc_LINC00507_FREM3_Markers))
Visium <- AddModuleScore(Visium, features=list(Exc_LINC00507_FREM3_Markers), assay="SpatialLogCounts", name = "Exc_LINC00507_FREM3_Score_LogCounts")
VlnPlot(Visium, features = "Exc_LINC00507_FREM3_Score_LogCounts1", group.by = "layer_guess_reordered_short", pt.size=0) 


Markers_Exc_RORB_ERBB4 <- FindMarkers(M0_RNA, ident.1 = "Exc_RORB_ERBB4")
Markers_Exc_RORB_ERBB4.bkp <- Markers_Exc_RORB_ERBB4
Markers_Exc_RORB_ERBB4$Gene <- rownames(Markers_Exc_RORB_ERBB4)

Markers_Exc_RORB_ERBB4 <- Markers_Exc_RORB_ERBB4[Markers_Exc_RORB_ERBB4$p_val_adj == 0,]
#Markers_Exc_RORB_ERBB4 <- Markers_Exc_RORB_ERBB4[Markers_Exc_RORB_ERBB4$avg_log2FC > 2.5,]
Markers_Exc_RORB_ERBB4$Gene <- rownames(Markers_Exc_RORB_ERBB4)

Exc_RORB_ERBB4_Markers <- features$ENSG[match(Markers_Exc_RORB_ERBB4$Gene, features$Symbol)]
table(is.na(Exc_RORB_ERBB4_Markers))
Visium <- AddModuleScore(Visium, features=list(Exc_RORB_ERBB4_Markers), assay="SpatialLogCounts", name = "Exc_RORB_ERBB4_Score_LogCounts")
VlnPlot(Visium, features = "Exc_RORB_ERBB4_Score_LogCounts1", group.by = "layer_guess_reordered_short", pt.size=0) 

Extratelencephalic_Markers <- fread("../../../Markers_Extratelencephalic_Neurons.txt", data.table = FALSE, header=FALSE)
Extratelencephalic_Markers <- Extratelencephalic_Markers$V1
Extratelencephalic_Markers_ENSG <- features$ENSG[match(Extratelencephalic_Markers, features$Symbol)]
table(is.na(Extratelencephalic_Markers))

Idents(M0_RNA) <- M0_RNA$WNN_L25
M0_RNA <- AddModuleScore(M0_RNA, features=list(Extratelencephalic_Markers), assay="RNA", name = "Extratelencephalic_Markers_Score_LogCounts")
VlnPlot(M0_RNA, features = "Extratelencephalic_Markers_Score_LogCounts1", group.by = "WNN_L4", pt.size=0, idents = c("Exc_FEZF2"))
VlnPlot(M0_RNA, features = "Extratelencephalic_Markers_Score_LogCounts1", group.by = "WNN_L4", pt.size=0) + NoLegend()
VlnPlot(M0_RNA, features = "Extratelencephalic_Markers_Score_LogCounts1", group.by = "WNN_L4", pt.size=1, idents = c("Exc_FEZF2"))



Visium <- AddModuleScore(Visium, features=list(Extratelencephalic_Markers_ENSG), assay="SpatialLogCounts", name = "Extratelencephalic_Markers_Score_LogCounts")
VlnPlot(Visium, features = "Exc_FEZF2_Score_LogCounts1", group.by = "layer_guess_reordered_short", pt.size=0) 


Idents(M0_RNA) <-  M0_RNA$WNN_L25
Markers_Exc_FEZF2 <- FindMarkers(M0_RNA, ident.1 = "Exc_FEZF2")
Markers_Exc_FEZF2.bkp <- Markers_Exc_FEZF2
Markers_Exc_FEZF2$Gene <- rownames(Markers_Exc_FEZF2)

Markers_Exc_FEZF2 <- Markers_Exc_FEZF2[Markers_Exc_FEZF2$p_val_adj == 0,]
Markers_Exc_FEZF2 <- Markers_Exc_FEZF2[Markers_Exc_FEZF2$avg_log2FC > 2.5,]
Markers_Exc_FEZF2$Gene <- rownames(Markers_Exc_FEZF2)

Exc_FEZF2_Markers <- features$ENSG[match(Markers_Exc_FEZF2$Gene, features$Symbol)]
table(is.na(Exc_FEZF2_Markers))
Visium <- AddModuleScore(Visium, features=list(Exc_FEZF2_Markers), assay="SpatialLogCounts", name = "Exc_FEZF2_Score_LogCounts")
VlnPlot(Visium, features = "Exc_FEZF2_Score_LogCounts1", group.by = "layer_guess_reordered_short", pt.size=0) 


Idents(M0_RNA) <- M0_RNA$WNN_L4
Markers_FEZF2_NTNG1 <- FindMarkers(M0_RNA, ident.1 = "Exc_FEZF2_NTNG1")
Markers_FEZF2_NTNG1.bkp <- Markers_FEZF2_NTNG1
Markers_FEZF2_NTNG1$Gene <- rownames(Markers_FEZF2_NTNG1)

Markers_FEZF2_NTNG1 <- Markers_FEZF2_NTNG1[Markers_FEZF2_NTNG1$p_val_adj == 0,]
Markers_FEZF2_NTNG1 <- Markers_FEZF2_NTNG1[Markers_FEZF2_NTNG1$avg_log2FC > 1,]
Markers_FEZF2_NTNG1$Gene <- rownames(Markers_FEZF2_NTNG1)

Exc_FEZF2_NTNG1_Markers <- features$ENSG[match(Markers_FEZF2_NTNG1$Gene, features$Symbol)]
table(is.na(Exc_FEZF2_NTNG1_Markers))
Visium <- AddModuleScore(Visium, features=list(Exc_FEZF2_NTNG1_Markers), assay="SpatialLogCounts", name = "Exc_FEZF1_NTNG1_Score_LogCounts")
VlnPlot(Visium, features = "Exc_FEZF1_NTNG1_Score_LogCounts1", group.by = "layer_guess_reordered_short", pt.size=0)

OPC_Markers <- features$ENSG[match(Markers_OPC$Gene, features$Symbol)]
table(is.na(OPC_Markers))
Visium <- AddModuleScore(Visium, features=list(OPC_Markers), name = "OPC_Score_Counts")
Visium <- AddModuleScore(Visium, features=list(OPC_Markers), assay="SpatialLogCounts", name = "OPC_Score_LogCounts")
Visium@meta.data %>% 
  group_by(sample_id, layer_guess_reordered) %>%
  summarize(Score=mean(Exc_FEZF1_NTNG1_Score_LogCounts1)) %>% 
  ggplot2::ggplot() + 
  aes(layer_guess_reordered, Score, fill=layer_guess_reordered) + 
  geom_boxplot() + 
  geom_jitter(width=0.15) 


VlnPlot(Visium, features = "OPC_Score_Counts1", group.by = "layer_guess_reordered_short", pt.size=0) + 
  VlnPlot(Visium, features = "OPC_Score_LogCounts1", group.by = "layer_guess_reordered_short", pt.size=0)

Visium@meta.data %>% 
  group_by(sample_id, layer_guess_reordered) %>%
  summarize(Score=mean(LINC00507_Score_LogCounts1)) %>% 
  ggplot2::ggplot() + 
  aes(layer_guess_reordered, Score, fill=layer_guess_reordered) + 
  geom_boxplot() + 
  geom_jitter(width=0.15)

Visium@meta.data %>% 
  group_by(sample_id, layer_guess_reordered) %>%
  summarize(Score=mean(ALS_GWAS_1_LogCounts_Score1)) %>% 
  ggplot2::ggplot() + 
  aes(layer_guess_reordered, Score, fill=layer_guess_reordered) + 
  geom_boxplot() + 
  geom_jitter(width=0.15)

Visium@meta.data %>% 
  group_by(sample_id, layer_guess_reordered) %>%
  summarize(Score=mean(ALS_GWAS_2_LogCounts_Score1)) %>% 
  ggplot2::ggplot() + 
  aes(layer_guess_reordered, Score, fill=layer_guess_reordered) + 
  geom_boxplot() + 
  geom_jitter(width=0.15)

Visium@meta.data %>% 
  group_by(sample_id, layer_guess_reordered) %>%
  summarize(Score=mean(Oligodendrocytes_Score_LogCounts1)) %>% 
  ggplot2::ggplot() + 
  aes(layer_guess_reordered, Score, fill=layer_guess_reordered) + 
  geom_boxplot() + 
  geom_jitter(width=0.15)  

Visium@meta.data %>% 
  group_by(sample_id, layer_guess_reordered) %>%
  summarize(Score=mean(Exc_FEZF2_Score_LogCounts1)) %>% 
  ggplot2::ggplot() + 
  aes(layer_guess_reordered, Score, fill=layer_guess_reordered) + 
  geom_boxplot() + 
  geom_jitter(width=0.15) 

all(Visium$Original_Barcode == colnames(spe))
all(Visium$Original_Barcode == rownames(spe@colData))
spe@colData$LINC00507_Score_LogCounts1 <- Visium$LINC00507_Score_LogCounts1
spe@colData$Oligodendrocytes_Score_LogCounts1 <- Visium$Oligodendrocytes_Score_LogCounts1
spe@colData$OPC_Score_LogCounts1 <- Visium$OPC_Score_LogCounts1
spe@colData$Exc_LINC00507_FREM3_Score_LogCounts1 <- Visium$Exc_LINC00507_FREM3_Score_LogCounts1
spe@colData$ALS_GWAS_1_LogCounts_Score1 <- Visium$ALS_GWAS_1_LogCounts_Score1
spe@colData$ALS_GWAS_2_LogCounts_Score1 <- Visium$ALS_GWAS_2_LogCounts_Score1
spe@colData$Exc_RORB_ERBB4_Score_LogCounts1 <- Visium$Exc_RORB_ERBB4_Score_LogCounts1
spe@colData$Exc_FEZF2_Score_LogCounts1 <- Visium$Exc_FEZF2_Score_LogCounts1
spe@colData$Extratelencephalic_Markers_Score_LogCounts1 <- Visium$Extratelencephalic_Markers_Score_LogCounts1
spe@colData$Exc_FEZF1_NTNG1_Score_LogCounts1 <- Visium$Exc_FEZF1_NTNG1_Score_LogCounts1

a <- vis_clus(
  spe = spe,
  clustervar = "layer_guess_reordered",
  sampleid = "151676",
  colors = libd_layer_colors,
  spatial = FALSE, 
  ... = " LIBD Layers"
) 

b <- vis_gene(
  spe = spe,
  sampleid = "151676",
  geneid = "LINC00507_Score_LogCounts1"
)

a + b




a <- vis_clus(
  spe = spe,
  clustervar = "layer_guess_reordered",
  sampleid = "151676",
  colors = libd_layer_colors,
  spatial = FALSE, 
  ... = " LIBD Layers"
) 

b <- vis_gene(
  spe = spe,
  sampleid = "151676",
  geneid = "Oligodendrocytes_Score_LogCounts1"
)

a + b

a <- vis_clus(
  spe = spe,
  clustervar = "layer_guess_reordered",
  sampleid = "151676",
  colors = libd_layer_colors,
  spatial = FALSE, 
  ... = " LIBD Layers"
) 

b <- vis_gene(
  spe = spe,
  sampleid = "151676",
  geneid = "OPC_Score_LogCounts1"
)

a + b

a <- vis_clus(
  spe = spe,
  clustervar = "layer_guess_reordered",
  sampleid = "151676",
  colors = libd_layer_colors,
  spatial = FALSE, 
  ... = " LIBD Layers"
) 

b <- vis_gene(
  spe = spe,
  sampleid = "151676",
  geneid = "Exc_LINC00507_FREM3_Score_LogCounts1"
)

a + b

a <- vis_clus(
  spe = spe,
  clustervar = "layer_guess_reordered",
  sampleid = "151676",
  colors = libd_layer_colors,
  spatial = FALSE, 
  ... = " LIBD Layers"
) 

b <- vis_gene(
  spe = spe,
  sampleid = "151676",
  geneid = "ALS_GWAS_1_LogCounts_Score1"
)

a + b


a <- vis_clus(
  spe = spe,
  clustervar = "layer_guess_reordered",
  sampleid = "151676",
  colors = libd_layer_colors,
  spatial = FALSE, 
  ... = " LIBD Layers"
) 

b <- vis_gene(
  spe = spe,
  sampleid = "151676",
  geneid = "ALS_GWAS_2_LogCounts_Score1"
)

a + b 

b <- vis_gene(
  spe = spe,
  sampleid = "151676",
  geneid = "Exc_RORB_ERBB4_Score_LogCounts1"
)

a + b

b <- vis_gene(
  spe = spe,
  sampleid = "151676",
  geneid = "Exc_FEZF2_Score_LogCounts1", 
  spatial=FALSE
)

a + b 

b <- vis_gene(
  spe = spe,
  sampleid = "151676",
  geneid = "Extratelencephalic_Markers_Score_LogCounts1", 
  spatial=FALSE
)

a + b 

b <- vis_gene(
  spe = spe,
  sampleid = "151676",
  geneid = "Exc_FEZF1_NTNG1_Score_LogCounts1", 
  spatial=FALSE
)

a + b 

spe$Exc_FEZF1_NTNG1_Score_LogCounts1_Thresholded <- spe$Exc_FEZF1_NTNG1_Score_LogCounts1
spe$Exc_FEZF1_NTNG1_Score_LogCounts1_Thresholded[spe$Exc_FEZF1_NTNG1_Score_LogCounts1_Thresholded<0.08] <- 0

b <- vis_gene(
  spe = spe,
  sampleid = "151676",
  geneid = "Exc_FEZF1_NTNG1_Score_LogCounts1_Thresholded", 
  spatial=FALSE
)

a + b

for (sample_id in unique(spe$sample_id)){

  a <- vis_clus(
    spe = spe,
    clustervar = "layer_guess_reordered",
    sampleid = sample_id,
    colors = libd_layer_colors,
    spatial = FALSE, 
    ... = " LIBD Layers"
  ) 
  
  b <- vis_gene(
    spe = spe,
    sampleid = sample_id,
    geneid = "Exc_FEZF1_NTNG1_Score_LogCounts1_Thresholded", 
    spatial=FALSE
  )
  
  plot(a + b) 
  readline("Next?")
}


for (sample_id in unique(spe$sample_id)){
  
  a <- vis_clus(
    spe = spe,
    clustervar = "layer_guess_reordered",
    sampleid = sample_id,
    colors = libd_layer_colors,
    spatial = FALSE, 
    ... = " LIBD Layers"
  ) 
  
  b <- vis_gene(
    spe = spe,
    sampleid = sample_id, 
    geneid = "Extratelencephalic_Markers_Score_LogCounts1", 
    spatial=FALSE
  )
  plot(a + b) 
  readline("Next?")
}

VAT1L = features$ENSG[which(features$Symbol=="VAT1L")]
vis_gene(
  spe = spe,
  sampleid = sample_id, 
  geneid = VAT1L, 
  spatial=FALSE
)

VlnPlot(
  M0_RNA, 
  features="VAT1L", 
  group.by = "WNN_L4", 
  pt.size=0
) + NoLegend()

VlnPlot(
  M0_RNA, 
  features="Extratelencephalic_Markers_Score_LogCounts1", 
  group.by = "WNN_L4", 
  pt.size=0
) + NoLegend()

M0_RNA <- AddModuleScore(M0_RNA, features=list(ALS_GWAS_1), assay="RNA", name="ALS_GWAS_1_Score")
M0_RNA <- AddModuleScore(M0_RNA, features=list(ALS_GWAS_2), assay="RNA", name="ALS_GWAS_2_Score")

VlnPlot(
  M0_RNA, 
  features=M0_RNA$
  group.by = "WNN_L4", 
  pt.size=0
) + NoLegend() 

VlnPlot(
  M0_RNA, 
  features="ALS_GWAS_1_Score1",
  group.by = "WNN_L4", 
  pt.size=0
) + NoLegend() 

VlnPlot(
  M0_RNA, 
  features="ALS_GWAS_2_Score1",
  group.by = "WNN_L4", 
  pt.size=0
) + NoLegend()

df <- data.frame(
  UMAP1=M0_RNA@reductions$wumap@cell.embeddings[,1], 
  UMAP2=M0_RNA@reductions$wumap@cell.embeddings[,2], 
  Case=M0_RNA$Case, 
  Case_Type=M0_RNA$Case_Type, 
  WNN_L25=M0_RNA$WNN_L25, 
  WNN_L4=M0_RNA$WNN_L4, 
  ALS_GWAS1_Score = M0_RNA$ALS_GWAS_1_Score1, 
  ALS_GWAS2_Score = M0_RNA$ALS_GWAS_2_Score1
)

ggplot(df[df$Case %in% c("HC", "ALS"),]) +
  geom_half_violin(
    aes(x = WNN_L25, y = ALS_GWAS1_Score, split = Case, fill = Case),
    position = "identity"
  )

ggplot(df[df$Case %in% c("HC", "ALS_FTD") & df$WNN_L25=="Exc_FEZF2",]) +
  geom_half_violin(
    aes(x = WNN_L4, y = ALS_GWAS2_Score, split = Case, fill = Case),
    position = "identity"
  ) + 
  theme(
    axis.text.x = element_text(angle=45, vjust=0)
  )

Results_FEZF2_NTNG1 <- qread("../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/12_SVs/DESeq_Results_12SVs_Index.qrds")
Res <- qread("../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/12_SVs/DESeq_Results_12SVs_ALS.qrds", nthr=nthr)
a <- Res[[91]]
a <- data.frame(a)
a <- a[!is.na(a$padj),]
a <- a[a$padj<0.05,]

Idents(M0_RNA) <- M0_RNA$Case
VlnPlot(M0_RNA, features="RBFOX3", group.by="WNN_L25", pt.size=1) 
VlnPlot(M0_RNA, features="RBFOX3", group.by="WNN_L25", idents = "HC", pt.size=1)  + NoLegend()
  VlnPlot(M0_RNA, features="RBFOX3", group.by="WNN_L25", idents = "ALS", pt.size=1)  + NoLegend()
  VlnPlot(M0_RNA, features="RBFOX3", group.by="WNN_L25", idents = "ALS_FTD", pt.size=1)  + NoLegend() 

M0_RNA@meta.data %>% 
  group_by(ID) %>% 
  summarize(Case=dplyr::first(Case), Cells=length(WNN_L4), CellNumber=sum(WNN_L4=="Exc_LINC00507_FREM3")) %>% 
  as.data.frame() %>% 
  ggplot() + 
  aes(Case, CellNumber/Cells, fill=Case) + 
  geom_boxplot() + geom_jitter(width=0.15)

M0_RNA@meta.data %>% 
  group_by(ID) %>% 
  summarize(Case=dplyr::first(Case), Cells=length(WNN_L4), CellNumber=sum(WNN_L4=="Exc_FEZF2_NTNG1")) %>% 
  as.data.frame() %>% 
  ggplot() + 
  aes(Case, CellNumber/Cells, fill=Case) + 
  geom_boxplot() + geom_jitter(width=0.15)

  
M0_RNA@meta.data %>% 
  group_by(ID) %>% 
  summarize(Case=dplyr::first(Case), Cells=length(WNN_L4), CellNumber=sum(WNN_L25=="Astrocytes")) %>% 
  as.data.frame() %>% 
  ggplot() + 
  aes(Case, CellNumber/Cells, fill=Case) + 
  geom_boxplot() + geom_jitter(width=0.15) 

M0_RNA@meta.data %>% 
  group_by(ID) %>% 
  summarize(Case=dplyr::first(Case), Cells=length(WNN_L4), CellNumber=sum(WNN_L25=="Oligodendrocytes")) %>% 
  as.data.frame() %>% 
  ggplot() + 
  aes(Case, CellNumber/Cells, fill=Case) + 
  geom_boxplot() + geom_jitter(width=0.15) 

M0_RNA@meta.data %>% 
  group_by(ID) %>% 
  summarize(Case=dplyr::first(Case), Cells=length(WNN_L4), CellNumber=sum(WNN_L25=="Microglia")) %>% 
  as.data.frame() %>% 
  ggplot() + 
  aes(Case, CellNumber/Cells, fill=Case) + 
  geom_boxplot() + geom_jitter(width=0.15)

M0_RNA <- AddModuleScore(M0_RNA, features=list(ALS_Genes), name="ALS_Genes_Score")
VlnPlot(M0_RNA, features="ALS_Genes_Score1", group.by="WNN_L4", idents = "HC", pt.size=0)  + NoLegend() 

################# 
vis_clus(
  spe = spe,
  clustervar = "layer_guess_reordered",
  sampleid = "151673",
  colors = libd_layer_colors,
  alpha=2, 
  spatial = TRUE, 
  point_size = 2.4,
  
  ... = " LIBD Layers"
)

JL <- qread("/home/genomics/Bioinfo_Data/RDS/SingleCell/Multiome_JL/spe_deconvoluted_WNNL25.qrds", nthr)
all(JL$key==spe$key)
spe@colData$CARD_Oligodendrocytes <-JL$Exc_RORB
vis_gene(
  spe = spe[,!is.na(spe$layer_guess_reordered)],
  sampleid = "151673",
  geneid = "CARD_Oligodendrocytes", 
  spatial = FALSE, 
  point_size = 1.2, 
  viridis = FALSE 
) + 
  scale_fill_gradientn(colours = c("lightblue", "lightyellow", "red"), na.value = "#00000000")
