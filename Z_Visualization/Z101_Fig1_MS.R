  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 

library(tidyverse)
library(data.table)
library(qs) 
library(Seurat)
library(Signac)
library(ggnewscale)
library(patchwork)
library(monocle3)
library(cicero) 
library(ggrastr)
library(scattermore)
library(pbapply)

qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)




  ### 1.0 Load data ------------------------------------------------------------

M0 <- qread(
  paste0(
    "../Data/SeuratObjects/",  
    "M0_Full",
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


ColDict_ATAC_L4 <- setNames(
  object = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "ATAC_L4")$Color, 
  nm = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "ATAC_L4")$ATAC_L4
)

setequal(names(ColDict_ATAC_L4), unique(M0$ATAC_L4))
any(duplicated(ColDict_ATAC_L4))
any(duplicated(names(ColDict_ATAC_L4)))


ColDict_RNA_L4 <- setNames(
  object = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "RNA_L4")$Color, 
  nm = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "RNA_L4")$RNA_L4
)

setequal(names(ColDict_RNA_L4), unique(M0$RNA_L4))
any(duplicated(ColDict_RNA_L4))
any(duplicated(names(ColDict_RNA_L4)))


ColDict_WNN_L2 <- setNames(
  object = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "WNN_L2")$Color, 
  nm = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "WNN_L2")$WNN_L2
)

setequal(names(ColDict_WNN_L2), unique(M0$WNN_L2))
setdiff(names(ColDict_WNN_L2), unique(M0$WNN_L2))
any(duplicated(ColDict_WNN_L2))
any(duplicated(names(ColDict_WNN_L2)))


ColDict_WNN_L4 <- setNames(
  object = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "WNN_L4")$Color, 
  nm = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "WNN_L4")$WNN_L4
)

setequal(names(ColDict_WNN_L4), unique(M0$WNN_L4))
setdiff(names(ColDict_WNN_L4), unique(M0$WNN_L4))
any(duplicated(ColDict_WNN_L4))
any(duplicated(names(ColDict_WNN_L4)))


ColDict_Case_Type <- setNames(
  object = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "Case_Type")$Color, 
  nm = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "Case_Type")$Case_Type
)




  ### 2.0 UMAP RNA -------------------------------------------------------------

RNA_Umap <- data.frame(
  CellId = M0$CellId, 
  UMAP1 = M0@reductions$umap@cell.embeddings[,1], 
  UMAP2 = M0@reductions$umap@cell.embeddings[,2], 
  RNA_L1 = as.character(M0$RNA_L1), 
  RNA_L15 = as.character(M0$RNA_L15), 
  RNA_L2 = as.character(M0$RNA_L2), 
  RNA_L25 = as.character(M0$RNA_L25), 
  RNA_L3 = as.character(M0$RNA_L3), 
  RNA_L4 = as.character(M0$RNA_L4) 
) 


ggplot(RNA_Umap) + 
  aes(UMAP1, UMAP2, col=RNA_L4) + 
  geom_point(size=0.8, alpha=0.6) + 
  scale_color_manual(values=ColDict_RNA_L4) + 
  theme_void() + 
  NoLegend() 

ggsave(
  filename = paste0(
    "../Data/Visualization/Figures/Fig1/", 
    "UMAP_RNA_L4", 
    ".tiff"
  ),
  ggplot(RNA_Umap) + 
    aes(UMAP1, UMAP2, col=RNA_L4) + 
    geom_point(size=0.1, alpha=0.6) + 
    scale_color_manual(values=ColDict_RNA_L4) + 
    theme_void() + 
    NoLegend(), 
  dpi=300, width = 700, height = 700, units = "px"
)

ggsave(
  filename = paste0(
    "../Data/Visualization/Figures/Fig1/", 
    "UMAP_RNA_L4", 
    ".jpg"
  ),
  ggplot(RNA_Umap) + 
    aes(UMAP1, UMAP2, col=RNA_L4) + 
    geom_point(size=0.1, alpha=0.6) + 
    scale_color_manual(values=ColDict_RNA_L4) + 
    theme_void() + 
    NoLegend(), 
  dpi=300, width = 700, height = 700, units = "px"
)




  ### 3.0 UMAP ATAC ------------------------------------------------------------

ATAC_Umap <- data.frame(
  CellId = M0$CellId, 
  UMAP1 = M0@reductions$aumap@cell.embeddings[,1], 
  UMAP2 = M0@reductions$aumap@cell.embeddings[,2], 
  ATAC_L1 = as.character(M0$ATAC_L1), 
  ATAC_L15 = as.character(M0$ATAC_L15), 
  ATAC_L2 = as.character(M0$ATAC_L2), 
  ATAC_L25 = as.character(M0$ATAC_L25), 
  ATAC_L3 = as.character(M0$ATAC_L3), 
  ATAC_L4 = as.character(M0$ATAC_L4) 
) 

ggplot(ATAC_Umap) + 
  aes(UMAP1, UMAP2, col=ATAC_L4) + 
  geom_point(size=0.8, alpha=0.6) + 
  scale_color_manual(values=ColDict_ATAC_L4) + 
  theme_void() + 
  NoLegend() 

ggsave(
  filename = paste0(
    "../Data/Visualization/Figures/Fig1/", 
    "UMAP_ATAC_L4", 
    ".tiff"
  ),
  ggplot(ATAC_Umap) + 
    aes(UMAP1, UMAP2, col=ATAC_L4) + 
    geom_point(size=0.1, alpha=0.6) + 
    scale_color_manual(values=ColDict_ATAC_L4) + 
    theme_void() + 
    NoLegend(), 
  dpi=300, width = 700, height = 700, units = "px"
  
)

ggsave(
  filename = paste0(
    "../Data/Visualization/Figures/Fig1/", 
    "UMAP_ATAC_L4", 
    ".jpg"
  ),
  ggplot(ATAC_Umap) + 
    aes(UMAP1, UMAP2, col=ATAC_L4) + 
    geom_point(size=0.1, alpha=0.6) + 
    scale_color_manual(values=ColDict_ATAC_L4) + 
    theme_void() + 
    NoLegend(), 
  dpi=300, width = 700, height = 700, units = "px"
  
)




  ### 4.0 UMAP WNN -------------------------------------------------------------

WNN_Umap <- data.frame(
  CellId = M0_RNA$CellId, 
  UMAP1 = M0_RNA@reductions$wumap@cell.embeddings[,1], 
  UMAP2 = M0_RNA@reductions$wumap@cell.embeddings[,2], 
  WNN_L1 = as.character(M0_RNA$WNN_L1), 
  WNN_L15 = as.character(M0_RNA$WNN_L15), 
  WNN_L2 = as.character(M0_RNA$WNN_L2), 
  WNN_L25 = as.character(M0_RNA$WNN_L25), 
  WNN_L3 = as.character(M0_RNA$WNN_L3), 
  WNN_L4 = as.character(M0_RNA$WNN_L4), 
) 


table(M0$WNN_L2)
table(as.character(M0$WNN_L4[M0$WNN_L2=="Astrocytes"]))
table(as.character(M0$WNN_L4[M0$WNN_L2=="Microglia"]))
table(as.character(M0$WNN_L4[M0$WNN_L2=="Oligodendrocytes"]))
table(as.character(M0$WNN_L4[M0$WNN_L2=="OPC"])) 
table(as.character(M0$WNN_L4[M0$WNN_L2=="Exc_Neurons"]))
table(as.character(M0$WNN_L4[M0$WNN_L2=="Inh_Neurons"]))




  ### 5.0 UMAP Case Type -------------------------------------------------------

WNN_Umap$Case_Type <- M0_RNA$Case_Type

ggplot(WNN_Umap[WNN_Umap$Case_Type=="ALS",]) + 
  aes(UMAP1, UMAP2, col=Case_Type) + 
  geom_point(size=0.8, alpha=0.6) + 
  scale_color_manual(values=ColDict_Case_Type) + 
  theme_void() + 
  NoLegend() + 
  
ggplot(WNN_Umap[WNN_Umap$Case_Type=="ALS_FTD",]) + 
  aes(UMAP1, UMAP2, col=Case_Type) + 
  geom_point(size=0.8, alpha=0.6) + 
  scale_color_manual(values=ColDict_Case_Type) + 
  theme_void() + 
  NoLegend() + 
  
ggplot(WNN_Umap[WNN_Umap$Case_Type=="C9_ALS_FTD",]) + 
  aes(UMAP1, UMAP2, col=Case_Type) + 
  geom_point(size=0.8, alpha=0.6) + 
  scale_color_manual(values=ColDict_Case_Type) + 
  theme_void() + 
  NoLegend() + 
  
ggplot(WNN_Umap[WNN_Umap$Case_Type=="HC",]) + 
  aes(UMAP1, UMAP2, col=Case_Type) + 
  geom_point(size=0.8, alpha=0.6) + 
  scale_color_manual(values=ColDict_Case_Type) + 
  theme_void() + 
  NoLegend() 


ggsave(
  filename = paste0(
    "../Data/Visualization/Figures/Fig1/", 
    "UMAP_WNN_Case-Type", 
    ".pdf"
  ), 
  dpi=300, width = 1400, height = 1400, units = "px"
)

ggsave(
  filename = paste0(
    "../Data/Visualization/Figures/Fig1/", 
    "UMAP_WNN_Case-Type", 
    ".tiff"
  ), 
  dpi=300, width = 1400, height = 1400, units = "px"
)

ggsave(
  filename = paste0(
    "../Data/Visualization/Figures/Fig1/", 
    "UMAP_WNN_Case-Type", 
    ".jpg"
  ), 
  dpi=300, width = 1400, height = 1400, units = "px"
)

rm(WNN_Umap)




  ### 6.0 WNN_L2 Clusters RNA Marker Heatmap -----------------------------------



    ## 6.1 Generate Metacell sample count matrices -----------------------------

CellTypes <- unique(M0_RNA$WNN_L2) 

WNN_Metacells_Dgx <- list() 
WNN_Metacells_Assignments <- list() 
WNN_Metacells_Samples <- list() 

for (CellType in CellTypes){
  message(paste0("Working on ", CellType, "... ")) 
  M0.tmp <- subset(M0_RNA, subset=WNN_L2==CellType) 
  umap.tmp <- data.frame(
    umap_coord1 = M0.tmp@reductions$wumap@cell.embeddings[,1], 
    umap_coord2 = M0.tmp@reductions$wumap@cell.embeddings[,2]
  ) 
  exp_mat.tmp <- M0.tmp@assays$RNA@counts
  cell_meta.tmp <- M0.tmp@meta.data 
  feature_meta.tmp <- M0.tmp@assays$RNA@meta.features 
  feature_meta.tmp$Gene <- rownames(feature_meta.tmp)
  feature_meta.tmp$gene_short_name <- feature_meta.tmp$Gene 
  feature_meta.tmp$site_name <- feature_meta.tmp$Gene
  
  cds.tmp <- new_cell_data_set(
    expression_data = exp_mat.tmp, 
    cell_metadata = cell_meta.tmp, 
    gene_metadata = feature_meta.tmp
  )
  rm(M0.tmp) 
  all(rownames(umap.tmp)==colnames(exp_mat.tmp))
  
  cicero_cds.tmp <- make_cicero_cds(
    cds.tmp, 
    reduced_coordinates = umap.tmp, 
    k=100, 
    return_agg_info = TRUE
  )
  
  rm(cds.tmp, umap.tmp, exp_mat.tmp, cell_meta.tmp, feature_meta.tmp)
  
  set.seed(311415)
  samp.tmp <- sample(
    x = unique(as.character(cicero_cds.tmp[[2]]$agg_cell)), 
    size=100, 
    replace = FALSE 
  )
  
  
  WNN_Metacells_Dgx[[CellType]] <- cicero_cds.tmp[[1]]@assays@data@listData$counts[,colnames(cicero_cds.tmp[[1]]@assays@data@listData$counts) %in% samp.tmp]
  WNN_Metacells_Assignments[[CellType]] <- cicero_cds.tmp[[2]] 
  WNN_Metacells_Samples[[CellType]] <- samp.tmp
  rm(samp.tmp, cicero_cds.tmp)
  
}
rm(CellType)



    ## 6.2 Plot Heatmap markers ---- 

for (i in 2:length(WNN_Metacells_Dgx)){
  message(paste0(names(WNN_Metacells_Dgx)[i], " features concordantly sorted: ", all(rownames(WNN_Metacells_Dgx[[1]])==rownames(WNN_Metacells_Dgx[[i]])))) 
  if(!all(rownames(WNN_Metacells_Dgx[[1]])==rownames(WNN_Metacells_Dgx[[i]]))){stop("Features not concordantly sorted!")} 
}
for (i in 1:length(WNN_Metacells_Dgx)){
  colnames(WNN_Metacells_Dgx[[i]]) <- paste0(
    names(WNN_Metacells_Dgx)[i], 
    "_", 
    colnames(WNN_Metacells_Dgx[[i]])
  )
}

dgx_WNN <- do.call(
  cbind, 
  WNN_Metacells_Dgx[c(1:5,7)]
)
dim(dgx_WNN)
colnames(dgx_WNN)

WNN_HM <- CreateSeuratObject(counts=dgx_WNN)
WNN_HM
WNN_HM <- NormalizeData(WNN_HM)
WNN_HM <- ScaleData(WNN_HM) 

WNN_HM$tmp.id <- as.character(WNN_HM$orig.ident)
class(WNN_HM$tmp.id)
WNN_HM$tmp.id <- str_replace_all(WNN_HM$tmp.id, "Astrocytes", "ASC")
WNN_HM$tmp.id <- str_replace_all(WNN_HM$tmp.id, "Exc", "EXC")
WNN_HM$tmp.id <- str_replace_all(WNN_HM$tmp.id, "Inh", "INH") 
WNN_HM$tmp.id <- str_replace_all(WNN_HM$tmp.id, "Microglia", "MG")  
WNN_HM$tmp.id <- str_replace_all(WNN_HM$tmp.id, "Oligodendrocytes", "OLG")   
WNN_HM$tmp.id <- str_replace_all(WNN_HM$tmp.id, "OPC", "OPC")  
WNN_HM$tmp.id <- factor(WNN_HM$tmp.id, levels=c(
  "INH", 
  "EXC", 
  "OPC", 
  "OLG", 
  "ASC", 
  "MG"
))
Idents(WNN_HM) <- WNN_HM$tmp.id

p1 <- rasterize(DoHeatmap(
  WNN_HM, 
  features=c(
    "P2RY12", 
    "TYROBP",
    "CSF1R", 
    "SLC1A2", 
    "AQP4", 
    "GFAP", 
    "MOG",
    "MBP", 
    "OPALIN",
    "CSPG4",
    "PDGFRA",  
    "SNAP25", 
    "SATB2", 
    "SLC17A7", 
    "ADARB2",  
    "GAD2", 
    "GAD1"
  ), 
  group.colors = c("ASC"="#6FD22E", "EXC"="#FF2C19", "INH"="#4A90E2", "MG"="#F1D91C", "OLG"="#BE3DDD", "OPC"="#F941D0")
) + 
  scale_fill_viridis_c() + 
  theme(
    axis.text.y = element_text(angle=135, vjust=1, color="black"), 
  )
)


p1
ggsave(
  filename = paste0(
    "../Data/Visualization/Figures/Fig1/", 
    "RNA_Marker_WNN_L2_Cluster_100_Metacells_Heatmap.pdf"
  ), 
  plot = p1, 
)



  ### 7.0 Plot Heatmap dendrograms ---------------------------------------------



    ## 7.1 RNA Assay -----------------------------------------------------------

Idents(M0_RNA) <- M0_RNA$WNN_L2

M0_RNA <- BuildClusterTree(
  M0_RNA, 
  assay="RNA", 
  features=Features(M0_RNA)
) 


PlotClusterTree(
  M0_RNA, direction="leftward"
) 

ggsave(
  filename = paste0(
    "../Data/Visualization/Figures/Fig1/", 
    "RNA_WNN_L2_Tree.pdf"
  ), 
  plot = PlotClusterTree(
    M0_RNA, direction="leftwards"
  ) 
) 



    ## 7.2 ATAC Assay ----------------------------------------------------------

M0_ATAC <- qread(
  paste0(
    "../Data/SeuratObjects/", 
    "M0_ATAC", 
    ".qrds"
  ), 
  nthr=nthr
)

Idents(M0_ATAC) <- M0_ATAC$WNN_L2
DefaultAssay(M0_ATAC) 

M0_ATAC <- BuildClusterTree(
  M0_ATAC, 
  assay="ATAC", 
  features=Features(M0_ATAC)
) 

PlotClusterTree(
  M0_ATAC, 
  direction = "leftwards"
) 

ggsave(
  filename = paste0(
    "../Data/Visualization/Figures/Fig1/", 
    "ATAC_WNN_L2_Tree.pdf"
  ), 
  plot = PlotClusterTree(
    M0_ATAC, 
    direction="leftwards"
  )
) 








