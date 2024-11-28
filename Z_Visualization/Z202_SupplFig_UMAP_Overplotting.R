

  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 

library(tidyverse)
library(data.table)
library(qs) 
library(Seurat)
library(ggnewscale)


qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)




  ### 1.0 Load data ------------------------------------------------------------

M0_RNA <- qread(
  paste0(
    "../Data/SeuratObjects/",  
    "M0_RNA",
    ".qrds"
  ),
  nthr=nthr 
)




  ### 2.0 UMAP WNN_L2 Density -------------------------------------------------- 

table(M0_RNA$WNN_L2)

UMAP.tmp <- data.frame(Embeddings(M0_RNA, reduction = "wumap"))
all(rownames(UMAP.tmp)==rownames(M0_RNA@meta.data))
UMAP.tmp$CellType <- "NA"
UMAP.tmp$CellType[M0_RNA@meta.data$WNN_L2 %in% c("Exc_Neurons")] <- "Exc_Neurons"
UMAP.tmp$CellType[M0_RNA@meta.data$WNN_L2 %in% c("Inh_Neurons")] <- "Inh_Neurons"
UMAP.tmp$CellType[M0_RNA@meta.data$WNN_L2 %in% c("Astrocytes", "Microglia", "Oligodendrocytes", "OPC")] <- "Glial_Cells"
table(UMAP.tmp$CellType) 

p1 <- ggplot(UMAP.tmp) + 
  aes(wUMAP_1, wUMAP_2) + 
  stat_density_2d(aes(fill = ..level..), geom = "polygon", n=100, h=1) + 
  scale_fill_gradientn(colors=c("#163B56", "#FFFFFF")) + 
  theme_void() 

p1

ggsave(
  filename = paste0(
    "../Data/Visualization/Figures/SupplFig_UMAP_Overplotting/", 
    "M0_UMAP_AllCells_Density", 
    ".jpg"
  ),
  p1, 
  dpi=300, width = 1800, height = 1400, units = "px"
)

ggsave(
  filename = paste0(
    "../Data/Visualization/Figures/SupplFig_UMAP_Overplotting/", 
    "M0_UMAP_AllCells_Density", 
    ".pdf"
  ),
  p1, 
  dpi=300, width = 1800, height = 1400, units = "px"
)

x.tmp <- layer_scales(p1)$x$range$range
y.tmp <- layer_scales(p1)$y$range$range




  ### 3.0 UMAP WNN_L2 HexBins -------------------------------------------------- 

p1 <- ggplot() + 
  geom_hex(aes(wUMAP_1, wUMAP_2), data=UMAP.tmp[UMAP.tmp$CellType=="Exc_Neurons",], bins=120, col="black", size=0.1, show.legend = TRUE) + 
  scale_fill_gradientn(name="Legend", colors=c("#AA0000", "#FFEEDD", "#FFFFFF"), limits=c(0,500)) + 
  
  new_scale_fill() + 
  geom_hex(aes(wUMAP_1, wUMAP_2), data=UMAP.tmp[UMAP.tmp$CellType=="Inh_Neurons",], bins=120, col="black", size=0.1, show.legend = TRUE) + 
  scale_fill_gradientn(name="Legend", colors=c("#0000AA", "#DDEEFF", "#FFFFFF"), limits=c(0,500) ) +  
  
  new_scale_fill() + 
  geom_hex(aes(wUMAP_1, wUMAP_2), data=UMAP.tmp[UMAP.tmp$CellType=="Glial_Cells",], bins=120, col="black", size=0.1) + 
  scale_fill_gradientn(name="Legend", colors=c("#00AA00", "#DDFFDD", "#FFFFFF"), limits=c(0,500)) + 
  scale_x_continuous(limits=x.tmp) + 
  scale_y_continuous(limits=x.tmp) + 
  theme_void() 
  
p1 

ggsave(
  filename = paste0(
    "../Data/Visualization/Figures/SupplFig_UMAP_Overplotting/", 
    "M0_UMAP_AllCells_HexBin", 
    ".jpg"
  ),
  p1, 
  dpi=300, width = 1800, height = 1400, units = "px"
)

ggsave(
  filename = paste0(
    "../Data/Visualization/Figures/SupplFig_UMAP_Overplotting/", 
    "M0_UMAP_AllCells_HexBin", 
    ".pdf"
  ),
  p1, 
  dpi=300, width = 1800, height = 1800, units = "px"
)
