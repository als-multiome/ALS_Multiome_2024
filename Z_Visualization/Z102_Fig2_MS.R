


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

ColDict_WNN_L2 <- setNames(
  object = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "WNN_L2")$Color, 
  nm = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "WNN_L2")$WNN_L2
)

setequal(names(ColDict_WNN_L2), unique(M0_RNA$WNN_L2))
setdiff(names(ColDict_WNN_L2), unique(M0_RNA$WNN_L2))
any(duplicated(ColDict_WNN_L2))
any(duplicated(names(ColDict_WNN_L2)))


ColDict_WNN_L4 <- setNames(
  object = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "WNN_L4")$Color, 
  nm = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "WNN_L4")$WNN_L4
)

setequal(names(ColDict_WNN_L4), unique(M0_RNA$WNN_L4))
setdiff(names(ColDict_WNN_L4), unique(M0_RNA$WNN_L4))
any(duplicated(ColDict_WNN_L4))
any(duplicated(names(ColDict_WNN_L4)))




  ### 2.0 UMAPs WNN L2 Cell Types ----------------------------------------------



    ## 2.1 Exc_Neurons ---------------------------------------------------------


      # 2.1.1 DimPlot ---- -----------------------------------------------------

p1 <- DimPlot(M0, cells.highlight = Cells(M0)[M0$WNN_L2 == "Exc_Neurons"]) + NoLegend() 
p1 
layer_scales(p1)$x$range$range
layer_scales(p1)$y$range$range


p1 <- DimPlot(
  M0, 
  cells = Cells(M0)[M0$WNN_L2=="Exc_Neurons"], 
  pt.size=0.5
  ) + 
  scale_color_manual(values=ColDict_WNN_L4) + 
  scale_x_continuous(limits=c(-26,20)) + 
  scale_y_continuous(limits=c(-27,17)) + 
  theme_void() + 
  NoLegend() 
p1

ggsave(
  filename = paste0(
    "../Data/Visualization/Figures/Fig2/", 
    "UMAP_Exc_WNN_L4", 
    ".tiff"
  ), 
  plot = p1, 
  dpi=300, 
  width = 2000, 
  height = 2000, 
  units = "px"
) 


      # 2.1.2 Highlight WNN_L3 Cell Types --------------------------------------

table(M0$WNN_L3)

DimPlot(
  M0, 
  cells = Cells(M0)[M0$WNN_L2=="Exc_Neurons"], 
  cells.highlight = Cells(M0)[M0$WNN_L3=="Exc_FEZF2"],
  pt.size=1
) +
  theme_void() + 
  NoLegend() 


DimPlot(
  M0, 
  cells = Cells(M0)[M0$WNN_L2=="Exc_Neurons"], 
  cells.highlight = Cells(M0)[M0$WNN_L3=="Exc_LINC00507"],
  pt.size=1
) +
  theme_void() + 
  NoLegend() 


DimPlot(
  M0, 
  cells = Cells(M0)[M0$WNN_L2=="Exc_Neurons"], 
  cells.highlight = Cells(M0)[M0$WNN_L3=="Exc_RORB"],
  pt.size=1
) +
  theme_void() + 
  NoLegend() 


DimPlot(
  M0, 
  cells = Cells(M0)[M0$WNN_L2=="Exc_Neurons"], 
  cells.highlight = Cells(M0)[M0$WNN_L3=="Exc_THEMIS"],
  pt.size=1
) +
  theme_void() + 
  NoLegend() 


      # 2.1.3 Exc_Neurons WNN_L4 Clusters Table --------------------------------

M0@meta.data[M0$WNN_L2=="Exc_Neurons",] %>% 
  group_by(WNN_L3) %>% 
    summarise(Clusters=unique(WNN_L4)) %>% 
      View()

M0@meta.data[M0$WNN_L2=="Exc_Neurons",] %>% 
  group_by(WNN_L3) %>% 
    summarise(Clusters=unique(WNN_L4)) %$% 
      WNN_L3 %>% 
        table() 



    ## 2.2 Inh_Neurons ---------------------------------------------------------


      # 2.2.1 DimPlot ----------------------------------------------------------

p1 <- DimPlot(M0, cells.highlight = Cells(M0)[M0$WNN_L2 == "Inh_Neurons"]) + NoLegend() 

p1 <- DimPlot(
  M0, 
  cells = Cells(M0)[M0$WNN_L2=="Inh_Neurons"], 
  pt.size=0.5
) + 
  scale_color_manual(values=ColDict_WNN_L4) + 
  scale_x_continuous(limits=c(-26,20)) + 
  scale_y_continuous(limits=c(-27,17)) + 
  theme_void() + 
  NoLegend() 
p1

ggsave(
  filename = paste0(
    "../Data/Visualization/Figures/Fig2/", 
    "UMAP_Inh_WNN_L4", 
    ".tiff"
  ), 
  plot = p1, 
  dpi=300, 
  width = 2000, 
  height = 2000, 
  units = "px"
) 


      # 2.2.2 Highlight WNN_L3 Cell Types --------------------------------------

table(M0$WNN_L3)

DimPlot(
  M0, 
  cells = Cells(M0)[M0$WNN_L2=="Inh_Neurons"], 
  cells.highlight = Cells(M0)[M0$WNN_L3=="Inh_LAMP5"],
  pt.size=1
) +
  theme_void() + 
  NoLegend() 

DimPlot(
  M0, 
  cells = Cells(M0)[M0$WNN_L2=="Inh_Neurons"], 
  cells.highlight = Cells(M0)[M0$WNN_L3=="Inh_PAX6"],
  pt.size=1
) +
  theme_void() + 
  NoLegend() 

DimPlot(
  M0, 
  cells = Cells(M0)[M0$WNN_L2=="Inh_Neurons"], 
  cells.highlight = Cells(M0)[M0$WNN_L3=="Inh_PVALB"],
  pt.size=1
) +
  theme_void() + 
  NoLegend() 
 

DimPlot(
  M0, 
  cells = Cells(M0)[M0$WNN_L2=="Inh_Neurons"], 
  cells.highlight = Cells(M0)[M0$WNN_L3=="Inh_SST"],
  pt.size=1
) +
  theme_void() + 
  NoLegend() 

DimPlot(
  M0, 
  cells = Cells(M0)[M0$WNN_L2=="Inh_Neurons"], 
  cells.highlight = Cells(M0)[M0$WNN_L3=="Inh_TAFA1"],
  pt.size=1
) +
  theme_void() + 
  NoLegend()  

DimPlot(
  M0, 
  cells = Cells(M0)[M0$WNN_L2=="Inh_Neurons"], 
  cells.highlight = Cells(M0)[M0$WNN_L3=="Inh_VIP"],
  pt.size=1
) +
  theme_void() + 
  NoLegend() 


      # 2.2.3 Inh_Neurons WNN_L4 Clusters Table -------------------------------- 

M0@meta.data[M0$WNN_L2=="Inh_Neurons",] %>% 
  group_by(WNN_L3) %>% 
    summarise(Clusters=unique(WNN_L4)) %>% 
      View()

M0@meta.data[M0$WNN_L2=="Inh_Neurons",] %>% 
  group_by(WNN_L3) %>% 
    summarise(Clusters=unique(WNN_L4)) %$% 
      WNN_L3 %>% 
        table() 




    ## 2.3 Glia cells ----------------------------------------------------------


      # 2.3.1 DimPlot ---- 

table(M0$WNN_L2)

p1 <- DimPlot(M0, cells.highlight = Cells(M0)[M0$WNN_L2 %in% c("Astrocytes", "Microglia", "Oligodendrocytes", "OPC")]) + NoLegend() 
p1 

p1 <- DimPlot(
  M0, 
  cells = Cells(M0)[M0$WNN_L2 %in% c("Astrocytes", "Microglia", "Oligodendrocytes", "OPC")], 
  pt.size=0.5
) + 
  scale_color_manual(values=ColDict_WNN_L4) + 
  scale_x_continuous(limits=c(-26,20)) + 
  scale_y_continuous(limits=c(-27,17)) + 
  theme_void() + 
  NoLegend() 
p1

ggsave(
  filename = paste0(
    "../Data/Visualization/Figures/Fig2/", 
    "UMAP_Glia_WNN_L4", 
    ".tiff"
  ), 
  plot = p1, 
  dpi=300, 
  width = 2000, 
  height = 2000, 
  units = "px"
) 


      # 2.3.2 Highlight WNN_L3 Cell Types --------------------------------------

table(M0$WNN_L3)

DimPlot(
  M0, 
  cells = Cells(M0)[M0$WNN_L2 %in% c("Astrocytes", "Microglia", "Oligodendrocytes", "OPC")], 
  cells.highlight = Cells(M0)[M0$WNN_L3=="Astrocytes"],
  pt.size=1
) +
  theme_void() + 
  NoLegend() 

DimPlot(
  M0, 
  cells = Cells(M0)[M0$WNN_L2 %in% c("Astrocytes", "Microglia", "Oligodendrocytes", "OPC")], 
  cells.highlight = Cells(M0)[M0$WNN_L3=="Microglia"],
  pt.size=1
) +
  theme_void() + 
  NoLegend() 

DimPlot(
  M0, 
  cells = Cells(M0)[M0$WNN_L2 %in% c("Astrocytes", "Microglia", "Oligodendrocytes", "OPC")], 
  cells.highlight = Cells(M0)[M0$WNN_L3=="Oligodendrocytes"],
  pt.size=1
) +
  theme_void() + 
  NoLegend() 


DimPlot(
  M0, 
  cells = Cells(M0)[M0$WNN_L2 %in% c("Astrocytes", "Microglia", "Oligodendrocytes", "OPC")], 
  cells.highlight = Cells(M0)[M0$WNN_L3=="OPC"],
  pt.size=1
) +
  theme_void() + 
  NoLegend() 


      # 5.3.3 Glial Cells WNN_L4 Clusters Table -------------------------------- 

M0@meta.data[M0$WNN_L2 %in% c("Astrocytes", "Microglia", "Oligodendrocytes", "OPC"),] %>% 
  group_by(WNN_L3) %>% 
  summarise(Clusters=unique(WNN_L4)) %>% 
  View()

M0@meta.data[M0$WNN_L2 %in% c("Astrocytes", "Microglia", "Oligodendrocytes", "OPC"),] %>% 
  group_by(WNN_L3) %>% 
  summarise(Clusters=unique(WNN_L4)) %$% 
  WNN_L3 %>% 
  table() 




  ### 5.0 UMAP WNN_L2 ----------------------------------------------------------

table(M0$WNN_L2)

UMAP.tmp <- data.frame(Embeddings(M0, reduction = "wumap"))
all(rownames(UMAP.tmp)==rownames(M0@meta.data))
UMAP.tmp$CellType <- "NA"
UMAP.tmp$CellType[M0@meta.data$WNN_L2 %in% c("Exc_Neurons")] <- "Exc_Neurons"
UMAP.tmp$CellType[M0@meta.data$WNN_L2 %in% c("Inh_Neurons")] <- "Inh_Neurons"
UMAP.tmp$CellType[M0@meta.data$WNN_L2 %in% c("Astrocytes", "Microglia", "Oligodendrocytes", "OPC")] <- "Glial_Cells"
table(UMAP.tmp$CellType) 

p1 <- ggplot(UMAP.tmp) + 
  aes(wUMAP_1, wUMAP_2) + 
  stat_density_2d(aes(fill = ..level..), geom = "polygon", n=100, h=1.5) + 
  theme_void()
p1

x.tmp <- layer_scales(p1)$x$range$range
y.tmp <- layer_scales(p1)$y$range$range


p1 <- ggplot() + 
  geom_hex(aes(wUMAP_1, wUMAP_2), data=UMAP.tmp[UMAP.tmp$CellType=="Exc_Neurons",], bins=120, col="black", size=0.1) + 
  scale_fill_gradientn(colors=c("#AA0000", "#FFEEDD", "#FFFFFF"), limits=c(0,500)) + 
  
  new_scale_fill() + 
  geom_hex(aes(wUMAP_1, wUMAP_2), data=UMAP.tmp[UMAP.tmp$CellType=="Inh_Neurons",], bins=120, col="black", size=0.1) + 
  scale_fill_gradientn(colors=c("#0000AA", "#DDEEFF", "#FFFFFF"), limits=c(0,500) ) +  
  
  new_scale_fill() + 
  geom_hex(aes(wUMAP_1, wUMAP_2), data=UMAP.tmp[UMAP.tmp$CellType=="Glial_Cells",], bins=120, col="black", size=0.1) + 
  scale_fill_gradientn(colors=c("#00AA00", "#DDFFDD", "#FFFFFF"), limits=c(0,500)) + 
  scale_x_continuous(limits=x.tmp) + 
  scale_y_continuous(limits=x.tmp) + 
  theme_void()
p1 + NoLegend()
