

  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 

library(qs)  
library(tidyverse)
library(ggraph)

library(data.table)
library(Seurat)
library(ggnewscale)
library(patchwork)
library(ggrastr)
library(pbapply)
library(DESeq2)
library(data.tree)
library(networkD3) 
library(igraph)



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


ColDict_WNN_L15 <- setNames(
  object = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "WNN_L15")$Color, 
  nm = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "WNN_L15")$WNN_L15
)

ColDict_WNN_L25 <- setNames(
  object = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "WNN_L25")$Color, 
  nm = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "WNN_L25")$WNN_L25
)

ColDict_WNN_L3 <- setNames(
  object = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "WNN_L3")$Color, 
  nm = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "WNN_L3")$WNN_L3
)

ColDict_WNN_L4 <- setNames(
  object = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "WNN_L4")$Color, 
  nm = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "WNN_L4")$WNN_L4
)




  ### 2.0 Generate hierarchical cell-type tree from curated hierarchical cell-type annotation ----   

Tree <- M0_RNA@meta.data[,colnames(M0_RNA@meta.data) %in% c(
  "AllCells", 
  "WNN_L15", 
  "WNN_L25", 
  "WNN_L3", 
  "WNN_L4", 
  "CellId"
),]

Tree <- Tree %>% 
          group_by(AllCells, WNN_L15, WNN_L25, WNN_L3, WNN_L4) %>% 
          summarise(
          Cells=length(CellId)
          ) 

Tree <-rowwise(Tree) %>% 
  mutate(pathString = paste0("AllCells/", 
                             WNN_L15, "/", 
                             #if (WNN_L25!=WNN_L15) WNN_L25, if (WNN_L25!=WNN_L15) "/",
                             WNN_L25, "/",
                             if (WNN_L3!=WNN_L25) {WNN_L3} else{ paste0(WNN_L25, "_2")}, "/",
                             if (WNN_L4!=WNN_L3) {WNN_L4} else{ paste0(WNN_L3, "_3")}, "/"
  )
) 




  ### 3.0 Plot the hierarchical cell-type tree ---------------------------------

ColDict_WNN_L3_2 <- ColDict_WNN_L3
names(ColDict_WNN_L3_2) <- paste0(names(ColDict_WNN_L3_2), "_", 2)
ColDict_WNN_L4_2 <- ColDict_WNN_L4
names(ColDict_WNN_L4_2) <- paste0(names(ColDict_WNN_L4_2), "_", 3)
ColDict_All <- c(
    ColDict_WNN_L15, 
    ColDict_WNN_L25, 
    ColDict_WNN_L3, 
    ColDict_WNN_L3_2, 
    ColDict_WNN_L4, 
    ColDict_WNN_L4_2
)


p1 <- ggraph(as.Node(Tree), "dendrogram") + 
  geom_edge_elbow(flipped=FALSE) + 
  geom_node_label(aes(label = name, fill=name), hjust = 0, repel = FALSE) + 
  scale_fill_manual(values=ColDict_All) + 
  coord_flip() + 
  scale_y_reverse(expand=c(0.5,0,1,0)) + 
  theme_void() + NoLegend()
p1 



ggsave(
  filename = paste0(
    "../Data/Visualization/Figures/SupplFig_Hierarchical_CellTypes_Trees/", 
    "Hierarchical_CellType_Annotation_Dendrogram", 
    ".pdf"
  ),
  p1,
  dpi=300, width = 3600, height = 4800, units = "px"
)




  ### 4.0 Plot stats -----------------------------------------------------------

round(prop.table(table(M0_RNA$WNN_L15))*100,0)
round(prop.table(table(M0_RNA$WNN_L25))*100,0)

round(prop.table(table(M0_RNA$WNN_L25[M0_RNA$WNN_L15=="Inh_Neurons"]))*100,0)
round(prop.table(table(M0_RNA$WNN_L3[M0_RNA$WNN_L25=="Inh_TAFA1_VIP"]))*100,0)
round(prop.table(table(M0_RNA$WNN_L3[M0_RNA$WNN_L25=="Inh_LAMP5_PAX6"]))*100,0)


round(prop.table(table(M0_RNA$WNN_L25[M0_RNA$WNN_L15=="Glia"]))*100,0)
round(prop.table(table(M0_RNA$WNN_L25[M0_RNA$WNN_L15=="Exc_Neurons"]))*100,0)

round(prop.table(table(as.character(M0_RNA$WNN_L4)))*100, 1)[grep("Inh_VIP", names(round(prop.table(table(as.character(M0_RNA$WNN_L4)))*100, 2)))]
round(prop.table(table(as.character(M0_RNA$WNN_L4[M0_RNA$WNN_L3=="Inh_VIP"])))*100, 0)

round(prop.table(table(as.character(M0_RNA$WNN_L4)))*100, 1)[grep("Inh_TAFA1", names(round(prop.table(table(as.character(M0_RNA$WNN_L4)))*100, 2)))]
round(prop.table(table(as.character(M0_RNA$WNN_L4[M0_RNA$WNN_L3=="Inh_TAFA1"])))*100, 0)

round(prop.table(table(as.character(M0_RNA$WNN_L4)))*100, 1)[grep("Inh_SST", names(round(prop.table(table(as.character(M0_RNA$WNN_L4)))*100, 2)))]
round(prop.table(table(as.character(M0_RNA$WNN_L4[M0_RNA$WNN_L3=="Inh_SST"])))*100, 0) 

round(prop.table(table(as.character(M0_RNA$WNN_L4)))*100, 1)[grep("Inh_PVALB", names(round(prop.table(table(as.character(M0_RNA$WNN_L4)))*100, 2)))]
round(prop.table(table(as.character(M0_RNA$WNN_L4[M0_RNA$WNN_L3=="Inh_PVALB"])))*100, 0) 

round(prop.table(table(as.character(M0_RNA$WNN_L4)))*100, 1)[grep("Inh_PAX6", names(round(prop.table(table(as.character(M0_RNA$WNN_L4)))*100, 2)))]
round(prop.table(table(as.character(M0_RNA$WNN_L4[M0_RNA$WNN_L3=="Inh_PAX6"])))*100, 0) 

round(prop.table(table(as.character(M0_RNA$WNN_L4)))*100, 1)[grep("Inh_LAMP5", names(round(prop.table(table(as.character(M0_RNA$WNN_L4)))*100, 2)))]
round(prop.table(table(as.character(M0_RNA$WNN_L4[M0_RNA$WNN_L3=="Inh_LAMP5"])))*100, 0) 

round(prop.table(table(as.character(M0_RNA$WNN_L4)))*100, 1)[grep("Oligodendrocytes", names(round(prop.table(table(as.character(M0_RNA$WNN_L4)))*100, 2)))]
round(prop.table(table(as.character(M0_RNA$WNN_L4[M0_RNA$WNN_L3=="Oligodendrocytes"])))*100, 0) 

round(prop.table(table(as.character(M0_RNA$WNN_L4)))*100, 1)[grep("OPC", names(round(prop.table(table(as.character(M0_RNA$WNN_L4)))*100, 2)))]
round(prop.table(table(as.character(M0_RNA$WNN_L4[M0_RNA$WNN_L3=="OPC"])))*100, 0) 

round(prop.table(table(as.character(M0_RNA$WNN_L4)))*100, 1)[grep("Microglia", names(round(prop.table(table(as.character(M0_RNA$WNN_L4)))*100, 2)))]
round(prop.table(table(as.character(M0_RNA$WNN_L4[M0_RNA$WNN_L3=="Microglia"])))*100, 0) 

round(prop.table(table(as.character(M0_RNA$WNN_L4)))*100, 1)[grep("Astrocytes", names(round(prop.table(table(as.character(M0_RNA$WNN_L4)))*100, 2)))]
round(prop.table(table(as.character(M0_RNA$WNN_L4[M0_RNA$WNN_L3=="Astrocytes"])))*100, 0) 

round(prop.table(table(as.character(M0_RNA$WNN_L4)))*100, 1)[grep("Exc_THEMIS", names(round(prop.table(table(as.character(M0_RNA$WNN_L4)))*100, 2)))]
round(prop.table(table(as.character(M0_RNA$WNN_L4[M0_RNA$WNN_L3=="Exc_THEMIS"])))*100, 0) 

round(prop.table(table(as.character(M0_RNA$WNN_L4)))*100, 1)[grep("Exc_RORB", names(round(prop.table(table(as.character(M0_RNA$WNN_L4)))*100, 2)))]
round(prop.table(table(as.character(M0_RNA$WNN_L4[M0_RNA$WNN_L3=="Exc_RORB"])))*100, 0) 

round(prop.table(table(as.character(M0_RNA$WNN_L4)))*100, 1)[grep("Exc_LINC00507", names(round(prop.table(table(as.character(M0_RNA$WNN_L4)))*100, 2)))]
round(prop.table(table(as.character(M0_RNA$WNN_L4[M0_RNA$WNN_L3=="Exc_LINC00507"])))*100, 0) 

round(prop.table(table(as.character(M0_RNA$WNN_L4)))*100, 1)[grep("Exc_FEZF2", names(round(prop.table(table(as.character(M0_RNA$WNN_L4)))*100, 2)))]
round(prop.table(table(as.character(M0_RNA$WNN_L4[M0_RNA$WNN_L3=="Exc_FEZF2"])))*100, 0) 

round(prop.table(table(M0_RNA$WNN_L3))*100,0)
round(prop.table(table(M0_RNA$WNN_L4))*100,0)
round(prop.table(table(M0_RNA$WNN_L15))*100,0)





