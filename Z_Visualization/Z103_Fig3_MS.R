

  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 

library(qs)
library(tidyverse)
library(data.table)
library(spatialLIBD)
library(Seurat)
library(EnsDb.Hsapiens.v86) 


qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)




  ### 1.0 Load data ------------------------------------------------------------

M0_WNN_L25_HC_ModuleScores <- qread(
  paste0(
    "../Data/Annotations/Module_Scores/RNA/HC_Derived/", 
    "M0_WNN_L25_HC_ModuleScores", 
    ".qrds"
  ), 
  nthr=nthr
)    

M0_WNN_L4_HC_ModuleScores <- qread(
  paste0(
    "../Data/Annotations/Module_Scores/RNA/HC_Derived/", 
    "M0_WNN_L4_HC_ModuleScores", 
    ".qrds"
  ), 
  nthr=nthr
)    

M0_ETNC_HC_ModuleScores <- qread(
  paste0(
    "../Data/Annotations/Module_Scores/RNA/HC_Derived/", 
    "M0_ETNC_HC_ModuleScores", 
    ".qrds"
  ), 
  nthr=nthr
)    


LIBD_Spatial_Seurat_WNN_L25_HC_ModuleScores <- qread(
  paste0(
    "../Data/Annotations/Module_Scores/RNA/HC_Derived/", 
    "LIBD_Spatial_Seurat_WNN_L25_HC_ModuleScores", 
    ".qrds"
  ), 
  nthr=nthr
)    

LIBD_Spatial_Seurat_WNN_L4_HC_ModuleScores <- qread(
  paste0(
    "../Data/Annotations/Module_Scores/RNA/HC_Derived/", 
    "LIBD_Spatial_Seurat_WNN_L4_HC_ModuleScores", 
    ".qrds"
  ), 
  nthr=nthr
)   

LIBD_Spatial_Seurat_ETNC_HC_ModuleScores <- qread(
  paste0(
    "../Data/Annotations/Module_Scores/RNA/HC_Derived/", 
    "LIBD_Spatial_Seurat_ETNC_HC_ModuleScores", 
    ".qrds"
  ), 
  nthr=nthr
)


spe <- qread( 
  paste0(
    "../Data/SpatialData/", 
    "LIBD_Spatial_Spe", 
    ".qrds"
  ), 
  nthr = nthr
)

all(colnames(spe)==LIBD_Spatial_Seurat_WNN_L25_HC_ModuleScores$Original_Barcode)
all(spe$key==LIBD_Spatial_Seurat_WNN_L25_HC_ModuleScores$key) 

all(colnames(spe)==LIBD_Spatial_Seurat_ETNC_HC_ModuleScores$Original_Barcode)
all(spe$key==LIBD_Spatial_Seurat_ETNC_HC_ModuleScores$key)  

Files  = setNames(
  read.csv("../Data/cfg/Files_List.txt", header=FALSE, sep="\t")[,2], 
  nm=read.csv("../Data/cfg/Files_List.txt", header=FALSE, sep="\t")[,1]
)


CARD <- qread(
  Files["CARD_WNN_L25"], 
  nthr=nthr
)


all(CARD$key==spe$key)

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

ColDict_WNN_L4 <- setNames(
  object = readxl::read_xlsx(paste0(
    "../Data/Visualization/", 
    "ALS_Brain_Multiome_ColDicts", 
    ".xlsx"
  ), 
  sheet = "WNN_L4"
  )$Color, 
  nm = readxl::read_xlsx(paste0(
    "../Data/Visualization/", 
    "ALS_Brain_Multiome_ColDicts", 
    ".xlsx"
  ), 
  sheet = "WNN_L4"
  )$WNN_L4
)


libd_layer_colors2 <- libd_layer_colors
libd_layer_colors2["WM"] <- "#A4A4A4"
libd_layer_colors3 <- libd_layer_colors2
libd_layer_colors3["WM"] <- "#515151" 




  ### 2.0 Plot combined spatial and single-cell data per WNN_L25 CellType ------

for (i in grep("500fts", colnames(LIBD_Spatial_Seurat_WNN_L25_HC_ModuleScores))){
  
  name=colnames(LIBD_Spatial_Seurat_WNN_L25_HC_ModuleScores)[i]
  spe@colData[[name]] <- base::scale(LIBD_Spatial_Seurat_WNN_L25_HC_ModuleScores[[name]])[,1]
  
  
  df = M0_WNN_L25_HC_ModuleScores
  CellType=str_split(str_split(name, "WNN_L25_", simplify=TRUE)[,2], "_ModuleScore_", simplify=TRUE)[,1]
  
  spe@colData[[paste0("CARD_", CellType)]] <- CARD[[CellType]]
  
  
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
  
  
  df$Score <- base::scale(M0_WNN_L25_HC_ModuleScores[[name]])[,1]
  fill_cols <- c(
    setNames(ColDict_WNN_L25, nm=str_replace_all(names(ColDict_WNN_L25), "_", " ")), 
    "Exc Neurons"="#AA0000", 
    "Other Exc Neurons"="#AA0000", 
    "Inh Neurons"="#0000AA", 
    "Other Inh Neurons"="#0000AA", 
    "Glia"="#00AA00", 
    "Other Glia"="#00AA00"
  )
  
  
  p1 <- ggplot(df) + 
    aes(tmp, Score, fill=tmp) + 
    geom_hline(yintercept = 0, size=0.8, col="#00000066", linetype=2) + 
    geom_boxplot(fatten=2, size=0.8, col="#000000", outlier.shape = NA, ) + 
    scale_fill_manual(values = fill_cols) + 
    ylab("Gene signature score\n") + 
    theme_classic()   + 
    theme(
      axis.title.x = element_blank(), 
      axis.title.y = element_text(color="#000000", face = "bold.italic", size=12), 
      axis.text.x = element_text(color="#000000", angle=45, hjust=1, face="bold.italic", size = 12 ),  
      axis.text.y = element_text(color="#000000", face = "bold.italic", size=10), 
      title = element_blank(), 
      legend.position = "Null", 
      axis.line = element_line(color="#000000", linewidth=0.6)
    ) 
  
  p2 <- vis_clus(
    spe = spe,
    clustervar = "layer_guess_reordered",
    sampleid = "151674",
    colors = libd_layer_colors3,
    spatial = TRUE, 
    point_size = 1.5,
    
    ... = " LIBD Layers"
  ) + 
    theme(
      title = element_blank(), 
      legend.position = "bottom"
    ) 
  
  p2$layers[[1]]$aes_params$colour="#000000"
  
  p3 <- vis_gene(
    spe = spe[,!is.na(spe$layer_guess_reordered)],
    sampleid = "151674",
    geneid = name, 
    spatial = TRUE, 
    point_size = 1.6
  ) + theme(
    title = element_blank(), 
    legend.position = "bottom"
  )
  
  p4 <- LIBD_Spatial_Seurat_WNN_L25_HC_ModuleScores %>% 
    mutate(Score=base::scale(.data[[name]])) %>% 
    group_by(layer_guess_reordered, sample_id) %>% 
    dplyr::filter(!is.na(layer_guess_reordered)) %>% 
    summarize(Score=mean(Score)) %>% 
    ggplot() + 
    aes(layer_guess_reordered, Score, fill=layer_guess_reordered) + 
    geom_hline(yintercept = 0, size=0.8, col="#00000066", linetype=2) + 
    geom_boxplot(fatten=1.6, size=0.6, color="#000000CC", outlier.shape=NA) + 
    geom_jitter(fill="#00000088", col="#000000FF",  shape=21, size=2.4, width=0.15) + 
    scale_fill_manual(values=libd_layer_colors2) + 
    ggtitle(name) + 
    ylab("Gene signature score\n") + 
    theme_classic() + 
    theme(
      axis.title.x = element_blank(), 
      axis.title.y = element_text(color="#000000", face = "bold.italic", size=12), 
      axis.text.x = element_text(color="#000000", angle=45, hjust=1, face="bold.italic", size=12),  
      axis.text.y = element_text(color="#000000", face = "bold.italic", size=10), 
      title = element_blank(), 
      legend.position = "Null", 
      axis.line = element_line(color="#000000", linewidth=0.6)
    )
  
  
  p5 <- vis_gene(
    spe = spe[,!is.na(spe$layer_guess_reordered)],
    sampleid = "151674",
    geneid = paste0("CARD_", CellType), 
    spatial = FALSE, 
    point_size = 1.6, 
    viridis = FALSE 
  ) + 
  scale_fill_gradientn(colours = c("lightblue", "lightyellow", "red"), na.value = "#00000000") + 
  theme(
    title = element_blank(), 
    legend.position = "bottom", 
    legend.text = element_text(color="#000000", angle=45, hjust=1)
  )
  
  
  
  p_all <- plot(p1 + p2 + p3 + p4 + p5) + 
    patchwork::plot_layout(ncol=5, nrow=1, , widths=c(1,2,2,2,2), heights = c(1,1,1,1,1)) 
  
  ggsave(
    plot = p_all, 
    filename = paste0(
      "../Data/Visualization/Figures/Fig3/Combined_Spatial_", 
      CellType, 
      ".pdf"
    ), 
    dpi=300, 
    width = 5400, 
    height = 1800, 
    units = "px"
  ) 
  
  rm(df)
}




  ### 3.0 Plot extraterencephalic neuron markers -------------------------------

spe@colData["ETNC_Score_500fts"] <- base::scale(LIBD_Spatial_Seurat_WNN_L4_HC_ModuleScores$WNN_L4_Exc_FEZF2_NTNG1_ModuleScore_500fts_1)[,1]


df = M0_ETNC_HC_ModuleScores

df$tmp <- df$WNN_L25
df$tmp[df$WNN_L4=="Exc_FEZF2_NTNG1"] <- "Exc_FEZF2_NTNG1"
df$tmp[df$WNN_L25=="Exc_FEZF2" & df$WNN_L4!="Exc_FEZF2_NTNG1"] <- "Other Exc FEZF2"
df$tmp[df$WNN_L15=="Exc_Neurons" & df$WNN_L25!="Exc_FEZF2"] <- "Other Exc"
df$tmp[df$WNN_L15=="Inh_Neurons"] <- "Inh"
df$tmp[df$WNN_L15=="Glia"] <- "Glia"
df$tmp <- factor(df$tmp, levels=c("Exc_FEZF2_NTNG1", "Other Exc FEZF2", "Other Exc", "Inh", "Glia"))

fill_cols <- c(
  setNames(ColDict_WNN_L25, nm=str_replace_all(names(ColDict_WNN_L25), "_", " ")), 
  
  "Exc_FEZF2_NTNG1" = as.character(ColDict_WNN_L4["Exc_FEZF2_NTNG1"]), 
  "Other Exc FEZF2" = as.character(ColDict_WNN_L25["Exc_FEZF2"]), 
  "Other Exc"="#AA0000", 
  "Inh"="#0000AA", 
  "Glia"="#00AA00"
)


p1 <- ggplot(df) + 
  aes(tmp, ETNC_ModuleScore_500fts_1, fill=tmp) + 
  geom_hline(yintercept = 0, size=0.8, col="#00000066", linetype=2) + 
  geom_boxplot(fatten=2, size=0.8, col="#000000", outlier.shape = NA, ) + 
  scale_fill_manual(values = fill_cols) + 
  ylab("Gene signature score\n") + 
  theme_classic()   + 
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_text(color="#000000", face = "bold.italic", size=12), 
    axis.text.x = element_text(color="#000000", angle=45, hjust=1, face="bold.italic", size = 12 ),  
    axis.text.y = element_text(color="#000000", face = "bold.italic", size=10), 
    title = element_blank(), 
    legend.position = "Null", 
    axis.line = element_line(color="#000000", linewidth=0.6)
  ) 

ggsave(
  plot = p1, 
  filename = paste0(
    "../Data/Visualization/Figures/Fig3/", 
    "Exc_FEZF2_M0_RNA_Score",
    ".pdf"
  ), 
  dpi=300, 
  width = 900, 
  height = 1800, 
  units = "px"
) 


p3 <- vis_gene(
  spe = spe[,!is.na(spe$layer_guess_reordered)],
  sampleid = "151674",
  geneid = "ETNC_Score_500fts", 
  spatial = TRUE, 
  point_size = 1.6
) + theme(
  title = element_blank(), 
  legend.position = "bottom"
)

ggsave(
  plot = p3, 
  filename = paste0(
    "../Data/Visualization/Figures/Fig3/Exc_FEZF2_Spatial",
    ".pdf"
  ), 
  dpi=300, 
  width = 900, 
  height = 1800, 
  units = "px"
) 



