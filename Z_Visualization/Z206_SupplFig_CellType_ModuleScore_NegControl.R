

  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 

library(qs)
library(tidyverse)
library(spatialLIBD)


qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)




  ### 1.0 Load data ------------------------------------------------------------

M0_HC_Random_ModuleScores <- qread(
  paste0(
    "../Data/Annotations/Module_Scores/RNA/HC_Derived/", 
    "M0_HC_Random_ModuleScores", 
    ".qrds"
  ), 
  nthr=nthr
)

M0_ATAC_HC_Random_ModuleScores <- qread(
  paste0(
    "../Data/Annotations/Module_Scores/ATAC/HC_Derived/", 
    "M0_HC_Random_ModuleScores", 
    ".qrds"
  ), 
  nthr=nthr
)

LIBD_Spatial_HC_Random_ModuleScores <- qread(
  paste0(
    "../Data/Annotations/Module_Scores/RNA/HC_Derived/", 
    "LIBD_Spatial_HC_Random_ModuleScores", 
    ".qrds"
  ), 
  nthr=nthr
)


M0_RNA <- qread(
  paste0(
    "../Data/SeuratObjects/", 
    "M0_RNA", 
    ".qrds"
  )
)

M0_RNA_Metadata <- M0_RNA@meta.data 
rm(M0_RNA)


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

names(ColDict_WNN_L25) <- str_replace_all(
  names(ColDict_WNN_L25), 
  pattern = "_", 
  replacement = " "
)

libd_layer_colors2 <- libd_layer_colors
libd_layer_colors2["WM"] <- "#A4A4A4" 




  ### 2.0 Visualize M0_RNA HC Random scores ------------------------------------

for (Comparison in names(M0_HC_Random_ModuleScores)){
  
  
  M0_HC_Random_ModuleScores[[Comparison]]$WNN_L25 <- as.character(M0_HC_Random_ModuleScores[[Comparison]]$WNN_L25) 
  M0_HC_Random_ModuleScores[[Comparison]]$WNN_L25 <- str_replace_all(M0_HC_Random_ModuleScores[[Comparison]]$WNN_L25, "_", " ")
  M0_HC_Random_ModuleScores[[Comparison]]$WNN_L25 <- factor(
    M0_HC_Random_ModuleScores[[Comparison]]$WNN_L25, 
    levels = c(
      "Exc LINC00507", 
      "Exc RORB", 
      "Exc THEMIS", 
      "Exc FEZF2", 
      "Inh LAMP5 PAX6", 
      "Inh TAFA1 VIP", 
      "Inh SST", 
      "Inh PVALB", 
      "OPC", 
      "Oligodendrocytes", 
      "Astrocytes", 
      "Microglia"
    )
  )
  
 p1 <- ggplot(M0_HC_Random_ModuleScores[[Comparison]]) + 
    aes(WNN_L25, ScoreScaledCombined, fill=WNN_L25) + 
    geom_hline(yintercept = 0, lwd=1, linetype=2, col="#AAAAAA")  +
    geom_hline(yintercept = c(-1, 1), lwd=1, linetype=2, col="#FAC8C3") + 
    geom_boxplot(fatten=1, col="#000000", outlier.shape=NA) + 
    stat_summary(
      geom="point",
      fun="mean",
      shape=21,
      fill="orangered2",
      size=5
    ) + 
    scale_y_continuous(limits=c(-3,3), breaks = seq(-3, 3, 1)) + 
    scale_fill_manual(values=ColDict_WNN_L25) + 
    ggtitle(Comparison) + 
    ylab("Random genes signature score \n") + 
    theme_classic() + 
    theme(
      axis.title.x = element_blank(), 
      axis.title.y = element_text(color="#000000", face = "bold.italic", size=18), 
      axis.text.x = element_text(color="#000000", angle=45, hjust=1, face="bold", size = 15 ),  
      axis.text.y = element_text(color="#000000", face = "bold.italic", size=15), 
      legend.position = "Null", 
      axis.line = element_line(color="#000000", linewidth=0.6)
    ) 
  
  ggsave(
    plot=p1, 
    filename=paste0(
      "../Data/Visualization/Figures/SupplFig_CellType_ModuleScore_NegControl/", 
      "M0_HC_Random_ModuleScores_", 
      Comparison, 
      ".pdf"
    ), 
    dpi=300, 
    width = 1800, 
    height = 1800, 
    units = "px"
  ) 
 
  rm(p1)
}

rm(M0_HC_Random_ModuleScores, Comparison)


 
  ### 3.0 Visualize M0_ATAC HC Random scores ------------------------------------

for (Comparison in names(M0_ATAC_HC_Random_ModuleScores)){
  
  
  M0_ATAC_HC_Random_ModuleScores[[Comparison]]$WNN_L25 <- as.character(M0_ATAC_HC_Random_ModuleScores[[Comparison]]$WNN_L25) 
  M0_ATAC_HC_Random_ModuleScores[[Comparison]]$WNN_L25 <- str_replace_all(M0_ATAC_HC_Random_ModuleScores[[Comparison]]$WNN_L25, "_", " ")
  M0_ATAC_HC_Random_ModuleScores[[Comparison]]$WNN_L25 <- factor(
    M0_ATAC_HC_Random_ModuleScores[[Comparison]]$WNN_L25, 
    levels = c(
      "Exc LINC00507", 
      "Exc RORB", 
      "Exc THEMIS", 
      "Exc FEZF2", 
      "Inh LAMP5 PAX6", 
      "Inh TAFA1 VIP", 
      "Inh SST", 
      "Inh PVALB", 
      "OPC", 
      "Oligodendrocytes", 
      "Astrocytes", 
      "Microglia"
    )
  )
  
  p1 <- ggplot(M0_ATAC_HC_Random_ModuleScores[[Comparison]]) + 
    aes(WNN_L25, ScoreScaledCombined, fill=WNN_L25) + 
    geom_hline(yintercept = 0, lwd=1, linetype=2, col="#AAAAAA")  +
    geom_hline(yintercept = c(-1, 1), lwd=1, linetype=2, col="#FAC8C3") + 
    geom_boxplot(fatten=1, col="#000000", outlier.shape=NA) + 
    stat_summary(
      geom="point",
      fun="mean",
      shape=21,
      fill="orangered2",
      size=5
    ) + 
    scale_y_continuous(limits=c(-3,3), breaks = seq(-3, 3, 1)) + 
    scale_fill_manual(values=ColDict_WNN_L25) + 
    ggtitle(Comparison) + 
    ylab("Random genes signature score \n") + 
    theme_classic() + 
    theme(
      axis.title.x = element_blank(), 
      axis.title.y = element_text(color="#000000", face = "bold.italic", size=18), 
      axis.text.x = element_text(color="#000000", angle=45, hjust=1, face="bold", size = 15 ),  
      axis.text.y = element_text(color="#000000", face = "bold.italic", size=15), 
      legend.position = "Null", 
      axis.line = element_line(color="#000000", linewidth=0.6)
    ) 
  
  ggsave(
    plot=p1, 
    filename=paste0(
      "../Data/Visualization/Figures/SupplFig_ATAC_CellType_ModuleScore_NegControl/", 
      "M0_ATAC_HC_Random_ModuleScores_", 
      Comparison, 
      ".pdf"
    ), 
    dpi=300, 
    width = 1800, 
    height = 1800, 
    units = "px"
  ) 
  
  rm(p1)
}

rm(M0_ATAC_HC_Random_ModuleScores, M0_RNA_Metadata, ColDict_WNN_L25, Comparison)



  ### 4.0 Visualize LIBD Spatial HC Random scores ------------------------------

for (Comparison in names(LIBD_Spatial_HC_Random_ModuleScores)){
  
  df <- LIBD_Spatial_HC_Random_ModuleScores[[Comparison]]
  df <- df[!is.na(df$layer_guess_reordered),]
  p1 <- ggplot(df) + 
    aes(layer_guess_reordered, ScoreScaledCombined, fill=layer_guess_reordered) + 
    geom_hline(yintercept = 0, lwd=1, linetype=2, col="#AAAAAA")  +
    geom_hline(yintercept = c(-1, 1), lwd=1, linetype=2, col="#FAC8C3") + 
    geom_boxplot(fatten=1, col="#000000", outlier.shape=NA) + 
    stat_summary(
      geom="point",
      fun="mean",
      shape=21,
      fill="orangered2",
      size=5
    ) + 
    scale_y_continuous(limits=c(-3,3), breaks = seq(-3, 3, 1)) + 
    scale_fill_manual(values=libd_layer_colors2) + 
    ggtitle(Comparison) + 
    ylab("Random genes signature score \n") + 
    theme_classic() + 
    theme(
      axis.title.x = element_blank(), 
      axis.title.y = element_text(color="#000000", face = "bold.italic", size=15), 
      axis.text.x = element_text(color="#000000", angle=45, hjust=1, face="bold", size = 12),  
      axis.text.y = element_text(color="#000000", face = "bold.italic", size=12), 
      legend.position = "Null", 
      axis.line = element_line(color="#000000", linewidth=0.6)
    ) 
  
  ggsave(
    plot=p1, 
    filename=paste0(
      "../Data/Visualization/Figures/SupplFig_CellType_ModuleScore_NegControl/", 
      "LIBD_Spatial_HC_Random_ModuleScores_", 
      Comparison, 
      ".pdf"
    ), 
    dpi=300, 
    width = 1800, 
    height = 1800, 
    units = "px"
  ) 
  
  rm(p1, df)
}

rm(LIBD_Spatial_HC_Random_ModuleScores, libd_layer_colors2, Comparison)


