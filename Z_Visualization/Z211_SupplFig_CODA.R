

  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 

library(qs)
library(tidyverse)
library(ggrepel)



qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)




  ### 1.0 Load data ------------------------------------------------------------

M0_RNA <- qread(
  paste0(
    "../Data/SeuratObjects/", 
    "M0_RNA", 
    ".qrds"
  ), 
  nthr = nthr
)

MASC_ALS_WNN_L25 <- qread(
  paste0(
    "../Data/Coda/MASC/", 
    "MASC_ALS_WNN_L25", 
    ".qrds"
  ), 
  nthr = nthr
)

MASC_ALSFTD_WNN_L25 <- qread(
  paste0(
    "../Data/Coda/MASC/", 
    "MASC_ALSFTD_WNN_L25", 
    ".qrds"
  ), 
  nthr = nthr
)


ColDict_Case <- setNames(
  object = xlsx::read.xlsx(
    "../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", 
    sheetName = "Case"
  )[["Color"]], 
  nm = xlsx::read.xlsx(
    "../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", 
    sheetName = "Case"
  )[["Case"]]
) 




  ### 2.0 Define AUX functions -------------------------------------------------

only_value <- function(x){
  if(is.null(x)) stop("Only_Value: Please define variable!")
  if(length(unique(x))==1){
    return(unique(x)[1])
  } else {
    stop("Only_Value: Unequivocal variable!")
  }
}




  ### 3.0 Cell-type composition boxplots ---------------------------------------

M0_RNA@meta.data %>% 
  filter(Case %in% c("ALS", "HC")) %>% 
  group_by(ID, WNN_L25) %>% 
  summarize(CellNumber = n(), Case = only_value(Case)) %>% 
  group_by(ID) %>% mutate(TotalCells = sum(CellNumber)) %>% 
  ungroup() %>% 
  mutate(Proportion = CellNumber/TotalCells) %>% 
  mutate(WNN_L25 = factor(WNN_L25, levels = c(
    "Exc_RORB", 
    "Exc_LINC00507", 
    "Exc_THEMIS", 
    "Exc_FEZF2", 
    "Inh_LAMP5_PAX6", 
    "Inh_TAFA1_VIP", 
    "Inh_SST", 
    "Inh_PVALB", 
    "OPC", 
    "Oligodendrocytes", 
    "Astrocytes", 
    "Microglia"
  ))) %>% 
  mutate(Case = as.character(Case)) %>% 
  mutate(Case = factor(Case, levels = c("HC", "ALS"))) %>% 
  
  ggplot() + 
  aes(WNN_L25, Proportion, fill=Case) + 
  geom_boxplot(position = position_dodge(0.64), width = 0.5, fatten=2, outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.05), alpha = 0.2, size = 2, color = "black") + 
  scale_y_continuous(limits=c(0, 0.6)) + 
  scale_fill_manual(values=ColDict_Case) + 
  xlab("") + 
  theme_classic() + 
  theme(
    axis.line = element_line(colour = "#000000", linewidth = 1), 
    axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1.0, colour = "#000000"), 
    axis.text.y = element_text(size = 10, face = "bold.italic", hjust = 1.0, colour = "#000000"), 
    axis.ticks = element_line(linewidth = 1, colour = "#000000"), 
    axis.ticks.length =  unit(0.2, "cm"), 
  )

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/SupplFig_Coda/", 
    "Coda_ALS_HC_WNN_L25_Boxplot", 
    ".pdf"
  ),
  width = 11, 
  height = 5.5, 
  units="in"
)  


M0_RNA@meta.data %>% 
  filter(Case %in% c("ALS_FTD", "HC")) %>% 
  group_by(ID, WNN_L25) %>% 
  summarize(CellNumber = n(), Case = only_value(Case)) %>% 
  group_by(ID) %>% mutate(TotalCells = sum(CellNumber)) %>% 
  ungroup() %>% 
  mutate(Proportion = CellNumber/TotalCells) %>% 
  mutate(WNN_L25 = factor(WNN_L25, levels = c(
    "Exc_RORB", 
    "Exc_LINC00507", 
    "Exc_THEMIS", 
    "Exc_FEZF2", 
    "Inh_LAMP5_PAX6", 
    "Inh_TAFA1_VIP", 
    "Inh_SST", 
    "Inh_PVALB", 
    "OPC", 
    "Oligodendrocytes", 
    "Astrocytes", 
    "Microglia"
  ))) %>% 
  mutate(Case = as.character(Case)) %>% 
  mutate(Case = factor(Case, levels = c("HC", "ALS_FTD"))) %>% 
  
  ggplot() + 
  aes(WNN_L25, Proportion, fill=Case) + 
  geom_boxplot(position = position_dodge(0.64), width = 0.5, fatten=2, outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.05), alpha = 0.2, size = 2, color = "black") + 
  scale_y_continuous(limits=c(0, 0.6)) + 
  scale_fill_manual(values=ColDict_Case) + 
  xlab("") + 
  theme_classic() + 
  theme(
    axis.line = element_line(colour = "#000000", linewidth = 1), 
    axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1.0, colour = "#000000"), 
    axis.text.y = element_text(size = 10, face = "bold.italic", hjust = 1.0, colour = "#000000"), 
    axis.ticks = element_line(linewidth = 1, colour = "#000000"), 
    axis.ticks.length =  unit(0.2, "cm"), 
  )

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/SupplFig_Coda/", 
    "Coda_ALS_FTD_HC_WNN_L25_Boxplot", 
    ".pdf"
  ),
  width = 11, 
  height = 5.5, 
  units="in"
)  




  ### 3.0 MASC Plots -----------------------------------------------------------

MASC_ALS_WNN_L25[[1]] %>% 
  as.data.frame() %>% 
  ggplot() + 
  aes(log2(CaseALS.OR), -log10(fdr), fill=fdr<0.05, label = str_replace_all(cluster, "cluster", "" )) + 
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = 2, col = "#BBBBBB") + 
  geom_hline(yintercept = -log10(0.05), linetype = 2, col = "#BBBBBB") + 
  geom_point(size=3, pch=21) + 
  geom_text_repel(max.overlaps = 12) + 
  scale_x_continuous(limits=c(-2,2)) + 
  scale_y_continuous(limits=c(0,2)) + 
  scale_fill_manual(values=c("FALSE"="#CCCCCC", "TRUE"="#EE1111")) + 
  theme_classic()

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/SupplFig_Coda/", 
    "Coda_ALS_HC_WNN_L25_MASC_Plot", 
    ".pdf"
  ),
  width = 5.2, 
  height = 4.2, 
  units="in"
)  

  
MASC_ALSFTD_WNN_L25[[1]] %>% 
  as.data.frame() %>% 
  ggplot() + 
  aes(log2(CaseALS_FTD.OR), -log10(fdr), fill=fdr<0.05, label = str_replace_all(cluster, "cluster", "" )) + 
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = 2, col = "#BBBBBB") + 
  geom_hline(yintercept = -log10(0.05), linetype = 2, col = "#BBBBBB") + 
  geom_point(size=3, pch=21) + 
  geom_text_repel(max.overlaps = 12) + 
  scale_x_continuous(limits=c(-2,2)) + 
  scale_y_continuous(limits=c(0,2)) + 
  scale_fill_manual(values=c("FALSE"="#CCCCCC", "TRUE"="#EE1111")) + 
  theme_classic()

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/SupplFig_Coda/", 
    "Coda_ALS_FTD_HC_WNN_L25_MASC_Plot", 
    ".pdf"
  ),
  width = 5.2, 
  height = 4.2, 
  units="in"
)  

  