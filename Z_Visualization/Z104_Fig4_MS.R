

### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 


library(qs) 
library(tidyverse) 
library(DESeq2)
library(magrittr)
library(pheatmap)
library(reshape2)
library(ComplexHeatmap)




### 1.0 Load data ------------------------------------------------------------ 

M0_RNA <- qread(
  paste0(
    "../Data/SeuratObjects/", 
    "M0_RNA", 
    ".qrds"
  ), 
  nthr=nthr
)

M0_Meta <- M0_RNA@meta.data 
rm(M0_RNA)

DE_Results_SVA <- qread(
  paste0(
    "../Data/DE/WNN/AllCase/SVA/12_SVs/", 
    "DESeq_Results_12SVs_Index", 
    ".qrds"
  )
)

ColDict_Case <- setNames(
  object = readxl::read_xlsx(paste0(
    "../Data/Visualization/", 
    "ALS_Brain_Multiome_ColDicts", 
    ".xlsx"
  ), 
  sheet = "Case"
  )$Color, 
  nm = readxl::read_xlsx(paste0(
    "../Data/Visualization/", 
    "ALS_Brain_Multiome_ColDicts", 
    ".xlsx"
  ), 
  sheet = "Case"
  )$Case
)


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

DDS_List <- qread(
  paste0(
    "../Data/DE/WNN/AllCase/SVA/12_SVs/", 
    "DDS_SVA_12sv_list",
    ".qrds"
  ), 
  nthr=nthr
)

DDS_List_Index <- qread(
  paste0(
    "../Data/DE/WNN/AllCase/SVA/12_SVs/", 
    "DESeq_Results_12SVs_Index",
    ".qrds"
  ), 
  nthr=nthr
)




### 2.0 Plot TotalReads and nFeatures ----------------------------------------



## 2.1 Summarize Data ------------------------------------------------------

Results_WNN_L25 <- DE_Results_SVA[
  DE_Results_SVA$CellTypeLevel=="WNN_L25" & 
    DE_Results_SVA$Comparison=="All_Cases"
  ,]


Results_WNN_L25$CellType = factor(
  Results_WNN_L25$CellType, 
  levels=c(
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
  )
) 

M0_Meta_ALSvsHC <- M0_Meta[M0_Meta$Case %in% c("HC", "ALS"),]
M0_Meta_ALSFTDvsHC <- M0_Meta[M0_Meta$Case %in% c("HC", "ALS_FTD"),]

Results_WNN_L25$TotalCells_ALSvsHC <- 
  table(M0_Meta_ALSvsHC$WNN_L25)[match(Results_WNN_L25$CellType, names(table(M0_Meta_ALSvsHC$WNN_L25)))]

Results_WNN_L25$TotalCells_ALSFTDvsHC <- 
  table(M0_Meta_ALSFTDvsHC$WNN_L25)[match(Results_WNN_L25$CellType, names(table(M0_Meta_ALSFTDvsHC$WNN_L25)))]

ind.tmp <- which(
  DDS_List_Index$CellTypeLevel=="WNN_L25" & 
    DDS_List_Index$Comparison=="All_Cases"
)

TotalReads_ALS <- c() 
TotalReads_ALSFTD <- c() 

for (i in 1:length(ind.tmp)){
  colSums <- colSums(DDS_List[[ind.tmp[i]]]@assays@data$counts)
  TotalReads_ALS[i] <- sum(colSums[which(names(colSums) %in% unique(M0_Meta_ALSvsHC$ID))]) 
  names(TotalReads_ALS)[i] <- DDS_List_Index$CellType[ind.tmp[i]]
  TotalReads_ALSFTD[i] <- sum(colSums[which(names(colSums) %in% unique(M0_Meta_ALSFTDvsHC$ID))]) 
  names(TotalReads_ALSFTD)[i] <- DDS_List_Index$CellType[ind.tmp[i]]
}
rm(ind.tmp,i)

setequal(names(TotalReads_ALS), Results_WNN_L25$CellType)
Results_WNN_L25$TotalReads_ALSvsHC <- TotalReads_ALS[match(Results_WNN_L25$CellType, names(TotalReads_ALS))]
Results_WNN_L25$TotalReads_ALSFTDvsHC <- TotalReads_ALSFTD[match(Results_WNN_L25$CellType, names(TotalReads_ALSFTD))]
rm(TotalReads_ALS, TotalReads_ALSFTD)

Results_WNN_L25 <- Results_WNN_L25[match(levels(Results_WNN_L25$CellType), Results_WNN_L25$CellType),]



## 2.2 Plot ALSvsHC --------------------------------------------------------

ggplot(
  Results_WNN_L25
) + 
  aes(x=CellType) + 
  geom_col(aes(y=TotalCells_ALSvsHC, fill=CellType), col="#000000") + 
  geom_line(aes(y=TotalReads_ALSvsHC/1e4, group="All"),col="#BBBBBB") + 
  geom_point(aes(y=TotalReads_ALSvsHC/1e4), pch=21, size=3, fill="#BBBBBB") + 
  scale_fill_manual(values=ColDict_WNN_L25) + 
  scale_y_continuous(
    limits=c(0,50000),
    sec.axis = sec_axis(~.*1e-2, name="Second Axis"), 
    expand=c(0,0,0.1,0)
  ) + 
  theme_classic() + 
  theme(
    legend.position = "NULL", 
    axis.text.x = element_text(angle=45, hjust=1)
  )

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig4/", 
    "DE_RNA_WNN_L25_nReads_nCells_ALSvsHC",
    ".pdf"
  ), 
  width = 6.92, 
  height = 4.21, 
  units="in"
)



## 2.3 Plot ALSFTDvsHC -----------------------------------------------------

ggplot(
  Results_WNN_L25
) + 
  aes(x=CellType) + 
  geom_col(aes(y=TotalCells_ALSFTDvsHC, fill=CellType), col="#000000") + 
  geom_line(aes(y=TotalReads_ALSFTDvsHC/1e4, group="All"),col="#BBBBBB") + 
  geom_point(aes(y=TotalReads_ALSFTDvsHC/1e4), pch=21, size=3, fill="#BBBBBB") + 
  scale_fill_manual(values=ColDict_WNN_L25) + 
  scale_y_continuous(
    limits=c(0,50000),
    sec.axis = sec_axis(~.*1e-2, name="Second Axis"), 
    expand=c(0,0,0.1,0)
  ) + 
  theme_classic() + 
  theme(
    legend.position = "NULL", 
    axis.text.x = element_text(angle=45, hjust=1)
  )

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig4/", 
    "DE_RNA_WNN_L25_nReads_nCells_ALSFTDvsHC",
    ".pdf"
  ), 
  width = 6.92, 
  height = 4.21, 
  units="in"
)




### 3.0 Plot DEG Numbers U/D-Regulated WNN_L25 -------------------------------



## 3.1 Sanity check: double-check comparison direction ---------------------

san.tmp <- DDS_List[[which(DDS_List_Index$CellTypeLevel=="AllCells" & DDS_List_Index$Comparison=="All_Cases")]]
san.tmp <- results(san.tmp, name = "Case_ALS_vs_HC", alpha=0.05)
san.tmp <- data.frame(san.tmp)
san.tmp <- san.tmp[!is.na(san.tmp$padj),]
san.tmp <- san.tmp[san.tmp$padj < 0.05,]
san.tmp <- san.tmp[order(san.tmp$log2FoldChange, decreasing = TRUE),]

best.hit <- rownames(san.tmp)[1]
View(san.tmp[1,])

DDS.tmp <- DDS_List[[which(DDS_List_Index$CellTypeLevel=="AllCells" & DDS_List_Index$Comparison=="All_Cases")]]
plotCounts(DDS.tmp, gene = best.hit, intgroup = "Case", returnData = TRUE, normalized = TRUE) %>% 
  ggplot() + aes(Case, count) + 
  geom_boxplot(aes(fill=Case), outlier.shape = NA) + geom_jitter(size=3,pch=21,width=0.15, fill="#000000", alpha=0.4) + 
  scale_y_continuous(trans="log10")


rm(san.tmp, best.hit, DDS.tmp)



## 3.2 ALSvsHC -------------------------------------------------------------

ggplot(
  Results_WNN_L25
) + 
  aes(CellType, Up_q_0.05_ALS, fill=CellType) + 
  geom_col(col="#000000", fill="#FF5744") + 
  scale_y_continuous(expand=c(0,0,0.2,0), limits=c(0,1800)) + 
  theme_classic() + 
  theme(
    axis.text.x = element_text(angle=43, hjust=1.0), 
    panel.grid.major.y = element_line(color="#DDDDDD")
  ) 

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig4/", 
    "DE_RNA_WNN_L25_nDEG_Up_ALS",
    ".pdf"
  ), 
  width = 6.92, 
  height = 2.833, 
  units="in"
)


ggplot(
  Results_WNN_L25
) + 
  aes(CellType, Down_q_0.05_ALS, fill=CellType) + 
  geom_col(fill="#3354E2", col="#000000") + 
  scale_y_reverse(expand=c(0.1,0,0,0), limits=c(1800,0)) + 
  theme_classic() + 
  theme(
    axis.text.x = element_text(angle=43, hjust=1.0), 
    panel.grid.major.y = element_line(color="#DDDDDD")
  ) 

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig4/", 
    "DE_RNA_WNN_L25_nDEG_Down_ALS", 
    ".pdf"
  ), 
  width = 6.92, 
  height = 2.833, 
  units="in"
)



## 3.3 ALSFTDvsHC ----------------------------------------------------------

ggplot(
  Results_WNN_L25
) + 
  aes(CellType, Up_q_0.05_ALSFTD, fill=CellType) + 
  geom_col(col="#000000", fill="#FF5744") + 
  scale_y_continuous(expand=c(0,0,0.2,0), limits=c(0,1800)) + 
  theme_classic() + 
  theme(
    axis.text.x = element_text(angle=43, hjust=1.0), 
    panel.grid.major.y = element_line(color="#DDDDDD")
  ) 

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig4/", 
    "DE_RNA_WNN_L25_nDEG_Up_ALSFTD",
    ".pdf"
  ), 
  width = 6.92, 
  height = 2.833, 
  units="in"
)


ggplot(
  Results_WNN_L25
) + 
  aes(CellType, Down_q_0.05_ALSFTD, fill=CellType) + 
  geom_col(fill="#3354E2", col="#000000") + 
  scale_y_reverse(expand=c(0.1,0,0,0), limits=c(1800,0)) + 
  theme_classic() + 
  theme(
    axis.text.x = element_text(angle=43, hjust=1.0), 
    panel.grid.major.y = element_line(color="#DDDDDD")
  ) 

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig4/", 
    "DE_RNA_WNN_L25_nDEG_Down_ALSFTD", 
    ".pdf"
  ), 
  width = 6.92, 
  height = 2.833, 
  units="in"
)




### 4.0 Plot CellType Color Rectangles ---------------------------------------

ggplot(
  Results_WNN_L25
) + 
  aes(x=CellType) + 
  geom_col(aes(y=100, fill=CellType), col="#000000") + 
  scale_fill_manual(values=ColDict_WNN_L25) + 
  theme(
    legend.position = "NULL", 
    axis.text.x = element_text(angle=45, hjust=1), 
    panel.background = element_blank(), 
    panel.grid = element_blank()
  )

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig4/", 
    "CellTypeColors_Rectangles", 
    ".pdf"
  )
)




### 5.0 Plot CellType composition --------------------------------------------



## 5.1 ALSvsHC -------------------------------------------------------------

Exc_RORB_ALS <- M0_Meta_ALSvsHC[which(M0_Meta_ALSvsHC$WNN_L25=="Exc_RORB"),]

Exc_RORB_ALS$WNN_L4 <- droplevels(Exc_RORB_ALS$WNN_L4)
table(Exc_RORB_ALS$WNN_L4, Exc_RORB_ALS$Case)

Exc_RORB_ALS %>% 
  group_by(WNN_L4, Case) %>% 
  summarize(nCells=length(CellId)) %>% 
  ggplot() + 
  aes(nCells, Case, fill=WNN_L4) + 
  geom_col(col="black") + 
  scale_x_continuous(expand=c(0,0,0.1,0)) + 
  theme_classic() 

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig4/", 
    "WNN_L25_CellType_Composition_Exc_RORB_ALS", 
    ".pdf"
  ), 
  width = 6.92, 
  height = 2.833, 
  units="in"
)



Exc_LINC00507_ALS <- M0_Meta_ALSvsHC[which(M0_Meta_ALSvsHC$WNN_L25=="Exc_LINC00507"),]

Exc_LINC00507_ALS$WNN_L4 <- droplevels(Exc_LINC00507_ALS$WNN_L4)
table(Exc_LINC00507_ALS$WNN_L4, Exc_LINC00507_ALS$Case)

Exc_LINC00507_ALS %>% 
  group_by(WNN_L4, Case) %>% 
  summarize(nCells=length(CellId)) %>% 
  ggplot() + 
  aes(nCells, Case, fill=WNN_L4) + 
  geom_col(col="black") + 
  scale_x_continuous(expand=c(0,0,0.1,0)) + 
  theme_classic() 

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig4/", 
    "WNN_L25_CellType_Composition_Exc_LINC00507_ALS", 
    ".pdf"
  ), 
  width = 6.92, 
  height = 2.833, 
  units="in"
)

rm(Exc_RORB_ALS, Exc_LINC00507_ALS)

## 5.1 ALSFTDvsHC ----------------------------------------------------------

Exc_RORB_ALSFTD <- M0_Meta_ALSFTDvsHC[which(M0_Meta_ALSFTDvsHC$WNN_L25=="Exc_RORB"),]

Exc_RORB_ALSFTD$WNN_L4 <- droplevels(Exc_RORB_ALSFTD$WNN_L4)
table(Exc_RORB_ALSFTD$WNN_L4, Exc_RORB_ALSFTD$Case)

Exc_RORB_ALSFTD %>% 
  group_by(WNN_L4, Case) %>% 
  summarize(nCells=length(CellId)) %>% 
  ggplot() + 
  aes(nCells, Case, fill=WNN_L4) + 
  geom_col(col="black") + 
  scale_x_continuous(expand=c(0,0,0.1,0)) + 
  theme_classic() 

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig4/", 
    "WNN_L25_CellType_Composition_Exc_RORB_ALSFTD", 
    ".pdf"
  ), 
  width = 6.92, 
  height = 2.833, 
  units="in"
)



Exc_LINC00507_ALSFTD <- M0_Meta_ALSFTDvsHC[which(M0_Meta_ALSFTDvsHC$WNN_L25=="Exc_LINC00507"),]

Exc_LINC00507_ALSFTD$WNN_L4 <- droplevels(Exc_LINC00507_ALSFTD$WNN_L4)
table(Exc_LINC00507_ALSFTD$WNN_L4, Exc_LINC00507_ALSFTD$Case)

Exc_LINC00507_ALSFTD %>% 
  group_by(WNN_L4, Case) %>% 
  summarize(nCells=length(CellId)) %>% 
  ggplot() + 
  aes(nCells, Case, fill=WNN_L4) + 
  geom_col(col="black") + 
  scale_x_continuous(expand=c(0,0,0.1,0)) + 
  theme_classic() 

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig4/", 
    "WNN_L25_CellType_Composition_Exc_LINC00507_ALSFTD", 
    ".pdf"
  ), 
  width = 6.92, 
  height = 2.833, 
  units="in"
)




### 6.0 Plot DEG Numbers U/D-Regulated WNN_L4 --------------------------------



## 6.1 Exc_RORB ------------------------------------------------------------

RORB_Subtypes <- M0_Meta$WNN_L4[M0_Meta$WNN_L25=="Exc_RORB"] |> 
  as.character() |> 
  unique()

Results_Exc_RORB <- DE_Results_SVA[
  DE_Results_SVA$CellTypeLevel=="WNN_L4" & 
    DE_Results_SVA$CellType %in% RORB_Subtypes & 
    DE_Results_SVA$Comparison=="All_Cases"
  ,]

df = data.frame(
  CellType=rep(Results_Exc_RORB$CellType, 2), 
  Type=rep(c("Up", "Down"), each=nrow(Results_Exc_RORB)), 
  nDEGs=c(Results_Exc_RORB$Up_q_0.05_ALS, Results_Exc_RORB$Down_q_0.05_ALS)
)

df$Type <- factor(df$Type, levels=c("Up", "Down"))
ggplot(
  df
) + 
  aes(str_replace(str_replace_all(CellType, "Exc_RORB_", ""), "-113535578-113537604", ""), nDEGs, fill=Type) + 
  geom_bar(col="#000000", position="stack", stat="identity") + 
  scale_y_continuous(expand=c(0,0,0.02,0), limits=c(0,620)) + 
  scale_fill_manual(values=c("Down"="blue", "Up"="red")) + 
  theme_classic() + 
  theme(
    axis.text.x = element_text(angle=43, hjust=1.0), 
    panel.grid.major.y = element_line(color="#DDDDDD")
  ) 

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig4/", 
    "DE_RNA_Exc_RORB_nDEG_Up_ALS",
    ".pdf"
  ), 
  width = 6.12, 
  height = 2.833, 
  units="in"
)

rm(df) 



df = data.frame(
  CellType=rep(Results_Exc_RORB$CellType, 2), 
  Type=rep(c("Up", "Down"), each=nrow(Results_Exc_RORB)), 
  nDEGs=c(Results_Exc_RORB$Up_q_0.05_ALSFTD, Results_Exc_RORB$Down_q_0.05_ALSFTD)
)

df$Type <- factor(df$Type, levels=c("Up", "Down"))
ggplot(
  df
) + 
  aes(str_replace(str_replace_all(CellType, "Exc_RORB_", ""), "-113535578-113537604", ""), nDEGs, fill=Type) + 
  geom_bar(col="#000000", position="stack", stat="identity") + 
  scale_y_continuous(expand=c(0,0,0.02,0), limits=c(0,620)) + 
  scale_fill_manual(values=c("Down"="blue", "Up"="red")) + 
  theme_classic() + 
  theme(
    axis.text.x = element_text(angle=43, hjust=1.0), 
    panel.grid.major.y = element_line(color="#DDDDDD")
  ) 

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig4/", 
    "DE_RNA_Exc_RORB_nDEG_Up_ALSFTD",
    ".pdf"
  ), 
  width = 6.12, 
  height = 2.833, 
  units="in"
)

rm(df) 
rm(RORB_Subtypes, Results_Exc_RORB)



## 6.2 Exc_LINC00507 -------------------------------------------------------

LINC00507_Subtypes <- M0_Meta$WNN_L4[M0_Meta$WNN_L25=="Exc_LINC00507"] |> 
  as.character() |> 
  unique()

Results_Exc_LINC00507 <- DE_Results_SVA[
  DE_Results_SVA$CellTypeLevel=="WNN_L4" & 
    DE_Results_SVA$CellType %in% LINC00507_Subtypes & 
    DE_Results_SVA$Comparison=="All_Cases"
  ,]

df = data.frame(
  CellType=rep(Results_Exc_LINC00507$CellType, 2), 
  Type=rep(c("Up", "Down"), each=nrow(Results_Exc_LINC00507)), 
  nDEGs=c(Results_Exc_LINC00507$Up_q_0.05_ALS, Results_Exc_LINC00507$Down_q_0.05_ALS)
)

df$Type <- factor(df$Type, levels=c("Up", "Down"))
ggplot(
  df
) + 
  aes(str_replace_all(CellType, "Exc_LINC00507_", ""), nDEGs, fill=Type) + 
  geom_bar(col="#000000", position="stack", stat="identity") + 
  scale_y_continuous(expand=c(0,0,0.1,0), limits=c(0,1500)) + 
  scale_fill_manual(values=c("Down"="blue", "Up"="red")) + 
  theme_classic() + 
  theme(
    axis.text.x = element_text(angle=43, hjust=1.0), 
    panel.grid.major.y = element_line(color="#DDDDDD")
  ) 

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig4/", 
    "DE_RNA_Exc_LINC00507_nDEG_Up_ALS",
    ".pdf"
  ), 
  width = 3.12, 
  height = 2.833, 
  units="in"
)

rm(df) 



df = data.frame(
  CellType=rep(Results_Exc_LINC00507$CellType, 2), 
  Type=rep(c("Up", "Down"), each=nrow(Results_Exc_LINC00507)), 
  nDEGs=c(Results_Exc_LINC00507$Up_q_0.05_ALSFTD, Results_Exc_LINC00507$Down_q_0.05_ALSFTD)
)

df$Type <- factor(df$Type, levels=c("Up", "Down"))
ggplot(
  df
) + 
  aes(str_replace_all(CellType, "Exc_LINC00507_", ""), nDEGs, fill=Type) + 
  geom_bar(col="#000000", position="stack", stat="identity") + 
  #  scale_y_continuous(expand=c(0,0,0.2,0), limits=c(0,1800)) + 
  scale_fill_manual(values=c("Down"="blue", "Up"="red")) + 
  theme_classic() + 
  theme(
    axis.text.x = element_text(angle=43, hjust=1.0), 
    panel.grid.major.y = element_line(color="#DDDDDD")
  ) 

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig4/", 
    "DE_RNA_Exc_LINC00507_nDEG_Up_ALSFTD",
    ".pdf"
  ), 
  width = 6.12, 
  height = 2.833, 
  units="in"
)

rm(df) 
rm(LINC00507_Subtypes, Results_Exc_LINC00507)
rm(Exc_RORB_ALSFTD, Exc_LINC00507_ALSFTD, M0_Meta_ALSvsHC, M0_Meta_ALSFTDvsHC, Results_WNN_L25)



### 7.0 All DEGs Heatmap WNN_L2.5   ------------------------------------------



## 7.1 Collect All DEGs ----------------------------------------------------

ind_WNN25 <- which(
  DDS_List_Index$CellTypeLevel=="WNN_L25" & 
    DDS_List_Index$Comparison=="All_Cases"
)

ind=ind_WNN25[12]

ALS_DEGs_WNN25 <- lapply(
  ind_WNN25, 
  FUN=function(ind){
    DDS_List[[ind]] |> 
      results(alpha=0.05, name="Case_ALS_vs_HC") |> 
      data.frame() %>% 
      filter(!is.na(padj)) |> 
      filter(padj<0.05) |> 
      rownames() 
  } 
) |> unlist() |> unique() 

ALSFTD_DEGs_WNN25 <- lapply(
  ind_WNN25, 
  FUN=function(ind){
    DDS_List[[ind]] |> 
      results(alpha=0.05, name="Case_ALS_FTD_vs_HC") |> 
      data.frame() %>% 
      filter(!is.na(padj)) |> 
      filter(padj<0.05) |> 
      rownames() 
  } 
) |> unlist() |> unique() 


table(ALS_DEGs_WNN25 %in% ALSFTD_DEGs_WNN25)
table(ALSFTD_DEGs_WNN25 %in% ALS_DEGs_WNN25)

All_DEGs_WNN_25 <- unique(
  c(
    ALS_DEGs_WNN25, 
    ALSFTD_DEGs_WNN25
  )
)


## Collect L2FoldChanges ---- 

All_DEGs_WNN_25_L2FC <- data.frame(
  Gene=sort(All_DEGs_WNN_25)
)

for(i in 1:length(ind_WNN25)){
  
  tmp <- DDS_List[[ind_WNN25[i]]] |> 
    results(alpha=0.05, name="Case_ALS_vs_HC") |> 
    data.frame() 
  
  All_DEGs_WNN_25_L2FC[[paste0("LFC_ALS_", DDS_List_Index$CellType[ind_WNN25[i]])]] <- tmp$log2FoldChange[match(All_DEGs_WNN_25_L2FC$Gene, rownames(tmp))]
  rm(tmp) 
  
  tmp <- DDS_List[[ind_WNN25[i]]] |> 
    results(alpha=0.05, name="Case_ALS_FTD_vs_HC") |> 
    data.frame() 
  
  All_DEGs_WNN_25_L2FC[[paste0("LFC_ALSFTD_", DDS_List_Index$CellType[ind_WNN25[i]])]] <- tmp$log2FoldChange[match(All_DEGs_WNN_25_L2FC$Gene, rownames(tmp))]
  rm(tmp) 
  
}

rownames(All_DEGs_WNN_25_L2FC) <- All_DEGs_WNN_25_L2FC$Gene
All_DEGs_WNN_25_L2FC %<>% select(!Gene)


heatmap(t(as.matrix(All_DEGs_WNN_25_L2FC)), scale = "col")
Heatmap(t(as.matrix(All_DEGs_WNN_25_L2FC)))

All_DEGs_WNN_25_L2FC_Scaled_Genewise <- (scale(t(All_DEGs_WNN_25_L2FC)))
Heatmap(All_DEGs_WNN_25_L2FC_Scaled_Genewise)


All_DEGs_WNN_25_PVal <- data.frame(
  Gene=sort(All_DEGs_WNN_25)
)

for(i in 1:length(ind_WNN25)){
  
  tmp <- DDS_List[[ind_WNN25[i]]] |> 
    results(alpha=0.05, name="Case_ALS_vs_HC") |> 
    data.frame() 
  
  All_DEGs_WNN_25_PVal[[paste0("LFC_ALS_", DDS_List_Index$CellType[ind_WNN25[i]])]] <- tmp$pvalue[match(All_DEGs_WNN_25_PVal$Gene, rownames(tmp))]
  rm(tmp) 
  
  tmp <- DDS_List[[ind_WNN25[i]]] |> 
    results(alpha=0.05, name="Case_ALS_FTD_vs_HC") |> 
    data.frame() 
  
  All_DEGs_WNN_25_PVal[[paste0("LFC_ALSFTD_", DDS_List_Index$CellType[ind_WNN25[i]])]] <- tmp$pvalue[match(All_DEGs_WNN_25_PVal$Gene, rownames(tmp))]
  rm(tmp) 
  
}


rownames(All_DEGs_WNN_25_PVal) <- All_DEGs_WNN_25_PVal$Gene
All_DEGs_WNN_25_PVal %<>% select(!Gene)


## Collect PAdjs 

All_DEGs_WNN_25_PAdj <- data.frame(
  Gene=sort(All_DEGs_WNN_25)
)

for(i in 1:length(ind_WNN25)){
  
  tmp <- DDS_List[[ind_WNN25[i]]] |> 
    results(alpha=0.05, name="Case_ALS_vs_HC") |> 
    data.frame() 
  
  All_DEGs_WNN_25_PAdj[[paste0("LFC_ALS_", DDS_List_Index$CellType[ind_WNN25[i]])]] <- tmp$padj[match(All_DEGs_WNN_25_PAdj$Gene, rownames(tmp))]
  rm(tmp) 
  
  tmp <- DDS_List[[ind_WNN25[i]]] |> 
    results(alpha=0.05, name="Case_ALS_FTD_vs_HC") |> 
    data.frame() 
  
  All_DEGs_WNN_25_PAdj[[paste0("LFC_ALSFTD_", DDS_List_Index$CellType[ind_WNN25[i]])]] <- tmp$padj[match(All_DEGs_WNN_25_PAdj$Gene, rownames(tmp))]
  rm(tmp) 
  
}

rownames(All_DEGs_WNN_25_PAdj) <- All_DEGs_WNN_25_PAdj$Gene
All_DEGs_WNN_25_PAdj %<>% select(!Gene)

##

All_DEGs_WNN_25_QVal <- as.matrix(All_DEGs_WNN_25_PVal) 
for (i in 1:ncol(All_DEGs_WNN_25_QVal)){
  All_DEGs_WNN_25_QVal[,i] <- p.adjust(All_DEGs_WNN_25_QVal[,i],method="BH")
}

All_DEGs_WNN_25_L2FC_Sign <- as.matrix(All_DEGs_WNN_25_L2FC)
All_DEGs_WNN_25_L2FC_PVal <- as.matrix(All_DEGs_WNN_25_PVal) 
All_DEGs_WNN_25_L2FC_PAdj <- as.matrix(All_DEGs_WNN_25_PAdj) 
All_DEGs_WNN_25_L2FC_Sign_QVal <- as.matrix(All_DEGs_WNN_25_L2FC)

rownames(All_DEGs_WNN_25_L2FC_QVal)
colnames(All_DEGs_WNN_25_L2FC_QVal)

all(rownames(All_DEGs_WNN_25_PVal)==rownames(All_DEGs_WNN_25_QVal))
all(colnames(All_DEGs_WNN_25_PVal)==colnames(All_DEGs_WNN_25_QVal))

plot(as.matrix(All_DEGs_WNN_25_PVal), All_DEGs_WNN_25_QVal, pch=".")


all(rownames(All_DEGs_WNN_25_L2FC)==rownames(All_DEGs_WNN_25_PVal))
all(colnames(All_DEGs_WNN_25_L2FC)==colnames(All_DEGs_WNN_25_PVal))

all(rownames(All_DEGs_WNN_25_L2FC)==rownames(All_DEGs_WNN_25_PAdj))
all(colnames(All_DEGs_WNN_25_L2FC)==colnames(All_DEGs_WNN_25_PAdj))

rownames(All_DEGs_WNN_25_PAdj)
colnames(All_DEGs_WNN_25_PAdj)

for(i in 1:24){
  message(colnames(All_DEGs_WNN_25_PAdj)[i])
  #table(is.na(All_DEGs_WNN_25_PAdj[,i]))
  print(table(All_DEGs_WNN_25_PAdj[,i]<0.05))
  readline("Next?")
}


table(All_DEGs_WNN_25_L2FC_Sign==0)
table(is.na(All_DEGs_WNN_25_L2FC_PAdj))
table(All_DEGs_WNN_25_L2FC_PAdj<0.05)

i=1
for(i in 1:24){
  All_DEGs_WNN_25_L2FC_Sign_QVal[which(is.na(All_DEGs_WNN_25_L2FC_PVal[,i])),i] <- 0 
  All_DEGs_WNN_25_L2FC_Sign_QVal[which(All_DEGs_WNN_25_L2FC_PVal[,i]>=0.05),i] <- 0
}
table(All_DEGs_WNN_25_L2FC_Sign_QVal==0)


for(i in 1:24){
  message(colnames(All_DEGs_WNN_25_L2FC_Sign)[i])
  #table(is.na(All_DEGs_WNN_25_PAdj[,i]))
  print(table(All_DEGs_WNN_25_L2FC_Sign[,i]==0))
  readline("Next?")
}

pheatmap(All_DEGs_WNN_25_L2FC_Sign)  
pheatmap(All_DEGs_WNN_25_L2FC_Sign_QVal) 
rownames(All_DEGs_WNN_25_L2FC_Sign)
All_DEGs_WNN_25_L2FC_Sign2 <- t(scale(t(All_DEGs_WNN_25_L2FC_Sign_QVal)))


gg.df <- melt(All_DEGs_WNN_25_L2FC_Sign_QVal, varnames=c("Gene","Comparison")) 
h <- hclust(dist(All_DEGs_WNN_25_L2FC_Sign_QVal))
gg.df$Gene <- factor(gg.df$Gene, levels=h$labels[h$order])
rm(h)

h <- hclust(dist(t(All_DEGs_WNN_25_L2FC_Sign_QVal)))
gg.df$Comparison <- factor(gg.df$Comparison, levels=h$labels[h$order])
rm(h)

ggplot(gg.df) + 
  aes(Gene, Comparison, fill=value) + 
  geom_tile() + 
  scale_fill_gradientn(colours = c("#0000FF", "white", "#FF0000"), limits=c(-1,1), na.value = "#FFFFFF") + 
  #  scale_alpha_continuous(limit=c(0,3)) + 
  theme_classic() + 
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.line.x = element_blank(), 
    axis.title.x = element_blank(), 
    panel.background = element_blank(), 
    axis.line.y = element_blank(), 
    axis.text.y = element_text(hjust=1), 
    axis.title.y = element_blank(), 
    panel.border = element_rect(color="#000000", size=1, fill = NA)
  )

Case=rep("ALS", ncol(All_DEGs_WNN_25_L2FC_QVal))
Case[grep("ALSFTD", colnames(All_DEGs_WNN_25_L2FC_QVal))] <- "ALS_FTD"
CellTypes <- str_replace_all(str_replace_all(colnames(All_DEGs_WNN_25_L2FC_QVal), "LFC_ALS_", ""), "LFC_ALSFTD_", "")
cols <- unique(Case)
cols <- setNames(rep(NA, length(cols)), nm = cols)
cols <- ColDict_Case[match(names(cols), names(ColDict_Case))]
cols
cols2 <- unique(CellTypes)
cols2 <- setNames(rep(NA, length(cols2)), nm = cols2)
cols2 <- ColDict_WNN_L25[match(names(cols2), names(ColDict_WNN_L25))]
cols2
cols3 <- c(cols, cols2)
cols3

library(circlize)
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

row_ha = rowAnncolsrow_ha = rowAnnotation(
  Case=Case, 
  CellType=CellTypes, 
  col=list(Case=cols, CellType=cols2), 
  gp = gpar(col="black")  
)

Heatmap(
  matrix=t(All_DEGs_WNN_25_L2FC_Sign_QVal), 
  name = "All_DEGs_Sign_QVal", 
  clustering_distance_row = "pearson", 
  clustering_method_rows="complete", 
  clustering_distance_columns  = "manhattan", 
  clustering_method_columns = "average", 
  show_column_dend = FALSE, 
  show_column_names = FALSE, 
  left_annotation = row_ha, 
  col=col_fun, 
  border_gp = gpar(col = "black", lty = 1), 
  use_raster = FALSE
  
)

class(a)
ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig4/", 
    "DE_RNA_WNN_L25_Heatmap_ALS_ALSFTD_L2FCs",
    ".pdf"
  ), 
  width = 14.85, 
  height = 4.82, 
  units="in"
)

ALS_DEGs_WNN25_list <- lapply(
  ind_WNN25, 
  FUN=function(ind){
    DDS_List[[ind]] |> 
      results(alpha=0.05, name="Case_ALS_vs_HC") |> 
      data.frame() %>% 
      filter(!is.na(padj)) |> 
      filter(padj<0.05) |> 
      rownames() 
  } 
) 

ALS_DEGs_WNN25_list_L2FCs <- lapply(
  ind_WNN25, 
  FUN=function(ind){
    DDS_List[[ind]] |> 
      results(alpha=0.05, name="Case_ALS_vs_HC") |> 
      data.frame() %>% 
      filter(!is.na(padj)) |> 
      filter(padj<0.05) |> 
      select(log2FoldChange) 
  } 
) 

ALSFTD_DEGs_WNN25_list <- lapply(
  ind_WNN25, 
  FUN=function(ind){
    DDS_List[[ind]] |> 
      results(alpha=0.05, name="Case_ALS_FTD_vs_HC") |> 
      data.frame() %>% 
      filter(!is.na(padj)) |> 
      filter(padj<0.05) |> 
      rownames() 
  } 
) 

ALSFTD_DEGs_WNN25_list_L2FCs <- lapply(
  ind_WNN25, 
  FUN=function(ind){
    DDS_List[[ind]] |> 
      results(alpha=0.05, name="Case_ALS_FTD_vs_HC") |> 
      data.frame() %>% 
      filter(!is.na(padj)) |> 
      filter(padj<0.05) |> 
      select(log2FoldChange)
  } 
) 

CellTypes_list <- lapply(
  ind_WNN25, 
  FUN=function(ind){
    DDS_List_Index$CellType[[ind]]
  }
) |> unlist()


names(ALS_DEGs_WNN25_list) <- CellTypes_list
names(ALSFTD_DEGs_WNN25_list) <- CellTypes_list 

names(ALS_DEGs_WNN25_list_L2FCs) <- CellTypes_list
names(ALSFTD_DEGs_WNN25_list_L2FCs) <- CellTypes_list 



All_WNN_L25_Overlap_Percent <- c() 
All_WNN_L25_Overlap_PVal <- c() 
All_WNN_L25_Overlap_PVal
NDEG <- c() 
NDEG_Concordant <- c() 
NDEG_NonConcordant <- c() 


for (i in 1:length(ALSFTD_DEGs_WNN25_list)){
  ALS <- ALS_DEGs_WNN25_list[[i]] 
  ALS_FTD <- ALSFTD_DEGs_WNN25_list[[i]] 
  all <- length(unique(c(ALS, ALS_FTD))) 
  overlap <- intersect(ALS, ALS_FTD) 
  tmp <- length(overlap)/all
  All_WNN_L25_Overlap_Percent <- c(All_WNN_L25_Overlap_Percent, tmp) 
  if(length(overlap)==0){
    NDEG <- c(NDEG, NA) 
    NDEG_Concordant <- c( NDEG_Concordant, NA) 
    NDEG_NonConcordant <- c(NDEG_NonConcordant, NA) 
  }else{
    ALS_L2FCs <- ALS_DEGs_WNN25_list_L2FCs[[i]]$log2FoldChange[match(overlap, rownames(ALS_DEGs_WNN25_list_L2FCs[[i]]))]
    ALS_FTD_L2FCs <- ALSFTD_DEGs_WNN25_list_L2FCs[[i]]$log2FoldChange[match(overlap, rownames(ALSFTD_DEGs_WNN25_list_L2FCs[[i]]))]
    NDEG <- c(NDEG, length(overlap)) 
    NDEG_Concordant <- c(NDEG_Concordant, sum(ALS_L2FCs*ALS_FTD_L2FCs > 0))
    NDEG_NonConcordant <- c(NDEG_NonConcordant, length(overlap)- sum(ALS_L2FCs*ALS_FTD_L2FCs > 0)) 
  }
  a <- table(ALS %in% ALS_FTD)["TRUE"]
  b <- table(ALS %in% ALS_FTD)["FALSE"]
  c <- table(ALS_FTD %in% ALS)["FALSE"]
  d <- DE_Results_SVA$NonzeroCounts_ALS[ind_WNN25[i]] 
  if(any(is.na(c(a,b,c,d)))){
    All_WNN_L25_Overlap_PVal <- c(All_WNN_L25_Overlap_PVal, NA)
  } else{
    All_WNN_L25_Overlap_PVal <- c(All_WNN_L25_Overlap_PVal, fisher.test(matrix(c(a,b,c,d),2,2))$p.value)
    
  }
  rm(ALS, ALS_FTD, all, overlap, tmp, a, b, c, d) 
  
}


All_WNN_L25_Overlap_Percent 
All_WNN_L25_Overlap_PVal 
NDEG 
NDEG_Concordant 
NDEG_NonConcordant  
NDEG_Concordant_Percent <- (NDEG_Concordant/(NDEG_Concordant+NDEG_NonConcordant))*All_WNN_L25_Overlap_Percent 
NDEG_NonConcordant_Percent <- All_WNN_L25_Overlap_Percent  - NDEG_Concordant_Percent
NDEG_Concordant_Percent
NDEG_NonConcordant_Percent

(All_WNN_L25_Overlap_Percent*100) |> round(1) -> All_WNN_L25_Overlap_Percent
names(All_WNN_L25_Overlap_Percent) <- names(ALSFTD_DEGs_WNN25_list)
names(All_WNN_L25_Overlap_PVal) <- names(ALSFTD_DEGs_WNN25_list)
names(NDEG_Concordant_Percent) <- names(ALSFTD_DEGs_WNN25_list)
names(NDEG_NonConcordant_Percent) <- names(ALSFTD_DEGs_WNN25_list)

All_WNN_L25_Overlap_Percent
All_WNN_L25_Overlap_PVal 
All_WNN_L25_Overlap_QVal <- p.adjust(All_WNN_L25_Overlap_PVal, method="fdr")
All_WNN_L25_Overlap_QVal |> round(3)

df <- data.frame (
  CellType=names(All_WNN_L25_Overlap_Percent), 
  OverlapPercent <- All_WNN_L25_Overlap_Percent, 
  OverlapFisherPVal <- All_WNN_L25_Overlap_PVal, 
  OverlapFisherQCal <- All_WNN_L25_Overlap_QVal,
)

all(names(All_WNN_L25_Overlap_Percent)==names(NDEG_Concordant_Percent))
all(names(All_WNN_L25_Overlap_Percent)==names(NDEG_NonConcordant_Percent))

df$CellType <- factor(
  df$CellType, 
  levels=rev(
    c(
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
    )
  )
)

ggplot(df) + 
  aes(OverlapPercent, CellType, fill=CellType) + 
  geom_bar(stat="identity", col="#000000") + 
  scale_x_continuous(expand=c(0,0,0.1,0)) + 
  scale_fill_manual(values=ColDict_WNN_L25) + 
  ylab("") + 
  theme_classic() + 
  theme(
    legend.position = "Null"
  )

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig4/", 
    "DE_RNA_WNN_L25_ALS_ALSFTD_DEG_Overlap_Percent",
    ".pdf"
  ), 
  width = 3.10, 
  height = 6.21, 
  units="in"
)

df2 <- data.frame(
  CellType=rep(names(All_WNN_L25_Overlap_Percent), 2), 
  NDEG_Type=rep(c("Concordant", "NonConcordant"), each=length(All_WNN_L25_Overlap_Percent)), 
  NDEG_Percent=c(NDEG_Concordant_Percent, NDEG_NonConcordant_Percent)
)

df2$CellType <- factor(
  df2$CellType, 
  levels=rev(
    c(
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
    )
  )
)


df2$NDEG_Type <- factor(df2$NDEG_Type, levels=rev(c("Concordant", "NonConcordant")))
df2$CellType_Concordance <- paste0(df2$CellType, "_", df2$NDEG_Type)
ColDict_Concordant <- ColDict_WNN_L25
names(ColDict_Concordant) <- paste0(names(ColDict_Concordant), "_", "Concordant")
ColDict_Concordant

ColDict_NonConcordant <- rep("#FFFFFF", length(ColDict_Concordant))
names(ColDict_NonConcordant) <- paste0(names(ColDict_WNN_L25), "_", "NonConcordant")
ColDict_NonConcordant
ColDict_WNN_L25_Concordance <- c(
  ColDict_Concordant, 
  ColDict_NonConcordant
)
ColDict_WNN_L25_Concordance

df2$CellType_Concordance <- factor(
  df2$CellType_Concordance, 
  levels=
    c(
      "Exc_RORB_NonConcordant", 
      "Exc_RORB_Concordant",
      "Exc_LINC00507_NonConcordant", 
      "Exc_LINC00507_Concordant", 
      "Exc_THEMIS_NonConcordant",
      "Exc_THEMIS_Concordant",
      "Exc_FEZF2_NonConcordant", 
      "Exc_FEZF2_Concordant", 
      "Inh_LAMP5_PAX6_NonConcordant", 
      "Inh_LAMP5_PAX6_Concordant", 
      "Inh_TAFA1_VIP_NonConcordant", 
      "Inh_TAFA1_VIP_Concordant", 
      "Inh_SST_NonConcordant", 
      "Inh_SST_Concordant", 
      "Inh_PVALB_NonConcordant", 
      "Inh_PVALB_Concordant", 
      "OPC_NonConcordant", 
      "OPC_Concordant", 
      "Oligodendrocytes_NonConcordant", 
      "Oligodendrocytes_Concordant", 
      "Astrocytes_NonConcordant",         
      "Astrocytes_Concordant",                                        
      "Microglia_NonConcordant",
      "Microglia_Concordant"
    )
  
)
library(ggpattern)
ggplot(df2) + 
  aes(NDEG_Percent*100, CellType, fill=CellType_Concordance, pattern_density=NDEG_Type) + 
  geom_bar_pattern(stat="identity", col="#000000", pattern="stripe") + 
  scale_x_continuous(expand=c(0,0,0.1,0)) + 
  scale_fill_manual(values=ColDict_WNN_L25_Concordance) + 
  scale_pattern_density_manual(values = c(Concordant = 0, NonConcordant=0.3)) + 
  ylab("") + 
  theme_classic() + 
  theme(
    legend.position = "Null"
  )

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig4/", 
    "DE_RNA_WNN_L25_ALS_ALSFTD_DEG_Overlap_Concordant_Percent",
    ".pdf"
  ), 
  width = 3.10, 
  height = 6.21, 
  units="in"
)

DDS_NoSva <- qread("../Data/DE/WNN/AllCase/DDS_DESEqed_list.qrds", nthreads = nthr)

lapply(ALS_DEGs_WNN25_list, FUN=function(x){grep("NRXN3", x, value=TRUE)}) |> unlist()
lapply(ALS_DEGs_WNN25_list, FUN=function(x){grep("EPH", x, value=TRUE)}) |> unlist()
lapply(ALS_DEGs_WNN25_list, FUN=function(x){grep("MT-", x, value=TRUE)}) |> unlist()
lapply(ALS_DEGs_WNN25_list, FUN=function(x){grep("ERBB4", x, value=TRUE)}) |> unlist()
lapply(ALS_DEGs_WNN25_list, FUN=function(x){grep("LINC00507", x, value=TRUE)}) |> unlist()
lapply(ALS_DEGs_WNN25_list, FUN=function(x){grep("THEMIS", x, value=TRUE)}) |> unlist()
lapply(ALS_DEGs_WNN25_list, FUN=function(x){grep("FREM3", x, value=TRUE)}) |> unlist()

DDS_NoSva <- qread("../Data/DE/WNN/AllCase/DDS_DESEqed_list.qrds", nthreads = nthr)

ALS_DEGs_WNN25_NoSVA <- lapply(
  ind_WNN25, 
  FUN=function(ind){
    DDS_NoSva[[ind]] |> 
      results(alpha=0.05, name="Case_ALS_vs_HC") |> 
      data.frame() %>% 
      #      filter(!is.na(padj)) |> 
      filter(pvalue<0.05) |> 
      rownames() 
  } 
) 

ALSFTD_DEGs_WNN25_NoSva <- lapply(
  ind_WNN25, 
  FUN=function(ind){
    DDS_NoSva[[ind]] |> 
      results(alpha=0.05, name="Case_ALS_FTD_vs_HC") |> 
      data.frame() %>% 
      #      filter(!is.na(padj)) |> 
      filter(pvalue<0.05) |> 
      rownames() 
  } 
) 

names(ALS_DEGs_WNN25_NoSVA) <- CellTypes_list
names(ALSFTD_DEGs_WNN25_NoSva) <- CellTypes_list  

lapply(ALS_DEGs_WNN25_NoSVA, FUN=function(x){grep("ERBB4", x, value=TRUE)}) |> unlist()
lapply(ALS_DEGs_WNN25_NoSVA, FUN=function(x){grep("FREM3", x, value=TRUE)}) |> unlist()
lapply(ALS_DEGs_WNN25_NoSVA, FUN=function(x){grep("NRXN3", x, value=TRUE)}) |> unlist()

lapply(ALS_DEGs_WNN25_list, FUN=function(x){grep("FRMD4B", x, value=TRUE)}) |> unlist()
lapply(ALS_DEGs_WNN25_NoSVA, FUN=function(x){grep("FRMD4B", x, value=TRUE)}) |> unlist()

lapply(ALSFTD_DEGs_WNN25_list, FUN=function(x){grep("LNX2", x, value=TRUE)}) |> unlist()
lapply(ALSFTD_DEGs_WNN25_NoSva, FUN=function(x){grep("LNX2", x, value=TRUE)}) |> unlist()

# Compare to Pseudobulk 

All_res_ALS <- results(DDS_List[[1]], alpha=0.05, name="Case_ALS_vs_HC")
All_res_ALS <- All_res_ALS[!is.na(All_res_ALS$padj),]        
All_res_ALS <- All_res_ALS[All_res_ALS$padj<0.05,]        
All_res_ALS <- rownames(All_res_ALS)        

table(ALS_DEGs_WNN25 %in% All_res_ALS)         
table(All_res_ALS %in% All_res_ALS)

All_res_ALS_FTD <- results(DDS_List[[1]], alpha=0.05, name="Case_ALS_FTD_vs_HC")
All_res_ALS_FTD <- All_res_ALS_FTD[!is.na(All_res_ALS_FTD$padj),]        
All_res_ALS_FTD <- All_res_ALS_FTD[All_res_ALS_FTD$padj<0.05,]        
All_res_ALS_FTD <- rownames(All_res_ALS_FTD)     

table(All_res_ALS_FTD %in% ALSFTD_DEGs_WNN25)
table(ALSFTD_DEGs_WNN25 %in% All_res_ALS_FTD)

All_res <- unique(c(All_res_ALS, All_res_ALS_FTD))
table(All_res %in% All_DEGs_WNN_25)
table(All_DEGs_WNN_25 %in% All_res)

Pineda_Sheets <- readxl::excel_sheets("/home/genomics/Bioinfo_Data/Datasets/Single_Cell/Single-cell/Pineda_ALS_scRNA_M1/Pineda_DEGs/6_1_deg_res_BA4_SALS.xlsx") 
Pineda_Results <- list() 
for (i in 1:length(Pineda_Sheets)){
  Pineda_Results[[Pineda_Sheets[i]]] <- readxl::read_excel("/home/genomics/Bioinfo_Data/Datasets/Single_Cell/Single-cell/Pineda_ALS_scRNA_M1/Pineda_DEGs/6_1_deg_res_BA4_SALS.xlsx", sheet = Pineda_Sheets[i])
}
Pineda_Res_Filtered <- lapply(
  Pineda_Results, 
  FUN=function(x){
    x %>% 
      filter(!is.na(padj)) %>% 
      filter(padj<0.05) %>% 
      select(Gene, logFC)
  }
)

Pineda_Res_NoNAs <- lapply(
  Pineda_Results, 
  FUN=function(x){
    x$logFC[match(Pineda_All_DEGs, x$Gene)]
  }
)

Pineda_DEGs <- lapply(
  Pineda_Res_Filtered, 
  FUN=function(x){
    x %>% 
      data.frame() %>% 
      select(Gene) %>% 
      unlist %>% 
      unname() 
  }
)
Pineda_DEGs

lapply(Pineda_DEGs, FUN=function(x){grep("ERBB4",x, value=TRUE)})

Pineda_All_DEGs <- unique(unlist(Pineda_DEGs))

Pineda_Common <- lapply(
  Pineda_Res_Filtered, 
  FUN=function(x){
    x$logFC[match(Pineda_All_DEGs, x$Gene)]
  }
)

Pineda_Common_Hm <- as.matrix(data.frame(Pineda_Common))
Pineda_Common_Hm[is.na(Pineda_Common_Hm)] <- 0
Heatmap(as.matrix(Pineda_Common_Hm))

Pineda_Res_NoNAs <- as.matrix(data.frame(Pineda_Res_NoNAs))
Pineda_Res_NoNAs[is.na(Pineda_Res_NoNAs)] <- 0
Heatmap(Pineda_Res_NoNAs)

plot(Pineda_Res_NoNAs[,1], Pineda_Res_NoNAs[,2])
plot(Pineda_Res_NoNAs[,2], Pineda_Res_NoNAs[,3])
plot(Pineda_Res_NoNAs[,3], Pineda_Res_NoNAs[,4])
cor.test(Pineda_Res_NoNAs[,3], Pineda_Res_NoNAs[,4])
plot(Pineda_Res_NoNAs[,4], Pineda_Res_NoNAs[,5])
abline(h=0, col="red")
abline(v=0, col="red")
table(Pineda_Res_NoNAs[,1]*Pineda_Res_NoNAs[,2]>0)
cor(Pineda_Res_NoNAs) |> mean()

table(Pineda_All_DEGs %in% All_DEGs_WNN_25)
table(All_DEGs_WNN_25 %in% Pineda_All_DEGs)

We_Res_ALS <- lapply(
  DDS_List[ind_WNN25], 
  FUN=function(x){
    results(x, alpha=0.05, name="Case_ALS_vs_HC") %>% 
      data.frame %>% 
      rownames_to_column(var="Gene") %>% 
      as_tibble()
  }
)
names(We_Res_ALS) <- DE_Results_SVA$CellType[ind_WNN25]

We_DEGs_ALS <- lapply(
  We_Res_ALS, 
  FUN=function(x){
    x %>% 
      filter(!is.na(padj)) %>% 
      filter(padj<0.05) %>% 
      select(Gene, log2FoldChange, padj)
  }
)
We_DEGs_ALS

We_Common_DEGs <- lapply(
  We_DEGs_ALS, 
  FUN=function(x){
    x %>% select(Gene)
  }
) %>% unlist() %>% unique

Rand_DEGs <- lapply(
  DDS_List[ind_WNN25], 
  FUN=function(x){
    results(x, alpha=0.05, name="Case_ALS_vs_HC") %>% data.frame() %>% rownames_to_column(var="Gene")
  }
)

Rand_DEGs_Genes <- lapply(
  Rand_DEGs, 
  FUN=function(x){
    x %>% select(Gene)
  }
)


Rand_DEGs_Genes <- Rand_DEGs_Genes %>% unlist() %>% unique() %>% sample(size=1315)

table(All_DEGs_WNN_25 %in% We_Common_DEGs)

We_Res <- lapply(
  We_Res_ALS, 
  FUN=function(x){
    x$log2FoldChange[match(ALS_DEGs_WNN25, x$Gene)]
  }
)

Rand_DEGs <- lapply(
  Rand_DEGs, 
  FUN=function(x){
    x %>% as.tibble
  }
)

We_Res_Rand <- lapply(
  Rand_DEGs, 
  FUN=function(x){
    x$log2FoldChange[match(Rand_DEGs_Genes, x$Gene)]
  }
)

We_Res <- data.frame(We_Res) %>% as.matrix()
We_Res[is.na(We_Res)] <- 0 

We_Res_Rand <- data.frame(We_Res_Rand) %>% as.matrix()
We_Res_Rand[is.na(We_Res_Rand)] <- 0 

plot(We_Res[,1], We_Res[,2])
plot(We_Res[,2], We_Res[,3])
plot(We_Res[,3], We_Res[,4])

cor(We_Res)
mean(cor(We_Res))

cor(We_Res_Rand)
mean(cor(We_Res_Rand))

DDS_List2 <- DDS_List 
ind_WNN252=ind_WNN25

Get_Rand_Multivar_Corr <- function(DDS_List=DDS_List2, ind_WNN25=ind_WNN252, seed=314159265){
  
  Rand_DEGs <- lapply(
    DDS_List[ind_WNN25], 
    FUN=function(x){
      results(x, alpha=0.05, name="Case_ALS_vs_HC") %>% data.frame() %>% rownames_to_column(var="Gene")
    }
  )
  
  Rand_DEGs_Genes <- lapply(
    Rand_DEGs, 
    FUN=function(x){
      x %>% select(Gene)
    }
  )
  set.seed(seed)
  Rand_DEGs_Genes <- Rand_DEGs_Genes %>% unlist() %>% unique() %>% sample(size=1315)
  
  Rand_DEGs <- lapply(
    Rand_DEGs, 
    FUN=function(x){
      x %>% as.tibble
    }
  )
  We_Res_Rand <- lapply(
    Rand_DEGs, 
    FUN=function(x){
      x$log2FoldChange[match(Rand_DEGs_Genes, x$Gene)]
    }
  )
  
  
  We_Res_Rand <- data.frame(We_Res_Rand) %>% as.matrix()
  We_Res_Rand[is.na(We_Res_Rand)] <- 0 
  
  return(
    mean(return_Upper_triangle_Values(cor(We_Res_Rand)))
  )
  
}

Ps <- lapply(
  c(1:10), 
  FUN=function(x){
    Get_Rand_Multivar_Corr(seed=x)
  }
)
Ps %<>% unlist()
summary(Ps)
hist(Ps,br=20)
plot(density(Ps))
t.test(Ps, mu = mean(return_Upper_triangle_Values(cor(We_Res))), alternative = "less")
plot(density(cor(We_Res)))

return_Upper_triangle_Values <- function(a){
  b <- c()
  for(i in 1:(nrow(a)-1)){
    b <- c(b, a[i,c((i+1):ncol(a))])
  }
  return(b)
}
mean(return_Upper_triangle_Values(cor(We_Res)))
mean(return_Upper_triangle_Values(cor(Pineda_Res_NoNAs)))
  