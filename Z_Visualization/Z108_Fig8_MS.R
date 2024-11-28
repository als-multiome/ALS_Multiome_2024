
source("~/ALS_Brain_Multiome.Rcfg") 

library(qs) 
library(tidyverse) 
library(Seurat)
library(reshape2)

qs::set_SeuratObjectqs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)

M0_RNA <- qread(
  "../Data/SeuratObjects/M0_RNA.qrds", nthr=nthr
)

ALS_ALSFTD_RNA_WNN_L25_Signature <- qread(
  "../Data/Annotations/Signatures/RNA/Signature_RNA_ALS_ALSFTD_WNN_L25.qrds", 
  nthr=nthr
)

M0_RNA <- AddModuleScore(
  M0_RNA, 
  features = list(ALS_ALSFTD_RNA_WNN_L25_Signature$ID_10X), 
  ctrl = 1000, 
  name = "ALS_ALSFTD_RNA_WNN_L25_Signature_Score"
)

L2FC_Shrink_Results_12SVs_RNA_ALS <- qread("../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/12_SVs/L2FC_Shrink_Results_12SVs_ALS.qrds")
L2FC_Shrink_Results_12SVs_RNA_ALSFTD <- qread("../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/12_SVs/L2FC_Shrink_Results_12SVs_ALSFTD.qrds")
L2FC_Shrink_Results_12SVs_RNA_Index <- qread("../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/SVA/12_SVs/DESeq_Results_12SVs_Index.qrds")
ind <- which(
  L2FC_Shrink_Results_12SVs_RNA_Index$CellTypeLevel=="WNN_L25" & 
    L2FC_Shrink_Results_12SVs_RNA_Index$Comparison=="All_Cases"
)

L2FCs_ALS <- lapply(
  L2FC_Shrink_Results_12SVs_RNA_ALS[ind], 
  function(x){
    
    y <- x 
    y$log2FoldChange[is.na(y$padj)] <- 0 
    y$log2FoldChange[y$padj>=0.05] <- 0 
    return(
      y$log2FoldChange[match(ALS_ALSFTD_RNA_WNN_L25_Signature$ID_10X, rownames(y))]
    )
  }
)

DEG_Logical_ALS <- lapply(
  L2FC_Shrink_Results_12SVs_RNA_ALS[ind], 
  function(x){
    
    y <- x 
    y$Sign <- FALSE
    y$Sign[y$padj<0.05] <- TRUE
    return(
      y$Sign[match(ALS_ALSFTD_RNA_WNN_L25_Signature$ID_10X, rownames(y))]
    )
  }
)
lapply(
  DEG_Logical_ALS, 
  function(x){
    return(
      table(x, useNA="always")
    )
  }
)


L2FCs_ALSFTD <- lapply(
  L2FC_Shrink_Results_12SVs_RNA_ALSFTD[ind], 
  function(x){
    
    y <- x 
    y$log2FoldChange[is.na(y$padj)] <- 0 
    y$log2FoldChange[y$padj>=0.05] <- 0 
    return(
      y$log2FoldChange[match(ALS_ALSFTD_RNA_WNN_L25_Signature$ID_10X, rownames(y))]
    )
  }
)


DEG_Logical_ALSFTD <- lapply(
  L2FC_Shrink_Results_12SVs_RNA_ALSFTD[ind], 
  function(x){
    
    y <- x 
    y$Sign <- FALSE
    y$Sign[y$padj<0.05] <- TRUE
    return(
      y$Sign[match(ALS_ALSFTD_RNA_WNN_L25_Signature$ID_10X, rownames(y))]
    )
  }
)
lapply(
  DEG_Logical_ALSFTD, 
  function(x){
    return(
      table(x, useNA="always")
    )
  }
) 

names.tmp  <- names(L2FCs_ALS[[1]])
sapply(
  c(L2FCs_ALS, L2FCs_ALSFTD),  
  function(x) all(names(x)==names.tmp)
) %>% all() 
sapply(
  c(DEG_Logical_ALS, DEG_Logical_ALSFTD),  
  function(x) all(names(x)==names.tmp)
) %>% all()

sapply(
  c(L2FCs_ALS, L2FCs_ALSFTD), 
  function(x){
    sum(is.na(x))
  }
)
sapply(
  c(DEG_Logical_ALS, DEG_Logical_ALSFTD), 
  function(x){
    sum(is.na(x))
  }
)

list.tmp <- lapply(  
  c(L2FCs_ALS,L2FCs_ALSFTD),  
  function(x){
    x[which(is.na(x))] <- 0 
    return(x)
  }
) 

list.tmp_ALS <- lapply(  
  c(L2FCs_ALS),  
  function(x){
    x[which(is.na(x))] <- 0 
    return(x)
  }
) 

list.tmp_ALSFTD <- lapply(  
  c(L2FCs_ALSFTD),  
  function(x){
    x[which(is.na(x))] <- 0 
    return(x)
  }
) 

sapply(
  list.tmp, 
  function(x){
    sum(is.na(x))
  }
)

DEG_In_N <- rep(0, length(DEG_Logical_ALS[[1]]))
DEG_In_N_ALS <- rep(0, length(DEG_Logical_ALS[[1]]))
DEG_In_N_ALSFTD <- rep(0, length(DEG_Logical_ALSFTD[[1]]))

sapply(
  DEG_Logical_ALS, 
  function(x){
    DEG_In_N <<- DEG_In_N + case_when(x==TRUE ~ 1, x==FALSE ~ 0) 
    DEG_In_N_ALS <<- DEG_In_N_ALS + case_when(x==TRUE ~ 1, x==FALSE ~ 0) 
  }
) 

sapply(
  DEG_Logical_ALSFTD, 
  function(x){
    DEG_In_N <<- DEG_In_N + case_when(x==TRUE ~ 1, x==FALSE ~ 0) 
    DEG_In_N_ALSFTD <<- DEG_In_N_ALSFTD + case_when(x==TRUE ~ 1, x==FALSE ~ 0) 
  }
) 

summary(DEG_In_N)
summary(DEG_In_N_ALS) 
summary(DEG_In_N_ALSFTD)

all(
  DEG_In_N==
    DEG_In_N_ALS + DEG_In_N_ALSFTD
)

Gene_Up <- rep(0, length(list.tmp[[1]]))
Gene_Down <- rep(0, length(list.tmp[[1]]))

sapply(
  list.tmp, 
  function(x){
    Gene_Up <<- Gene_Up + case_when(x>0 ~ x, x <=0 ~ 0) 
    Gene_Down <<- Gene_Down + case_when(x >= 0 ~ 0, x <0 ~ x) 
  }
) 
summary(Gene_Up)
summary(Gene_Down)

barplot(sort(Gene_Up))
barplot(sort(Gene_Down))



Gene_Up_ALS <- rep(0, length(list.tmp_ALS[[1]]))
Gene_Down_ALS <- rep(0, length(list.tmp_ALS[[1]]))

sapply(
  list.tmp_ALS, 
  function(x){
    Gene_Up_ALS <<- Gene_Up_ALS + case_when(x>0 ~ x, x <=0 ~ 0) 
    Gene_Down_ALS <<- Gene_Down_ALS + case_when(x >= 0 ~ 0, x <0 ~ x) 
  }
) 
summary(Gene_Up_ALS)
summary(Gene_Down_ALS)

barplot(sort(Gene_Up_ALS))
barplot(sort(Gene_Down_ALS))



Gene_Up_ALSFTD <- rep(0, length(list.tmp_ALSFTD[[1]]))
Gene_Down_ALSFTD <- rep(0, length(list.tmp_ALSFTD[[1]]))

sapply(
  list.tmp_ALSFTD, 
  function(x){
    Gene_Up_ALSFTD <<- Gene_Up_ALSFTD + case_when(x>0 ~ x, x <=0 ~ 0) 
    Gene_Down_ALSFTD <<- Gene_Down_ALSFTD + case_when(x >= 0 ~ 0, x <0 ~ x) 
  }
) 
summary(Gene_Up_ALSFTD)
summary(Gene_Down_ALSFTD)

barplot(sort(Gene_Up_ALSFTD))
barplot(sort(Gene_Down_ALSFTD))


df <- data.frame(
  Gene = names.tmp, 
  Signif_In = ALS_ALSFTD_RNA_WNN_L25_Signature$Signif_In[match(names.tmp, ALS_ALSFTD_RNA_WNN_L25_Signature$ID_10X)], 
  Sign_In_N_All = DEG_In_N, 
  Sign_In_N_ALS = DEG_In_N_ALS, 
  Sign_In_N_ALSFTD = DEG_In_N_ALSFTD,
  Gene_Up_All = Gene_Up, 
  Gene_Down_All = Gene_Down, 
  Gene_Up_ALS = Gene_Up_ALS, 
  Gene_Down_ALS = Gene_Down_ALS, 
  Gene_Up_ALSFTD = Gene_Up_ALSFTD, 
  Gene_Down_ALSFTD = Gene_Down_ALSFTD
)

Gene_order <- df$Gene[order(df$Gene_Up_All + df$Gene_Down_All)]
df <- df[order(df$Gene_Up_All + df$Gene_Down_All),]
df$Gene <- factor(
  df$Gene, 
  levels= Gene_order
)

table(
  round(df$Gene_Up_All,1)==round(df$Gene_Up_ALS + df$Gene_Up_ALSFTD,1)
) 

table(
  round(df$Gene_Down_All,1)==round(df$Gene_Down_ALS + df$Gene_Down_ALSFTD,1)
)



df %>% 
  filter(Signif_In=="Both") %>% 
  dplyr::select(
    Gene, 
    Gene_Up_ALS, 
    Gene_Up_ALSFTD, 
    Gene_Down_ALS, 
    Gene_Down_ALSFTD
  ) %>% 
  melt() %>% 
  mutate(Disease = str_split(variable, "_",simplify=TRUE)[,3]) %>% 
  mutate(Direction = str_split(variable, "_", simplify=TRUE)[,2]) %>% 
  ggplot() + 
    aes(
      x = Gene, 
      y = value, 
      fill = interaction(Disease, Direction)
    ) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    values = c(
      "ALS.Down" = "#97cef7", 
      "ALS.Up" = "#f7ba97",
      "ALSFTD.Down" = "#598cb3", 
      "ALSFTD.Up" = "#db621d"
    )
  ) +
  geom_hline(yintercept = 0, color = "black") +
  scale_y_continuous(limits=c(-5,5), breaks=c(-2.5, 0, 2.5)) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title = element_blank(), 
    legend.title = element_blank(), 
    axis.line = element_blank(), 
    panel.grid.major.y = element_line(color = "grey80"),  # Major y-axis grid lines
    panel.grid.minor.y = element_blank(), 
    panel.grid.major.x = element_blank(),  # Remove x-axis major grid lines (optional)
    panel.grid.minor.x = element_blank()   # Remove x-axis minor grid lines (optional)
  ) 

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig8/", 
    "Genes_Barplot_ALS_ALSFTD_Both", 
    ".pdf"
  ), 
  width = 9, 
  height = 9,
  units = "in"
)

summary(df$Sign_In_N_All)
summary(df$Sign_In_N_All[df$Signif_In=="Both"])
df %>% 
  dplyr::select(Gene, Signif_In, Sign_In_N_All) %>% 
  pivot_longer(cols=starts_with("Sign_In"), names_to = "Variable", values_to="Value" ) %>% 
  filter(Signif_In == "Both") %>% 
  ggplot() + 
  aes(x = Gene, y = Variable, fill = factor(Value, levels=1:12)) +
  geom_tile() +                          
  geom_hline(yintercept = 0, color = "black") +
  scale_fill_viridis_d(option="magma", direction = -1, n=12) + 
  theme_classic() +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title = element_blank(), 
    legend.title = element_blank(), 
    axis.line = element_blank(), 
    panel.grid.major.y = element_line(color = "grey80"), 
    panel.grid.minor.y = element_blank(), 
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank()   
  ) 




FANS <- qread("../Data/DE/FANS/FANS_DE_TDP43_MAST.qrds")
ATAC_Peaks <- qread("../Data/Annotations/Genomic_Regions/M0_ATAC_Peak_Links.qrds")
ATAC <- qread("../Data/Annotations/Signatures/ATAC/Signature_ATAC_ALS_ALSFTD_WNN_L25.qrds")
ATAC_Genes <- unique(ATAC_Peaks$gene[match(ATAC$ID2, ATAC_Peaks$peak)])

df %>% 
  dplyr::select(Gene) %>% 
  mutate(In_ATAC = Gene %in% ATAC_Genes) %>%
  ggplot() + 
  aes(x = Gene, y = 1, fill = factor(In_ATAC)) +
  geom_tile() +                          
  geom_hline(yintercept = 0, color = "black") +
  scale_fill_manual(values=c("TRUE"="#000000", "FALSE"="#FFFFFF")) + 
  theme_classic() +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title = element_blank(), 
    legend.title = element_blank(), 
    axis.line = element_blank(), 
    panel.grid.major.y = element_line(color = "grey80"), 
    panel.grid.minor.y = element_blank(), 
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank()   
  ) 



df %>% 
  filter(Signif_In %in% c("ALS", "Both")) %>% 
  dplyr::select(
    Gene, 
    Gene_Up_ALS, 
    Gene_Up_ALSFTD, 
    Gene_Down_ALS, 
    Gene_Down_ALSFTD
  ) %>% 
  melt() %>% 
  mutate(Disease = str_split(variable, "_",simplify=TRUE)[,3]) %>% 
  mutate(Direction = str_split(variable, "_", simplify=TRUE)[,2]) %>% 
  ggplot() + 
  aes(
    x = Gene, 
    y = value, 
    fill = interaction(Disease, Direction)
  ) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    values = c(
      "ALS.Down" = "#97cef7", 
      "ALS.Up" = "#f7ba97",
      "ALSFTD.Down" = "#598cb3", 
      "ALSFTD.Up" = "#db621d"
    )
  ) +
  geom_hline(yintercept = 0, color = "black") +
  scale_y_continuous(limits=c(-5,5), breaks=c(-2.5, 0, 2.5)) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title = element_blank(), 
    legend.title = element_blank(), 
    axis.line = element_blank(), 
    panel.grid.major.y = element_line(color = "grey80"),  # Major y-axis grid lines
    panel.grid.minor.y = element_blank(), 
    panel.grid.major.x = element_blank(),  # Remove x-axis major grid lines (optional)
    panel.grid.minor.x = element_blank()   # Remove x-axis minor grid lines (optional)
  ) 

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig8/", 
    "Genes_Barplot_ALS_Only", 
    ".pdf"
  ), 
  width = 9, 
  height = 9,
  units = "in"
)

df %>% 
  #filter(Gene %in% c("NRXN3", "EPHA4", "OPTN", "TBK1", "TARDBP", "SOD1", "NEFL", "FUS", "C0orf72", "KIF5C")) %>% 
  filter(Gene %in% ATAC_Genes[which(! ATAC_Genes %in% FANS$gene[FANS$fdr<0.05])]) %>% 
  dplyr::select(
    Gene, 
    Gene_Up_ALS, 
    Gene_Up_ALSFTD, 
    Gene_Down_ALS, 
    Gene_Down_ALSFTD
  ) %>% 
  melt() %>% 
  mutate(Disease = str_split(variable, "_",simplify=TRUE)[,3]) %>% 
  mutate(Direction = str_split(variable, "_", simplify=TRUE)[,2]) %>% 
  ggplot() + 
  aes(
    x = Gene, 
    y = value, 
    fill = interaction(Disease, Direction)
  ) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    values = c(
      "ALS.Down" = "#97cef7", 
      "ALS.Up" = "#f7ba97",
      "ALSFTD.Down" = "#598cb3", 
      "ALSFTD.Up" = "#db621d"
    )
  ) +
  geom_hline(yintercept = 0, color = "black") +
  scale_y_continuous(limits=c(-5,5), breaks=c(-2.5, 0, 2.5)) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title = element_blank(), 
    legend.title = element_blank(), 
    axis.line = element_blank(), 
    panel.grid.major.y = element_line(color = "grey80"),  # Major y-axis grid lines
    panel.grid.minor.y = element_blank(), 
    panel.grid.major.x = element_blank(),  # Remove x-axis major grid lines (optional)
    panel.grid.minor.x = element_blank()   # Remove x-axis minor grid lines (optional)
  ) 

summary(df$Sign_In_N_ALS[df$Signif_In=="ALS"])
df %>% 
  filter(Signif_In=="ALS") %>% 
  dplyr::select(
    Gene, Sign_In_N_ALS
  ) %>% 
ggplot() + 
  aes(
    x = Gene, 
    y = 1, 
    fill = factor(Sign_In_N_ALS)
  ) + 
  geom_col() + 
  geom_hline(yintercept = 0, color = "black") +
  scale_fill_viridis_d(option="magma", direction = -1, n=12) + 
  theme_classic() +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title = element_blank(), 
    legend.title = element_blank(), 
    axis.line = element_blank(), 
    panel.grid.major.y = element_line(color = "grey80"),  # Major y-axis grid lines
    panel.grid.minor.y = element_blank(), 
    panel.grid.major.x = element_blank(),  # Remove x-axis major grid lines (optional)
    panel.grid.minor.x = element_blank()   # Remove x-axis minor grid lines (optional)
  ) 

df %>% 
  filter(Signif_In %in% c("ALSFTD")) %>% 
  dplyr::select(
    Gene, 
    Gene_Up_ALS, 
    Gene_Up_ALSFTD, 
    Gene_Down_ALS, 
    Gene_Down_ALSFTD
  ) %>% 
  melt() %>% 
  mutate(Disease = str_split(variable, "_",simplify=TRUE)[,3]) %>% 
  mutate(Direction = str_split(variable, "_", simplify=TRUE)[,2]) %>% 
  ggplot() + 
  aes(
    x = Gene, 
    y = value, 
    fill = interaction(Disease, Direction)
  ) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    values = c(
      "ALS.Down" = "#97cef7", 
      "ALS.Up" = "#f7ba97",
      "ALSFTD.Down" = "#598cb3", 
      "ALSFTD.Up" = "#db621d"
    )
  ) +
  geom_hline(yintercept = 0, color = "black") +
  scale_y_continuous(limits=c(-5,5), breaks=c(-2.5, 0, 2.5)) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title = element_blank(), 
    legend.title = element_blank(), 
    axis.line = element_blank(), 
    panel.grid.major.y = element_line(color = "grey80"),  # Major y-axis grid lines
    panel.grid.minor.y = element_blank(), 
    panel.grid.major.x = element_blank(),  # Remove x-axis major grid lines (optional)
    panel.grid.minor.x = element_blank()   # Remove x-axis minor grid lines (optional)
  ) 

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig8/", 
    "Genes_Barplot_ALSFTD_Only", 
    ".pdf"
  ), 
  width = 9, 
  height = 9,
  units = "in"
)

df %>% 
  filter(Signif_In=="ALSFTD") %>% 
  dplyr::select(
    Gene, Sign_In_N_ALSFTD
  ) %>% 
ggplot() + 
  aes(
    x = Gene, 
    y = 1, 
    fill = Sign_In_N_ALSFTD
  ) + 
  geom_col() + 
  scale_fill_continuous(
  ) +
  geom_hline(yintercept = 0, color = "black") +
  scale_y_continuous(limits=c(-5,5), breaks=c(-2.5, 0, 2.5)) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title = element_blank(), 
    legend.title = element_blank(), 
    axis.line = element_blank(), 
    panel.grid.major.y = element_line(color = "grey80"),  # Major y-axis grid lines
    panel.grid.minor.y = element_blank(), 
    panel.grid.major.x = element_blank(),  # Remove x-axis major grid lines (optional)
    panel.grid.minor.x = element_blank()   # Remove x-axis minor grid lines (optional)
  ) 

library(clusterProfiler)
library(org.Hs.eg.db)
GEX_Features <- qread("../Data/Annotations/Genes/GEX_Features_ENSEMBL_ENTREZ.qrds")


Universe_AllCase_ALS_ALSFTD_List <- qread("../Data/Annotations/Genes/Universe_AllCase_ALS_ALSFTD_List.qrds")
universe <- unlist(Universe_AllCase_ALS_ALSFTD_List[2:13]) %>% 
  unique() 

universe <- GEX_Features$ENTREZID[match(universe, GEX_Features$ID_10X)]
universe <- universe[!is.na(universe)]

df %>% 
  #filter(Signif_In=="Both" ) %>% 
  filter(Gene %in% ATAC_Genes[which(! ATAC_Genes %in% FANS$gene[FANS$fdr<0.05])]) %>% 
  dplyr::select(Gene, Gene_Up_All, Gene_Down_All) %>% 
  mutate(valueTog = Gene_Up_All + Gene_Down_All) -> df2 
table(df2$valueTog==0)
table(df2$valueTog<0)


df2 <- df2[order(df2$valueTog, decreasing = TRUE),]

genelist <- df2$valueTog
names(genelist) <- df2$Gene 
names(genelist) <- GEX_Features$ENTREZID[match(names(genelist), GEX_Features$ID_10X)]
genelist <- genelist[!is.na(names(genelist))]
genelist2 = names(genelist)

Enrich_KEGG_ATAC_Up_And_Down_0.05 <- enrichKEGG(
  gene = genelist2, 
  organism = "hsa", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05, 
  universe = universe
)

genelist <- df2$valueTog
names(genelist) <- df2$Gene 
names(genelist) <- GEX_Features$ENTREZID[match(names(genelist), GEX_Features$ID_10X)]
genelist <- genelist[!is.na(names(genelist))]


GSEA_KEGG_ATAC_Up_And_Down_0.05 <- gseKEGG(
  gene = genelist, 
  pvalueCutoff = 0.05
)

genelist <- df2$valueTog
names(genelist) <- df2$Gene 
genelist <- genelist[genelist<=0]
names(genelist) <- GEX_Features$ENTREZID[match(names(genelist), GEX_Features$ID_10X)]
genelist <- genelist[!is.na(names(genelist))]
genelist2 = names(genelist)

Enrich_KEGG_ATAC_Down_0.05 <- enrichKEGG(
  gene = genelist2, 
  organism = "hsa", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05, 
  universe = universe
)

genelist <- df2$valueTog
names(genelist) <- df2$Gene 
genelist <- genelist[genelist>0]
names(genelist) <- GEX_Features$ENTREZID[match(names(genelist), GEX_Features$ID_10X)]
genelist <- genelist[!is.na(names(genelist))]
genelist2 = names(genelist)

Enrich_KEGG_ATAC_Up_0.05 <- enrichKEGG(
  gene = genelist2, 
  organism = "hsa", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05, 
  universe = universe
)




genelist3 = FANS$gene[FANS$fdr < 0.05][which(! FANS$gene[FANS$fdr<0.05] %in% ATAC_Genes)]
genelist3 = genelist3[which(genelist3 %in% df$Gene[df$Signif_In=="Both"])]
genelist3 = GEX_Features$ENTREZID[match(genelist3, GEX_Features$ID_10X)]
genelist3 <- genelist3[!is.na(genelist3)]
Enrich_KEGG_FANS <-enrichKEGG(
  gene = genelist3, 
  organism = "hsa", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.01, 
  universe = universe
)

## For Figure 
df <- data.frame(
  Gene = names.tmp, 
  Signif_In = ALS_ALSFTD_RNA_WNN_L25_Signature$Signif_In[match(names.tmp, ALS_ALSFTD_RNA_WNN_L25_Signature$ID_10X)], 
  Sign_In_N_All = DEG_In_N, 
  Sign_In_N_ALS = DEG_In_N_ALS, 
  Sign_In_N_ALSFTD = DEG_In_N_ALSFTD,
  Gene_Up_All = Gene_Up, 
  Gene_Down_All = Gene_Down, 
  Gene_Up_ALS = Gene_Up_ALS, 
  Gene_Down_ALS = Gene_Down_ALS, 
  Gene_Up_ALSFTD = Gene_Up_ALSFTD, 
  Gene_Down_ALSFTD = Gene_Down_ALSFTD
)

Gene_order <- df$Gene[order(df$Gene_Up_All + df$Gene_Down_All)]
df <- df[order(df$Gene_Up_All + df$Gene_Down_All),]
df$Gene <- factor(
  df$Gene, 
  levels= Gene_order
)

table(
  round(df$Gene_Up_All,1)==round(df$Gene_Up_ALS + df$Gene_Up_ALSFTD,1)
) 

table(
  round(df$Gene_Down_All,1)==round(df$Gene_Down_ALS + df$Gene_Down_ALSFTD,1)
)

df %>% 
  filter(Signif_In=="ALS") %>% 
  dplyr::select(
    Gene, 
    Gene_Up_ALS, 
    Gene_Down_ALS
  ) %>% 
  melt() %>% 
  mutate(Disease = str_split(variable, "_",simplify=TRUE)[,3]) %>% 
  mutate(Direction = str_split(variable, "_", simplify=TRUE)[,2]) %>% 
  ggplot() + 
  aes(
    x = Gene, 
    y = value, 
    fill = interaction(Disease, Direction)
  ) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    values = c(
      "ALS.Down" = "#97cef7", 
      "ALS.Up" = "#f7ba97",
      "ALSFTD.Down" = "#598cb3", 
      "ALSFTD.Up" = "#db621d"
    )
  ) +
  geom_hline(yintercept = 0, color = "black") +
  scale_y_continuous(limits=c(-5,5), breaks=c(-2.5, 0, 2.5)) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title = element_blank(), 
    legend.title = element_blank(), 
    axis.line = element_blank(), 
    panel.grid.major.y = element_line(color = "grey80"),  # Major y-axis grid lines
    panel.grid.minor.y = element_blank(), 
    panel.grid.major.x = element_blank(),  # Remove x-axis major grid lines (optional)
    panel.grid.minor.x = element_blank()   # Remove x-axis minor grid lines (optional)
  ) 

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig8/", 
    "Genes_Barplot_ALS", 
    ".pdf"
  ), 
  width = 9, 
  height = 9,
  units = "in"
)


df %>% 
  filter(Signif_In=="ALSFTD") %>% 
  dplyr::select(
    Gene, 
    Gene_Up_ALSFTD, 
    Gene_Down_ALSFTD
  ) %>% 
  melt() %>% 
  mutate(Disease = str_split(variable, "_",simplify=TRUE)[,3]) %>% 
  mutate(Direction = str_split(variable, "_", simplify=TRUE)[,2]) %>% 
  ggplot() + 
  aes(
    x = Gene, 
    y = value, 
    fill = interaction(Disease, Direction)
  ) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    values = c(
      "ALS.Down" = "#97cef7", 
      "ALS.Up" = "#f7ba97",
      "ALSFTD.Down" = "#598cb3", 
      "ALSFTD.Up" = "#db621d"
    )
  ) +
  geom_hline(yintercept = 0, color = "black") +
  scale_y_continuous(limits=c(-5,5), breaks=c(-2.5, 0, 2.5)) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title = element_blank(), 
    legend.title = element_blank(), 
    axis.line = element_blank(), 
    panel.grid.major.y = element_line(color = "grey80"),  # Major y-axis grid lines
    panel.grid.minor.y = element_blank(), 
    panel.grid.major.x = element_blank(),  # Remove x-axis major grid lines (optional)
    panel.grid.minor.x = element_blank()   # Remove x-axis minor grid lines (optional)
  ) 

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig8/", 
    "Genes_Barplot_ALSFTD", 
    ".pdf"
  ), 
  width = 9, 
  height = 9,
  units = "in"
)


df %>% 
  filter(Signif_In=="Both") %>% 
  dplyr::select(
    Gene, 
    Gene_Up_ALS, 
    Gene_Down_ALS,
    Gene_Up_ALSFTD, 
    Gene_Down_ALSFTD
  ) %>% 
  melt() %>% 
  mutate(Disease = str_split(variable, "_",simplify=TRUE)[,3]) %>% 
  mutate(Direction = str_split(variable, "_", simplify=TRUE)[,2]) %>% 
  ggplot() + 
  aes(
    x = Gene, 
    y = value, 
    fill = interaction(Disease, Direction)
  ) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    values = c(
      "ALS.Down" = "#97cef7", 
      "ALS.Up" = "#f7ba97",
      "ALSFTD.Down" = "#598cb3", 
      "ALSFTD.Up" = "#db621d"
    )
  ) +
  geom_hline(yintercept = 0, color = "black") +
  scale_y_continuous(limits=c(-5,5), breaks=c(-2.5, 0, 2.5)) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title = element_blank(), 
    legend.title = element_blank(), 
    axis.line = element_blank(), 
    panel.grid.major.y = element_line(color = "grey80"),  # Major y-axis grid lines
    panel.grid.minor.y = element_blank(), 
    panel.grid.major.x = element_blank(),  # Remove x-axis major grid lines (optional)
    panel.grid.minor.x = element_blank()   # Remove x-axis minor grid lines (optional)
  ) 

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig8/", 
    "Genes_Barplot_ALS_ALSFTD", 
    ".pdf"
  ), 
  width = 9, 
  height = 9,
  units = "in"
)

df <- df %>% 
  mutate(In_FANS_TDP43 = Gene %in% FANS$gene[FANS$fdr<0.05]) 

df <- df %>% 
  mutate(In_ATAC = Gene %in% ATAC_Genes) 

df %>% 
  dplyr::select(Gene, Signif_In, Sign_In_N_All) %>% 
  pivot_longer(cols=starts_with("Sign_In"), names_to = "Variable", values_to="Value" ) %>% 
  filter(Signif_In == "Both") %>% 
  ggplot() + 
  aes(x = Gene, y = Variable, fill = factor(Value, levels=1:12)) +
  geom_tile() +                          
  geom_hline(yintercept = 0, color = "black") +
  scale_fill_viridis_d(option="magma", direction = -1, n=12) + 
  theme_classic() +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title = element_blank(), 
    legend.title = element_blank(), 
    axis.line = element_blank(), 
    panel.grid.major.y = element_line(color = "grey80"), 
    panel.grid.minor.y = element_blank(), 
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank()   
  ) 

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig8/", 
    "Genes_Signif_In_ALS_ALSFTD", 
    ".pdf"
  ), 
  width = 9, 
  height = 9,
  units = "in"
)

df %>% 
  dplyr::select(Gene, Signif_In, In_FANS_TDP43) %>% 
  filter(Signif_In == "Both") %>% 
  ggplot() + 
  aes(x = Gene, y = 1, fill = factor(In_FANS_TDP43)) +
  geom_tile() +                         
  scale_fill_manual(values=c("TRUE"="#00135e", "FALSE"="#FFFFFF")) + 
  theme_classic() +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title = element_blank(), 
    legend.title = element_blank(), 
    axis.line = element_blank(), 
    panel.grid.major.y = element_line(color = "grey80"), 
    panel.grid.minor.y = element_blank(), 
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank()   
  ) 

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig8/", 
    "Genes_FANS_TDP43_ALS_ALSFTD", 
    ".pdf"
  ), 
  width = 9, 
  height = 9,
  units = "in"
)

df %>% 
  dplyr::select(Gene, Signif_In, In_FANS_TDP43) %>% 
  filter(Signif_In == "Both") %>% 
  ggplot() + 
  aes(x = Gene, y = 1, fill = factor(In_FANS_TDP43)) +
  geom_tile() +                         
  scale_fill_manual(values=c("TRUE"="#00778c", "FALSE"="#FFFFFF")) + 
  theme_classic() +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title = element_blank(), 
    legend.title = element_blank(), 
    axis.line = element_blank(), 
    panel.grid.major.y = element_line(color = "grey80"), 
    panel.grid.minor.y = element_blank(), 
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank()   
  ) 

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig8/", 
    "Genes_FANS_TDP43_ALS_ALSFTD", 
    ".pdf"
  ), 
  width = 9, 
  height = 9,
  units = "in"
)


df %>% 
  dplyr::select(Gene, Signif_In, In_ATAC) %>% 
  filter(Signif_In == "Both") %>% 
  ggplot() + 
  aes(x = Gene, y = 1, fill = factor(In_ATAC)) +
  geom_tile() +                         
  scale_fill_manual(values=c("TRUE"="#491273", "FALSE"="#FFFFFF")) + 
  theme_classic() +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title = element_blank(), 
    legend.title = element_blank(), 
    axis.line = element_blank(), 
    panel.grid.major.y = element_line(color = "grey80"), 
    panel.grid.minor.y = element_blank(), 
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank()   
  ) 

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig8/", 
    "Genes_ATAC_ALS_ALSFTD", 
    ".pdf"
  ), 
  width = 9, 
  height = 9,
  units = "in"
)

df %>% 
  #filter(Signif_In=="Both") %>% 
  filter(In_FANS_TDP43==FALSE & In_ATAC==TRUE) %>% 
  dplyr::select(
    Gene, 
    Gene_Up_ALS, 
    Gene_Down_ALS,
    Gene_Up_ALSFTD, 
    Gene_Down_ALSFTD
  ) %>% 
  melt() %>% 
  mutate(Disease = str_split(variable, "_",simplify=TRUE)[,3]) %>% 
  mutate(Direction = str_split(variable, "_", simplify=TRUE)[,2]) %>% 
  ggplot() + 
  aes(
    x = Gene, 
    y = value, 
    fill = interaction(Disease, Direction)
  ) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    values = c(
      "ALS.Down" = "#97cef7", 
      "ALS.Up" = "#f7ba97",
      "ALSFTD.Down" = "#598cb3", 
      "ALSFTD.Up" = "#db621d"
    )
  ) +
  geom_hline(yintercept = 0, color = "black") +
  scale_y_continuous(limits=c(-5,5), breaks=c(-2.5, 0, 2.5)) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title = element_blank(), 
    legend.title = element_blank(), 
    axis.line = element_blank(), 
    panel.grid.major.y = element_line(color = "grey80"),  # Major y-axis grid lines
    panel.grid.minor.y = element_blank(), 
    panel.grid.major.x = element_blank(),  # Remove x-axis major grid lines (optional)
    panel.grid.minor.x = element_blank()   # Remove x-axis minor grid lines (optional)
  ) 

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig8/", 
    "Genes_Barplot_ATAC_Only_ALS_ALSFTD", 
    ".pdf"
  ), 
  width = 9, 
  height = 9,
  units = "in"
)

df %>% 
  dplyr::select(Gene, Signif_In, In_FANS_TDP43, In_ATAC) %>% 
  #filter(Signif_In=="Both") %>% 
  filter(In_FANS_TDP43==FALSE & In_ATAC==TRUE) %>% 
  ggplot() + 
  aes(x = Gene, y = 1, fill = factor(In_FANS_TDP43)) +
  geom_tile() +                         
  scale_fill_manual(values=c("TRUE"="#00778c", "FALSE"="#FFFFFF")) + 
  theme_classic() +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title = element_blank(), 
    legend.title = element_blank(), 
    axis.line = element_blank(), 
    panel.grid.major.y = element_line(color = "grey80"), 
    panel.grid.minor.y = element_blank(), 
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank()   
  ) 

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig8/", 
    "Genes_ATAC_TDP43_ALS_ALSFTD", 
    ".pdf"
  ), 
  width = 9, 
  height = 9,
  units = "in"
)


df %>% 
  dplyr::select(Gene, Signif_In, In_FANS_TDP43, In_ATAC) %>% 
  #filter(Signif_In=="Both") %>% 
  filter(In_FANS_TDP43==FALSE & In_ATAC==TRUE) %>% 
  ggplot() + 
  aes(x = Gene, y = 1, fill = factor(In_ATAC)) +
  geom_tile() +                         
  scale_fill_manual(values=c("TRUE"="#491273", "FALSE"="#FFFFFF")) + 
  theme_classic() +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title = element_blank(), 
    legend.title = element_blank(), 
    axis.line = element_blank(), 
    panel.grid.major.y = element_line(color = "grey80"), 
    panel.grid.minor.y = element_blank(), 
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank()   
  ) 

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig8/", 
    "Genes_ATAC_ATAC_ALS_ALSFTD", 
    ".pdf"
  ), 
  width = 9, 
  height = 9,
  units = "in"
)


df$FANS_TDP43_ATAC <- "None"
df$FANS_TDP43_ATAC[df$In_FANS_TDP43 & !df$In_ATAC] <- "FANS"
df$FANS_TDP43_ATAC[!df$In_FANS_TDP43 & df$In_ATAC] <- "ATAC"
df$FANS_TDP43_ATAC[df$In_FANS_TDP43 & df$In_ATAC] <- "FANS_ATAC"
df$FANS_TDP43_ATAC <- factor(df$FANS_TDP43_ATAC)
table(df$FANS_TDP43_ATAC)

df %>% 
  dplyr::select(Gene, Signif_In, FANS_TDP43_ATAC) %>% 
  filter(Signif_In=="Both") %>% 
  ggplot() + 
  aes(x = Gene, y = 1, fill=FANS_TDP43_ATAC) +
  geom_tile() +                         
  scale_fill_manual(
    values=c(
      "ATAC"="#6FB000", 
      "FANS"="#B32100", 
      "FANS_ATAC"="#FAD900", 
      "None"="#FFFFFF"
    )
  ) + 
  theme_classic() +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title = element_blank(), 
    legend.title = element_blank(), 
    axis.line = element_blank(), 
    panel.grid.major.y = element_line(color = "grey80"), 
    panel.grid.minor.y = element_blank(), 
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank()   
  ) 

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig8/", 
    "Genes_Heatmap_FANS_TDP43_ATAC_RGY", 
    ".pdf"
  ), 
  width = 9, 
  height = 9,
  units = "in"
)


#### 
library(igraph)
library(tidygraph)
library(dendextend)
library(igraph)
library(ggraph)
library(circlize)

plot_Enrich_as_Tree <- function(GSEA_Results, first=NULL, save=TRUE, name, width, height){
  
  geneLists <- lapply(GSEA_Results$geneID, function(x) unlist(strsplit(as.character(x), "/")))
  if(!is.null(first)){
    geneLists <- geneLists[1:first]
  }
  names(geneLists) <- GSEA_Results$Description[1:first]
  geneLists <- geneLists[!is.na(geneLists)]
  
  n <- length(geneLists)
  w <- matrix(NA, nrow = n, ncol = n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      u <- unlist(geneLists[i])
      v <- unlist(geneLists[j])
      w[i, j] <- length(intersect(u, v)) / length(unique(c(u, v)))
    }
  }
  
  
  Terms <- names(geneLists)
  
  rownames(w) <- Terms
  colnames(w) <- Terms
  
  
  hclust <- as.dist(1 - w) %>% hclust(method = "average")
  dend <- as.dendrogram(hclust) 
  
  
  # Get the vertex IDs of the leaves (bottom level nodes)
  leaf_ids <- which(V(graph_dend)$leaf)
  
  # Assign cluster labels to the leaf nodes
  V(graph_dend)$cluster[leaf_ids] <- clusters  # Add the cluster variable to the leaf nodes
  
  library(tidygraph)
  tbl_graph(dend)
  a <- ggraph(dend, 'dendrogram', height = height) + 
    geom_edge_elbow() + 
    geom_node_point(aes(filter = leaf, color = label))  # Color only the leaves
  
  a$data$NES <- GSEA_Kegg_ByL2FC@result$NES[match(a$data$label, GSEA_Kegg_ByL2FC@result$Description)] 
  
  count_GSEA_CoreEnrichment_Genes <- function(x){
    if(is.na(x)) {return(NA)} else{
      return(sum(str_split(x, "/", simplify=TRUE)[1,] != ""))
    }
  }
  a$data$nGenes <- sapply(
    GSEA_Kegg_ByL2FC@result$core_enrichment[match(a$data$label, GSEA_Kegg_ByL2FC@result$Description)], 
    FUN = count_GSEA_CoreEnrichment_Genes
  )
  colFun=colorRamp2(
    c(
      -max(abs(a$data$NES), na.rm=T), 
      0, 
      max(abs(a$data$NES), na.rm=T)
    ), 
    c(
      "blue", 
      "#EEEEEE", 
      "red"
    )
  )
  
  library(ggrepel)
  ggraph(a$data, 'dendrogram', height = height) + 
    geom_edge_elbow(strength = 1, flipped = FALSE, width=1) + 
    geom_node_point(aes(filter = leaf, fill = NES, size = nGenes), pch=21, col="#000000") + 
    geom_node_text(aes(filter = leaf, label = label), nudge_y = 0.1, hjust = 0) + 
    scale_size(
      limits=c(
        min(a$data$nGenes, na.rm=TRUE)-0.25*(max(a$data$nGenes, na.rm=TRUE)-min(a$data$nGenes, na.rm=TRUE)), 
        max(a$data$nGenes, na.rm=TRUE)
      )
    ) + 
    scale_fill_gradientn(
      colours = c(
        "#0000FF", 
        "#05f7ff", 
        "white", 
        "#f5dc00",
        "#FF0000"
      ), 
      limits=c(
        -max(abs(a$data$NES), na.rm=TRUE), 
        max(abs(a$data$NES), na.rm=TRUE)
      )
      
    ) + coord_flip(xlim = c(0, 8), ylim=c(1, -3)) + scale_y_reverse() + 
    theme(
      panel.background = element_blank(), 
      
    ) 
  
  
  
  GSEA_Results_Direction <- GSEA_Results$NES<0
  
  mycolors <- sort(rainbow(20))[c(1, 20, 10, 11, 2, 19, 3, 12, 4, 13, 5, 14, 6, 15, 7, 16, 8, 17, 9, 18)]
  
  leafType <- as.factor(gsub(" .*", "", GSEA_Results_Direction[ix]))
  
  if (max(nchar(GSEA_Results[ix])) >= 1) { # if "Up regulated or Downregulated"; not "A", "B"
    # leafColors = c("green","red")  else  # mycolors # k-Means
    leafColors <- mycolors[1:2]
  } else { # convert c("B","D","E") to c(2, 4, 5)
    # leafType= as.factor( gsub(" .*","", enrichedTerms$Direction[ix] ) )
    leafType <- match(gsub(" .*", "", GSEA_Results$Direction[ix]), toupper(letters))
    
    
    
  }
  leafColors <- c("#FF0000", "#0000FF")
  
  leafSize <- as.numeric(abs(GSEA_Results$NES[ix])) # leaf size represent NES values
  leafSize <- .9 * (leafSize - min(leafSize)) / (max(leafSize) - min(leafSize) + 1e-50) + .1 # scale more aggressively
  
  
  print(
    dend %>%
      as.dendrogram(hang = -1) %>%
      set("leaves_pch", 19) %>% # type of marker
      set("leaves_cex", leafSize*3) %>% # Size
      set("leaves_col", leafColors[leafType]) %>% # up or down genes
      plot(horiz = TRUE, axes = FALSE)
  ) 
  
  if(save){
    pdf(
      name, 
      width = width, 
      height = height
    ) 
    dend %>%
      as.dendrogram(hang = -1) %>%
      set("leaves_pch", 19) %>% # type of marker
      set("leaves_cex", leafSize*3) %>% # Size
      set("leaves_col", leafColors[leafType]) %>% # up or down genes #leafColors
      plot(horiz = TRUE, axes = FALSE)
    dev.off()
  }
  
}

Return_GSEA_as_Tree <- function(GSEA_Results, name, save=TRUE, width, height){
  geneLists <- lapply(GSEA_Results$core_enrichment, function(x) unlist(strsplit(as.character(x), "/")))
  names(geneLists) <- GSEA_Results$Description
  geneLists <- geneLists[!is.na(geneLists)]
  
  n <- length(geneLists)
  w <- matrix(NA, nrow = n, ncol = n)
  
  for (i in 1:n) {
    for (j in i:n) {
      u <- unlist(geneLists[i])
      v <- unlist(geneLists[j])
      w[i, j] <- length(intersect(u, v)) / length(unique(c(u, v)))
    }
  }
  
  for (i in 1:n) {
    for (j in 1:(i - 1)) {
      w[i, j] <- w[j, i]
    }
  }
  
  if (0) {
    total_elements <- 30000
    n <- length(geneLists)
    w <- matrix(rep(0, n * n), nrow = n, ncol = n)
    
    for (i in 1:n) {
      for (j in (i + 1):n) {
        u <- unlist(geneLists[i])
        v <- unlist(geneLists[j])
        xx <- length(intersect(u, v))
        if (xx == 0) {
          next
        }
        mm <- length(u)
        nn <- total_elements - mm
        kk <- length(v)
        w[i, j] <- -sqrt(-phyper(xx - 1, mm, nn, kk, lower.tail = FALSE, log.p = TRUE))
      }
    }
    
    
    for (i in 1:n) {
      for (j in 1:(i - 1)) {
        w[i, j] <- w[j, i]
      }
    }
  }
  
  Terms <- paste(
    names(geneLists)
  )
  rownames(w) <- Terms
  colnames(w) <- Terms
  rightMargin=33
  par(mar = c(0, 0, 1, rightMargin)) 
  
  dend <- as.dist(1 - w) %>%
    hclust(method = "average")
  ix <- dend$order # permutated order of leaves
  
  GSEA_Results_Direction <- GSEA_Results$NES<0
  
  mycolors <- sort(rainbow(20))[c(1, 20, 10, 11, 2, 19, 3, 12, 4, 13, 5, 14, 6, 15, 7, 16, 8, 17, 9, 18)]
  
  leafType <- as.factor(gsub(" .*", "", GSEA_Results_Direction[ix]))
  
  if (max(nchar(GSEA_Results[ix])) >= 1) { # if "Up regulated or Downregulated"; not "A", "B"
    # leafColors = c("green","red")  else  # mycolors # k-Means
    leafColors <- mycolors[1:2]
  } else { # convert c("B","D","E") to c(2, 4, 5)
    # leafType= as.factor( gsub(" .*","", enrichedTerms$Direction[ix] ) )
    leafType <- match(gsub(" .*", "", GSEA_Results$Direction[ix]), toupper(letters))
    
    
    
  }
  leafColors <- c("#FF0000", "#0000FF")
  
  leafSize <- as.numeric(abs(GSEA_Results$NES[ix])) # leaf size represent NES values
  leafSize <- .9 * (leafSize - min(leafSize)) / (max(leafSize) - min(leafSize) + 1e-50) + .1 # scale more aggressively
  
  
  dend2 <- dend %>%
    as.dendrogram(hang = -1) %>%
    set("leaves_pch", 19) %>% # type of marker
    set("leaves_cex", leafSize*3) %>% # Size
    set("leaves_col", leafColors[leafType]) %>% # up or down genes
    
    return(dend2)
}

####
Enrich_KEGG_ATAC_Up_And_Down_0.05_First10 <- Enrich_KEGG_ATAC_Up_And_Down_0.05@result[1:10,]
GSEA_Results = Enrich_KEGG_ATAC_Up_And_Down_0.05_First10
first = NULL
  
geneLists <- lapply(GSEA_Results$geneID, function(x) unlist(strsplit(as.character(x), "/")))
if(!is.null(first)){
  geneLists <- geneLists[1:first]
}
names(geneLists) <- GSEA_Results$Description
geneLists <- geneLists[!is.na(geneLists)]
  
  n <- length(geneLists)
  w <- matrix(NA, nrow = n, ncol = n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      u <- unlist(geneLists[i])
      v <- unlist(geneLists[j])
      w[i, j] <- length(intersect(u, v)) / length(unique(c(u, v)))
    }
  }
  
  
  Terms <- names(geneLists)
  
  rownames(w) <- Terms
  colnames(w) <- Terms
  
  
  hclust <- as.dist(1 - w) %>% hclust(method = "average")
  dend <- as.dendrogram(hclust) 

  
  a <- ggraph(dend, 'dendrogram', height = height) + 
  geom_edge_elbow() + 
  geom_node_point(aes(filter = leaf, color = label))  # Color only the leaves

a$data$FoldEnrichment <- Enrich_KEGG_ATAC_Up_And_Down_0.05_First10$FoldEnrichment[match(a$data$label, Enrich_KEGG_ATAC_Up_And_Down_0.05_First10$Description)] 

count_GSEA_CoreEnrichment_Genes <- function(x){
  if(is.na(x)) {return(NA)} else{
    return(sum(str_split(x, "/", simplify=TRUE)[1,] != ""))
  }
}
a$data$nGenes <- sapply(
  Enrich_KEGG_ATAC_Up_And_Down_0.05_First10$geneID[match(a$data$label, Enrich_KEGG_ATAC_Up_And_Down_0.05_First10$Description)], 
  FUN = count_GSEA_CoreEnrichment_Genes
)
colFun=colorRamp2(
  c(
    -max(abs(a$data$FoldEnrichment), na.rm=T), 
    0, 
    max(abs(a$data$FoldEnrichment), na.rm=T)
  ), 
  c(
    "blue", 
    "#EEEEEE", 
    "red"
  )
)

library(ggrepel)
ggraph(a$data, 'dendrogram', height = height) + 
  geom_edge_elbow(strength = 1, flipped = FALSE, width=1) + 
  geom_node_point(aes(filter = leaf, fill = FoldEnrichment, size = nGenes), pch=21, col="#000000") + 
  geom_node_text(aes(filter = leaf, label = label), nudge_y = 0.1, hjust = 0) + 
  scale_size(
    limits=c(
      min(a$data$nGenes, na.rm=TRUE)-0.25*(max(a$data$nGenes, na.rm=TRUE)-min(a$data$nGenes, na.rm=TRUE)), 
      max(a$data$nGenes, na.rm=TRUE)
    )
  ) + 
  scale_fill_gradientn(
    colours = c(
      "white", 
      "#f5dc00",
      "#FF0000"
    ), 
    limits=c(
      0, 
      max(abs(a$data$FoldEnrichment), na.rm=TRUE)
    )
    
  ) + coord_flip(xlim = c(0, 10), ylim=c(1, -3)) + scale_y_reverse() + 
  theme(
    panel.background = element_blank(), 
    
  ) 

ATAC_genes_from_top10_KEGGs <- unlist(geneLists)
sort(table(ATAC_genes_from_top10_KEGGs))
table(table(ATAC_genes_from_top10_KEGGs)>1)
ATAC_genes_from_top10_KEGGs <- names(table(ATAC_genes_from_top10_KEGGs)[table(ATAC_genes_from_top10_KEGGs)>1])


ATAC_genes_from_top10_KEGGs <- GEX_Features$SYMBOL[match(ATAC_genes_from_top10_KEGGs, GEX_Features$ENTREZID)]
ATAC_genes_from_top10_KEGGs <- sort(ATAC_genes_from_top10_KEGGs)


df %>% 
  filter(Gene %in% ATAC_genes_from_top10_KEGGs) %>% 
  dplyr::select(
    Gene, 
    Gene_Up_ALS, 
    Gene_Down_ALS,
    Gene_Up_ALSFTD, 
    Gene_Down_ALSFTD
  ) %>% 
  melt() %>% 
  mutate(Disease = str_split(variable, "_",simplify=TRUE)[,3]) %>% 
  mutate(Direction = str_split(variable, "_", simplify=TRUE)[,2]) %>% 
  ggplot() + 
  aes(
    x = Gene, 
    y = value, 
    fill = interaction(Disease, Direction)
  ) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    values = c(
      "ALS.Down" = "#97cef7", 
      "ALS.Up" = "#f7ba97",
      "ALSFTD.Down" = "#598cb3", 
      "ALSFTD.Up" = "#db621d"
    )
  ) +
  geom_hline(yintercept = 0, color = "black") +
  scale_y_continuous(limits=c(-5,5), breaks=c(-2.5, 0, 2.5)) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title = element_blank(), 
    legend.title = element_blank(), 
    axis.line = element_blank(), 
    panel.grid.major.y = element_line(color = "grey80"),  # Major y-axis grid lines
    panel.grid.minor.y = element_blank(), 
    panel.grid.major.x = element_blank(),  # Remove x-axis major grid lines (optional)
    panel.grid.minor.x = element_blank()   # Remove x-axis minor grid lines (optional)
  ) + 
  coord_flip()

ggsave(
  paste0(
    "../Data/Visualization/Figures/Fig8/", 
    "Genes_Barplot_ATAC_Only_Top10_KEGGs", 
    ".pdf"
  ), 
  width = 9, 
  height = 9,
  units = "in"
)

write.csv(Enrich_KEGG_ATAC_Up_And_Down_0.05@result, "../Data/Visualization/Tables/Suppl.Table_Enrich_KEGG_ATAC_ALS_ALSFTD_0.05.csv", row.names = TRUE, col.names = TRUE, quote = FALSE)
geneLists3 <- lapply(Enrich_KEGG_ATAC_Up_And_Down_0.05@result$geneID, function(x) unlist(strsplit(as.character(x), "/")))
AllGenes <- unique(unlist(geneLists3))
AllGenes <- GEX_Features$SYMBOL[match(AllGenes, GEX_Features$ENTREZID)]
AllGenes <- AllGenes[!is.na(AllGenes)]
AllGenes <- unique(AllGenes)
AllGenes <- sort(AllGenes)
AllGenes
write.csv(df2$Gene, "../Data/Visualization/Tables/ALS_ALSFTD_ATAC_Only_Genes.csv", row.names=FALSE, col.names=TRUE, quote=FALSE)
