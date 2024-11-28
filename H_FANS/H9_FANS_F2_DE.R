

  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 

library(qs) 
library(Seurat)
library(tidyverse) 
library(SeuratObject) 
library(MAST)
library(zinbwave)
library(patchwork)




  ### 1.0 Load data ------------------------------------------------------------

F2 <- qread(
  paste0(
    "../Data/SeuratObjects/", 
    "F2", 
    ".qrds"
  ), 
  nthr=nthr
)




  ### 2.0 DGE Analysis with Seurat's FindMarkers ------------------------------- 

FANS_DE_TDP43_Wilcox <- Seurat::FindMarkers(
  F2, 
  ident.1 = "TDP43_Low", 
  ident.2 = "TDP43_High", 
  test.use="MAST", 
  latent.vars = "Sex"
  )

FANS_DE_TDP43_Wilcox$Gene <- rownames(FANS_DE_TDP43_Wilcox)




  ### 3.0 DGE Analysis with MAST ----------------------------------------------- 

options(mc.cores=46)

counts <- GetAssayData(F2[["RNA"]], layer="counts")
counts_norm <- counts
counts_norm@x <- log2((10^6*counts@x/rep.int(colSums(counts), diff(counts@p)))+1)
  
data_cell_meta <- F2@meta.data
data_feature_meta <- tibble(gene_name=row.names(counts_norm))
  
scaRaw <- FromMatrix(
  exprsArray = Matrix::as.matrix(counts_norm), 
  cData = data_cell_meta, 
  fData = data_feature_meta
)
  
sca <- scaRaw[which(freq(scaRaw)>0),]
  
cdr <- colSums(assay(sca)>0)
colData(sca)$nFeature_SCT <- cdr
  
colData(sca)$cn_genes_on <- scale(colData(sca)$nFeature_SCT)
colData(sca)$cn_nCount_SCT <- scale(colData(sca)$nCount_RNA)
colData(sca)$cn_percent_mt <- scale(colData(sca)$PctMt)
expressed_genes <- freq(sca) > 0.05 

sca <- sca[expressed_genes,] 
  
cond <- factor(colData(sca)$TDP43)
table(cond)
cond <- relevel(cond, "TDP43_High")
colData(sca)$TDP43 <- cond

zlmCond <- zlm(
  ~TDP43 + cn_genes_on + cn_nCount_SCT + cn_percent_mt + (1 | Sex) + (1|orig.ident), 
  sca, 
  method="glmer", 
  ebayes=FALSE, 
  strictConvergence=FALSE, 
  fitArgsD = list(nAGQ=0)
)
  
  
zlmCond_conv <- zlmCond@converged %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  as_tibble() %>%
  dplyr::rename(conv_C = C, conv_D = D)
  
lrt_term <- paste0("TDP43", "TDP43_Low")
  
summary_cond <- summary(zlmCond, doLRT = lrt_term, fitArgsD = list(nAGQ = 0))
  
summary_Dt <- summary_cond$datatable
  
fcHurdle <- merge(summary_Dt[contrast==lrt_term & component=='H',.(primerid, `Pr(>Chisq)`)], 
                    summary_Dt[contrast==lrt_term & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') 
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  
fcHurdle_df <- fcHurdle %>%
  as_tibble() %>%
    arrange(fdr) %>%
      dplyr::rename(
        gene = primerid,
        p_value = `Pr(>Chisq)`,
        model_log2FC = coef
        )
  
data_cpm <- t(assay(sca)) %>%
  as_tibble() %>%
    mutate(TDP43 = colData(sca)$TDP43, ID = colData(sca)$orig.ident) %>%
      pivot_longer(names_to = "gene", values_to = "log2_CPM", cols = -c(TDP43, ID))

data_fc <- data_cpm %>%
    group_by(gene, TDP43) %>%
      summarise(mean_log2_CPM = mean(log2_CPM)) %>%
        ungroup() %>%
          mutate(
            disease_group = case_when(
              TDP43 == "TDP43_Low" ~ "avg_CPM_disease_1",
              TDP43 == "TDP43_High" ~ "avg_CPM_disease_2"
              )
            ) %>%
            dplyr::select(-"TDP43") %>%
              pivot_wider(names_from = "disease_group", values_from = "mean_log2_CPM") %>%
                mutate(avg_log2FC = avg_CPM_disease_1 - avg_CPM_disease_2)
  
final_res <- fcHurdle_df %>%
  left_join(zlmCond_conv) %>%
    left_join(data_fc)
  

FANS_DE_TDP43_MAST <- data.frame(final_res) 




  ### 4.0 Save data ------------------------------------------------------------

qsave(
  FANS_DE_TDP43_Wilcox, 
  paste0(
    "../Data/DE/FANS/", 
    "FANS_DE_TDP43_Wilcox", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  FANS_DE_TDP43_MAST, 
  paste0(
    "../Data/DE/FANS/", 
    "FANS_DE_TDP43_MAST", 
    ".qrds"
  ), 
  nthr=nthr
)


