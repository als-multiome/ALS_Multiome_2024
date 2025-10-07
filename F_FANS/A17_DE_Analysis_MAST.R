# A17_DE_Analysis_MAST.R 



  ### 0.0 Load libraries ------------------------------------------------------- 

library(qs2)
library(tidyverse)
library(Seurat)
library(zinbwave) 
library(MAST)
library(Matrix)



  
  ### 1.0 Load data ------------------------------------------------------------

FANS <- qs_read(
    "../Data/FANS/SeuratObjects/FANS.Unsubsetted.qs2"
)

FANS.SubsetCells <- qs_read(
  "../Data/FANS/SeuratObjects/FANS.SubsetCells.qs2"
)

FANS.SubsetCounts <- qs_read(
  "../Data/FANS/SeuratObjects/FANS.SubsetCounts.qs2"
)

Samples <- qs_read(
  "../Data/FANS/Samples/Samples.qs2"
)




  ### 2.0 Define MAST Workflow -------------------------------------------------

MAST_DE <- function(
    SeuratObject 
 ){
  
  GroupingVar <- "TDP43"
  StudyCase <- "Low"
  RefCase <- "High"
  
  SeuratObject[["GroupingVar"]] <- SeuratObject[["TDP43"]]
  counts <- GetAssayData(SeuratObject[["RNA"]], layer="counts")
  counts_norm <- counts
  counts_norm@x <- log2((10^6*counts@x/rep.int(colSums(counts), diff(counts@p)))+1)
  
  data_cell_meta <- SeuratObject@meta.data
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
  colData(sca)$cn_percent_mt <- scale(colData(sca)$percent.mt)
  
  expressed_genes <- freq(sca) > 0.05 
  sca <- sca[expressed_genes,] 
  
  cond <- factor(colData(sca)$GroupingVar)
  table(cond)
  cond <- relevel(cond, RefCase)
  colData(sca)$GroupingVar <- cond
  
  zlmCond <- zlm(
    ~GroupingVar + cn_genes_on + cn_nCount_SCT + cn_percent_mt + (1|Sample_donor), 
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
  
  lrt_term <- paste0("GroupingVar", StudyCase)
  
  summary_cond <- summary(zlmCond, doLRT = lrt_term, fitArgsD = list(nAGQ = 0))
  
  summary_Dt <- summary_cond$datatable
  
  fcHurdle <- merge(summary_Dt[contrast==lrt_term & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summary_Dt[contrast==lrt_term & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  
  fcHurdle_df <- fcHurdle %>%
    as_tibble() %>%
    arrange(fdr) %>%
    dplyr::rename(gene = primerid,
                  p_value = `Pr(>Chisq)`,
                  model_log2FC = coef)
  
  # compute avg log2FC by hand
  data_cpm <- t(assay(sca)) %>%
    as_tibble() %>%
    mutate(GroupingVar = colData(sca)$GroupingVar, ID = colData(sca)$ID) %>%
    pivot_longer(names_to = "gene", values_to = "log2_CPM", cols = -c(GroupingVar, ID))
  data_fc <- data_cpm %>%
    group_by(gene, GroupingVar) %>%
    summarise(mean_log2_CPM = mean(log2_CPM)) %>%
    ungroup() %>%
    mutate(disease_group = case_when(GroupingVar == StudyCase ~ "avg_CPM_disease_1",
                                     GroupingVar == RefCase ~ "avg_CPM_disease_2")) %>%
    dplyr::select(-GroupingVar) %>%
    pivot_wider(names_from = "disease_group", values_from = "mean_log2_CPM") %>%
    mutate(avg_log2FC = avg_CPM_disease_1 - avg_CPM_disease_2)
  
  final_res <- fcHurdle_df %>%
    left_join(zlmCond_conv) %>%
    left_join(data_fc)
  
  return(final_res)
  
}



MAST_DE_WNN_L4 <- function(
    SeuratObject 
){
  
  GroupingVar <- "TDP43"
  StudyCase <- "Low"
  RefCase <- "High"
  
  SeuratObject[["GroupingVar"]] <- SeuratObject[["TDP43"]]
  counts <- GetAssayData(SeuratObject[["RNA"]], layer="counts")
  counts_norm <- counts
  counts_norm@x <- log2((10^6*counts@x/rep.int(colSums(counts), diff(counts@p)))+1)
  
  data_cell_meta <- SeuratObject@meta.data
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
  colData(sca)$cn_percent_mt <- scale(colData(sca)$percent.mt)
  
  expressed_genes <- freq(sca) > 0.05 
  sca <- sca[expressed_genes,] 
  
  cond <- factor(colData(sca)$GroupingVar)
  table(cond)
  cond <- relevel(cond, RefCase)
  colData(sca)$GroupingVar <- cond
  
  zlmCond <- zlm(
    ~GroupingVar + cn_genes_on + cn_nCount_SCT + cn_percent_mt + Sex + (1|ID_WNN_L4_Predicted), 
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
  
  lrt_term <- paste0("GroupingVar", StudyCase)
  
  summary_cond <- summary(zlmCond, doLRT = lrt_term, fitArgsD = list(nAGQ = 0))
  
  summary_Dt <- summary_cond$datatable
  
  fcHurdle <- merge(summary_Dt[contrast==lrt_term & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summary_Dt[contrast==lrt_term & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  
  fcHurdle_df <- fcHurdle %>%
    as_tibble() %>%
    arrange(fdr) %>%
    dplyr::rename(gene = primerid,
                  p_value = `Pr(>Chisq)`,
                  model_log2FC = coef)
  
  # compute avg log2FC by hand
  data_cpm <- t(assay(sca)) %>%
    as_tibble() %>%
    mutate(GroupingVar = colData(sca)$GroupingVar, ID = colData(sca)$ID) %>%
    pivot_longer(names_to = "gene", values_to = "log2_CPM", cols = -c(GroupingVar, ID))
  data_fc <- data_cpm %>%
    group_by(gene, GroupingVar) %>%
    summarise(mean_log2_CPM = mean(log2_CPM)) %>%
    ungroup() %>%
    mutate(disease_group = case_when(GroupingVar == StudyCase ~ "avg_CPM_disease_1",
                                     GroupingVar == RefCase ~ "avg_CPM_disease_2")) %>%
    dplyr::select(-GroupingVar) %>%
    pivot_wider(names_from = "disease_group", values_from = "mean_log2_CPM") %>%
    mutate(avg_log2FC = avg_CPM_disease_1 - avg_CPM_disease_2)
  
  final_res <- fcHurdle_df %>%
    left_join(zlmCond_conv) %>%
    left_join(data_fc)
  
  return(final_res)
  
}


  ### 3.0 DE MAST analysis -----------------------------------------------------



    ## 3.1 All Cells -----------------------------------------------------------

FANS_DE_MAST_Unsubsetted <- MAST_DE(
  FANS
)

FANS_DE_MAST_SubsetCells <- MAST_DE(
  FANS.SubsetCells
) 


FANS_DE_MAST_SubsetCounts <- MAST_DE(
  FANS.SubsetCounts
)


table(FANS_DE_MAST_Unsubsetted$fdr < 0.05)
summary(FANS_DE_MAST_Unsubsetted$model_log2FC)
summary(FANS_DE_MAST_Unsubsetted$avg_log2FC)


table(FANS_DE_MAST_SubsetCells$fdr < 0.05)
summary(FANS_DE_MAST_SubsetCells$model_log2FC)
summary(FANS_DE_MAST_SubsetCells$avg_log2FC)


table(FANS_DE_MAST_SubsetCounts$fdr < 0.05)
summary(FANS_DE_MAST_SubsetCounts$model_log2FC)
summary(FANS_DE_MAST_SubsetCounts$avg_log2FC)

plot(FANS_DE_MAST_Unsubsetted$avg_log2FC, FANS_DE_MAST_Unsubsetted$model_log2FC)
abline(0,1,col="red")

plot(FANS_DE_MAST_SubsetCells$avg_log2FC, FANS_DE_MAST_SubsetCells$model_log2FC)
abline(0,1,col="red")

plot(FANS_DE_MAST_SubsetCounts$avg_log2FC, FANS_DE_MAST_SubsetCounts$model_log2FC)
abline(0,1,col="red")



    ## 3.2 Exc_LINC00507 Cells -------------------------------------------------

table(
  FANS$ID_WNN_L25_Predicted, 
  FANS$TDP43
)

table(
  FANS.SubsetCells$ID_WNN_L25_Predicted, 
  FANS.SubsetCells$TDP43
)

table(
  FANS.SubsetCounts$ID_WNN_L25_Predicted, 
  FANS.SubsetCounts$TDP43
)

LINC00507.Unsubsetted <- subset(FANS, subset = ID_WNN_L25_Predicted == "Exc_LINC00507")
LINC00507.SubsetCells <- subset(FANS.SubsetCells, subset = ID_WNN_L25_Predicted == "Exc_LINC00507")
LINC00507.SubsetCounts <- subset(FANS.SubsetCounts, subset = ID_WNN_L25_Predicted == "Exc_LINC00507")

LINC00507_DE_MAST_Unsubsetted <- MAST_DE(
  LINC00507.Unsubsetted
)

LINC00507_DE_MAST_Unsubsetted_WNN_L4Cov <- MAST_DE_WNN_L4(
  LINC00507.Unsubsetted
)


qs_save(
  LINC00507_DE_MAST_Unsubsetted, 
  "../Data/FANS/DE/LINC00507_DE_MAST_Unsubsetted.qs2"
)

qs_save(
  LINC00507_DE_MAST_Unsubsetted_WNN_L4Cov, 
  "../Data/FANS/DE/MAST/WNN_L25/LINC00507_DE_MAST_Unsubsetted_WNN_L4Cov.qs2"
)

table(LINC00507_DE_MAST_Unsubsetted_WNN_L4Cov$fdr<0.05, useNA="always")


LINC00507_DE_MAST_SubsetCells <- MAST_DE(
  LINC00507.SubsetCells
)
qs_save(
  LINC00507_DE_MAST_SubsetCells, 
  "../Data/FANS/DE/LINC00507_DE_MAST_SubsetCells.qs2"
)

LINC00507_DE_MAST_SubsetCounts <- MAST_DE(
  LINC00507.SubsetCounts
)
qs_save(
  LINC00507_DE_MAST_SubsetCounts, 
  "../Data/FANS/DE/LINC00507_DE_MAST_SubsetCounts.qs2"
)



    ## 3.3 Exc_THEMIS Cells ----------------------------------------------------

THEMIS.Unsubsetted <- subset(FANS, subset = ID_WNN_L25_Predicted == "Exc_THEMIS")
THEMIS.SubsetCells <- subset(FANS.SubsetCells, subset = ID_WNN_L25_Predicted == "Exc_THEMIS")
THEMIS.SubsetCounts <- subset(FANS.SubsetCounts, subset = ID_WNN_L25_Predicted == "Exc_THEMIS")

THEMIS_DE_MAST_Unsubsetted <- MAST_DE(
  THEMIS.Unsubsetted
)
qs_save(
  THEMIS_DE_MAST_Unsubsetted, 
  "../Data/FANS/DE/THEMIS_DE_MAST_Unsubsetted.qs2"
)

THEMIS_DE_MAST_SubsetCells <- MAST_DE(
  THEMIS.SubsetCells
)
qs_save(
  THEMIS_DE_MAST_SubsetCells, 
  "../Data/FANS/DE/THEMIS_DE_MAST_SubsetCells.qs2"
)

THEMIS_DE_MAST_SubsetCounts <- MAST_DE(
  THEMIS.SubsetCounts
)
qs_save(
  THEMIS_DE_MAST_SubsetCounts, 
  "../Data/FANS/DE/THEMIS_DE_MAST_SubsetCounts.qs2"
)


THEMIS_DE_MAST_Unsubsetted_WNN_L4Cov <- MAST_DE_WNN_L4(
  THEMIS.Unsubsetted
)

qs_save(
  THEMIS_DE_MAST_Unsubsetted_WNN_L4Cov, 
  "../Data/FANS/DE/MAST/WNN_L25/THEMIS_DE_MAST_Unsubsetted_WNN_L4Cov.qs2"
)

    ## 3.4 Exc_RORB Cells ------------------------------------------------------

RORB.Unsubsetted <- subset(FANS, subset = ID_WNN_L25_Predicted == "Exc_RORB")
RORB.SubsetCells <- subset(FANS.SubsetCells, subset = ID_WNN_L25_Predicted == "Exc_RORB")
RORB.SubsetCounts <- subset(FANS.SubsetCounts, subset = ID_WNN_L25_Predicted == "Exc_RORB")

RORB_DE_MAST_Unsubsetted <- MAST_DE(
  RORB.Unsubsetted
)
qs_save(
  RORB_DE_MAST_Unsubsetted, 
  "../Data/FANS/DE/RORB_DE_MAST_Unsubsetted.qs2"
)

RORB_DE_MAST_SubsetCells <- MAST_DE(
  RORB.SubsetCells
)
qs_save(
  RORB_DE_MAST_SubsetCells, 
  "../Data/FANS/DE/RORB_DE_MAST_SubsetCells.qs2"
)

RORB_DE_MAST_SubsetCounts <- MAST_DE(
  RORB.SubsetCounts
)
qs_save(
  RORB_DE_MAST_SubsetCounts, 
  "../Data/FANS/DE/RORB_DE_MAST_SubsetCounts.qs2"
)

RORB_DE_MAST_Unsubsetted_WNN_L4Cov <- MAST_DE_WNN_L4(
  RORB.Unsubsetted
)

qs_save(
  RORB_DE_MAST_Unsubsetted_WNN_L4Cov, 
  "../Data/FANS/DE/MAST/WNN_L25/RORB_DE_MAST_Unsubsetted_WNN_L4Cov.qs2"
)


    ## 3.4 Exc_FEZF2 Cells -----------------------------------------------------

FEZF2.Unsubsetted <- subset(FANS, subset = ID_WNN_L25_Predicted == "Exc_FEZF2")
FEZF2.SubsetCells <- subset(FANS.SubsetCells, subset = ID_WNN_L25_Predicted == "Exc_FEZF2")
FEZF2.SubsetCounts <- subset(FANS.SubsetCounts, subset = ID_WNN_L25_Predicted == "Exc_FEZF2")

FEZF2_DE_MAST_Unsubsetted <- MAST_DE(
  FEZF2.Unsubsetted
) 

FEZF2_DE_MAST_Unsubsetted_WNN_L4Cov <- MAST_DE_WNN_L4(
  FEZF2.Unsubsetted
)

qs_save(
  FEZF2_DE_MAST_Unsubsetted, 
  "../Data/FANS/DE/FEZF2_DE_MAST_Unsubsetted.qs2"
)

qs_save(
  FEZF2_DE_MAST_Unsubsetted_WNN_L4Cov, 
  "../Data/FANS/DE/MAST/WNN_L25/FEZF2_DE_MAST_Unsubsetted_WNN_L4Cov.qs2"
)

###
table(FEZF2_DE_MAST_Unsubsetted$fdr<0.05, useNA="always")
table(FEZF2_DE_MAST_Unsubsetted_WNN_L4Cov$fdr<0.05, useNA="always")

FEZF2_DE_DESeq2_Unsubsetted <- qs_read("../Data/FANS/DE/DESeq2/WNN_L25/FEZF2_DE_DESeq2_Unsubsetted.qs2")
table(FEZF2_DE_DESeq2_Unsubsetted$res_L2FC_Shrink$padj<0.05, useNA="always")

F.MAST <- FEZF2_DE_MAST_Unsubsetted$gene[FEZF2_DE_MAST_Unsubsetted$fdr<0.05]
F.MAST_WNNL4 <- FEZF2_DE_MAST_Unsubsetted_WNN_L4Cov$gene[FEZF2_DE_MAST_Unsubsetted_WNN_L4Cov$fdr<0.05]
F.DESeq2 <- rownames(FEZF2_DE_DESeq2_Unsubsetted$res_L2FC_Shrink)[which(FEZF2_DE_DESeq2_Unsubsetted$res_L2FC_Shrink$padj<0.05)]

table(F.MAST %in% F.MAST_WNNL4, useNA="always")
table(F.MAST_WNNL4 %in% F.MAST, useNA="always")

table(F.MAST %in% F.DESeq2, useNA="always")
table(F.DESeq2 %in% F.MAST, useNA="always")

table(F.MAST_WNNL4 %in% F.DESeq2, useNA="always")
table(F.DESeq2 %in% F.MAST_WNNL4, useNA="always")

write.csv(FEZF2_DE_MAST_Unsubsetted_WNN_L4Cov, "/home/veselin/NAS/Labmembers/Grozdanov_Veselin/080925_FEZF2_MAST_WNNL4Cov.csv", quote=FALSE)

###
FEZF2_DE_MAST_SubsetCells <- MAST_DE(
  FEZF2.SubsetCells
)
qs_save(
  FEZF2_DE_MAST_SubsetCells, 
  "../Data/FANS/DE/FEZF2_DE_MAST_SubsetCells.qs2"
)

FEZF2_DE_MAST_SubsetCounts <- MAST_DE(
  FEZF2.SubsetCounts
)
qs_save(
  FEZF2_DE_MAST_SubsetCounts, 
  "../Data/FANS/DE/FEZF2_DE_MAST_SubsetCounts.qs2"
)


    ## 3.5 Oligodendrocytes ------------------------------------------------------

Oligodendrocytes.Unsubsetted <- subset(FANS, subset = ID_WNN_L25_Predicted == "Oligodendrocytes")
Oligodendrocytes.SubsetCells <- subset(FANS.SubsetCells, subset = ID_WNN_L25_Predicted == "Oligodendrocytes")
Oligodendrocytes.SubsetCounts <- subset(FANS.SubsetCounts, subset = ID_WNN_L25_Predicted == "Oligodendrocytes")

Oligodendrocytes_DE_MAST_Unsubsetted <- MAST_DE(
  Oligodendrocytes.Unsubsetted
)
qs_save(
  Oligodendrocytes_DE_MAST_Unsubsetted, 
  "../Data/FANS/DE/Oligodendrocytes_DE_MAST_Unsubsetted.qs2"
)

Oligodendrocytes_DE_MAST_SubsetCells <- MAST_DE(
  Oligodendrocytes.SubsetCells
)
qs_save(
  Oligodendrocytes_DE_MAST_SubsetCells, 
  "../Data/FANS/DE/Oligodendrocytes_DE_MAST_SubsetCells.qs2"
)

Oligodendrocytes_DE_MAST_SubsetCounts <- MAST_DE(
  Oligodendrocytes.SubsetCounts
)
qs_save(
  Oligodendrocytes_DE_MAST_SubsetCounts, 
  "../Data/FANS/DE/Oligodendrocytes_DE_MAST_SubsetCounts.qs2"
)



    ## 3.6 Inh SST neurons -----------------------------------------------------

SST.Unsubsetted <- subset(FANS, subset = ID_WNN_L25_Predicted == "Inh_SST")
SST.SubsetCells <- subset(FANS.SubsetCells, subset = ID_WNN_L25_Predicted == "Inh_SST")
SST.SubsetCounts <- subset(FANS.SubsetCounts, subset = ID_WNN_L25_Predicted == "Inh_SST")

SST_DE_MAST_Unsubsetted <- MAST_DE(
  SST.Unsubsetted
)
qs_save(
  SST_DE_MAST_Unsubsetted, 
  "../Data/FANS/DE/SST_DE_MAST_Unsubsetted.qs2"
)

SST_DE_MAST_SubsetCells <- MAST_DE(
  SST.SubsetCells
)
qs_save(
  SST_DE_MAST_SubsetCells, 
  "../Data/FANS/DE/SST_DE_MAST_SubsetCells.qs2"
)

SST_DE_MAST_SubsetCounts <- MAST_DE(
  SST.SubsetCounts
)
qs_save(
  SST_DE_MAST_SubsetCounts, 
  "../Data/FANS/DE/SST_DE_MAST_SubsetCounts.qs2"
) 




    ## 3.7 LINC00507 FREM3 neurons ---------------------------------------------

LINC00507_FREM3.Unsubsetted <- subset(FANS, subset = ID_WNN_L4_Predicted == "Exc_LINC00507_FREM3")
LINC00507_FREM3.SubsetCells <- subset(FANS.SubsetCells, subset = ID_WNN_L4_Predicted == "Exc_LINC00507_FREM3")
LINC00507_FREM3.SubsetCounts <- subset(FANS.SubsetCounts, subset = ID_WNN_L4_Predicted == "Exc_LINC00507_FREM3")

LINC00507_FREM3_DE_MAST_Unsubsetted <- MAST_DE(
  LINC00507_FREM3.Unsubsetted
)
qs_save(
  LINC00507_FREM3_DE_MAST_Unsubsetted, 
  "../Data/FANS/DE/LINC00507_FREM3_DE_MAST_Unsubsetted.qs2"
)

LINC00507_FREM3_DE_MAST_SubsetCells <- MAST_DE(
  LINC00507_FREM3.SubsetCells
)
qs_save(
  LINC00507_FREM3_DE_MAST_SubsetCells, 
  "../Data/FANS/DE/LINC00507_FREM3_DE_MAST_SubsetCells.qs2"
)

LINC00507_FREM3_DE_MAST_SubsetCounts <- MAST_DE(
  LINC00507_FREM3.SubsetCounts
)
qs_save(
  LINC00507_FREM3_DE_MAST_SubsetCounts, 
  "../Data/FANS/DE/LINC00507_FREM3_DE_MAST_SubsetCounts.qs2"
) 




    ## 3.8 FEZF2 NTNG1 neurons ---------------------------------------------

FEZF2_NTNG1.Unsubsetted <- subset(FANS, subset = ID_WNN_L4_Predicted == "Exc_FEZF2_NTNG1")

FEZF2_NTNG1_DE_MAST_Unsubsetted <- MAST_DE(
  FEZF2_NTNG1.Unsubsetted
)


qs_save(
  FEZF2_NTNG1_DE_MAST_Unsubsetted, 
  "../Data/FANS/DE/MAST/WNN_L4/FEZF2_NTNG1_DE_MAST_Unsubsetted.qs2"
)

    ## 3.9 THEMIS LINC00343 neurons ---------------------------------------------

THEMIS_LINC00343.Unsubsetted <- subset(FANS, subset = ID_WNN_L4_Predicted == "Exc_THEMIS_LINC00343")
THEMIS_LINC00343.SubsetCells <- subset(FANS.SubsetCells, subset = ID_WNN_L4_Predicted == "Exc_THEMIS_LINC00343")
THEMIS_LINC00343.SubsetCounts <- subset(FANS.SubsetCounts, subset = ID_WNN_L4_Predicted == "Exc_THEMIS_LINC00343")

THEMIS_LINC00343_DE_MAST_Unsubsetted <- MAST_DE(
  THEMIS_LINC00343.Unsubsetted
)
qs_save(
  THEMIS_LINC00343_DE_MAST_Unsubsetted, 
  "../Data/FANS/DE/THEMIS_LINC00343_DE_MAST_Unsubsetted.qs2"
)

THEMIS_LINC00343_DE_MAST_SubsetCells <- MAST_DE(
  THEMIS_LINC00343.SubsetCells
)
qs_save(
  THEMIS_LINC00343_DE_MAST_SubsetCells, 
  "../Data/FANS/DE/THEMIS_LINC00343_DE_MAST_SubsetCells.qs2"
)

THEMIS_LINC00343_DE_MAST_SubsetCounts <- MAST_DE(
  THEMIS_LINC00343.SubsetCounts
)
qs_save(
  THEMIS_LINC00343_DE_MAST_SubsetCounts, 
  "../Data/FANS/DE/THEMIS_LINC00343_DE_MAST_SubsetCounts.qs2"
) 



 
    ## 3.10 RORB ADGRL4 neurons ---------------------------------------------

RORB_ADGRL4.Unsubsetted <- subset(FANS, subset = ID_WNN_L4_Predicted == "Exc_RORB_ADGRL4")
RORB_ADGRL4.SubsetCells <- subset(FANS.SubsetCells, subset = ID_WNN_L4_Predicted == "Exc_RORB_ADGRL4")
RORB_ADGRL4.SubsetCounts <- subset(FANS.SubsetCounts, subset = ID_WNN_L4_Predicted == "Exc_RORB_ADGRL4")

RORB_ADGRL4_DE_MAST_Unsubsetted <- MAST_DE(
  RORB_ADGRL4.Unsubsetted
)
qs_save(
  RORB_ADGRL4_DE_MAST_Unsubsetted, 
  "../Data/FANS/DE/RORB_ADGRL4_DE_MAST_Unsubsetted.qs2"
)

RORB_ADGRL4_DE_MAST_SubsetCells <- MAST_DE(
  RORB_ADGRL4.SubsetCells
)
qs_save(
  RORB_ADGRL4_DE_MAST_SubsetCells, 
  "../Data/FANS/DE/RORB_ADGRL4_DE_MAST_SubsetCells.qs2"
)

RORB_ADGRL4_DE_MAST_SubsetCounts <- MAST_DE(
  RORB_ADGRL4.SubsetCounts
)
qs_save(
  RORB_ADGRL4_DE_MAST_SubsetCounts, 
  "../Data/FANS/DE/RORB_ADGRL4_DE_MAST_SubsetCounts.qs2"
) 



 
    ## 3.11 RORB LNX2 neurons ---------------------------------------------

RORB_LNX2.Unsubsetted <- subset(FANS, subset = ID_WNN_L4_Predicted == "Exc_RORB_LNX2")
RORB_LNX2.SubsetCells <- subset(FANS.SubsetCells, subset = ID_WNN_L4_Predicted == "Exc_RORB_LNX2")
RORB_LNX2.SubsetCounts <- subset(FANS.SubsetCounts, subset = ID_WNN_L4_Predicted == "Exc_RORB_LNX2")

RORB_LNX2_DE_MAST_Unsubsetted <- MAST_DE(
  RORB_LNX2.Unsubsetted
)
qs_save(
  RORB_LNX2_DE_MAST_Unsubsetted, 
  "../Data/FANS/DE/RORB_LNX2_DE_MAST_Unsubsetted.qs2"
)

RORB_LNX2_DE_MAST_SubsetCells <- MAST_DE(
  RORB_LNX2.SubsetCells
)
qs_save(
  RORB_LNX2_DE_MAST_SubsetCells, 
  "../Data/FANS/DE/RORB_LNX2_DE_MAST_SubsetCells.qs2"
)

RORB_LNX2_DE_MAST_SubsetCounts <- MAST_DE(
  RORB_LNX2.SubsetCounts
)
qs_save(
  RORB_LNX2_DE_MAST_SubsetCounts, 
  "../Data/FANS/DE/RORB_LNX2_DE_MAST_SubsetCounts.qs2"
) 




  ### 4.0 Save data ------------------------------------------------------------



    ## 4.1 All cells     -------------------------------------------------------

qs_save(
  FANS_DE_MAST_Unsubsetted, 
  "../Data/FANS/DE/FANS_DE_MAST_Unsubsetted.qs2"
)

qs_save(
  FANS_DE_MAST_SubsetCells, 
  "../Data/FANS/DE/FANS_DE_MAST_SubsetCells.qs2"
)

qs_save(
  FANS_DE_MAST_SubsetCounts, 
  "../Data/FANS/DE/FANS_DE_MAST_SubsetCounts.qs2"
)
