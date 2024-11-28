

  ### 0.0 Load libraries -------------------------------------------------------

source("~/ALS_Brain_Multiome.Rcfg") 

library(qs) 
library(Seurat)
library(Signac)
library(tidyverse)




  ### 1.0 Load data ------------------------------------------------------------


      # M0 Seurat single-cell data object --------------------------------------

M0 <- qread(
  paste0(
    "../Data/SeuratObjects/", 
    "M0", 
    ".qrds"
  ), 
  nthr = nthr
)




  ### 2.0 Initiate a hierarchical cell-type annotation df ----------------------

M0_hclust <- data.frame(
  gex_barcode = M0$gex_barcode, 
  atac_barcode = M0$atac_barcode, 
  UMAP1 = M0@reductions$umap@cell.embeddings[,1], 
  UMAP2 = M0@reductions$umap@cell.embeddings[,2], 
  ATAC_UMAP1 = M0@reductions$aumap@cell.embeddings[,1], 
  ATAC_UMAP2 = M0@reductions$aumap@cell.embeddings[,2], 
  WNN_UMAP1 = M0@reductions$wumap@cell.embeddings[,1], 
  WNN_UMAP2 = M0@reductions$wumap@cell.embeddings[,2]
)
rownames(M0_hclust) <- rownames(M0@meta.data)




  ### 3.0 Add RNA hierarchical-level cell-type annotation ----------------------



    ## 3.1 RNA CellType Level 4 (finest clustering) ----------------------------

M0_hclust$RNA_L4 <- M0$CellType_Level5
table(M0_hclust$RNA_L4) |> length()
table(M0_hclust$RNA_L4) 
table(M0_hclust$RNA_L4) |> prop.table() %>% "*"(100) |> round(1) |> sort(decreasing=TRUE)
table(M0_hclust$RNA_L4) |> names() |> sort() 
table(is.na(M0_hclust$RNA_L4))
M0$RNA_L4 <- M0_hclust$RNA_L4



    ## 3.2 RNA CellType Level 3 ------------------------------------------------

M0_hclust$RNA_L3 <- "NA"

M0_hclust$RNA_L3[M0_hclust$RNA_L4 %in% c(
  "Astrocytes_CHI3L1", 
  "Astrocytes_ELMO1", 
  "Astrocytes_HPSE2", 
  "Astrocytes_LHFPL3", 
  "Astrocytes_TNC"
)] <- "Astrocytes"

M0_hclust$RNA_L3[M0_hclust$RNA_L4 %in% c(
  "Choroid_Plexus_Epithelium", 
  "Endothelial_Cells", 
  "Pericytes",
  "VLMC"
)] <- "Epithelial"

M0_hclust$RNA_L3[M0_hclust$RNA_L4 %in% c(
  "Exc_FEZF2_NTNG1", 
  "Exc_FEZF2_PKD2L1", 
  "Exc_FEZF2_TRPM3",  
  "Exc_FEZF2_ZFHX3" 
)] <- "Exc_FEZF2"


M0_hclust$RNA_L3[M0_hclust$RNA_L4 %in% c(
  "Exc_LINC00507_AP001977.1", 
  "Exc_LINC00507_FREM3", 
  "Exc_LINC00507_FRMD4B" 
)] <- "Exc_LINC00507" 

M0_hclust$RNA_L3[M0_hclust$RNA_L4 %in% c(
  "Exc_RORB_ADGRL4", 
  "Exc_RORB_CUX2", 
  "Exc_RORB_ERBB4", 
  "Exc_RORB_LNX2", 
  "Exc_RORB_QKI", 
  "Exc_RORB_RTKN2", 
  "Exc_RORB_TNNT2"
)] <- "Exc_RORB" 

M0_hclust$RNA_L3[M0_hclust$RNA_L4 %in% c(
  "Exc_THEMIS_LINC00343", 
  "Exc_THEMIS_SMYD1"
)] <- "Exc_THEMIS"

M0_hclust$RNA_L3[M0_hclust$RNA_L4 %in% c(
  "Inh_LAMP5_CHST9", 
  "Inh_LAMP5_NMBR", 
  "Inh_LAMP5_RAB11FIP1"
)] <- "Inh_LAMP5" 

M0_hclust$RNA_L3[M0_hclust$RNA_L4 %in% c(
  "Inh_PAX6_NXPH1"
)] <- "Inh_PAX6"

M0_hclust$RNA_L3[M0_hclust$RNA_L4 %in% c(
  "Inh_PVALB_COL15A1", 
  "Inh_PVALB_KCNIP2", 
  "Inh_PVALB_LRIG3", 
  "Inh_PVALB_ZFPM2-AS1"
)] <- "Inh_PVALB" 

M0_hclust$RNA_L3[M0_hclust$RNA_L4 %in% c(
  "Inh_SST_BRINP3", 
  "Inh_SST_EPB41L4A", 
  "Inh_SST_GPC5", 
  "Inh_SST_NPY", 
  "Inh_SST_PAWR", 
  "Inh_SST_SGCZ"
)] <- "Inh_SST"   

M0_hclust$RNA_L3[M0_hclust$RNA_L4 %in% c(
  "Inh_TAFA1_ALCAM", 
  "Inh_TAFA1_FSTL5", 
  "Inh_TAFA1_MEIS2", 
  "Inh_TAFA1_TENM3"
)] <- "Inh_TAFA1"

M0_hclust$RNA_L3[M0_hclust$RNA_L4 %in% c(
  "Inh_VIP_AC006305.1", 
  "Inh_VIP_BSPRY", 
  "Inh_VIP_EXPH5",  
  "Inh_VIP_GRM8", 
  "Inh_VIP_HTR3A", 
  "Inh_VIP_RALYL", 
  "Inh_VIP_SOX11", 
  "Inh_VIP_TAC3"
)] <- "Inh_VIP"   

M0_hclust$RNA_L3[M0_hclust$RNA_L4 %in% c(
  "Lymphocytes"
)] <- "Lymphocytes"  

M0_hclust$RNA_L3[M0_hclust$RNA_L4 %in% c(
  "Microglia"
)] <- "Microglia" 

M0_hclust$RNA_L3[M0_hclust$RNA_L4 %in% c(
  "Oligodendrocytes"   
)] <- "Oligodendrocytes"  

M0_hclust$RNA_L3[M0_hclust$RNA_L4 %in% c(
  "OPC_ADGRV1", 
  "OPC_PMP2", 
  "OPC_S100B", 
  "Premyelinating_Oligodendrocytes"
)] <- "OPC" 


table(M0_hclust$RNA_L3)[order(names(table(M0_hclust$RNA_L3)))] 
table(is.na(M0_hclust$RNA_L3))
M0$RNA_L3 <- M0_hclust$RNA_L3



    ## 3.3 RNA CellType Level 2.5 ----------------------------------------------

M0_hclust$RNA_L25 <- M0_hclust$RNA_L3

M0_hclust$RNA_L25[M0_hclust$RNA_L3 %in% c(
  "Inh_TAFA1", 
  "Inh_VIP"
)] <- "Inh_TAFA1_VIP"

M0_hclust$RNA_L25[M0_hclust$RNA_L3 %in% c(
  "Inh_LAMP5", 
  "Inh_PAX6" 
)] <- "Inh_LAMP5_PAX6"  

table(M0_hclust$RNA_L25)[order(names(table(M0_hclust$RNA_L25)))]
M0$RNA_L25 <- M0_hclust$RNA_L25



    ## 3.4 RNA CellType Level 2 ------------------------------------------------

M0_hclust$RNA_L2 <- M0_hclust$RNA_L3

M0_hclust$RNA_L2[M0_hclust$RNA_L3 %in% c(
  "Exc_FEZF2", 
  "Exc_LINC00507", 
  "Exc_RORB", 
  "Exc_THEMIS" 
)] <- "Exc_Neurons"

M0_hclust$RNA_L2[M0_hclust$RNA_L3 %in% c(
  "Inh_LAMP5", 
  "Inh_PAX6", 
  "Inh_PVALB", 
  "Inh_SST", 
  "Inh_TAFA1", 
  "Inh_VIP"
)] <- "Inh_Neurons"

M0_hclust$RNA_L2[M0_hclust$RNA_L3 %in% c(
  "Epithelial", 
  "Lymphocytes"
)] <- "Other_Cells"  

M0_hclust$RNA_L2[M0_hclust$RNA_L3 %in% c(
  "OPC"
)] <- "OPC"  


table(M0_hclust$RNA_L2)[order(names(table(M0_hclust$RNA_L2)))]
M0$RNA_L2 <- M0_hclust$RNA_L2



  ## 3.5 RNA CellType Level 1.5 ------------------------------------------------

M0_hclust$RNA_L15 <- M0_hclust$RNA_L2


M0_hclust$RNA_L15[M0_hclust$RNA_L2 %in% c(
  "Astrocytes", 
  "Microglia", 
  "Oligodendrocytes", 
  "OPC", 
  "Other_Cells"
)] <- "Glia"


table(M0_hclust$RNA_L15)[order(names(table(M0_hclust$RNA_L15)))]

M0$RNA_L15 <- M0_hclust$RNA_L15



  ## 3.6 RNA CellType Level 1 --------------------------------------------------

M0_hclust$RNA_L1 <- M0_hclust$RNA_L2

M0_hclust$RNA_L1[M0_hclust$RNA_L2 %in% c(
  "Exc_Neurons", 
  "Inh_Neurons"
)] <- "Neuronal"

M0_hclust$RNA_L1[M0_hclust$RNA_L2 %in% c(
  "Astrocytes", 
  "Microglia", 
  "Oligodendrocytes", 
  "OPC", 
  "Other_Cells"
)] <- "Non-neuronal"

table(M0_hclust$RNA_L1)[order(names(table(M0_hclust$RNA_L1)))]
M0$RNA_L1 <- M0_hclust$RNA_L1




  ### 4.0 Add ATAC hierarchical-level cell-type annotation ---------------------



    ## 4.1 ATAC CellType Level 4 (finest clustering) ---------------------------

M0_hclust$ATAC_L4 <- M0$CellType_Level5_ATAC
table(M0_hclust$ATAC_L4) |> length()
table(M0_hclust$ATAC_L4) 
table(M0_hclust$ATAC_L4) |> prop.table() %>% "*"(100) |> round(1) |> sort(decreasing=TRUE)
table(M0_hclust$ATAC_L4) |> names() |> sort()
M0$ATAC_L4 <- M0_hclust$ATAC_L4



    ## 4.2 ATAC CellType Level 3 -----------------------------------------------

M0_hclust$ATAC_L3 <- "NA"

M0_hclust$ATAC_L3[M0_hclust$ATAC_L4 %in% c(
  "Astrocytes_CHI3L1", 
  "Astrocytes_ELMO1", 
  "Astrocytes_HPSE2", 
  "Astrocytes_ROBO2", 
  "Astrocytes_TNC" 
)] <- "Astrocytes"

M0_hclust$ATAC_L3[M0_hclust$ATAC_L4 %in% c(
  "Choroid_Plexus_Epithelium", 
  "Endothelial_Cells",
  "Pericytes", 
  "VLMC"
)] <- "Epithelial"


M0_hclust$ATAC_L3[M0_hclust$ATAC_L4 %in% c(
  "Exc_FEZF2_NTNG1", 
  "Exc_FEZF2_PKD2L1", 
  "Exc_FEZF2_TRPM3", 
  "Exc_FEZF2_ZFHX3" 
)] <- "Exc_FEZF2"


M0_hclust$ATAC_L3[M0_hclust$ATAC_L4 %in% c(
  "Exc_LINC00507_AP001977.1", 
  "Exc_LINC00507_FREM3", 
  "Exc_LINC00507_FRMD4B", 
  "Exc_LINC00507_ZNF385D"   
)] <- "Exc_LINC00507" 


M0_hclust$ATAC_L3[M0_hclust$ATAC_L4 %in% c(
  "Exc_RORB_ADGRL4", 
  "Exc_RORB_CUX2", 
  "Exc_RORB_ERBB4", 
  "Exc_RORB_LNX2", 
  "Exc_RORB_PDE1C", 
  "Exc_RORB_RTKN2", 
  "Exc_RORB_TNNT2"
)] <- "Exc_RORB" 


M0_hclust$ATAC_L3[M0_hclust$ATAC_L4 %in% c(
  "Exc_THEMIS_LINC00343", 
  "Exc_THEMIS_SMYD1"        
)] <- "Exc_THEMIS"


M0_hclust$ATAC_L3[M0_hclust$ATAC_L4 %in% c(
  "Inh_LAMP5_CHST9", 
  "Inh_LAMP5_NMBR", 
  "Inh_LAMP5_PIP4K2A", 
  "Inh_LAMP5_RAB11FIP1", 
  "Inh_LAMP5_SPOCK3"       
)] <- "Inh_LAMP5" 


M0_hclust$ATAC_L3[M0_hclust$ATAC_L4 %in% c(
  "Inh_PAX6_NXPH1"
)] <- "Inh_PAX6"


M0_hclust$ATAC_L3[M0_hclust$ATAC_L4 %in% c(
  "Inh_PVALB_COL15A1", 
  "Inh_PVALB_KCNIP2", 
  "Inh_PVALB_LRIG3", 
  "Inh_PVALB_NCKAP5", 
  "Inh_PVALB_ZFPM2-AS1"
)] <- "Inh_PVALB" 


M0_hclust$ATAC_L3[M0_hclust$ATAC_L4 %in% c(
  "Inh_SST_BRINP3", 
  "Inh_SST_CTNNA3", 
  "Inh_SST_EPB41L4A", 
  "Inh_SST_NPY", 
  "Inh_SST_PAWR", 
  "Inh_SST_SGCZ"    
)] <- "Inh_SST"   


M0_hclust$ATAC_L3[M0_hclust$ATAC_L4 %in% c(
  "Inh_TAFA1_ALCAM", 
  "Inh_TAFA1_FSTL5", 
  "Inh_TAFA1_PIP4K2A", 
  "Inh_TAFA1_TENM3"
)] <- "Inh_TAFA1"


M0_hclust$ATAC_L3[M0_hclust$ATAC_L4 %in% c(
  "Inh_VIP_AC006305.1", 
  "Inh_VIP_BSPRY", 
  "Inh_VIP_EXPH5", 
  "Inh_VIP_GRM8", 
  "Inh_VIP_RALYL", 
  "Inh_VIP_RNF220", 
  "Inh_VIP_SLC1A2", 
  "Inh_VIP_SOX11", 
  "Inh_VIP_TAC3"
)] <- "Inh_VIP"   


M0_hclust$ATAC_L3[M0_hclust$ATAC_L4 %in% c(
  "Lymphocytes"
)] <- "Lymphocytes"  


M0_hclust$ATAC_L3[M0_hclust$ATAC_L4 %in% c(
  "Microglia"
)] <- "Microglia" 


M0_hclust$ATAC_L3[M0_hclust$ATAC_L4 %in% c(
  "Oligodendrocytes" 
)] <- "Oligodendrocytes"  


M0_hclust$ATAC_L3[M0_hclust$ATAC_L4 %in% c(
  "OPC_ADGRV1", 
  "OPC_PLP1", 
  "OPC_PMP2",
  "OPC_S100B", 
  "Premyelinating_Oligodendrocytes" 
)] <- "OPC" 


table(M0_hclust$ATAC_L3)[order(names(table(M0_hclust$ATAC_L3)))]
M0$ATAC_L3 <- M0_hclust$ATAC_L3



    ## 4.3 ATAC CellType Level 2.5 ---------------------------------------------

M0_hclust$ATAC_L25 <- M0_hclust$ATAC_L3 

M0_hclust$ATAC_L25[M0_hclust$ATAC_L3 %in% c(
  "Inh_LAMP5", 
  "Inh_PAX6"
)] <- "Inh_LAMP5_PAX6"

M0_hclust$ATAC_L25[M0_hclust$ATAC_L3 %in% c(
  "Inh_TAFA1", 
  "Inh_VIP"
)] <- "Inh_TAFA1_VIP"

table(M0_hclust$ATAC_L25)[order(names(table(M0_hclust$ATAC_L25)))]
M0$ATAC_L25 <- M0_hclust$ATAC_L25 



    ## 4.4 ATAC CellType Level 2 -----------------------------------------------

M0_hclust$ATAC_L2 <- M0_hclust$ATAC_L3 

M0_hclust$ATAC_L2[M0_hclust$ATAC_L3 %in% c(
  "Exc_FEZF2", 
  "Exc_LINC00507", 
  "Exc_RORB", 
  "Exc_THEMIS"       
)] <- "Exc_Neurons"

M0_hclust$ATAC_L2[M0_hclust$ATAC_L3 %in% c(
  "Inh_LAMP5", 
  "Inh_PAX6", 
  "Inh_PVALB", 
  "Inh_SST", 
  "Inh_TAFA1", 
  "Inh_VIP"
)] <- "Inh_Neurons"

M0_hclust$ATAC_L2[M0_hclust$ATAC_L3 %in% c(
  "Epithelial", 
  "Lymphocytes"
)] <- "Other_Cells"  

M0_hclust$ATAC_L2[M0_hclust$ATAC_L3 %in% c(
  "OPC"
)] <- "OPC"  


table(M0_hclust$ATAC_L2)[order(names(table(M0_hclust$ATAC_L2)))]
M0$ATAC_L2 <- M0_hclust$ATAC_L2 



    ## 4.5 ATAC CellType Level 1.5 ---------------------------------------------

M0_hclust$ATAC_L15 <- M0_hclust$ATAC_L2 

M0_hclust$ATAC_L15[M0_hclust$ATAC_L2 %in% c(
  "Astrocytes", 
  "Microglia", 
  "Neurons", 
  "Oligodendrocytes", 
  "OPC", 
  "Other_Cells" 
)] <- "Glia"


table(M0_hclust$ATAC_L15)[order(names(table(M0_hclust$ATAC_L15)))]
M0$ATAC_L15 <- M0_hclust$ATAC_L15 



    ## 4.6 ATAC CellType Level 1 -----------------------------------------------

M0_hclust$ATAC_L1 <- M0_hclust$ATAC_L2 

M0_hclust$ATAC_L1[M0_hclust$ATAC_L2 %in% c(
  "Exc_Neurons", 
  "Inh_Neurons"
)] <- "Neuronal"


M0_hclust$ATAC_L1[M0_hclust$ATAC_L2 %in% c(
  "Astrocytes", 
  "Microglia", 
  "Neurons", 
  "Oligodendrocytes", 
  "OPC", 
  "Other_Cells" 
)] <- "Non-neuronal"


table(M0_hclust$ATAC_L1)[order(names(table(M0_hclust$ATAC_L1)))]
M0$ATAC_L1 <- M0_hclust$ATAC_L1 




  ### 5.0 Add WNN hierarchical-level cell-type annotation ----------------------



    ## 5.1 WNN CellType Level 4 (finest clustering) ----------------------------

M0_hclust$WNN_L4 <- M0$CellType_Level5_Weighted
table(M0_hclust$WNN_L4) |> length()
table(M0_hclust$WNN_L4) 
table(M0_hclust$WNN_L4) |> prop.table() %>% "*"(100) |> round(1) |> sort(decreasing=TRUE)
table(M0_hclust$WNN_L4) |> names() |> sort()
M0$WNN_L4 <- M0_hclust$WNN_L4



    ## 5.2 WNN CellType Level 3 ------------------------------------------------

M0_hclust$WNN_L3 <- "NA"


M0_hclust$WNN_L3[M0_hclust$WNN_L4 %in% c(
  "Astrocytes_CHI3L1", 
  "Astrocytes_ELMO1", 
  "Astrocytes_HPSE2",  
  "Astrocytes_LHFPL3", 
  "Astrocytes_PLSCR4", 
  "Astrocytes_ROBO2", 
  "Astrocytes_TNC"   
)] <- "Astrocytes"


M0_hclust$WNN_L3[M0_hclust$WNN_L4 %in% c(
  "Choroid_Plexus_Epithelium", 
  "Endothelial_Cells", 
  "Pericytes", 
  "VLMC"
)] <- "Epithelial"


M0_hclust$WNN_L3[M0_hclust$WNN_L4 %in% c(
  "Exc_FEZF2_NTNG1", 
  "Exc_FEZF2_PKD2L1", 
  "Exc_FEZF2_TRHDE", 
  "Exc_FEZF2_TRPM3", 
  "Exc_FEZF2_ZFHX3"            
)] <- "Exc_FEZF2"


M0_hclust$WNN_L3[M0_hclust$WNN_L4 %in% c(
  "Exc_LINC00507_AP001977.1", 
  "Exc_LINC00507_FREM3", 
  "Exc_LINC00507_FRMD4B", 
  "Exc_LINC00507_ZNF385D"   
)] <- "Exc_LINC00507" 


M0_hclust$WNN_L3[M0_hclust$WNN_L4 %in% c(
  "Exc_RORB_ADGRL4", 
  "Exc_RORB_chr9-113535578-113537604", 
  "Exc_RORB_CUX2", 
  "Exc_RORB_ERBB4", 
  "Exc_RORB_LHFPL3", 
  "Exc_RORB_LNX2", 
  "Exc_RORB_PDE1C", 
  "Exc_RORB_QKI", 
  "Exc_RORB_RTKN2", 
  "Exc_RORB_TNNT2", 
  "Exc_RORB_TSHZ2"
)] <- "Exc_RORB" 


M0_hclust$WNN_L3[M0_hclust$WNN_L4 %in% c(
  "Exc_THEMIS_LINC00343", 
  "Exc_THEMIS_MBP", 
  "Exc_THEMIS_SMYD1"
)] <- "Exc_THEMIS"


M0_hclust$WNN_L3[M0_hclust$WNN_L4 %in% c(
  "Inh_LAMP5_CHST9", 
  "Inh_LAMP5_NMBR", 
  "Inh_LAMP5_PIP4K2A", 
  "Inh_LAMP5_RAB11FIP1", 
  "Inh_LAMP5_SPOCK3" 
)] <- "Inh_LAMP5" 


M0_hclust$WNN_L3[M0_hclust$WNN_L4 %in% c(
  "Inh_PAX6_NXPH1" 
)] <- "Inh_PAX6"


M0_hclust$WNN_L3[M0_hclust$WNN_L4 %in% c(
  "Inh_PVALB_COL15A1", 
  "Inh_PVALB_KCNIP2", 
  "Inh_PVALB_LRIG3", 
  "Inh_PVALB_NCKAP5", 
  "Inh_PVALB_PLXDC2", 
  "Inh_PVALB_ZFPM2-AS1"      
)] <- "Inh_PVALB" 


M0_hclust$WNN_L3[M0_hclust$WNN_L4 %in% c(
  "Inh_SST_BRINP3", 
  "Inh_SST_CTNNA3", 
  "Inh_SST_EPB41L4A", 
  "Inh_SST_GPC5", 
  "Inh_SST_NPY", 
  "Inh_SST_PAWR", 
  "Inh_SST_SGCZ"  
)] <- "Inh_SST"   


M0_hclust$WNN_L3[M0_hclust$WNN_L4 %in% c(
  "Inh_TAFA1_ALCAM", 
  "Inh_TAFA1_FSTL5", 
  "Inh_TAFA1_GRIA4", 
  "Inh_TAFA1_PIP4K2A", 
  "Inh_TAFA1_TENM3"    
)] <- "Inh_TAFA1"


M0_hclust$WNN_L3[M0_hclust$WNN_L4 %in% c(
  "Inh_VIP_AC006305.1", 
  "Inh_VIP_BSPRY", 
  "Inh_VIP_EXPH5", 
  "Inh_VIP_GRM8", 
  "Inh_VIP_HTR3A", 
  "Inh_VIP_RALYL", 
  "Inh_VIP_RNF220", 
  "Inh_VIP_SOX11",
  "Inh_VIP_TAC3" 
)] <- "Inh_VIP"   


M0_hclust$WNN_L3[M0_hclust$WNN_L4 %in% c(
  "Lymphocytes" 
)] <- "Lymphocytes"  


M0_hclust$WNN_L3[M0_hclust$WNN_L4 %in% c(
  "Microglia"  
)] <- "Microglia" 


M0_hclust$WNN_L3[M0_hclust$WNN_L4 %in% c(
  "Oligodendrocytes" 
)] <- "Oligodendrocytes"  


M0_hclust$WNN_L3[M0_hclust$WNN_L4 %in% c(
  "OPC_ADGRV1", 
  "OPC_PLP1", 
  "OPC_PMP2", 
  "OPC_S100B", 
  "Premyelinating_Oligodendrocytes" 
)] <- "OPC" 


table(M0_hclust$WNN_L3)[order(names(table(M0_hclust$WNN_L3)))]
M0$WNN_L3 <- M0_hclust$WNN_L3



    ## 5.3 WNN CellType Level 2.5 ----------------------------------------------

M0_hclust$WNN_L25 <- M0_hclust$WNN_L3

M0_hclust$WNN_L25[M0_hclust$WNN_L3 %in% c(
  "Inh_LAMP5", 
  "Inh_PAX6" 
)] <- "Inh_LAMP5_PAX6" 

M0_hclust$WNN_L25[M0_hclust$WNN_L3 %in% c(
  "Inh_TAFA1", 
  "Inh_VIP"
)] <- "Inh_TAFA1_VIP" 


table(M0_hclust$WNN_L25)[order(names(table(M0_hclust$WNN_L25)))]
M0$WNN_L25 <- M0_hclust$WNN_L25



    ## 5.4 WNN CellType Level 2 ------------------------------------------------

M0_hclust$WNN_L2 <- M0_hclust$WNN_L3

M0_hclust$WNN_L2[M0_hclust$WNN_L3 %in% c(
  "Exc_FEZF2", 
  "Exc_LINC00507", 
  "Exc_RORB", 
  "Exc_THEMIS"       
)] <- "Exc_Neurons"

M0_hclust$WNN_L2[M0_hclust$WNN_L3 %in% c(
  "Inh_LAMP5", 
  "Inh_PAX6", 
  "Inh_PVALB", 
  "Inh_SST", 
  "Inh_TAFA1", 
  "Inh_VIP"
)] <- "Inh_Neurons"

M0_hclust$WNN_L2[M0_hclust$WNN_L3 %in% c(
  "Epithelial", 
  "Lymphocytes"
)] <- "Other_Cells"  

M0_hclust$WNN_L2[M0_hclust$WNN_L3 %in% c(
  "OPC" 
)] <- "OPC"  


table(M0_hclust$WNN_L2)[order(names(table(M0_hclust$WNN_L2)))]
M0$WNN_L2 <- M0_hclust$WNN_L2



    ## 5.5 WNN CellType Level 1.5 ----------------------------------------------

M0_hclust$WNN_L15 <- M0_hclust$WNN_L2 

M0_hclust$WNN_L15[M0_hclust$WNN_L2 %in% c(
  "Astrocytes", 
  "Microglia", 
  "Neurons", 
  "Oligodendrocytes", 
  "OPC", 
  "Other_Cells" 
)] <- "Glia" 



table(M0_hclust$WNN_L15)[order(names(table(M0_hclust$WNN_L15)))]
M0$WNN_L15 <- M0_hclust$WNN_L15



    ## 5.6 WNN CellType Level 1 ------------------------------------------------

M0_hclust$WNN_L1 <- M0_hclust$WNN_L2 

M0_hclust$WNN_L1[M0_hclust$WNN_L2 %in% c(
  "Exc_Neurons", 
  "Inh_Neurons"
)] <- "Neuronal"


M0_hclust$WNN_L1[M0_hclust$WNN_L2 %in% c(
  "Astrocytes", 
  "Microglia", 
  "Neurons", 
  "Oligodendrocytes", 
  "OPC", 
  "Other_Cells" 
)] <- "Non-neuronal" 



table(M0_hclust$WNN_L1)[order(names(table(M0_hclust$WNN_L1)))]
M0$WNN_L1 <- M0_hclust$WNN_L1




  ### 6.0 Add "AllCells" Label ------------------------------------------------- 

M0$AllCells <- "AllCells"




  ### 7.0 Visualize overlap of hierarchical cell-type annotations --------------

table(M0$ATAC_L1, M0$RNA_L1)
table(M0$ATAC_L15, M0$RNA_L15)
table(M0$ATAC_L15, M0$RNA_L15) %>% 
  "/"(rowSums(table(M0$ATAC_L15, M0$RNA_L15))) %>% 
    "*"(100) %>% 
      round(1)

table(M0$ATAC_L2, M0$RNA_L2)

table(M0$RNA_L1, M0$WNN_L1)
table(M0$RNA_L15, M0$WNN_L15)
table(M0$RNA_L15, M0$WNN_L15) %>% 
  "/"(rowSums(table(M0$RNA_L15, M0$WNN_L15))) %>% 
    "*"(100) %>% 
      round(1)

table(M0$RNA_L2, M0$WNN_L2)

table(M0$ATAC_L1, M0$WNN_L1)
table(M0$ATAC_L15, M0$WNN_L15)
table(M0$ATAC_L15, M0$WNN_L15) %>% 
  "/"(rowSums(table(M0$ATAC_L15, M0$WNN_L15))) %>% 
  "*"(100) %>% 
  round(1)

table(M0$ATAC_L2, M0$WNN_L2)




  ### 8.0 Export data ----------------------------------------------------------

qsave(
  M0_hclust,
  paste0(
    "../Data/Annotations/", 
    "M0_hclust", 
    ".qrds"
  ), 
  nthr = nthr
)


qsave(
  M0,
  paste0(
    "../Data/SeuratObjects/", 
    "M0", 
    ".qrds"
  ), 
  nthr = nthr
)



