
  ###  0.0 Load libraries ------------------------------------------------------ 
source("~/ALS_Brain_Multiome.Rcfg") 

library(qs)
library(tidyverse)
library(clusterProfiler)

  
qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)

 


  ### 1.0 Load data ------------------------------------------------------------

ALS_ALSFTD_Sign <- qread(
  "../Data/Annotations/Signatures/RNA/Signature_RNA_ALS_ALSFTD_WNN_L25.qrds"
)

GEX_Features <- qread(
  paste0(
    "../Data/Annotations/Genes/", 
    "GEX_Features", 
    ".qrds"
  ), 
  nthr
)

GEX_Features_ENSEMBL_ENTREZ <- qread(
  paste0(
    "../Data/Annotations/Genes/", 
    "GEX_Features_ENSEMBL_ENTREZ", 
    ".qrds"
  ), 
  nthr
)




  ### 2.0 KEGG Enrichment ------------------------------------------------------

Genes = ALS_ALSFTD_Sign$Gene
Genes = data.frame(Gene=Genes)
Genes$ENTREZID <- GEX_Features_ENSEMBL_ENTREZ$ENTREZID[match(Genes$Gene, GEX_Features_ENSEMBL_ENTREZ$ID_10X)]

genelist = Genes$ENTREZID
Enrich_Kegg_ALL_AND_WNN_L25 <- enrichKEGG(
  gene = genelist, 
  organism = "hsa", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05 
  #universe = GEX
)
View(as.data.frame(Enrich_Kegg_ALL_AND_WNN_L25))


  ### 3.0 Export table ---------------------------------------------------------

write.csv(
  as.data.frame(Enrich_Kegg_ALL_AND_WNN_L25), 
  "../Data/Visualization/Tables/SupplTab_ALS_ALSFTD_WNN_L25_Signature_KEGG_Enrichment.csv", 
  quote = FALSE, 
  row.names = FALSE
)
