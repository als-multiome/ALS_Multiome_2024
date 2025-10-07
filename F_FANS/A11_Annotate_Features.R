# A8_Annotate_Genes.R 

  ### 0.0 Load libraries ------------------------------------------------------- 

library(qs2)
library(tidyverse)
library(Seurat)
library(org.Hs.eg.db)  



      
  ### 1.0 Load data ------------------------------------------------------------
  
Samples <- qs_read(
  "../Data/FANS/Samples/Samples.qs2"
)

FANS <- qs_read(
  "../Data/FANS/SeuratObjects/FANS_Merged_QCFiltered.qs2"
)




  ### 2.0 Define AUX functions -------------------------------------------------

tab_intersect <- function(a,b, plot = TRUE){
  print(c("a"=as.numeric(table(a %in% b)["FALSE"]), "  a&b  "=as.numeric(table(a %in% b)["TRUE"]), "b"=as.numeric(table(b %in% a)["FALSE"])))
  if(plot){
    ggplot(
      data.frame(
        a=as.numeric(table(a %in% b)["FALSE"]), 
        ab=as.numeric(table(a %in% b)["TRUE"])), 
      b=as.numeric(table(b %in% a)["FALSE"])  
    )
  }
}




  ### 3.0 Feature annotation ---------------------------------------------------

any(duplicated(rownames(FANS)) | is.na(rownames(FANS)))

FANS_RNA_Features <- data.frame(
  ID_10X = rownames(FANS)
)



    ## 3.1 Collect Gene Symbols ------------------------------------------------

all_OrgHs_Symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")

tab_intersect(a=FANS_RNA_Features$ID_10X, b=all_OrgHs_Symbols)



    ## 3.2 Replace ALIASes in RNA features with gene SYMBOLS -------------------

grep("MT-", FANS_RNA_Features$ID_10X, value=TRUE)

tmp <- select(
  org.Hs.eg.db, 
  keys=str_replace_all(
    FANS_RNA_Features$ID_10X[!FANS_RNA_Features$ID_10X %in% all_OrgHs_Symbols], 
    "MT-", 
    "MT"
  ), 
  columns = "SYMBOL", 
  keytype = "ALIAS", 
  multiVals='list'
)

tmp <- tmp[!is.na(tmp$SYMBOL),]
any(is.na(tmp$ALIAS))
tmp <- unique(tmp)

table(duplicated(tmp$ALIAS))
table(duplicated(tmp$SYMBOL))

View(tmp[tmp$ALIAS %in% tmp$ALIAS[duplicated(tmp$ALIAS)],])

tmp$ALIAS[duplicated(tmp$ALIAS)]

FANS_RNA_Features$SYMBOL <- str_replace_all(
  FANS_RNA_Features$ID_10X, 
  "MT-", 
  "MT"
)
FANS_RNA_Features$SYMBOL[FANS_RNA_Features$SYMBOL %in% tmp$ALIAS] <- tmp$SYMBOL[match(
  FANS_RNA_Features$SYMBOL[FANS_RNA_Features$SYMBOL %in% tmp$ALIAS], tmp$ALIAS
)]

FANS_RNA_Features$SYMBOL[which(!FANS_RNA_Features$SYMBOL %in% keys(org.Hs.eg.db, keytype="SYMBOL"))] <- NA
table(is.na(FANS_RNA_Features$SYMBOL))





    ## 3.3 Mark GenBank accessions ---------------------------------------------

FANS_RNA_Features$Accession <- "NotAccession"
FANS_RNA_Features$Accession[grep("\\.", FANS_RNA_Features$ID_10X)] <- "Accession"

View(FANS_RNA_Features[FANS_RNA_Features$Accession=="Accession" & !is.na(FANS_RNA_Features$SYMBOL),])
FANS_RNA_Features$Accession[FANS_RNA_Features$Accession=="Accession" & !is.na(FANS_RNA_Features$SYMBOL)] <- "NotAccession"



table(
  is.na(FANS_RNA_Features$SYMBOL), 
  FANS_RNA_Features$Accession
) 


View(FANS_RNA_Features[FANS_RNA_Features$Accession=="NotAccession" & is.na(FANS_RNA_Features$SYMBOL),])



    ## 3.4 Annotate GENETYPEs --------------------------------------------------

GENETYPES <- AnnotationDbi::select(
  org.Hs.eg.db, 
  keys=FANS_RNA_Features$SYMBOL, 
  keytype="SYMBOL", 
  columns="GENETYPE"
)


GENETYPES <- unique(GENETYPES)
table(GENETYPES$GENETYPE)

table(duplicated(GENETYPES$SYMBOL))
View(GENETYPES[GENETYPES$SYMBOL %in% GENETYPES$SYMBOL[duplicated(GENETYPES$SYMBOL)],])
GENETYPES <- GENETYPES[-which(GENETYPES$SYMBOL %in% GENETYPES$SYMBOL[duplicated(GENETYPES$SYMBOL)] & GENETYPES$GENETYPE=="unknown"),]
table(duplicated(GENETYPES$SYMBOL))
table(GENETYPES$GENETYPE)

FANS_RNA_Features$GENETYPE <- NA 
FANS_RNA_Features$GENETYPE <- GENETYPES$GENETYPE[match(FANS_RNA_Features$SYMBOL, GENETYPES$SYMBOL)]

table(is.na(FANS_RNA_Features$GENETYPE), FANS_RNA_Features$Accession)
FANS_RNA_Features$GENETYPE[grep("MT-", FANS_RNA_Features$ID_10X)] <- "mitochondrial-encoded"
View(FANS_RNA_Features[FANS_RNA_Features$Accession!="Accession" & is.na(FANS_RNA_Features$GENETYPE),])


FANS_RNA_Features$GENETYPE2 <-FANS_RNA_Features$GENETYPE
FANS_RNA_Features$GENETYPE2[FANS_RNA_Features$Accession=="Accession"] <- "Accession"
FANS_RNA_Features$GENETYPE2[FANS_RNA_Features$Accession!="Accession" & is.na(FANS_RNA_Features$GENETYPE2)] <- "unknown"

tmp2 <- FANS_RNA_Features$GENETYPE
FANS_RNA_Features$GENETYPE <- FANS_RNA_Features$GENETYPE2
FANS_RNA_Features$GENETYPE2 <- tmp2
rm(tmp2)




  ### 4.0 Annotate ENSEMBL IDs -------------------------------------------------

ENSEMBLS <- AnnotationDbi::select(
  org.Hs.eg.db, 
  keys=FANS_RNA_Features$SYMBOL, 
  keytype="SYMBOL", 
  columns="ENSEMBL"
)


ENSEMBLS <- unique(ENSEMBLS)


FANS_RNA_Features_ENSEMBL_ENTREZ <- full_join(
  FANS_RNA_Features, 
  ENSEMBLS, 
  by=c("SYMBOL"="SYMBOL"), 
  relationship = "many-to-many"
)




  ### 5.0 Annotate ENSEMBL IDs and ENTREZIDs -----------------------------------

ENTREZIDS <- AnnotationDbi::select(
  org.Hs.eg.db, 
  keys=FANS_RNA_Features_ENSEMBL_ENTREZ$ENSEMBL, 
  keytype="ENSEMBL", 
  columns="ENTREZID"
)


ENTREZIDS <- unique(ENTREZIDS)
FANS_RNA_Features_ENSEMBL_ENTREZ <- full_join(
  FANS_RNA_Features_ENSEMBL_ENTREZ, 
  ENTREZIDS, 
  by=c("ENSEMBL"="ENSEMBL")
)




  ### 6.0 Annotate Chromosomes -------------------------------------------------

CHRs <- AnnotationDbi::select(
  org.Hs.eg.db, 
  keys=FANS_RNA_Features$SYMBOL, 
  keytype="SYMBOL", 
  columns="CHR"
)


FANS_RNA_Features$CHR <- CHRs$CHR[
  match(
    FANS_RNA_Features$SYMBOL, 
    CHRs$SYMBOL
  )
]




  ### 7.0 Save data ------------------------------------------------------------

qs_save(
  FANS_RNA_Features, 
  "../Data/FANS/Annotations/FANS_RNA_Features.qs2"
)

qs_save(
  FANS_RNA_Features_ENSEMBL_ENTREZ, 
  "../Data/FANS/Annotations/FANS_RNA_Features_ENSEMBL_ENTREZ.qs2"
)




