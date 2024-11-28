

  ### 0.0 Load libraries -------------------------------------------------------

source("~/ALS_Brain_Multiome.Rcfg") 

library(qs)
library(tidyverse)
library(org.Hs.eg.db)


qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)




  ### 1.0 Load data ------------------------------------------------------------

M0_RNA <- qread(
  paste0(
    "../Data/SeuratObjects/", 
    "M0_RNA", 
    ".qrds"
  ), 
  nthr=nthr
)




  ### 2.0 Define AUX functions -------------------------------------------------

tab_intersect <- function(a,b, plot=TRUE){
  print(c("a"=as.numeric(table(a %in% b)["FALSE"]), "  a&b  "=as.numeric(table(a %in% b)["TRUE"]), "b"=as.numeric(table(b %in% a)["FALSE"])))
  if(plot){
    ggplot(
      data.frame(
        a=as.numeric(table(a %in% b)["FALSE"]), 
        ab=))
  }
}


tab_intersect(a, b)




  ### 2.0 Annotate gene features -----------------------------------------------



    ## 2.1 Collect M0 GEX features ---------------------------------------------

any(duplicated(rownames(M0_RNA)) | is.na(rownames(M0_RNA)))

GEX_Features <- data.frame(
  ID_10X = rownames(M0_RNA)
)



    ## 2.2 Collect Gene Symbols ------------------------------------------------

all_OrgHs_Symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")

tab_intersect(GEX_Features$ID_10X, all_OrgHs_Symbols)



    ## 2.3 Replace ALIASes in GEX features with gene SYMBOLS -------------------

grep("MT-", GEX_Features$ID_10X, value=TRUE)

tmp <- select(
  org.Hs.eg.db, 
  keys=str_replace_all(
    GEX_Features$ID_10X[!GEX_Features$ID_10X %in% all_OrgHs_Symbols], 
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

GEX_Features$SYMBOL <- str_replace_all(
  GEX_Features$ID_10X, 
  "MT-", 
  "MT"
)
GEX_Features$SYMBOL[GEX_Features$SYMBOL %in% tmp$ALIAS] <- tmp$SYMBOL[match(GEX_Features$SYMBOL[GEX_Features$SYMBOL %in% tmp$ALIAS], tmp$ALIAS)]

GEX_Features$SYMBOL[which(!GEX_Features$SYMBOL %in% keys(org.Hs.eg.db, keytype="SYMBOL"))] <- NA
table(is.na(GEX_Features$SYMBOL))





    ## 2.4 Mark GenBank accessions ---------------------------------------------

GEX_Features$Accession <- "NotAccession"
GEX_Features$Accession[grep("\\.", GEX_Features$ID_10X)] <- "Accession"

View(GEX_Features[GEX_Features$Accession=="Accession" & !is.na(GEX_Features$SYMBOL),])
GEX_Features$Accession[GEX_Features$Accession=="Accession" & !is.na(GEX_Features$SYMBOL)] <- "NotAccession"



table(
  is.na(GEX_Features$SYMBOL), 
  GEX_Features$Accession
) 


View(GEX_Features[GEX_Features$Accession=="NotAccession" & is.na(GEX_Features$SYMBOL),])



    ## 2.5 Annotate GENETYPEs --------------------------------------------------

GENETYPES <- AnnotationDbi::select(
  org.Hs.eg.db, 
  keys=GEX_Features$SYMBOL, 
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

GEX_Features$GENETYPE <- NA 
GEX_Features$GENETYPE <- GENETYPES$GENETYPE[match(GEX_Features$SYMBOL, GENETYPES$SYMBOL)]

table(is.na(GEX_Features$GENETYPE), GEX_Features$Accession)
GEX_Features$GENETYPE[grep("MT-", GEX_Features$ID_10X)] <- "mitochondrial-encoded"
View(GEX_Features[GEX_Features$Accession!="Accession" & is.na(GEX_Features$GENETYPE),])


GEX_Features$GENETYPE2 <- GEX_Features$GENETYPE
GEX_Features$GENETYPE2[GEX_Features$Accession=="Accession"] <- "Accession"
GEX_Features$GENETYPE2[GEX_Features$Accession!="Accession" & is.na(GEX_Features$GENETYPE2)] <- "unknown"

tmp2 <- GEX_Features$GENETYPE
GEX_Features$GENETYPE <- GEX_Features$GENETYPE2
GEX_Features$GENETYPE2 <- tmp2
rm(tmp2)




  ### 3.0 Annotate ENSEMBL IDs -------------------------------------------------

ENSEMBLS <- AnnotationDbi::select(
  org.Hs.eg.db, 
  keys=GEX_Features$SYMBOL, 
  keytype="SYMBOL", 
  columns="ENSEMBL"
)


ENSEMBLS <- unique(ENSEMBLS)


GEX_Features_ENSEMBL_ENTREZ <- full_join(
  GEX_Features, 
  ENSEMBLS, 
  by=c("SYMBOL"="SYMBOL"), 
  relationship = "many-to-many"
)




  ### 4.0 Annotate ENSEMBL IDs and ENTREZIDs -----------------------------------

ENTREZIDS <- AnnotationDbi::select(
  org.Hs.eg.db, 
  keys=GEX_Features_ENSEMBL_ENTREZ$ENSEMBL, 
  keytype="ENSEMBL", 
  columns="ENTREZID"
)


ENTREZIDS <- unique(ENTREZIDS)
GEX_Features_ENSEMBL_ENTREZ <- full_join(
  GEX_Features_ENSEMBL_ENTREZ, 
  ENTREZIDS, 
  by=c("ENSEMBL"="ENSEMBL")
)




  ### 4.0 Save data ------------------------------------------------------------

qsave(
  GEX_Features, 
  paste0(
    "../Data/Annotations/Genes/", 
    "GEX_Features", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  GEX_Features_ENSEMBL_ENTREZ, 
  paste0(
    "../Data/Annotations/Genes/", 
    "GEX_Features_ENSEMBL_ENTREZ", 
    ".qrds"
  ), 
  nthr=nthr
)
