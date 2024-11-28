
  ### 0.0 Load libraries -------------------------------------------------------

source("~/ALS_Brain_Multiome.Rcfg")

library(data.table)
library(tidyverse)
library(GenomicRanges)
library(qs)




  ### 1.0 Load data ------------------------------------------------------------

APA_res <- fread(
  paste0(
    "../Data/APA/DaPars1_TDP43Pos_vs_TDP43Neg/", 
    "DaPars1_TDP43Pos_vs_TDP43Neg_All_Prediction_Results", 
    ".txt"
  ), 
  data.table=FALSE
)




  ### 2.0 Format features ------------------------------------------------------

APA_res$Transcript <- str_split(APA_res$Gene, "\\|", simplify=TRUE)[,1]
APA_res$Symbol <- str_split(APA_res$Gene, "\\|", simplify=TRUE)[,2]
APA_res$Chr <- str_split(APA_res$Gene, "\\|", simplify=TRUE)[,3]
APA_res$Strand <- str_split(APA_res$Gene, "\\|", simplify=TRUE)[,4]

table(is.na(APA_res$Strand))
table(APA_res$Strand)

APA_res$Start <- str_split(str_split(APA_res$Loci, ":", simplify=TRUE)[,2], "-", simplify=TRUE)[,1]
APA_res$End <- str_split(str_split(APA_res$Loci, ":", simplify=TRUE)[,2], "-", simplify=TRUE)[,2]

table(is.na(APA_res$Start))
table(is.na(APA_res$End))

table(APA_res$Start==APA_res$End)
table(APA_res$End > APA_res$Start)

APA.tmp <- data.frame(
  Gene=APA_res$Gene, 
  PDUI=APA_res$PDUI_Group_diff, 
  Region1=paste0(
    APA_res$Start, 
    "-", 
    APA_res$Predicted_Proximal_APA
  ), 
  Region2=paste0(
    APA_res$Predicted_Proximal_APA+1, 
    "-", 
    APA_res$End
  ), 
  Strand=APA_res$Strand, 
  Short=NA, 
  Long=NA
)

APA.tmp$Short[APA.tmp$Strand=="+"] <- "Region1"
APA.tmp$Short[APA.tmp$Strand=="-"] <- "Region2"

APA.tmp$Long[APA.tmp$Strand=="+"] <- "Region2"
APA.tmp$Long[APA.tmp$Strand=="-"] <- "Region1"

table(APA.tmp$Short)
table(APA.tmp$Long)
table(APA.tmp$Short==APA.tmp$Long)


APA_res$Region1 <- APA.tmp$Region1[match(APA_res$Gene, APA.tmp$Gene)]
APA_res$Region2 <- APA.tmp$Region2[match(APA_res$Gene, APA.tmp$Gene)]
APA_res$Short <- APA.tmp$Short[match(APA_res$Gene, APA.tmp$Gene)]
APA_res$Long <- APA.tmp$Long[match(APA_res$Gene, APA.tmp$Gene)]

rm(APA.tmp)

APA_res$TDP_Neg <- "NA"
APA_res$TDP_Neg[APA_res$PDUI_Group_diff == 0] <- "None"
APA_res$TDP_Neg[APA_res$PDUI_Group_diff < 0] <- "Long"
APA_res$TDP_Neg[APA_res$PDUI_Group_diff > 0] <- "Short"
table(APA_res$TDP_Neg)

APA_res$Region1_Class <- "NA"
APA_res$Region2_Class <- "NA"

APA_res$Region1_Class[(APA_res$Short=="Region1" & APA_res$TDP_Neg=="Short") | (APA_res$Long=="Region1" & APA_res$TDP_Neg=="Long")] <- "TDP43_Neg"
APA_res$Region1_Class[(APA_res$Short=="Region1" & APA_res$TDP_Neg=="Long") | (APA_res$Long=="Region1" & APA_res$TDP_Neg=="Short")] <- "TDP43_Pos"

table(APA_res$Region1_Class)


APA_res$Region2_Class[(APA_res$Short=="Region2" & APA_res$TDP_Neg=="Short") | (APA_res$Long=="Region2" & APA_res$TDP_Neg=="Long")] <- "TDP43_Neg"
APA_res$Region2_Class[(APA_res$Short=="Region2" & APA_res$TDP_Neg=="Long") | (APA_res$Long=="Region2" & APA_res$TDP_Neg=="Short")] <- "TDP43_Pos"

table(APA_res$Region2_Class)

table(APA_res$Region1_Class==APA_res$Region2_Class)
table(APA_res$Region1_Class[!is.na(APA_res$P_val) & APA_res$PDUI_Group_diff!=0]==APA_res$Region2_Class[!is.na(APA_res$P_val) & APA_res$PDUI_Group_diff!=0])




  ### 3.0 Filter significant hits ----------------------------------------------

APA_Sign <- APA_res[!is.na(APA_res$adjusted.P_val),]
APA_Sign <- APA_Sign[APA_Sign$adjusted.P_val<0.01,]
APA_Sign <- APA_Sign[abs(APA_Sign$PDUI_Group_diff)>=0.2,]




  ### 4.0 Generate a feature table ---------------------------------------------

APA_TDP43_Features <- data.frame(
  Gene=c(APA_Sign$Symbol, APA_Sign$Symbol),
  Strand=c(APA_Sign$Strand, APA_Sign$Strand),
  Chr=c(APA_Sign$Chr, APA_Sign$Chr), 
  Region=c(APA_Sign$Region1, APA_Sign$Region2), 
  Class=c(APA_Sign$Region1_Class, APA_Sign$Region2_Class)
)

APA_TDP43_Features <- APA_TDP43_Features[order(APA_TDP43_Features$Gene),]
APA_TDP43_Features$GeneID <- paste0(
  APA_TDP43_Features$Gene,
  "_",
  APA_TDP43_Features$Class
)
APA_TDP43_Features$Start <- as.numeric(str_split(APA_TDP43_Features$Region, "-", simplify = TRUE)[,1])
APA_TDP43_Features$End <- as.numeric(str_split(APA_TDP43_Features$Region, "-", simplify = TRUE)[,2])
table(APA_TDP43_Features$End > APA_TDP43_Features$Start)




  ### 5.0 Export data ----------------------------------------------------------

qsave(
  APA_res, 
  paste0(
    "../Data/APA/DaPars1_TDP43Pos_vs_TDP43Neg/", 
    "APA_Results", 
    ".qrds"
  ), 
  nthr=nthr
)

write.csv(
  APA_Sign, 
  paste0(
    "../Data/APA/DaPars1_TDP43Pos_vs_TDP43Neg/", 
    "APA_Sign_Results", 
    ".csv"
  ), 
  quote=FALSE
)

write.table(
  APA_TDP43_Features[c("GeneID", "Chr", "Start", "End", "Strand")],  
  paste0(
    "../Data/APA/DaPars1_TDP43Pos_vs_TDP43Neg/", 
    "APA_TDP43Neg_Features", 
    ".saf"
  ), 
  sep="\t",
  quote=FALSE, 
  row.names = FALSE
)





