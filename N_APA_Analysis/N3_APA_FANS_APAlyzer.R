

  ### 0.0 Load libraries -------------------------------------------------------

source("~/ALS_Brain_Multiome.Rcfg")

library(data.table)
library(tidyverse)
library(GenomicRanges)
library(qs)
library(APAlyzer)
library(Rsamtools)
library(repmis)




  ### 1.0 List Bam files -------------------------------------------------------

flsall <- c(
  "TDP_Pos"="../Data/FANS/Bam/FANS_TDP43_Pos.bam",
  "TDP_Neg"="../Data/FANS/Bam/FANS_TDP43_Neg.bam"
)
flsall




  ### 2.0 Get genomic reference ------------------------------------------------

URL="https://github.com/RJWANGbioinfo/PAS_reference_RData/blob/master/"
file="hg38_REF.RData"
source_data(paste0(URL,file,"?raw=True"))




  ### 3.0 Build PAS Reference Regions ------------------------------------------

PASREF=REF4PAS(refUTRraw_hg38,dfIPA_hg38,dfLE_hg38)
UTRdbraw=PASREF$UTRdbraw
dfIPA=PASREF$dfIPA
dfLE=PASREF$dfLE    




  ### 4.0 Build aUTR and cUTR References ---------------------------------------

UTRdbraw=REF3UTR(refUTRraw_hg38)
DFUTRraw=PASEXP_3UTR(UTRdbraw, flsall, Strandtype="forward")

head(DFUTRraw,2)

sampleTable2 = data.frame(samplename = c("TDP_Pos","TDP_Neg"),
                          condition = c("Phys","Path")) 
sampleTable2

test_3UTRsing=APAdiff(sampleTable2,DFUTRraw, 
                      conKET='Phys',
                      trtKEY='Path',
                      PAS='3UTR',
                      CUTreads=0,
                      p_adjust_methods="fdr")
qsave(
  DFUTRraw, 
  paste0(
    "../Data/APA/APAlyzer_TFP43Neg_vs_TDP43Pos/", 
    "DFUTRraw_TDPNeg_vs_TDPPos", 
    ".qrds"
  ), 
  nthr=nthr
)

qsave(
  test_3UTRsing, 
  paste0(
    "../Data/APA/APAlyzer_TFP43Neg_vs_TDP43Pos/", 
    "Res_3UTR_TDPNeg_vs_TDPPos", 
    ".qrds"
  ), 
  nthr=nthr
)



IPA_OUTraw=PASEXP_IPA(dfIPA, dfLE, flsall, Strandtype="forward", nts=24) 

head(IPA_OUTraw,2)

test_IPAsing=APAdiff(sampleTable2,
                     IPA_OUTraw, 
                     conKET='Phys',
                     trtKEY='Path',
                     PAS='IPA',
                     CUTreads=0,
                     p_adjust_methods="fdr") 
head(test_IPAsing,2)

qsave(
  test_IPAsing, 
  paste0(
    "../Data/APA/APAlyzer_TFP43Neg_vs_TDP43Pos/", 
    "Res_IPA_TDPNeg_vs_TDPPos", 
    ".qrds"
  ), 
  nthr=nthr
)


qsave(
  IPA_OUTraw, 
  paste0(
    "../Data/APA/APAlyzer_TFP43Neg_vs_TDP43Pos/", 
    "IPA_OUTraw_TDPNeg_vs_TDPPos", 
    ".qrds"
  ), 
  nthr=nthr
)



