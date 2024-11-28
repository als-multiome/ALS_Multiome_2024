############# Set Library Paths ################################################

.libPaths(c("/home/genomics/R/x86_64-pc-linux-gnu-library/4.4/", .libPaths()))
library(qs) 
library(DESeq2) 
library(BiocParallel) 

qs::set_trust_promises(TRUE) 
qs::set_trust_promises(TRUE) 
qs::set_trust_promises(TRUE)



############# Read Data ########################################################

qreadm("./B5_3_DE_12SVA_ALSFTD_Pre.qrdata",nthr=40)



############# Register Parallelization Plan ####################################

register(MulticoreParam(40)) 


############# Main DE Loop #####################################################

DDS_SVA_12sv_list <- list() 
for (i in 1:length(DDS_list)){
  DDS_SVA_12sv_list[[i]] <- NA
}

DESeq_SVA_12sv_Results_Index <- DESeq_Results_Index[,c(1:9)] 
rm(DESeq_Results_Index)

message(
  paste0(
    "Length of DDS_SVA_12sv_list same as Nrow DESeq_SVA_12sv_Results_Index: ", 
    length(DDS_SVA_12sv_list)==nrow(DESeq_SVA_12sv_Results_Index)
  )
)

DESeq_SVA_12sv_Results_Index$SVA <- NA  

for (i in 1:nrow(DESeq_SVA_12sv_Results_Index)){ 
  tryCatch({
    DDS.tmp <- DDS_list[[i]] 
    
    if(DESeq_SVA_12sv_Results_Index$Comparison[i]=="ALSFTDvsHC"){
      design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + Case
    }else{
      design(DDS.tmp) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + SV9 + SV10 + SV11 + SV12 + Rand 
    }
    message(paste0("Analysing DDS object #", i, "; Design: ", design(DDS.tmp)))
    DDS_SVA_12sv_list[[i]] <- DESeq(DDS.tmp, test="Wald", quiet=FALSE, parallel=TRUE) 
    DESeq_SVA_12sv_Results_Index$SVA[i] <- TRUE 
    message("########################################################")
    message(Sys.time())
    message(paste0("DESeq analysis #", i, " sucessfull!...")) 
    message("########################################################")
  }, error=function(e){
    DESeq_SVA_12sv_Results_Index$SVA[i] <<- FALSE 
    message("########################################################")
    message(Sys.time()) 
    message(paste0("DESeq analysis #", i, " failed :( "))
    message("########################################################")
  })
  tryCatch({rm(DDS.tmp)}, error=function(e){})
}
  
  message(Sys.time())
  message("Saving data....")
  
  do.call(qsavem,c(lapply(ls(), as.symbol), file="./B5_3_DE_12SVA_ALSFTD_Post.qrdata ", nthr=40))
  print("Job done!...")
  
print("Job done!...")

