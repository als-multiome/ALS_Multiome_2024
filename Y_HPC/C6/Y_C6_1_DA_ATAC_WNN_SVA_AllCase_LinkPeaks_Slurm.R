############# Set Library Paths ################################################

.libPaths(c("/home/genomics/R/x86_64-pc-linux-gnu-library/4.4/", .libPaths()))
library(qs) 
library(DESeq2) 
library(BiocParallel) 



############# Read Data ########################################################

qreadm("./C6_1_DA_ATAC_WNN_SVA_AllCase_LinkPeaks_Pre.qrdata",nthr=40)


############# Register Parallelization Plan ####################################

register(MulticoreParam(40)) 


############# Main DE Loop #####################################################


for (i in 1:length(DDS_SVA_12sv_list)){
  
  ######### Loop into new celltype 
  
  message("    ")
  message("#################################################################")
  
  print(Sys.time())
  message(paste0("Iteration #", i, "... \n"))
  message(paste0(
    "Working on ", 
    DESeq_SVA_12sv_Results_Index$CellTypeLevel[i], 
    " with cell type ", 
    DESeq_SVA_12sv_Results_Index$CellType[i], 
    "...")
  )
  message("#################################################################")  
  message("    ")
  
  
  if(i%%2==1){
    message(paste0("Comparing HC, ALS and ALS-FTD and asssigning the results under index nr. ", i, "... ")) 
  }else{
    message(paste0("Comparing Rand1, Rand2 and Rand3, asssigning the results under index nr. ", i, "... ")) 
  }
  
  tryCatch({
    DDS_SVA_12sv_list[[i]] <- DESeq(DDS_SVA_12sv_list[[i]], test="Wald", quiet = FALSE,  parallel=TRUE) 
    message(
      paste0(
        "Succesfully DESeqed Data! \n Design: ",
        design(DDS_SVA_12sv_list[[i]])
      )
    )
  }, error=function(e){
    tryCatch({
      DDS_SVA_12sv_list[[i]] <<- "NA"
      DESeq_SVA_12sv_Results_Index$Remarks[i] <<- "Results could not be calculated"
      message(
        paste0(
          "Could not DESeq Data... \nSetting DDS to 'NA'..."
        )
      )
    }, error=function(e){message("TOTAL ERROR!")})
  })
  
}  
  
message(Sys.time())
message("Saving data....")
do.call(qsavem,c(lapply(ls(), as.symbol), file="./C6_1_DA_ATAC_WNN_SVA_AllCase_LinkPeaks_Post.qrdata", nthr=40))
message(Sys.time())
print("Job done!...")

