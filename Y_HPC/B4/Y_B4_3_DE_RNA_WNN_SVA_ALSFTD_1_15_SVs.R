############# Set Library Paths ################################################

.libPaths(c("/home/genomics/R/x86_64-pc-linux-gnu-library/4.4/", .libPaths()))
library(qs) 
library(DESeq2) 
library(BiocParallel) 

qs::set_trust_promises(TRUE) 
qs::set_trust_promises(TRUE) 
qs::set_trust_promises(TRUE)


for (nSVs in 1:15){
  
  message(paste0("###### Working with ", nSVs, " SVs.... ######### "))
  message(Sys.time())
  
  
  ############# Read Data ########################################################
  
  DDS_list <- qread(paste0("./SVA_ALSFTD/DDS_", nSVs, "SVs_list.qrds"),nthr=40)
  
  
  nSVs=nSVs 
  DDS_SVA_list <- list() 
  for (i in 1:length(DDS_list)){
    DDS_SVA_list[[i]] <- NA 
  }
  rm(i)
  
  ############# Register Parallelization Plan ####################################
  
  register(MulticoreParam(40)) 
  
  
  ############# Main DE Loop #####################################################
  
  
  for (i in 1:length(DDS_list)){
    
    ######### Loop into new celltype 
    
    message("    ")
    message("#################################################################")
    
    print(Sys.time())
    message(paste0("Iteration #", i, "... \n"))
    
    message("#################################################################")  
    message("    ")
    
    
    if(i%%2==1){
      message(paste0("Comparing HC, ALS-FTD and assigning the results under index nr. ", i, "... ")) 
    }else{
      message(paste0("Comparing Rand1, Rand2 and assigning the results under index nr. ", i, "... ")) 
    }
    
    tryCatch({
      DDS_SVA_list[[i]] <- DESeq(DDS_list[[i]], test="Wald", quiet = FALSE,  parallel=TRUE) 
      message(
        paste0(
          "Succesfully DESeqed Data! \n Design: ",
          design(DDS_SVA_list[[i]])
        )
      )
    }, error=function(e){
      tryCatch({
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
  qsave(
    DDS_SVA_list, 
    file=paste0(
      "./SVA_ALSFTD/DDS_", 
      nSVs, 
      "SVs_list_DESeqed", 
      ".qrds"
    ), 
    nthr=40
  )
  
  message(Sys.time())
  
  tryCatch({rm(DDS_list, DDS_SVA_list)}, error=function(e){})
}

print("Job done!...")

