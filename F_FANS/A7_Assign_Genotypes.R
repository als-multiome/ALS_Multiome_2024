# A6_Assign_Genotypes.R




  ### 0.0 Load libraries -------------------------------------------------------

library(qs2)
library(Seurat)
library(tidyverse) 
library(data.table)



  
  ### 1.0 Load data ------------------------------------------------------------

Samples <- qs_read(
    "../Data/FANS/Samples/Samples.qs2"
)  
  
FANS <- qs_read(
  "../Data/FANS/SeuratObjects/FANS_Single_SeuratObjects_Filtered_scDblFinder_Marked.qs2"
)

Data_Dirs <- read.table("../Data/cfg/Files_list.txt", header = FALSE)
Data_Dirs$V2 <- str_replace_all(Data_Dirs$V2, "\\$HOME", "~")
FANS_Souporcell_Output_Dir <- Data_Dirs$V2[Data_Dirs$V1=="FANS_Souporcell_Output_Dir"]
Souporcell_Output_Dirs <- list.files(FANS_Souporcell_Output_Dir)

Souporcell_Clusters <- lapply(
  setNames(
    Souporcell_Output_Dirs, 
    nm = Souporcell_Output_Dirs
  ), 
  function(x){
    return(
      fread(
        paste0(
          FANS_Souporcell_Output_Dir, 
          "/", 
          x, 
          "/clusters.tsv"
        ), 
        data.table=FALSE
      )
    )
  }
) 

Souporcell_Clusters <- Souporcell_Clusters[
  match(
    names(FANS), 
    names(Souporcell_Clusters))[!is.na(match(names(FANS), names(Souporcell_Clusters)))]
] 

Souporcell_AmbientRNA <- lapply(
  setNames(
    Souporcell_Output_Dirs, 
    nm = Souporcell_Output_Dirs
  ), 
  function(x){
    return(
      fread(
        paste0(
          FANS_Souporcell_Output_Dir, 
          "/", 
          x, 
          "/ambient_rna.txt"
        ), 
        data.table=FALSE
      )
    )
  }
) 
Souporcell_AmbientRNA <- sapply(
  Souporcell_AmbientRNA, 
  function(x){
    return(
      str_replace_all(
        colnames(x)[5], 
        "%", 
        ""
      ) %>% as.numeric()
    )
  }
)




  ### 2.0 Add Souporcell Data --------------------------------------------------

sapply(
  setNames(
    names(Souporcell_Clusters), 
    nm = names(Souporcell_Clusters)
  ), 
  function(x){
    message(x)
    tmp <- Souporcell_Clusters[[x]]
    all_found <- all(
      rownames(FANS[[x]]@meta.data) %in% 
        tmp$barcode 
    )
    FANS[[x]]$Souporcell_Cluster <<- tmp$assignment[match(rownames(FANS[[x]]@meta.data), tmp$barcode)] 
    FANS[[x]]$Souporcell_Status <<- tmp$status[match(rownames(FANS[[x]]@meta.data), tmp$barcode)] 
    FANS[[x]]$DoubleCheck.tmp <<- tmp$barcode[match(rownames(FANS[[x]]@meta.data), tmp$barcode)]  
    return(all_found)
  }
) |> all()

      # Check correct assignment 
sapply(
  FANS, 
  function(x){
    if(any(colnames(x@meta.data)=="DoubleCheck.tmp")){
      return(all(rownames(x@meta.data) == x$DoubleCheck.tmp))
    } else{return(TRUE)}
  }
) |> all()

FANS <- lapply(
  FANS, 
  function(x){
    if("DoubleCheck.tmp" %in% colnames(x@meta.data)){
      x@meta.data <- x@meta.data %>% 
        select(-DoubleCheck.tmp)  
    } else {
      x$Souporcell_Cluster <- NA
      x$Souporcell_Status <- NA
    }
    
    return(
      x
    )
  }
)

Samples$Souporcell_PctAmbientRNA <- Souporcell_AmbientRNA[match(Samples$ID, names(Souporcell_AmbientRNA))]
Samples$Souporcell_PctAmbientRNA[Samples$Samples<2] <- NA 




  ### 3.0 Export data ----------------------------------------------------------

qs_save(
  FANS, 
  "../Data/FANS/SeuratObjects/FANS_Single_SeuratObjects_Filtered_scDblFinder_Marked_Genotypes_Assigned.qs2"
)

qs_save(
  Samples, 
  "../Data/FANS/Samples/Samples.qs2"
)



