# A4_Generate_Single_Seurat_Objects

  ### 0.0 Load libraries ------------------------------------------------------- 

library(Seurat)
library(data.table)
library(qs2)
library(tidyverse)

  
  
  
  ### 1.0 Load data ------------------------------------------------------------
  
Samples <- qs_read(
  "../Data/FANS/Samples/Samples.qs2"
)

Data_Dirs <- read.table("../Data/cfg/Files_list.txt", header = FALSE)
Data_Dirs$V2 <- str_replace_all(Data_Dirs$V2, "\\$HOME", "~")
FANS_Cellranger_Output_Dir <- Data_Dirs$V2[Data_Dirs$V1=="FANS_Cellranger_Output_Dir"]




  ### 2.0 Generate Seurat objects ----------------------------------------------



    # 2.1 Filter-out Sample FC5 due to low QC 

Samples <- Samples %>% 
  filter(ID != "FC5")

FANS_Seurats <- lapply(
  setNames(Samples$ID, nm = Samples$ID),  
  function(Sample){
    
    message(
      paste0(
        "Sample ", 
        Sample, 
        ": ..."
      )
    )
    
    filtered_data <- Read10X(
      file.path(
        FANS_Cellranger_Output_Dir, 
        sprintf("%s/outs/filtered_feature_bc_matrix/", Samples$Name[Samples$ID==Sample])
      ) 
    )  
    
    corrected_all <- Read10X_h5(
      paste0(
        "../Data/FANS/CellBender_Output/", 
        Sample, 
        "/", 
        Sample, 
        "_corrected_output_Seurat.h5"
      )
    )
    
    corrected <- corrected_all[
      ,which(colnames(corrected_all) %in% colnames(filtered_data))
    ]
    rm(corrected_all)
    
    if(all(colnames(corrected) == colnames(filtered_data))){
      message("Corrected and uncorrected data sorted concordantly!")
    } else {
      message("Corrected and uncorrected data NOT sorted concordantly! Sorting corrected data... ")
      corrected <- corrected[,match(colnames(filtered_data), colnames(corrected))]
    }
    
    CellBender_Cells <- fread(
      paste0(
        "../Data/FANS/CellBender_Output/",  
        Sample, 
        "/", 
        Sample, 
        "_corrected_output_cell_barcodes.csv"
      ), 
      data.table = FALSE, 
      header = FALSE
    )
    
    Seurat <- CreateSeuratObject(
      counts = corrected, 
      assay = "RNA", 
      project = Sample, 
      min.cells = 0, 
      min.features = 0
    )
    
    Seurat[["RNA_Uncorrected"]] <- CreateAssayObject(
      counts = filtered_data, 
      min.cells = 0, 
      min.features = 0
    ) 
    
    if(all(rownames(Seurat@assays$RNA) == rownames(Seurat@assays$RNA_Uncorrected)) & 
      all(colnames(Seurat@assays$RNA) == colnames(Seurat@assays$RNA_Uncorrected))){
      message("'RNA' and 'RNA_Uncorrected' Assays concordant! ")
    }else{
      warning("Warning! 'RNA' and 'RNA_Uncorrected' Assays not concordant! ")
    }
    
    Seurat$CellBender_IsCell <- rownames(Seurat@meta.data) %in% CellBender_Cells$V1 
    message(
      paste0(
        "Seurat: ", 
        dim(Seurat)[2], 
        " cells, ", 
        table(Seurat$CellBender_IsCell)["TRUE"], 
        "(", 
        round(table(Seurat$CellBender_IsCell)["TRUE"]*100/dim(Seurat)[2],0), 
        "%) called by CellBender too... "
      )
    )
    
    message(
      paste0(
        "CellBender called ", 
        if_else(
          (nrow(CellBender_Cells) - table(Seurat$CellBender_IsCell)["TRUE"]) <= 0, 
          0, 
          nrow(CellBender_Cells) - table(Seurat$CellBender_IsCell)["TRUE"]
        ), 
        " cells more... "
      )
    )
    
    message("") 
    
    return(Seurat)
    
  }
)




  ### 3.0 Calculate CellBender Ambient RNA Pct ---------------------------------

CellBender_Pct_Ambient <- sapply(
  FANS_Seurats, 
  function(x){
    (1-(sum(x@assays$RNA$counts)/sum(x@assays$RNA_Uncorrected$counts)))*100
  }
)
all(
  names(CellBender_Pct_Ambient) == 
    Samples$ID
)

Samples$CellBender_Pct_Ambient <- CellBender_Pct_Ambient




  ### 4.0 Add number of called cells to Samples data ---------------------------

all(names(FANS_Seurats) == Samples$ID) 
nCells <- sapply(
  FANS_Seurats, 
  function(x){
    return(
      dim(x)[2]
    )
  }
)
Samples$nCells_CR <- nCells[match(Samples$ID, names(nCells))]




  ### 4.0 Save data ------------------------------------------------------------

qs_save(
  FANS_Seurats, 
  "../Data/FANS/SeuratObjects/FANS_Single_SeuratObjects_Raw.qs2"
)

qs_save(
  Samples, 
  "../Data/FANS/Samples/Samples.qs2"
)







