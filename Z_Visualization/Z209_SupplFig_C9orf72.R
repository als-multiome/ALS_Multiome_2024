

  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 

library(qs)
library(tidyverse)
library(DESeq2)
library(sva)
library(ggrepel)



qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)




  ### 1.0 Load data ------------------------------------------------------------



    ## 1.1 RNA -----------------------------------------------------------------

M0_RNA <- qread(
  paste0(
    "../Data/SeuratObjects/", 
    "M0_RNA", 
    ".qrds"
  ), 
  nthr = nthr
)

DDS_list_RNA <- list() 

DDS_list_RNA[["NoCov"]] <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9_NonC9/", 
    "DDS_DESEqed_list", 
    ".qrds"
  ), 
  nthr = nthr
)

for(i in 1:14){
  
  DDS_list_RNA[[paste0("DDS_", i, "_SVs")]] <- qread(
    paste0(
      "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9_NonC9/SVA/1_14_SVs/", 
      paste0(
        "DDS_", i, "SVs_list_DESeqed", 
        ".qrds"
      )
    ), 
    nthr = nthr
  ) 
  rm(i)
}

Index_NoCov <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9_NonC9/", 
    "DESeq_Results_Index", 
    ".qrds"
  ), 
  nthr = nthr
)

Index_SVA_RNA <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/ALSFTD_C9_NonC9/SVA/1_14_SVs/", 
    "SVAs_14SVs_Index", 
    ".qrds"
  ), 
  nthr = nthr
)

all(
  Index_NoCov$CellTypeLevel[1:76] == 
    Index_SVA_RNA$CellTypeLevel
)

all(
  Index_NoCov$CellType[1:76] == 
    Index_SVA_RNA$CellType
) 

all(
  Index_NoCov$Comparison[1:76] == 
    Index_SVA_RNA$Comparison
)  

rm(Index_NoCov)
DDS_list_RNA[["NoCov"]] <- DDS_list_RNA[["NoCov"]][1:76]

dds_All_Cases_RNA <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/", 
    "DDS_list_RNA", 
    ".qrds"
  ), 
  nthr = nthr
)

Index_dds_All_Cases_RNA <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/", 
    "DESeq_Results_Index", 
    ".qrds"
  ), 
  nthr = nthr
)



    ## 1.2 ATAC ----------------------------------------------------------------

DDS_list_ATAC.tmp <- qread(
  paste0(
    "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/", 
    "DDS_DESeqed_list", 
    ".qrds"
  ), 
  nthr = nthr
) 

Index_ATAC.tmp <- qread(
  paste0(
    "../Data/DA/WNN/Pseudobulk/DESeq2/AllCase_LinkPeaks/", 
    "DESeq_Results_Index", 
    ".qrds"
  ), 
  nthr = nthr
) 



    ## 1.3 Compositional data --------------------------------------------------

MASC_C9_ALSFTD_WNN_L25 <- qread(
  paste0(
    "../Data/Coda/MASC/", 
    "MASC_C9_ALSFTD_WNN_L25", 
    ".qrds"
  ), 
  nthr = nthr
)

  

    ## 1.4 Color scales --------------------------------------------------------

ColDict_Case_Type <- setNames(
  object = xlsx::read.xlsx(
    "../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", 
    sheetName = "Case_Type"
  )[["Color"]], 
  nm = xlsx::read.xlsx(
    "../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", 
    sheetName = "Case_Type"
  )[["Case_Type"]]
) 




  ### 2.0 Define AUX functions -------------------------------------------------

get_N_DEGs <- function(x, alpha = 0.05, name){
  res <- as.data.frame(
    results(x, alpha = alpha, name = name) 
  )
  res <- res[!is.na(res$padj), ]
  return(
    sum(
      res$padj < alpha
    )
  )
} 


get_N_DEGs2 <- function(x, alpha = 0.05, name, direction=NULL){
  res <- as.data.frame(
    results(x, alpha = alpha, name = name) 
  )
  res <- res[!is.na(res$padj), ]
  stopifnot(!is.null(direction))
  if(direction == "up"){
    return(
      sum(
        res$padj < alpha & res$log2FoldChange > 0
      )
    )
  }
  if(direction == "down"){
    return(
      sum(
        res$padj < alpha & res$log2FoldChange < 0
      )
    )
  }
}



only_value <- function(x){
  if(is.null(x)) stop("Only_Value: Please define variable!")
  if(length(unique(x))==1){
    return(unique(x)[1])
  } else {
    stop("Only_Value: Unequivocal variable!")
  }
}




  ### 3.0 Compare DGEs ---------------------------------------------------------



    ## 3.1 All Cells -----------------------------------------------------------


      # 3.1.1 Double-check results names ---------------------------------------

ind <- which(
  Index_SVA_RNA$CellTypeLevel == "AllCells"
)

for (i in ind){

  res.name = tail(
    resultsNames(DDS_list_RNA[[1]][[i]]), 
    n = 1
  )
  
  print(res.name)
  print(
    all(
      sapply(
        DDS_list_RNA, 
        FUN = function(x){
          return(
            tail(
              resultsNames(x[[i]]),
              n = 1
            ) == res.name
          )
        }
      )
    )
  )
  
  
  rm(res.name, i)
  
}
rm(ind)



      # 3.1.2 Case w/o Sex as covariate ----------------------------------------

Case_woSex <- sapply(
  DDS_list_RNA, 
  FUN = function(x){
    get_N_DEGs(x[[1]], alpha = 0.05, name = "Case_Type_ALS_FTD_vs_C9_ALS_FTD")
  }
)


      # 3.1.3 Case with Sex as covariate ---------------------------------------

Case_withSex <- c() 
  
for (i in 1:length(DDS_list_RNA)){
  y = DDS_list_RNA[[i]][[1]]
  design(y) <- as.formula(
    paste0(
      as.character(design(y))[1], 
      " Sex + ", 
      as.character(design(y))[2]
    ) 
  ) 
  message(
    paste0(
      "Design: ", 
      design(y)
    )
  )
  y <- DESeq(y, test = "Wald", parallel = TRUE) 
  Case_withSex[i] <- get_N_DEGs(y) 
  names(Case_withSex)[i] <- names(DDS_list_RNA[i]) 
  rm(i)
}


      # 3.1.4 Rand w/o Sex as covariate -----------------------------------------

Rand_woSex <- sapply(
  DDS_list_RNA, 
  FUN = function(x){
    get_N_DEGs(x[[2]], alpha = 0.05, name = "Rand_R2_vs_R1")
  }
)


      # 3.1.5 C9_ALS_FTD vs. HC ------------------------------------------------


        # Select matching HCs --------------------------------------------------

dds <- dds_All_Cases_RNA[[1]]
sample_data <- colData(
  dds
) %>% as.data.frame 

table(
  sample_data$Case_Type
) 

table(
  sample_data$Sex[sample_data$Case_Type == "ALS_FTD"]
)

table(
  sample_data$Case_Type[sample_data$Case == "ALS_FTD"], 
  sample_data$Sex[sample_data$Case == "ALS_FTD"]
)

HCs_for_C9_ALSFTD <- c(
  "HC1", "HC2", "HC10", "HC12", "HC13", "HC15", "HC16", 
  "HC11", "HC14", "HC17"
) 

counts <- counts(dds, normalized = FALSE) 
coldata <- colData(dds)
all(
  colnames(counts) == rownames(coldata)
)

counts <- counts[
  , colnames(counts) %in% c(
    coldata$ID[coldata$Case_Type=="C9_ALS_FTD"], 
    HCs_for_C9_ALSFTD
  )
]

coldata <- coldata[
  rownames(coldata) %in% c(
    coldata$ID[coldata$Case_Type=="C9_ALS_FTD"], 
    HCs_for_C9_ALSFTD    
  ), 
]
all(
  colnames(counts) == rownames(coldata)
)

dds <- DESeqDataSetFromMatrix(
  countData = counts, 
  colData = coldata, 
  design = ~ Case
)
rm(counts, coldata, HCs_for_C9_ALSFTD)

table(
  dds$Case_Type, 
  dds$Sex
)


        # Add SVA --------------------------------------------------------------

dds <- estimateSizeFactors(dds)
dat  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1 
dat  <- dat[idx, ] 
mod  <- model.matrix(~ Case, colData(dds))
mod0 <- model.matrix(~ 1, colData(dds))

svseq <- svaseq(dat, mod, mod0, n.sv=14)
colnames(svseq$sv) <- paste0("SV", 1:ncol(svseq$sv))
colData(dds) <- cbind(colData(dds), svseq$sv)


DDS_Case_Control_woSex_list <- list()
dds$Sex <- factor(dds$Sex)

design(dds) <- ~ Case
register(MulticoreParam(nthr))  
DDS_Case_Control_woSex_list[["NoCov"]] <- DESeq(
  dds, 
  test = "Wald", 
  parallel = TRUE
)

for (i in 1:ncol(svseq$sv)){
  
  design(dds) <- as.formula(
    paste0(
      "~ ", 
      paste0(
        "SV", 
        1:i, 
        collapse = " + "
      ), 
      " + Case"
    )
  ) 
  
  message(
    design(dds)
  ) 
  
  DDS_Case_Control_woSex_list[[
    paste0(
      "DDS_", 
      i, 
      "SVs"
    )
  ]] <- DESeq(
    dds, 
    test = "Wald", 
    parallel = TRUE
  )
  
  message(
    get_N_DEGs(
      DDS_Case_Control_woSex_list[[
        paste0(
          "DDS_", 
          i, 
          "SVs"
        )
      ]]
    )
  )
  
}



DDS_Case_Control_withSex_list <- list()
dds$Sex <- factor(dds$Sex)

design(dds) <- ~ Sex + Case
register(MulticoreParam(nthr))  
DDS_Case_Control_withSex_list[["NoCov"]] <- DESeq(
  dds, 
  test = "Wald", 
  parallel = TRUE
)

for (i in 1:ncol(svseq$sv)){
  
  design(dds) <- as.formula(
    paste0(
      "~ Sex + ", 
      paste0(
        "SV", 
        1:i, 
        collapse = " + "
      ), 
      " + Case"
    )
  ) 
  
  message(
    design(dds)
  ) 
  
  DDS_Case_Control_withSex_list[[
    paste0(
      "DDS_", 
      i, 
      "SVs"
    )
  ]] <- DESeq(
    dds, 
    test = "Wald", 
    parallel = TRUE
  )
  
  message(
    get_N_DEGs(
      DDS_Case_Control_withSex_list[[
        paste0(
          "DDS_", 
          i, 
          "SVs"
        )
      ]]
    )
  )
  
}


rm(dds, DDS, dat, idx, mod, mod0, svseq, i)


C9_HC_Control_woSex <- sapply(
  DDS_Case_Control_woSex_list, 
  FUN = function(x){
    get_N_DEGs(x, alpha = 0.05, name = "Case_ALS_FTD_vs_HC")
  }
)

C9_HC_Control_withSex <- sapply(
  DDS_Case_Control_withSex_list, 
  FUN = function(x){
    get_N_DEGs(x, alpha = 0.05, name = "Case_ALS_FTD_vs_HC")
  }
)

rm(DDS_Case_Control_withSex_list, DDS_Case_Control_woSex_list, test_cases, sample_data)


      # 3.1.6 NonC9_ALS_FTD vs. HC ------------------------------------------------


        # Select matching HCs --------------------------------------------------

dds <- dds_All_Cases_RNA[[1]]
sample_data <- colData(
  dds
) %>% as.data.frame 

table(
  sample_data$Case_Type
) 

table(
  sample_data$Sex[sample_data$Case_Type == "C9_ALS_FTD"]
)

table(
  sample_data$Case_Type[sample_data$Case == "ALS_FTD"], 
  sample_data$Sex[sample_data$Case == "ALS_FTD"]
)

HCs_for_NonC9_ALSFTD <- c(
  c(
    sample(
      sample_data$ID[sample_data$Case == "HC" & sample_data$Sex == "f"], 
      size = 4, 
      replace = FALSE
    ), 
    sample(
      sample_data$ID[sample_data$Case == "HC" & sample_data$Sex == "m"], 
      size = 3, 
      replace = FALSE
    )
  )
) 

counts <- counts(dds, normalized = FALSE) 
coldata <- colData(dds)
all(
  colnames(counts) == rownames(coldata)
)

counts <- counts[
  , colnames(counts) %in% c(
    coldata$ID[coldata$Case_Type=="ALS_FTD"], 
    HCs_for_NonC9_ALSFTD
  )
]

coldata <- coldata[
  rownames(coldata) %in% c(
    coldata$ID[coldata$Case_Type=="ALS_FTD"], 
    HCs_for_NonC9_ALSFTD    
  ), 
]
all(
  colnames(counts) == rownames(coldata)
)

dds <- DESeqDataSetFromMatrix(
  countData = counts, 
  colData = coldata, 
  design = ~ Case
)
rm(counts, coldata, HCs_for_NonC9_ALSFTD)

table(
  dds$Case_Type, 
  dds$Sex
)


        # Add SVA --------------------------------------------------------------

dds <- estimateSizeFactors(dds)
dat  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1 
dat  <- dat[idx, ] 
mod  <- model.matrix(~ Case, colData(dds))
mod0 <- model.matrix(~ 1, colData(dds))

svseq <- svaseq(dat, mod, mod0, n.sv=14)
colnames(svseq$sv) <- paste0("SV", 1:ncol(svseq$sv))
colData(dds) <- cbind(colData(dds), svseq$sv)


DDS_Case_Control_woSex_list <- list()
dds$Sex <- factor(dds$Sex)

design(dds) <- ~ Case
register(MulticoreParam(nthr))  
DDS_Case_Control_woSex_list[["NoCov"]] <- DESeq(
  dds, 
  test = "Wald", 
  parallel = TRUE
)

for (i in 1:ncol(svseq$sv)){
  
  design(dds) <- as.formula(
    paste0(
      "~ ", 
      paste0(
        "SV", 
        1:i, 
        collapse = " + "
      ), 
      " + Case"
    )
  ) 
  
  message(
    design(dds)
  ) 
  
  DDS_Case_Control_woSex_list[[
    paste0(
      "DDS_", 
      i, 
      "SVs"
    )
  ]] <- DESeq(
    dds, 
    test = "Wald", 
    parallel = TRUE
  )
  
  message(
    get_N_DEGs(
      DDS_Case_Control_woSex_list[[
        paste0(
          "DDS_", 
          i, 
          "SVs"
        )
      ]]
    )
  )
  
}



DDS_Case_Control_withSex_list <- list()
dds$Sex <- factor(dds$Sex)

design(dds) <- ~ Sex + Case
register(MulticoreParam(nthr))  
DDS_Case_Control_withSex_list[["NoCov"]] <- DESeq(
  dds, 
  test = "Wald", 
  parallel = TRUE
)

for (i in 1:ncol(svseq$sv)){
  
  design(dds) <- as.formula(
    paste0(
      "~ Sex + ", 
      paste0(
        "SV", 
        1:i, 
        collapse = " + "
      ), 
      " + Case"
    )
  ) 
  
  message(
    design(dds)
  ) 
  
  DDS_Case_Control_withSex_list[[
    paste0(
      "DDS_", 
      i, 
      "SVs"
    )
  ]] <- DESeq(
    dds, 
    test = "Wald", 
    parallel = TRUE
  )
  
  message(
    get_N_DEGs(
      DDS_Case_Control_withSex_list[[
        paste0(
          "DDS_", 
          i, 
          "SVs"
        )
      ]]
    )
  )
  
}


rm(dds, DDS, dat, idx, mod, mod0, svseq, i)


NonC9_HC_Control_woSex <- sapply(
  DDS_Case_Control_woSex_list, 
  FUN = function(x){
    get_N_DEGs(x, alpha = 0.05, name = "Case_ALS_FTD_vs_HC")
  }
)

NonC9_HC_Control_withSex <- sapply(
  DDS_Case_Control_withSex_list, 
  FUN = function(x){
    get_N_DEGs(x, alpha = 0.05, name = "Case_ALS_FTD_vs_HC")
  }
)

rm(DDS_Case_Control_withSex_list, DDS_Case_Control_woSex_list, test_cases)





df <- data.frame(
  Case = c(rep(c("Case", "C9_Control", "Non_C9_Control"), each = 28), rep("Rand", 14))
)
df$CovCov = rep(0:13, 7)
df$Sex = c(rep(c("With", "Wo", "With", "Wo", "With", "Wo", "With"), each=14)) 
df$DEGs = c(
  Case_withSex[1:14], 
  Case_woSex[1:14], 
  C9_HC_Control_withSex[1:14], 
  C9_HC_Control_woSex[1:14], 
  NonC9_HC_Control_withSex[1:14], 
  NonC9_HC_Control_woSex[1:14], 
  Rand_woSex[1:14]
)
df$Tmp <- paste0(df$Case, "_", df$Sex)
ggplot(df) + 
  aes(CovCov, y = DEGs, col=Tmp) + 
  geom_line() + 
  geom_point(pch=21) 


      # 3.1.7 3/3, 3/3, 3/3 Comparison -----------------------------------------


        # Select matching HCs --------------------------------------------------

dds <- dds_All_Cases_RNA[[1]]
sample_data <- colData(
  dds
) %>% as.data.frame 

table(
  sample_data$Case_Type[sample_data$Case == "ALS_FTD"], 
  sample_data$Sex[sample_data$Case == "ALS_FTD"]
)

set.seed(seed = 42)
HCs_for_3x3 <- sample(
      sample_data$ID[sample_data$Case == "HC" & sample_data$Sex == "f"], 
      size = 3, 
      replace = FALSE
    )
set.seed(seed = 42)
HCs_for_3x3 <- c(
  HCs_for_3x3, 
  sample(
      sample_data$ID[sample_data$Case == "HC" & sample_data$Sex == "m"], 
      size = 3, 
      replace = FALSE
  )
)

set.seed(seed = 42)
C9s_for_3x3 <- sample(
  sample_data$ID[sample_data$Case_Type == "C9_ALS_FTD" & sample_data$Sex == "f"], 
  size = 3, 
  replace = FALSE
) 

set.seed(seed = 42)
C9s_for_3x3 <- c(
  C9s_for_3x3, 
  sample(
    sample_data$ID[sample_data$Case_Type == "C9_ALS_FTD" & sample_data$Sex == "m"], 
    size = 3, 
    replace = FALSE
  ) 
)

set.seed(seed = 42)
NonC9s_for_3x3 <- sample(
  sample_data$ID[sample_data$Case_Type == "ALS_FTD" & sample_data$Sex == "f"], 
  size = 3, 
  replace = FALSE
) 

set.seed(seed = 42)
NonC9s_for_3x3 <- c(
  NonC9s_for_3x3, 
  sample(
    sample_data$ID[sample_data$Case_Type == "ALS_FTD" & sample_data$Sex == "m"], 
    size = 3, 
    replace = FALSE
  ) 
)

any(
  NonC9s_for_3x3 %in% 
    C9s_for_3x3
)

Cases_3x3 <- unique(
  c(
    HCs_for_3x3, 
    NonC9s_for_3x3, 
    C9s_for_3x3
  )
)
rm(HCs_for_3x3, NonC9s_for_3x3, C9s_for_3x3)

counts <- counts(dds, normalized = FALSE) 
coldata <- colData(dds)
all(
  colnames(counts) == rownames(coldata)
)

counts <- counts[
  , colnames(counts) %in% Cases_3x3
]

coldata <- coldata[
  rownames(coldata) %in% Cases_3x3, 
]

all(
  colnames(counts) == rownames(coldata)
)

dds <- DESeqDataSetFromMatrix(
  countData = counts, 
  colData = coldata, 
  design = ~ Case_Type
)
rm(counts, coldata)

table(
  dds$Case_Type, 
  dds$Sex
)


        # Add SVA --------------------------------------------------------------

dds <- estimateSizeFactors(dds)
dat  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1 
dat  <- dat[idx, ] 
mod  <- model.matrix(~ Case_Type, colData(dds))
mod0 <- model.matrix(~ 1, colData(dds))

svseq <- svaseq(dat, mod, mod0, n.sv=14)
colnames(svseq$sv) <- paste0("SV", 1:ncol(svseq$sv))
colData(dds) <- cbind(colData(dds), svseq$sv)


DDS_Case_Control_woSex_list <- list()
dds$Sex <- factor(dds$Sex)

design(dds) <- ~ Case_Type
register(MulticoreParam(nthr))  
DDS_Case_Control_woSex_list[["NoCov"]] <- DESeq(
  dds, 
  test = "Wald", 
  parallel = TRUE
)

for (i in 1:ncol(svseq$sv)){
  
  design(dds) <- as.formula(
    paste0(
      "~ ", 
      paste0(
        "SV", 
        1:i, 
        collapse = " + "
      ), 
      " + Case_Type"
    )
  ) 
  
  message(
    design(dds)
  ) 
  
  DDS_Case_Control_woSex_list[[
    paste0(
      "DDS_", 
      i, 
      "SVs"
    )
  ]] <- DESeq(
    dds, 
    test = "Wald", 
    parallel = TRUE
  )
  
  print("N DEGs C9: ")
  message(
    get_N_DEGs(
      DDS_Case_Control_woSex_list[[
        paste0(
          "DDS_", 
          i, 
          "SVs"
        )
      ]], 
      name = "Case_Type_C9_ALS_FTD_vs_HC"
    )
  )

  print("N DEGs NonC9: ")
  message(
    get_N_DEGs(
      DDS_Case_Control_woSex_list[[
        paste0(
          "DDS_", 
          i, 
          "SVs"
        )
      ]], 
      name = "Case_Type_ALS_FTD_vs_HC"
    )
  )
  
}

C9_3x3_woSex <- sapply(
  DDS_Case_Control_woSex_list, 
  FUN = function(x){
    get_N_DEGs(x, alpha = 0.05, name = "Case_Type_C9_ALS_FTD_vs_HC")
  }
)

NonC9_3x3_woSex <- sapply(
  DDS_Case_Control_woSex_list, 
  FUN = function(x){
    get_N_DEGs(x, alpha = 0.05, name = "Case_Type_ALS_FTD_vs_HC"   )
  }
)


DDS_Case_Control_withSex_list <- list()
dds$Sex <- factor(dds$Sex)

design(dds) <- ~ Sex + Case_Type
register(MulticoreParam(nthr))  
DDS_Case_Control_withSex_list[["NoCov"]] <- DESeq(
  dds, 
  test = "Wald", 
  parallel = TRUE
)

for (i in 1:ncol(svseq$sv)){
  
  design(dds) <- as.formula(
    paste0(
      "~ Sex + ", 
      paste0(
        "SV", 
        1:i, 
        collapse = " + "
      ), 
      " + Case_Type"
    )
  ) 
  
  message(
    design(dds)
  ) 
  
  DDS_Case_Control_withSex_list[[
    paste0(
      "DDS_", 
      i, 
      "SVs"
    )
  ]] <- DESeq(
    dds, 
    test = "Wald", 
    parallel = TRUE
  )
  
  print("N DEGs C9: ")
  message(
    get_N_DEGs(
      DDS_Case_Control_withSex_list[[
        paste0(
          "DDS_", 
          i, 
          "SVs"
        )
      ]], 
      name = "Case_Type_C9_ALS_FTD_vs_HC"
    )
  )
  
  print("N DEGs NonC9: ")
  message(
    get_N_DEGs(
      DDS_Case_Control_withSex_list[[
        paste0(
          "DDS_", 
          i, 
          "SVs"
        )
      ]], 
      name = "Case_Type_ALS_FTD_vs_HC"   
    )
  )
  
}


rm(dds, DDS, dat, idx, mod, mod0, svseq, i)


C9_3x3_withSex <- sapply(
  DDS_Case_Control_withSex_list, 
  FUN = function(x){
    get_N_DEGs(x, alpha = 0.05, name = "Case_Type_C9_ALS_FTD_vs_HC")
  }
)

NonC9_3x3_withSex <- sapply(
  DDS_Case_Control_withSex_list, 
  FUN = function(x){
    get_N_DEGs(x, alpha = 0.05, name = "Case_Type_ALS_FTD_vs_HC")
  }
)

rm(DDS_Case_Control_withSex_list, DDS_Case_Control_woSex_list, Cases_3x3)


      # 3.1.8 Visualize Data ---------------------------------------------------

df <- data.frame(
  Case = c(rep(c("Case") , each = 28), rep("Rand", 14))
)
df$Cov = rep(0:13, 3)
df$Sex = c(rep(c("Wo", "With", "Wo"), each=14)) 
df$DEGs = c(
  Case_woSex[1:14], 
  Case_withSex[1:14], 
  Rand_woSex[1:14]
)
df$Tmp <- paste0(df$Case, "_", df$Sex)
ggplot(df) + 
  aes(x = Cov, y = DEGs, shape = Sex) + 
  geom_line(aes(col = Case)) + 
  geom_point(aes(fill = Case), size=3) + 
  scale_shape_manual(values = c("With"=24, "Wo"=21)) + 
  scale_y_continuous(limits = c(0, 7000)) + 
  scale_fill_manual(values = c("Case"="#126C94", "Rand"="#CCCCCC")) + 
  scale_color_manual(values = c("Case"="#126C94", "Rand"="#CCCCCC")) + 
  xlab("\n# of SVs") + 
  ylab("# DEGs \n") + 
  theme_classic() + 
  theme(
    axis.line = element_line(linewidth = 1, color = "#000000"), 
    axis.title = element_text(face = "bold.italic", size=14, color="#000000"), 
    axis.text = element_text(face="bold.italic", size=12, color = "#000000"), 
    axis.ticks = element_line(linewidth = 1, color = "#000000"), 
    axis.ticks.length = unit(0.2, "cm"), 
    legend.text = element_text(face = "bold.italic", size = 12, color = "#000000"), 
    legend.title = element_text(face = "bold.italic", size = 12, color = "#000000")
  )

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/SupplFig_C9/", 
    "DEGs_SVA_14SVs_C9vsNonC9_Rand", 
    ".pdf"
  ),
  width = 150, 
  height = 110, 
  units="mm"
)  
rm(df)


df <- data.frame(
  Case = c(rep(c("C9", "NonC9") , each = 28))
)
df$Cov = rep(0:13, 4)
df$Sex = c(rep(c("Wo", "With"), each=14)) 
df$DEGs = c(
  C9_HC_Control_woSex[1:14], 
  C9_HC_Control_withSex[1:14],
  NonC9_HC_Control_woSex[1:14], 
  NonC9_HC_Control_withSex[1:14]
)
df$Tmp <- paste0(df$Case, "_", df$Sex)
ggplot(df) + 
  aes(x = Cov, y = DEGs, shape = Sex) + 
  geom_line(aes(col = Case)) + 
  geom_point(aes(fill = Case), size=3) + 
  scale_shape_manual(values = c("With"=24, "Wo"=21)) + 
  scale_y_continuous(limits = c(0, 7000)) + 
  scale_fill_manual(values = c("C9"="#319DD8", "NonC9"="#126C94")) + 
  scale_color_manual(values = c("C9"="#319DD8", "NonC9"="#126C94")) + 
  xlab("\n# of SVs") + 
  ylab("# DEGs \n") + 
  theme_classic() + 
  theme(
    axis.line = element_line(linewidth = 1, color = "#000000"), 
    axis.title = element_text(face = "bold.italic", size=14, color="#000000"), 
    axis.text = element_text(face="bold.italic", size=12, color = "#000000"), 
    axis.ticks = element_line(linewidth = 1, color = "#000000"), 
    axis.ticks.length = unit(0.2, "cm"), 
    legend.text = element_text(face = "bold.italic", size = 12, color = "#000000"), 
    legend.title = element_text(face = "bold.italic", size = 12, color = "#000000")
  )

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/SupplFig_C9/", 
    "DEGs_SVA_14SVs_C9vsHC_NonC9_vs_HC", 
    ".pdf"
  ),
  width = 150, 
  height = 110, 
  units="mm"
)  
rm(df)


df <- data.frame(
  Case = c(rep(c("C9", "NonC9") , each = 28))
)
df$Cov = rep(0:13, 4)
df$Sex = c(rep(c("Wo", "With"), each=14)) 
df$DEGs = c(
  C9_3x3_woSex[1:14], 
  C9_3x3_withSex[1:14],
  NonC9_3x3_woSex[1:14], 
  NonC9_3x3_withSex[1:14]
)
df$Tmp <- paste0(df$Case, "_", df$Sex)
ggplot(df) + 
  aes(x = Cov, y = DEGs, shape = Sex) + 
  geom_line(aes(col = Case)) + 
  geom_point(aes(fill = Case), size=3) + 
  scale_shape_manual(values = c("With"=24, "Wo"=21)) + 
  scale_y_continuous(limits = c(0, 7000)) + 
  scale_fill_manual(values = c("C9"="#319DD8", "NonC9"="#126C94")) + 
  scale_color_manual(values = c("C9"="#319DD8", "NonC9"="#126C94")) + 
  xlab("\n# of SVs") + 
  ylab("# DEGs \n") + 
  theme_classic() + 
  theme(
    axis.line = element_line(linewidth = 1, color = "#000000"), 
    axis.title = element_text(face = "bold.italic", size=14, color="#000000"), 
    axis.text = element_text(face="bold.italic", size=12, color = "#000000"), 
    axis.ticks = element_line(linewidth = 1, color = "#000000"), 
    axis.ticks.length = unit(0.2, "cm"), 
    legend.text = element_text(face = "bold.italic", size = 12, color = "#000000"), 
    legend.title = element_text(face = "bold.italic", size = 12, color = "#000000")
  )


ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/SupplFig_C9/", 
    "DEGs_SVA_14SVs_C9vsHC_NonC9_vs_HC_3x3", 
    ".pdf"
  ),
  width = 150, 
  height = 110, 
  units="mm"
)  
rm(df)



    ## 3.2 WNN_L25 -------------------------------------------------------------


      # 3.2.1 Double-check results names ---------------------------------------

ind <- which(
  Index_SVA_RNA$CellTypeLevel == "WNN_L25"
)

res.name_odd <- resultsNames(
  DDS_list_RNA[["DDS_5_SVs"]][ind[ind %% 2 == 1][1]][[1]]
) %>% tail(n = 1)  
print(res.name_odd)

sapply(
  DDS_list_RNA[["DDS_5_SVs"]][ind[ind %% 2 == 1]], 
  FUN = function(x){
    return(
      tail(
        resultsNames(x), 
        n = 1
      ) == res.name_odd
    )
  }
) %>% all 

res.name_even <- resultsNames(
  DDS_list_RNA[["DDS_5_SVs"]][ind[ind %% 2 == 0][1]][[1]]
) %>% tail(n = 1)  
print(res.name_even)

sapply(
  DDS_list_RNA[["DDS_5_SVs"]][ind[ind %% 2 == 0]], 
  FUN = function(x){
    return(
      tail(
        resultsNames(x), 
        n = 1
      ) == res.name_even
    )
  }
) %>% all 




      # 3.2.2 Collect WNN_L25 # DEGs -------------------------------------------

WNN_L25_nDEGs_C9vsNonC9 <- sapply(
  DDS_list_RNA[["DDS_5_SVs"]][ind[ind %% 2 == 1]], 
  FUN = function(x){
    get_N_DEGs(
      x, 
      alpha = 0.05, 
      name = res.name_odd
    )
  }
)

names(WNN_L25_nDEGs_C9vsNonC9) <- Index_SVA_RNA$CellType[ind[ind %% 2 == 1]]
  

WNN_L25_nDEGs_C9vsNonC9_Up <- sapply(
  DDS_list_RNA[["DDS_5_SVs"]][ind[ind %% 2 == 1]], 
  FUN = function(x){
    get_N_DEGs2(
      x, 
      alpha = 0.05, 
      name = res.name_odd, 
      direction = "up"
    )
  }
)

names(WNN_L25_nDEGs_C9vsNonC9_Up) <- Index_SVA_RNA$CellType[ind[ind %% 2 == 1]]


WNN_L25_nDEGs_C9vsNonC9_Down <- sapply(
  DDS_list_RNA[["DDS_5_SVs"]][ind[ind %% 2 == 1]], 
  FUN = function(x){
    get_N_DEGs2(
      x, 
      alpha = 0.05, 
      name = res.name_odd, 
      direction = "down"
    )
  }
)

names(WNN_L25_nDEGs_C9vsNonC9_Down) <- Index_SVA_RNA$CellType[ind[ind %% 2 == 1]]

WNN_L25_nDEGs_C9vsNonC9_Up + WNN_L25_nDEGs_C9vsNonC9_Down == WNN_L25_nDEGs_C9vsNonC9

WNN_L25_nDEGs_R2vsR1 <- sapply(
  DDS_list_RNA[["DDS_5_SVs"]][ind[ind %% 2 == 0]], 
  FUN = function(x){
    get_N_DEGs(
      x, 
      alpha = 0.05, 
      name = res.name_even
    )
  }
)

names(WNN_L25_nDEGs_R2vsR1) <- Index_SVA_RNA$CellType[ind[ind %% 2 == 0]]

WNN_L25_nDEGs_C9vsNonC9
WNN_L25_nDEGs_R2vsR1 

rm(res.name_even, res.name_odd)


      # Double-check direction -------------------------------------------------

DDS.tmp <- DDS_list_RNA[["DDS_5_SVs"]][[which(
  Index_SVA_RNA$CellTypeLevel == "WNN_L25" & 
    Index_SVA_RNA$CellType == "Exc_RORB" & 
      Index_SVA_RNA$Comparison == "ALSFTD_C9vsALSFTD_NonC9"
)]]
resultsNames(DDS.tmp)
res.tmp <- results(
  DDS.tmp, 
  alpha = 0.05
)
summary(res.tmp)

plotCounts(DDS.tmp, gene = rownames(res.tmp)[res.tmp$log2FoldChange>0][which.min(res.tmp$padj[res.tmp$log2FoldChange>0])], intgroup = "Case_Type", normalized = TRUE)
res.tmp[rownames(res.tmp)=="AC122138.1",]

plotCounts(DDS.tmp, gene = rownames(res.tmp)[res.tmp$log2FoldChange<0][which.min(res.tmp$padj[res.tmp$log2FoldChange<0])], intgroup = "Case_Type", normalized = TRUE)
res.tmp[rownames(res.tmp)=="AC108868.1",]

rm(DDS.tmp, res.tmp)


Results_WNN_L25 <- data.frame(
  CellType = names(WNN_L25_nDEGs_C9vsNonC9)
)
Results_WNN_L25$Up_q_0.05_C9 <- WNN_L25_nDEGs_C9vsNonC9_Down[match(Results_WNN_L25$CellType, names(WNN_L25_nDEGs_C9vsNonC9_Down))] 
Results_WNN_L25$Down_q_0.05_C9 <- WNN_L25_nDEGs_C9vsNonC9_Up[match(Results_WNN_L25$CellType, names(WNN_L25_nDEGs_C9vsNonC9_Up))] 

ggplot(
  Results_WNN_L25
) + 
  aes(CellType, Up_q_0.05_C9, fill=CellType) + 
  geom_col(col="#000000", fill="#FF5744") + 
  scale_y_continuous(expand=c(0,0,0.2,0), limits=c(0,1800)) + 
  theme_classic() + 
  theme(
    axis.text.x = element_text(angle=43, hjust=1.0), 
    panel.grid.major.y = element_line(color="#DDDDDD")
  ) 

ggsave(
  paste0(
    "../Data/Visualization/Figures/SupplFig_C9/", 
    "DE_RNA_WNN_L25_nDEG_Up_C9",
    ".pdf"
  ), 
  width = 6.92, 
  height = 2.833, 
  units="in"
)


ggplot(
  Results_WNN_L25
) + 
  aes(CellType, Down_q_0.05_C9, fill=CellType) + 
  geom_col(fill="#3354E2", col="#000000") + 
  scale_y_reverse(expand=c(0.1,0,0,0), limits=c(1800,0)) + 
  theme_classic() + 
  theme(
    axis.text.x = element_text(angle=43, hjust=1.0), 
    panel.grid.major.y = element_line(color="#DDDDDD")
  ) 

ggsave(
  paste0(
    "../Data/Visualization/Figures/SupplFig_C9/", 
    "DE_RNA_WNN_L25_nDEG_Down_C9", 
    ".pdf"
  ), 
  width = 6.92, 
  height = 2.833, 
  units="in"
)

rm(
  list = ls()[
    which(
      as.character(ls()) %in% c("get_N_DEGs", "get_N_DEGs2", "nthr")
    )
  ]
)

rm(dat, dds, dds_All_Cases, DDS_list, Index_dds_All_Cases, Index_SVA, mod, mod0, Results_WNN_L25, 
   sample_data, svseq, C9_3x3_withSex, C9_3x3_woSex, C9_HC_Control_withSex, C9_HC_Control_woSex, 
   Case_withSex, Case_woSex, idx, ind, NonC9_3x3_withSex, NonC9_3x3_woSex, NonC9_HC_Control_withSex, 
   NonC9_HC_Control_woSex, Rand_woSex, res.name_even, res.name_odd, WNN_L25_nDEGs_C9vsNonC9, 
   WNN_L25_nDEGs_C9vsNonC9_Down, WNN_L25_nDEGs_C9vsNonC9_Up, WNN_L25_nDEGs_R2vsR1)




  ### 4.0 Compare DARs ---------------------------------------------------------



    ## 4.1 All Cells -----------------------------------------------------------

DDS_ATAC_AllCells <- DDS_list_ATAC.tmp[[
  which(
    Index_ATAC.tmp$CellTypeLevel == "AllCells" & 
      Index_ATAC.tmp$Comparison == "All_Cases"
  )
]]
sample_data_ATAC <- as.data.frame(
  colData(
    DDS_ATAC_AllCells
  )
) 


      # 4.2 Generate base DDS objects ------------------------------------------

sample_data_Case <- sample_data_ATAC[
  sample_data_ATAC$Case_Type %in% c(
    "ALS_FTD", "C9_ALS_FTD"
  ), 
]

sample_data_Case$Case_Type <- as.character(
  sample_data_Case$Case_Type
)
sample_data_Case$Sex <- as.character(
  sample_data_Case$Sex
)

counts.tmp <- counts(
  DDS_ATAC_AllCells, 
  normalized = FALSE
)
counts.tmp <- counts.tmp[
  , colnames(counts.tmp) %in% sample_data_Case$ID
]

sample_data_Case$Case_Type <- factor(
  sample_data_Case$Case_Type
)
sample_data_Case$Sex <- factor(
  sample_data_Case$Sex
)

DDS_ATAC_AllCells_Case <- DESeqDataSetFromMatrix(
  countData = counts.tmp, 
  colData = sample_data_Case, 
  design =  ~ Case_Type
)

rm(counts.tmp, sample_data_Case)


     
dds.tmp <- estimateSizeFactors(DDS_ATAC_AllCells_Case)
dat  <- counts(dds.tmp, normalized = TRUE)
idx  <- rowMeans(dat) > 1 
dat  <- dat[idx, ] 
mod  <- model.matrix(~ Case_Type, colData(dds.tmp))
mod0 <- model.matrix(~ 1, colData(dds.tmp))

svseq <- svaseq(dat, mod, mod0, n.sv=14)

colnames(svseq$sv) <- paste0("SV", 1:ncol(svseq$sv))
colData(dds.tmp) <- cbind(colData(dds.tmp), svseq$sv)
DDS_ATAC_AllCells_Case <- dds.tmp
rm(dds.tmp, dat, idx, mod, mod0, svseq)

DDS_ATAC_AllCells_Case$Case_Type


    ## 4.3 Case w/o Sex as covariate -------------------------------------------

ATAC_Case_AllCell_woSex_DDS_List <- list() 

DDS.tmp <- DDS_ATAC_AllCells_Case
design(DDS.tmp) <- ~ Case_Type

ATAC_Case_AllCell_woSex_DDS_List[["NoCov"]] <- DDS.tmp
rm(DDS.tmp)

for (i in 1:13){
  DDS.tmp <- DDS_ATAC_AllCells_Case
  design(DDS.tmp) <- as.formula(
    paste0(
      as.character(design(DDS.tmp))[1], 
      paste0(paste0("SV", 1:i), collapse = " + "), 
      " + ", 
      as.character(design(DDS.tmp))[2]
    ) 
  ) 
  message(
    paste0(
      "Design: ", 
      design(DDS.tmp)
    )
  )
  ATAC_Case_AllCell_woSex_DDS_List[[
    paste0(
      "DDS_", 
      i, 
      "_SVs"
    )
  ]] <- DDS.tmp 
  rm(i, DDS.tmp)
}

register(MulticoreParam(nthr))  

for (i in 1:length(ATAC_Case_AllCell_woSex_DDS_List)){
  
  message(
    paste0(
      "Comparison: ", 
      names(ATAC_Case_AllCell_woSex_DDS_List)[i]
    )
  )
  message(
    paste0(
      "Design: ", 
      design(ATAC_Case_AllCell_woSex_DDS_List[[i]])
    )
  )
  message(
    "Dimensions: ", 
    dim(ATAC_Case_AllCell_woSex_DDS_List[[i]])
  )
  
  ATAC_Case_AllCell_woSex_DDS_List[[i]] <- DESeq(
    ATAC_Case_AllCell_woSex_DDS_List[[i]], 
    test = "Wald", 
    parallel = TRUE
  ) 
  
  rm(i)
}



    ## 4.4 Case with Sex as covariate ------------------------------------------

ATAC_Case_AllCell_withSex_DDS_List <- list() 

DDS.tmp <- DDS_ATAC_AllCells_Case
design(DDS.tmp) <- ~ Sex + Case_Type

ATAC_Case_AllCell_withSex_DDS_List[["NoCov"]] <- DDS.tmp
rm(DDS.tmp)

for (i in 1:13){
  DDS.tmp <- DDS_ATAC_AllCells_Case
  design(DDS.tmp) <- as.formula(
    paste0(
      as.character(design(DDS.tmp))[1], 
      " Sex + ", 
      paste0(paste0("SV", 1:i), collapse = " + "), 
      " + ", 
      as.character(design(DDS.tmp))[2]
    ) 
  ) 
  message(
    paste0(
      "Design: ", 
      design(DDS.tmp)
    )
  )
  ATAC_Case_AllCell_withSex_DDS_List[[
    paste0(
      "DDS_", 
      i, 
      "_SVs"
    )
  ]] <- DDS.tmp 
  rm(i, DDS.tmp)
}


register(MulticoreParam(nthr))  

for (i in 1:length(ATAC_Case_AllCell_withSex_DDS_List)){
  
  message(
    paste0(
      "Comparison: ", 
      names(ATAC_Case_AllCell_withSex_DDS_List)[i]
    )
  )
  message(
    paste0(
      "Design: ", 
      design(ATAC_Case_AllCell_withSex_DDS_List[[i]])
    )
  )
  message(
    "Dimensions: ", 
    dim(ATAC_Case_AllCell_withSex_DDS_List[[i]])
  )
  
  ATAC_Case_AllCell_withSex_DDS_List[[i]] <- DESeq(
    ATAC_Case_AllCell_withSex_DDS_List[[i]], 
    test = "Wald", 
    parallel = TRUE
  ) 
  
  rm(i)
}


    ## 4.5 Rand with Sex as covariate ------------------------------------------ 

ATAC_Rand_AllCell_withSex_DDS_List <- list() 

DDS.tmp <- DDS_ATAC_AllCells_Case
DDS.tmp$Rand2 <-as.character(DDS.tmp$Rand)
DDS.tmp$Rand2[DDS.tmp$Rand2 %in% c("R2", "R3")] <- "R2"
DDS.tmp$Rand2 <-as.factor(DDS.tmp$Rand2)
table(DDS.tmp$Rand2, DDS.tmp$Case_Type)
table(DDS.tmp$Case_Type, DDS.tmp$Sex)
table(DDS.tmp$Rand2, DDS.tmp$Sex)

design(DDS.tmp) <- ~ Rand2
ATAC_Rand_AllCell_withSex_DDS_List[["NoCov"]] <- DDS.tmp
DDS.tmp2 <- DDS.tmp
rm(DDS.tmp)

for (i in 1:13){
  
  DDS.tmp <- DDS.tmp2
  design(DDS.tmp) <- as.formula(
    paste0(
      as.character(design(DDS.tmp))[1], 
      " Sex + ", 
      paste0(paste0("SV", 1:i), collapse = " + "), 
      " + ", 
      as.character(design(DDS.tmp))[2]
    ) 
  ) 
  message(
    paste0(
      "Design: ", 
      design(DDS.tmp)
    )
  )
  ATAC_Rand_AllCell_withSex_DDS_List[[
    paste0(
      "DDS_", 
      i, 
      "_SVs"
    )
  ]] <- DDS.tmp 
  rm(i, DDS.tmp)
}


register(MulticoreParam(nthr))  

for (i in 1:length(ATAC_Rand_AllCell_withSex_DDS_List)){
  
  message(
    paste0(
      "Comparison: ", 
      names(ATAC_Rand_AllCell_withSex_DDS_List)[i]
    )
  )
  message(
    paste0(
      "Design: ", 
      design(ATAC_Rand_AllCell_withSex_DDS_List[[i]])
    )
  )
  message(
    "Dimensions: ", 
    dim(ATAC_Rand_AllCell_withSex_DDS_List[[i]])
  )
  
  ATAC_Rand_AllCell_withSex_DDS_List[[i]] <- DESeq(
    ATAC_Rand_AllCell_withSex_DDS_List[[i]], 
    test = "Wald", 
    parallel = TRUE
  ) 
  
  rm(i)
}
 
rm(DDS.tmp2)



    ## 4.6 C9_ALS_FTD vs. HC ---------------------------------------------------


      # Select matching HCs ----------------------------------------------------

table(
  sample_data_ATAC$Case_Type
) 

table(
  sample_data_ATAC$Sex[sample_data_ATAC$Case_Type == "ALS_FTD"]
)

table(
  sample_data_ATAC$Case_Type[sample_data_ATAC$Case == "ALS_FTD"], 
  sample_data_ATAC$Sex[sample_data_ATAC$Case == "ALS_FTD"]
)

set.seed(seed=42)
HCs_for_C9_ALSFTD <- sample(
  sample_data_ATAC$ID[sample_data_ATAC$Case == "HC" & sample_data_ATAC$Sex == "f"], 
  size = 3, 
  replace = FALSE
) 

set.seed(seed=42)
HCs_for_C9_ALSFTD <- c(
  HCs_for_C9_ALSFTD, 
  sample(
    sample_data_ATAC$ID[sample_data_ATAC$Case == "HC" & sample_data_ATAC$Sex == "m"], 
    size = 7, 
    replace = FALSE
  ) 
) 

table(
  sample_data_ATAC$Case[sample_data_ATAC$Case_Type %in% c("ALS_FTD")], 
  sample_data_ATAC$Sex[sample_data_ATAC$Case_Type %in% c("ALS_FTD")]
)

table(
  sample_data_ATAC$Case[sample_data_ATAC$ID %in% HCs_for_C9_ALSFTD], 
  sample_data_ATAC$Sex[sample_data_ATAC$ID %in% HCs_for_C9_ALSFTD]
)

counts <- counts(DDS_ATAC_AllCells, normalized = FALSE) 
coldata <- colData(DDS_ATAC_AllCells)
all(
  colnames(counts) == rownames(coldata)
)

counts <- counts[
  , colnames(counts) %in% c(
    coldata$ID[coldata$Case_Type=="C9_ALS_FTD"], 
    HCs_for_C9_ALSFTD
  )
]

coldata <- coldata[
  rownames(coldata) %in% c(
    coldata$ID[coldata$Case_Type=="C9_ALS_FTD"], 
    HCs_for_C9_ALSFTD    
  ), 
]
all(
  colnames(counts) == rownames(coldata)
)

dds <- DESeqDataSetFromMatrix(
  countData = counts, 
  colData = coldata, 
  design = ~ Case
)
rm(counts, coldata, HCs_for_C9_ALSFTD)

table(
  dds$Case_Type, 
  dds$Sex
)


        # Add SVA --------------------------------------------------------------

dds <- estimateSizeFactors(dds)
dat  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1 
dat  <- dat[idx, ] 
mod  <- model.matrix(~ Case, colData(dds))
mod0 <- model.matrix(~ 1, colData(dds))

svseq <- svaseq(dat, mod, mod0, n.sv=14)
colnames(svseq$sv) <- paste0("SV", 1:ncol(svseq$sv))
colData(dds) <- cbind(colData(dds), svseq$sv)


ATAC_HC_C9_AllCell_woSex_list <- list()
dds$Sex <- factor(dds$Sex)

design(dds) <- ~ Case

register(MulticoreParam(nthr))  
ATAC_HC_C9_AllCell_woSex_list[["NoCov"]] <- DESeq(
  dds, 
  test = "Wald", 
  parallel = TRUE
)

for (i in 1:ncol(svseq$sv)){
  
  design(dds) <- as.formula(
    paste0(
      "~ ", 
      paste0(
        "SV", 
        1:i, 
        collapse = " + "
      ), 
      " + Case"
    )
  ) 
  
  message(
    design(dds)
  ) 
  
  ATAC_HC_C9_AllCell_woSex_list[[
    paste0(
      "DDS_", 
      i, 
      "SVs"
    )
  ]] <- DESeq(
    dds, 
    test = "Wald", 
    parallel = TRUE
  )
  
  message(
    get_N_DEGs(
      ATAC_HC_C9_AllCell_woSex_list[[
        paste0(
          "DDS_", 
          i, 
          "SVs"
        )
      ]]
    )
  )
  
}

ATAC_HC_C9_AllCell_withSex_list <- ATAC_HC_C9_AllCell_woSex_list

for (i in 1:length(ATAC_HC_C9_AllCell_withSex_list)){
  design(ATAC_HC_C9_AllCell_withSex_list[[i]]) <- 
    as.formula(
      paste0(
        as.character(design(ATAC_HC_C9_AllCell_withSex_list[[i]]))[1], 
        " Sex + ", 
        as.character(design(ATAC_HC_C9_AllCell_withSex_list[[i]]))[2] 
      )
    )
  print(design(ATAC_HC_C9_AllCell_withSex_list[[i]]))
}


for (i in 1:length(ATAC_HC_C9_AllCell_withSex_list)){
  
  message(design(ATAC_HC_C9_AllCell_withSex_list[[i]]))
  ATAC_HC_C9_AllCell_withSex_list[[i]] <- DESeq(
    ATAC_HC_C9_AllCell_withSex_list[[i]], 
    test="Wald", 
    parallel = TRUE
  )
  
}


    ## 4.7 NonC9_ALS_FTD vs. HC ------------------------------------------------

      # Select matching HCs ----------------------------------------------------

table(
  sample_data_ATAC$Case_Type
) 

table(
  sample_data_ATAC$Sex[sample_data_ATAC$Case_Type == "C9_ALS_FTD"]
)

table(
  sample_data_ATAC$Case_Type[sample_data_ATAC$Case == "ALS_FTD"], 
  sample_data_ATAC$Sex[sample_data_ATAC$Case == "ALS_FTD"]
)

set.seed(seed=42)
HCs_for_NonC9_ALSFTD <- sample(
  sample_data_ATAC$ID[sample_data_ATAC$Case == "HC" & sample_data_ATAC$Sex == "f"], 
  size = 4, 
  replace = FALSE
) 

set.seed(seed=42)
HCs_for_NonC9_ALSFTD <- c(
  HCs_for_NonC9_ALSFTD, 
  sample(
    sample_data_ATAC$ID[sample_data_ATAC$Case == "HC" & sample_data_ATAC$Sex == "m"], 
    size = 3, 
    replace = FALSE
  ) 
) 

table(
  sample_data_ATAC$Case[sample_data_ATAC$Case_Type %in% c("C9_ALS_FTD")], 
  sample_data_ATAC$Sex[sample_data_ATAC$Case_Type %in% c("C9_ALS_FTD")]
)

table(
  sample_data_ATAC$Case[sample_data_ATAC$ID %in% HCs_for_NonC9_ALSFTD], 
  sample_data_ATAC$Sex[sample_data_ATAC$ID %in% HCs_for_NonC9_ALSFTD]
)

counts <- counts(DDS_ATAC_AllCells, normalized = FALSE) 
coldata <- colData(DDS_ATAC_AllCells)
all(
  colnames(counts) == rownames(coldata)
)

counts <- counts[
  , colnames(counts) %in% c(
    coldata$ID[coldata$Case_Type=="ALS_FTD"], 
    HCs_for_NonC9_ALSFTD
  )
]

coldata <- coldata[
  rownames(coldata) %in% c(
    coldata$ID[coldata$Case_Type=="ALS_FTD"], 
    HCs_for_NonC9_ALSFTD    
  ), 
]
all(
  colnames(counts) == rownames(coldata)
)

dds <- DESeqDataSetFromMatrix(
  countData = counts, 
  colData = coldata, 
  design = ~ Case
)
rm(counts, coldata, HCs_for_NonC9_ALSFTD)

table(
  dds$Case_Type, 
  dds$Sex
)


      # Add SVA --------------------------------------------------------------

dds <- estimateSizeFactors(dds)
dat  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1 
dat  <- dat[idx, ] 
mod  <- model.matrix(~ Case, colData(dds))
mod0 <- model.matrix(~ 1, colData(dds))

svseq <- svaseq(dat, mod, mod0, n.sv=14)
colnames(svseq$sv) <- paste0("SV", 1:ncol(svseq$sv))
colData(dds) <- cbind(colData(dds), svseq$sv)


ATAC_HC_NonC9_AllCell_woSex_list <- list()
dds$Sex <- factor(dds$Sex)

design(dds) <- ~ Case

register(MulticoreParam(nthr))  
ATAC_HC_NonC9_AllCell_woSex_list[["NoCov"]] <- DESeq(
  dds, 
  test = "Wald", 
  parallel = TRUE
)

for (i in 1:ncol(svseq$sv)){
  
  design(dds) <- as.formula(
    paste0(
      "~ ", 
      paste0(
        "SV", 
        1:i, 
        collapse = " + "
      ), 
      " + Case"
    )
  ) 
  
  message(
    design(dds)
  ) 
  
  ATAC_HC_NonC9_AllCell_woSex_list[[
    paste0(
      "DDS_", 
      i, 
      "SVs"
    )
  ]] <- DESeq(
    dds, 
    test = "Wald", 
    parallel = TRUE
  )
  
  message(
    get_N_DEGs(
      ATAC_HC_NonC9_AllCell_woSex_list[[
        paste0(
          "DDS_", 
          i, 
          "SVs"
        )
      ]]
    )
  )
  
}

ATAC_HC_NonC9_AllCell_withSex_list <- ATAC_HC_NonC9_AllCell_woSex_list

for (i in 1:length(ATAC_HC_NonC9_AllCell_withSex_list)){
  design(ATAC_HC_NonC9_AllCell_withSex_list[[i]]) <- 
    as.formula(
      paste0(
        as.character(design(ATAC_HC_NonC9_AllCell_withSex_list[[i]]))[1], 
        " Sex + ", 
        as.character(design(ATAC_HC_NonC9_AllCell_withSex_list[[i]]))[2] 
      )
    )
  print(design(ATAC_HC_NonC9_AllCell_withSex_list[[i]]))
}


for (i in 1:length(ATAC_HC_NonC9_AllCell_withSex_list)){
  
  message(design(ATAC_HC_NonC9_AllCell_withSex_list[[i]]))
  ATAC_HC_NonC9_AllCell_withSex_list[[i]] <- DESeq(
    ATAC_HC_NonC9_AllCell_withSex_list[[i]], 
    test="Wald", 
    parallel = TRUE
  )
  
}


    ## 4.8 3/3, 3/3, 3/3 Comparison --------------------------------------------


        # Select matching HCs --------------------------------------------------

table(
  sample_data_ATAC$Case_Type[sample_data_ATAC$Case == "ALS_FTD"], 
  sample_data_ATAC$Sex[sample_data_ATAC$Case == "ALS_FTD"]
)

set.seed(seed = 42)
HCs_for_3x3 <- sample(
  sample_data_ATAC$ID[sample_data_ATAC$Case == "HC" & sample_data_ATAC$Sex == "f"], 
  size = 3, 
  replace = FALSE
)
set.seed(seed = 42)
HCs_for_3x3 <- c(
  HCs_for_3x3, 
  sample(
    sample_data_ATAC$ID[sample_data_ATAC$Case == "HC" & sample_data_ATAC$Sex == "m"], 
    size = 3, 
    replace = FALSE
  )
)

set.seed(seed = 42)
C9s_for_3x3 <- sample(
  sample_data_ATAC$ID[sample_data_ATAC$Case_Type == "C9_ALS_FTD" & sample_data_ATAC$Sex == "f"], 
  size = 3, 
  replace = FALSE
) 

set.seed(seed = 42)
C9s_for_3x3 <- c(
  C9s_for_3x3, 
  sample(
    sample_data_ATAC$ID[sample_data_ATAC$Case_Type == "C9_ALS_FTD" & sample_data_ATAC$Sex == "m"], 
    size = 3, 
    replace = FALSE
  ) 
)

set.seed(seed = 42)
NonC9s_for_3x3 <- sample(
  sample_data_ATAC$ID[sample_data_ATAC$Case_Type == "ALS_FTD" & sample_data_ATAC$Sex == "f"], 
  size = 3, 
  replace = FALSE
) 

set.seed(seed = 42)
NonC9s_for_3x3 <- c(
  NonC9s_for_3x3, 
  sample(
    sample_data_ATAC$ID[sample_data_ATAC$Case_Type == "ALS_FTD" & sample_data_ATAC$Sex == "m"], 
    size = 3, 
    replace = FALSE
  ) 
)

any(
  NonC9s_for_3x3 %in% 
    C9s_for_3x3
)

Cases_3x3 <- unique(
  c(
    HCs_for_3x3, 
    NonC9s_for_3x3, 
    C9s_for_3x3
  )
)
rm(HCs_for_3x3, NonC9s_for_3x3, C9s_for_3x3)


counts <- counts(DDS_ATAC_AllCells, normalized = FALSE) 
coldata <- colData(DDS_ATAC_AllCells)
all(
  colnames(counts) == rownames(coldata)
)

counts <- counts[
  , colnames(counts) %in% Cases_3x3
]

coldata <- coldata[
  rownames(coldata) %in% Cases_3x3, 
]
all(
  colnames(counts) == rownames(coldata)
)

dds <- DESeqDataSetFromMatrix(
  countData = counts, 
  colData = coldata, 
  design = ~ Case_Type
)
rm(counts, coldata, Cases_3x3)

table(
  dds$Case_Type, 
  dds$Sex
)


        # Add SVA --------------------------------------------------------------

dds <- estimateSizeFactors(dds)
dat  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1 
dat  <- dat[idx, ] 
mod  <- model.matrix(~ Case_Type, colData(dds))
mod0 <- model.matrix(~ 1, colData(dds))

svseq <- svaseq(dat, mod, mod0, n.sv=14)
colnames(svseq$sv) <- paste0("SV", 1:ncol(svseq$sv))
colData(dds) <- cbind(colData(dds), svseq$sv)


DDS_Case_Control_woSex_list <- list()
dds$Sex <- factor(dds$Sex)

design(dds) <- ~ Case_Type
register(MulticoreParam(nthr))  
DDS_Case_Control_woSex_list[["NoCov"]] <- DESeq(
  dds, 
  test = "Wald", 
  parallel = TRUE
)

for (i in 1:ncol(svseq$sv)){
  
  design(dds) <- as.formula(
    paste0(
      "~ ", 
      paste0(
        "SV", 
        1:i, 
        collapse = " + "
      ), 
      " + Case_Type"
    )
  ) 
  
  message(
    design(dds)
  ) 
  
  DDS_Case_Control_woSex_list[[
    paste0(
      "DDS_", 
      i, 
      "SVs"
    )
  ]] <- DESeq(
    dds, 
    test = "Wald", 
    parallel = TRUE
  )
  
  print("N DEGs C9: ")
  message(
    get_N_DEGs(
      DDS_Case_Control_woSex_list[[
        paste0(
          "DDS_", 
          i, 
          "SVs"
        )
      ]], 
      name = "Case_Type_C9_ALS_FTD_vs_HC"
    )
  )
  
  print("N DEGs NonC9: ")
  message(
    get_N_DEGs(
      DDS_Case_Control_woSex_list[[
        paste0(
          "DDS_", 
          i, 
          "SVs"
        )
      ]], 
      name = "Case_Type_ALS_FTD_vs_HC"
    )
  )
  
}


DDS_Case_Control_withSex_list <- list()
dds$Sex <- factor(dds$Sex)

design(dds) <- ~ Sex + Case_Type
register(MulticoreParam(nthr))  
DDS_Case_Control_withSex_list[["NoCov"]] <- DESeq(
  dds, 
  test = "Wald", 
  parallel = TRUE
)

for (i in 1:ncol(svseq$sv)){
  
  design(dds) <- as.formula(
    paste0(
      "~ Sex + ", 
      paste0(
        "SV", 
        1:i, 
        collapse = " + "
      ), 
      " + Case_Type"
    )
  ) 
  
  message(
    design(dds)
  ) 
  
  DDS_Case_Control_withSex_list[[
    paste0(
      "DDS_", 
      i, 
      "SVs"
    )
  ]] <- DESeq(
    dds, 
    test = "Wald", 
    parallel = TRUE
  )
  
  print("N DEGs C9: ")
  message(
    get_N_DEGs(
      DDS_Case_Control_withSex_list[[
        paste0(
          "DDS_", 
          i, 
          "SVs"
        )
      ]], 
      name = "Case_Type_C9_ALS_FTD_vs_HC"
    )
  )
  
  print("N DEGs NonC9: ")
  message(
    get_N_DEGs(
      DDS_Case_Control_withSex_list[[
        paste0(
          "DDS_", 
          i, 
          "SVs"
        )
      ]], 
      name = "Case_Type_ALS_FTD_vs_HC"   
    )
  )
  
}


rm(dds, DDS, dat, idx, mod, mod0, svseq, i)


ATAC_3x3_withSex <- DDS_Case_Control_withSex_list
ATAC_3x3_woSex <- DDS_Case_Control_woSex_list


rm(DDS_Case_Control_withSex_list, DDS_Case_Control_woSex_list, Cases_3x3)



    ## 4.9 Visualize Data ---------------------------------------------------

df <- data.frame(
  Case = c(rep(c("Case") , each = 28), rep("Rand", 14))
)
df$Cov = rep(0:13, 3)
df$Sex = c(rep(c("Wo", "With", "Wo"), each=14)) 
df$DEGs = c(
  sapply(ATAC_Case_AllCell_woSex_DDS_List[1:14], get_N_DEGs, name = "Case_Type_C9_ALS_FTD_vs_ALS_FTD"), 
  sapply(ATAC_Case_AllCell_withSex_DDS_List[1:14], get_N_DEGs, name = "Case_Type_C9_ALS_FTD_vs_ALS_FTD"), 
  sapply(ATAC_Rand_AllCell_withSex_DDS_List[1:14], get_N_DEGs, name = "Rand2_R2_vs_R1") 
)
df$Tmp <- paste0(df$Case, "_", df$Sex)
ggplot(df) + 
  aes(x = Cov, y = DEGs, shape = Sex) + 
  geom_line(aes(col = Case)) + 
  geom_point(aes(fill = Case), size=3) + 
  scale_shape_manual(values = c("With"=24, "Wo"=21)) + 
  scale_y_continuous(limits = c(0, 12000)) + 
  scale_fill_manual(values = c("Case"="#126C94", "Rand"="#CCCCCC")) + 
  scale_color_manual(values = c("Case"="#126C94", "Rand"="#CCCCCC")) + 
  xlab("\n# of SVs") + 
  ylab("# DEGs \n") + 
  theme_classic() + 
  theme(
    axis.line = element_line(linewidth = 1, color = "#000000"), 
    axis.title = element_text(face = "bold.italic", size=14, color="#000000"), 
    axis.text = element_text(face="bold.italic", size=12, color = "#000000"), 
    axis.ticks = element_line(linewidth = 1, color = "#000000"), 
    axis.ticks.length = unit(0.2, "cm"), 
    legend.text = element_text(face = "bold.italic", size = 12, color = "#000000"), 
    legend.title = element_text(face = "bold.italic", size = 12, color = "#000000")
  )

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/SupplFig_C9/", 
    "DARs_SVA_14SVs_C9vsNonC9_Rand", 
    ".pdf"
  ),
  width = 150, 
  height = 110, 
  units="mm"
)  
rm(df)


df <- data.frame(
  Case = c(rep(c("C9", "NonC9") , each = 28))
)
df$Cov = rep(0:13, 4)
df$Sex = c(rep(c("Wo", "With"), each=14)) 
df$DEGs = c(
  sapply(ATAC_HC_C9_AllCell_woSex_list[1:14], get_N_DEGs, name = "Case_ALS_FTD_vs_HC"),  
  sapply(ATAC_HC_C9_AllCell_withSex_list[1:14], get_N_DEGs, name = "Case_ALS_FTD_vs_HC"),  
  sapply(ATAC_HC_NonC9_AllCell_woSex_list[1:14], get_N_DEGs, name = "Case_ALS_FTD_vs_HC"),  
  sapply(ATAC_HC_NonC9_AllCell_withSex_list[1:14], get_N_DEGs, name = "Case_ALS_FTD_vs_HC")  
)
df$Tmp <- paste0(df$Case, "_", df$Sex)
ggplot(df) + 
  aes(x = Cov, y = DEGs, shape = Sex) + 
  geom_line(aes(col = Case)) + 
  geom_point(aes(fill = Case), size=3) + 
  scale_shape_manual(values = c("With"=24, "Wo"=21)) + 
  scale_y_continuous(limits = c(0, 12000)) + 
  scale_fill_manual(values = c("C9"="#319DD8", "NonC9"="#126C94")) + 
  scale_color_manual(values = c("C9"="#319DD8", "NonC9"="#126C94")) + 
  xlab("\n# of SVs") + 
  ylab("# DEGs \n") + 
  theme_classic() + 
  theme(
    axis.line = element_line(linewidth = 1, color = "#000000"), 
    axis.title = element_text(face = "bold.italic", size=14, color="#000000"), 
    axis.text = element_text(face="bold.italic", size=12, color = "#000000"), 
    axis.ticks = element_line(linewidth = 1, color = "#000000"), 
    axis.ticks.length = unit(0.2, "cm"), 
    legend.text = element_text(face = "bold.italic", size = 12, color = "#000000"), 
    legend.title = element_text(face = "bold.italic", size = 12, color = "#000000")
  )

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/SupplFig_C9/", 
    "DARs_SVA_14SVs_C9vsHC_NonC9_vs_HC", 
    ".pdf"
  ),
  width = 150, 
  height = 110, 
  units="mm"
)  
rm(df)


df <- data.frame(
  Case = c(rep(c("C9", "NonC9") , each = 28))
)
df$Cov = rep(0:13, 4)
df$Sex = c(rep(c("Wo", "With"), each=14)) 
df$DEGs = c(
  sapply(ATAC_3x3_woSex[1:14], get_N_DEGs, name = "Case_Type_C9_ALS_FTD_vs_HC"), 
  sapply(ATAC_3x3_withSex[1:14], get_N_DEGs, name = "Case_Type_C9_ALS_FTD_vs_HC"),
  sapply(ATAC_3x3_woSex[1:14], get_N_DEGs, name = "Case_Type_ALS_FTD_vs_HC"), 
  sapply(ATAC_3x3_withSex[1:14], get_N_DEGs, name = "Case_Type_ALS_FTD_vs_HC")
)
df$Tmp <- paste0(df$Case, "_", df$Sex)
ggplot(df) + 
  aes(x = Cov, y = DEGs, shape = Sex) + 
  geom_line(aes(col = Case)) + 
  geom_point(aes(fill = Case), size=3) + 
  scale_shape_manual(values = c("With"=24, "Wo"=21)) + 
  scale_y_continuous(limits = c(0, 12000)) + 
  scale_fill_manual(values = c("C9"="#319DD8", "NonC9"="#126C94")) + 
  scale_color_manual(values = c("C9"="#319DD8", "NonC9"="#126C94")) + 
  xlab("\n# of SVs") + 
  ylab("# DEGs \n") + 
  theme_classic() + 
  theme(
    axis.line = element_line(linewidth = 1, color = "#000000"), 
    axis.title = element_text(face = "bold.italic", size=14, color="#000000"), 
    axis.text = element_text(face="bold.italic", size=12, color = "#000000"), 
    axis.ticks = element_line(linewidth = 1, color = "#000000"), 
    axis.ticks.length = unit(0.2, "cm"), 
    legend.text = element_text(face = "bold.italic", size = 12, color = "#000000"), 
    legend.title = element_text(face = "bold.italic", size = 12, color = "#000000")
  )


ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/SupplFig_C9/", 
    "DARs_SVA_14SVs_C9vsHC_NonC9_vs_HC_3x3", 
    ".pdf"
  ),
  width = 150, 
  height = 110, 
  units="mm"
)  
rm(df)



    ## 4.10 WNN_L25 -------------------------------------------------------------

ind <- which(
  Index_ATAC.tmp$CellTypeLevel == "WNN_L25" & 
    Index_ATAC.tmp$Comparison == "All_Cases"
)

DDS_ATAC_WNN_L25 <- DDS_list_ATAC.tmp[ind]
names(DDS_ATAC_WNN_L25) <- Index_ATAC.tmp$CellType[ind]

register(MulticoreParam(nthr))  
ATAC_WNN_L25_List2 <- list()

for(i in 1:length(DDS_ATAC_WNN_L25)){
  
  dds <- DDS_ATAC_WNN_L25[[i]]
  dds <- dds[,colnames(dds) %in% dds$ID[dds$Case %in% c("ALS_FTD")]]
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds$Case_Type <- as.character(dds$Case_Type)
  dds$Case_Type <- factor(dds$Case_Type)
  dat  <- counts(dds, normalized = TRUE)
  idx  <- rowMeans(dat) > 1 
  dat  <- dat[idx, ] 
  mod  <- model.matrix(~ Case_Type, colData(dds))
  mod0 <- model.matrix(~ 1, colData(dds))
  
  svseq <- svaseq(dat, mod, mod0, n.sv=5)
  colnames(svseq$sv) <- paste0("SV", 1:ncol(svseq$sv))
  colData(dds) <- cbind(colData(dds), svseq$sv)
  
  
  dds$Sex <- factor(dds$Sex)
  
  design(dds) <- as.formula(
    paste0(
      "~ ", 
      paste0(
        "SV", 
        1:5, 
        collapse = " + "
      ), 
      " + Case_Type"
    )
  ) 
  
  message(
    design(dds)
  ) 
  
  
  ATAC_WNN_L25_List2[[names(DDS_ATAC_WNN_L25)[i]]] <- DESeq(
    dds, 
    test = "Wald", 
    parallel = TRUE
  )
  
}

WNN_L25_nDARs_C9vsNonC9 <- sapply(
  ATAC_WNN_L25_List2,  
  get_N_DEGs, 
  alpha = 0.05, 
  name = "Case_Type_C9_ALS_FTD_vs_ALS_FTD"
)


names(WNN_L25_nDARs_C9vsNonC9) <- names(ATAC_WNN_L25_List2)


WNN_L25_nDARs_C9vsNonC9_Up <- sapply(
  ATAC_WNN_L25_List2,  
  FUN = function(x){
    get_N_DEGs2(
      x, 
      alpha = 0.05, 
      name = "Case_Type_C9_ALS_FTD_vs_ALS_FTD", 
      direction = "up"
    )
  }
)

names(WNN_L25_nDARs_C9vsNonC9_Up) <- names(ATAC_WNN_L25_List2)


WNN_L25_nDARs_C9vsNonC9_Down <- sapply(
  ATAC_WNN_L25_List2, 
  FUN = function(x){
    get_N_DEGs2(
      x, 
      alpha = 0.05, 
      name = "Case_Type_C9_ALS_FTD_vs_ALS_FTD", 
      direction = "down"
    )
  }
)

names(WNN_L25_nDARs_C9vsNonC9_Down) <- names(ATAC_WNN_L25_List2)

WNN_L25_nDARs_C9vsNonC9_Up + WNN_L25_nDARs_C9vsNonC9_Down == WNN_L25_nDARs_C9vsNonC9


      # Double-check direction -------------------------------------------------

resultsNames(ATAC_WNN_L25_List2[[4]])
DDS.tmp <- ATAC_WNN_L25_List2[[4]]
res.tmp <- results(
  DDS.tmp, 
  alpha = 0.05, 
  name = "Case_Type_C9_ALS_FTD_vs_ALS_FTD"
)
summary(res.tmp)

plotCounts(DDS.tmp, gene = rownames(res.tmp)[res.tmp$log2FoldChange>0][which.min(res.tmp$padj[res.tmp$log2FoldChange>0])], intgroup = "Case_Type", normalized = TRUE)
res.tmp[rownames(res.tmp)=="chr8-143158839-143160771",]

plotCounts(DDS.tmp, gene = rownames(res.tmp)[res.tmp$log2FoldChange<0][which.min(res.tmp$padj[res.tmp$log2FoldChange<0])], intgroup = "Case_Type", normalized = TRUE)
res.tmp[rownames(res.tmp)=="chr1-16620254-16621379",]

rm(DDS.tmp, res.tmp)


Results_WNN_L25 <- data.frame(
  CellType = names(WNN_L25_nDARs_C9vsNonC9)
)
Results_WNN_L25$Up_q_0.05_C9 <- WNN_L25_nDARs_C9vsNonC9_Down[match(Results_WNN_L25$CellType, names(WNN_L25_nDARs_C9vsNonC9_Down))] 
Results_WNN_L25$Down_q_0.05_C9 <- WNN_L25_nDARs_C9vsNonC9_Up[match(Results_WNN_L25$CellType, names(WNN_L25_nDARs_C9vsNonC9_Up))] 

ggplot(
  Results_WNN_L25
) + 
  aes(CellType, Up_q_0.05_C9, fill=CellType) + 
  geom_col(col="#000000", fill="#FF5744") + 
  scale_y_continuous(expand=c(0,0,0.2,0), limits=c(0,1800)) + 
  theme_classic() + 
  theme(
    axis.text.x = element_text(angle=43, hjust=1.0), 
    panel.grid.major.y = element_line(color="#DDDDDD")
  ) 

ggsave(
  paste0(
    "../Data/Visualization/Figures/SupplFig_C9/", 
    "DA_ATAC_WNN_L25_nDEG_Up_C9",
    ".pdf"
  ), 
  width = 6.92, 
  height = 2.833, 
  units="in"
)


ggplot(
  Results_WNN_L25
) + 
  aes(CellType, Down_q_0.05_C9, fill=CellType) + 
  geom_col(fill="#3354E2", col="#000000") + 
  scale_y_reverse(expand=c(0.1,0,0,0), limits=c(1800,0)) + 
  theme_classic() + 
  theme(
    axis.text.x = element_text(angle=43, hjust=1.0), 
    panel.grid.major.y = element_line(color="#DDDDDD")
  ) 

ggsave(
  paste0(
    "../Data/Visualization/Figures/SupplFig_C9/", 
    "DA_ATAC_WNN_L25_nDEG_Down_C9", 
    ".pdf"
  ), 
  width = 6.92, 
  height = 2.833, 
  units="in"
)

rm(
  list = ls()[
    which(
      as.character(ls()) %in% c("get_N_DEGs", "get_N_DEGs2", "nthr")
    )
  ]
)

rm(dat, dds, dds_All_Cases, DDS_list, Index_dds_All_Cases, Index_SVA, mod, mod0, Results_WNN_L25, 
   sample_data, svseq, C9_3x3_withSex, C9_3x3_woSex, C9_HC_Control_withSex, C9_HC_Control_woSex, 
   Case_withSex, Case_woSex, idx, ind, NonC9_3x3_withSex, NonC9_3x3_woSex, NonC9_HC_Control_withSex, 
   NonC9_HC_Control_woSex, Rand_woSex, res.name_even, res.name_odd, WNN_L25_nDEGs_C9vsNonC9, 
   WNN_L25_nDEGs_C9vsNonC9_Down, WNN_L25_nDEGs_C9vsNonC9_Up, WNN_L25_nDEGs_R2vsR1)




  ### 5.0 Cell-type composition plots ------------------------------------------

M0_RNA@meta.data %>% 
  filter(Case_Type %in% c("C9_ALS_FTD", "ALS_FTD")) %>% 
  group_by(ID, WNN_L25) %>% 
  summarize(CellNumber = n(), Case_Type = only_value(Case_Type)) %>% 
  group_by(ID) %>% mutate(TotalCells = sum(CellNumber)) %>% 
  ungroup() %>% 
  mutate(Proportion = CellNumber/TotalCells) %>% 
  mutate(WNN_L25 = factor(WNN_L25, levels = c(
    "Exc_RORB", 
    "Exc_LINC00507", 
    "Exc_THEMIS", 
    "Exc_FEZF2", 
    "Inh_LAMP5_PAX6", 
    "Inh_TAFA1_VIP", 
    "Inh_SST", 
    "Inh_PVALB", 
    "OPC", 
    "Oligodendrocytes", 
    "Astrocytes", 
    "Microglia"
  ))) %>% 
  
  ggplot() + 
  aes(WNN_L25, Proportion, fill=Case_Type) + 
  geom_boxplot(position = position_dodge(0.64), width = 0.5, fatten=2, outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.05), alpha = 0.4, size = 2, color = "black") + 
  scale_y_continuous(limits=c(0, 0.6)) + 
  scale_fill_manual(values=ColDict_Case_Type) + 
  xlab("") + 
  theme_classic() + 
  theme(
    axis.line = element_line(colour = "#000000", linewidth = 1), 
    axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1.0, colour = "#000000"), 
    axis.text.y = element_text(size = 10, face = "bold.italic", hjust = 1.0, colour = "#000000"), 
    axis.ticks = element_line(linewidth = 1, colour = "#000000"), 
    axis.ticks.length =  unit(0.2, "cm"), 
  )

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/SupplFig_C9/", 
    "Coda_C9ALSFTD_ALSFTD_WNN_L25_Boxplot", 
    ".pdf"
  ),
  width = 11, 
  height = 5.5, 
  units="in"
)  


MASC_C9_ALSFTD_WNN_L25[[1]] %>% 
  as.data.frame() %>% 
  ggplot() + 
  aes(log2(Case_TypeC9_ALS_FTD.OR), -log10(fdr), fill=fdr<0.05, label = str_replace_all(cluster, "cluster", "" )) + 
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = 2, col = "#BBBBBB") + 
  geom_hline(yintercept = -log10(0.05), linetype = 2, col = "#BBBBBB") + 
  geom_point(size=3, pch=21) + 
  geom_text_repel(max.overlaps = 1) + 
  scale_x_continuous(limits=c(-2,2)) + 
  scale_y_continuous(limits=c(0,2)) + 
  scale_fill_manual(values=c("FALSE"="#CCCCCC", "TRUE"="#EE1111")) + 
  theme_classic()

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/SupplFig_C9/", 
    "Coda_C9ALSFTD_ALSFTD_WNN_L25_MASC_Plot", 
    ".pdf"
  ),
  width = 5.2, 
  height = 4.2, 
  units="in"
)  

  
  
  