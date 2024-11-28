

  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 

library(qs)
library(tidyverse)
library(DESeq2)
library(sva)


qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)




  ### 1.0 Load data ------------------------------------------------------------

DDS_list_All <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/", 
    "DDS_DESEqed_list", 
    ".qrds"
  ), 
  nthr = nthr
)

Index_All <- qread(
  paste0(
    "../Data/DE/WNN/Pseudobulk/DESeq2/AllCase/", 
    "DESeq_Results_Index", 
    ".qrds"
  ), 
  nthr = nthr
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




  ### 3.0 AllCells -------------------------------------------------------------

Sampled_Lists <- ls() 

for (seed_i in 47:51){
  ind <- which(
    Index_All$CellTypeLevel == "AllCells" & 
      Index_All$Comparison == "All_Cases"
  )
  
  DDS_AllCells <- DDS_list_All[[ind]] 
  rm(ind) 
  
  table(DDS_AllCells$Case, DDS_AllCells$Sex)
  
  seed <- seed_i 
  Cases <- c(
    DDS_AllCells$ID[DDS_AllCells$Case == "ALS_FTD"]  
  )
  
  set.seed(seed)
  Cases <- c(
    Cases, 
    sample(
      DDS_AllCells$ID[
        DDS_AllCells$Case == "HC" & 
          DDS_AllCells$Sex == "f"
      ], 
      size = 7, 
      replace = FALSE
    )
  )
  
  set.seed(seed)
  Cases <- c(
    Cases, 
    sample(
      DDS_AllCells$ID[
        DDS_AllCells$Case == "HC" & 
          DDS_AllCells$Sex == "m"
      ], 
      size = 10, 
      replace = FALSE
    )
  )
  
  set.seed(seed)
  Cases <- c(
    Cases, 
    sample(
      DDS_AllCells$ID[
        DDS_AllCells$Case == "ALS" & 
          DDS_AllCells$Sex == "f"
      ], 
      size = 7, 
      replace = FALSE
    )
  )
  
  set.seed(seed)
  Cases <- c(
    Cases, 
    sample(
      DDS_AllCells$ID[
        DDS_AllCells$Case == "ALS" & 
          DDS_AllCells$Sex == "m"
      ], 
      size = 10, 
      replace = FALSE
    )
  )
  
  counts.tmp <- counts(DDS_AllCells, normalized = FALSE)
  counts.tmp <- counts.tmp[, colnames(counts.tmp) %in% Cases]
  
  colData.tmp <- colData(DDS_AllCells)
  colData.tmp <- colData.tmp[rownames(colData.tmp) %in% Cases ,]
  
  all(
    colnames(counts.tmp) == 
      rownames(colData.tmp)
  ) 
  colData.tmp$Sex <- factor(colData.tmp$Sex)
  
  
  DDS.Sampled <- DESeqDataSetFromMatrix(
    countData = counts.tmp, 
    colData = colData.tmp, 
    design =  ~ Case
  )
  
  
  
  dds.tmp <- estimateSizeFactors(DDS.Sampled)
  dat  <- counts(dds.tmp, normalized = TRUE)
  idx  <- rowMeans(dat) > 1 
  dat  <- dat[idx, ] 
  mod  <- model.matrix(~ Case, colData(dds.tmp))
  mod0 <- model.matrix(~ 1, colData(dds.tmp))
  
  svseq <- svaseq(dat, mod, mod0, n.sv=15)
  
  colnames(svseq$sv) <- paste0("SV", 1:ncol(svseq$sv))
  colData(dds.tmp) <- cbind(colData(dds.tmp), svseq$sv)
  DDS.Sampled <- dds.tmp
  rm(dds.tmp, dat, idx, mod, mod0, svseq)
  
  AllCells_DDS_Sampled_List <- list() 
  
  DDS.tmp <- DDS.Sampled
  design(DDS.tmp) <- ~ Case
  
  AllCells_DDS_Sampled_List[["NoCov"]] <- DDS.tmp
  rm(DDS.tmp)
  
  for (i in 1:15){
    DDS.tmp <- DDS.Sampled
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
    AllCells_DDS_Sampled_List[[
      paste0(
        "DDS_", 
        i, 
        "_SVs"
      )
    ]] <- DDS.tmp 
    rm(i, DDS.tmp)
  }
  
  register(MulticoreParam(nthr))  
  
  for (i in 1:length(AllCells_DDS_Sampled_List)){
    
    message(
      paste0(
        "Comparison: ", 
        names(AllCells_DDS_Sampled_List)[i]
      )
    )
    message(
      paste0(
        "Design: ", 
        design(AllCells_DDS_Sampled_List[[i]])
      )
    )
    message(
      "Dimensions: ", 
      dim(AllCells_DDS_Sampled_List[[i]])
    )
    
    AllCells_DDS_Sampled_List[[i]] <- DESeq(
      AllCells_DDS_Sampled_List[[i]], 
      test = "Wald", 
      parallel = TRUE
    ) 
    
    rm(i)
  } 
  
  Sampled_Lists[[as.character(seed_i)]] <- list(AllCells_DDS_Sampled_List)
  rm(AllCells_DDS_Sampled_List) 
  rm(seed_i)
}


# Use samples with comparable library sizes 
df.tmp <- as.data.frame(colData(DDS_AllCells))
all(sizeFactors(DDS_AllCells)==df.tmp$sizeFactor)

ggplot(df.tmp) + 
  aes(Case_Type, sizeFactor) + 
  geom_boxplot() + 
  geom_jitter(size=3, pch=21, width = 0.15)



cases <- df.tmp$ID[df.tmp$Case_Type=="ALS_FTD"]
table(df.tmp$Sex[df.tmp$Case_Type=="ALS_FTD"])

a <- df.tmp$ID[df.tmp$Case=="ALS" & df.tmp$Sex=="f"]
b <- df.tmp$sizeFactor[match(a, df.tmp$ID)] 

a <- a[order(b, decreasing = TRUE)]
cases <- c(
  cases, 
  a[1:3]
)
a <- df.tmp$ID[df.tmp$Case=="ALS" & df.tmp$Sex=="m"]
b <- df.tmp$sizeFactor[match(a, df.tmp$ID)] 

a <- a[order(b, decreasing = TRUE)]
cases <- c(
  cases, 
  a[1:7]
)

a <- df.tmp$ID[df.tmp$Case=="HC" & df.tmp$Sex=="f" & df.tmp$sizeFactor < 3]
b <- df.tmp$sizeFactor[match(a, df.tmp$ID)] 

a <- a[order(b, decreasing = TRUE)]
cases <- c(
  cases, 
  a[1:3]
)
a <- df.tmp$ID[df.tmp$Case=="HC" & df.tmp$Sex=="m" & df.tmp$sizeFactor < 3]
b <- df.tmp$sizeFactor[match(a, df.tmp$ID)] 

a <- a[order(b, decreasing = TRUE)]
cases <- c(
  cases, 
  a[1:7]
)
rm(a, b)

df.tmp <- df.tmp[df.tmp$ID %in% cases, ]
ggplot(df.tmp) + 
  aes(Case, sizeFactor) + 
  geom_boxplot() + 
  geom_jitter(size=3, pch=21, width = 0.15)

c2 <- counts(DDS_AllCells, normalized = FALSE)
co2 <- colData(DDS_AllCells)

dds2 <- DESeqDataSetFromMatrix(
  countData = c2, 
  colData = co2, 
  design = ~ Case
)

DDS2 <- DESeq(
  dds2, 
  test = "Wald", 
  parallel = TRUE
)

resultsNames(DDS2)
summary(
  results(
    DDS2, 
    alpha = 0.05, 
    name = "Case_ALS_vs_HC"
  ), 
  alpha = 0.05
) 

summary(
  results(
    DDS2, 
    alpha = 0.05, 
    name = "Case_ALS_FTD_vs_HC"
  ), 
  alpha = 0.05
)
dds.tmp <- estimateSizeFactors(DDS.Sampled)
dat  <- counts(dds.tmp, normalized = TRUE)
idx  <- rowMeans(dat) > 1 
dat  <- dat[idx, ] 
mod  <- model.matrix(~ Case, colData(dds.tmp))
mod0 <- model.matrix(~ 1, colData(dds.tmp))

svseq <- svaseq(dat, mod, mod0, n.sv=15)

colnames(svseq$sv) <- paste0("SV", 1:ncol(svseq$sv))
colData(dds.tmp) <- cbind(colData(dds.tmp), svseq$sv)
DDS.Sampled <- dds.tmp
rm(dds.tmp, dat, idx, mod, mod0, svseq)

AllCells_DDS_Sampled_List <- list() 

DDS.tmp <- DDS.Sampled
design(DDS.tmp) <- ~ Case

AllCells_DDS_Sampled_List[["NoCov"]] <- DDS.tmp
rm(DDS.tmp)

for (i in 1:15){
  DDS.tmp <- DDS.Sampled
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
  AllCells_DDS_Sampled_List[[
    paste0(
      "DDS_", 
      i, 
      "_SVs"
    )
  ]] <- DDS.tmp 
  rm(i, DDS.tmp)
}

##seed
1+1
rm(DDS_AllCells, seed, Cases, counts.tmp, colData.tmp, DDS.Sampled, )
