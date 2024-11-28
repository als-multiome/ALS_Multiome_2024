
  ### 0.0 Load libraries -------------------------------------------------------

source("~/ALS_Brain_Multiome.Rcfg") 

library(qs) 
library(Seurat)
library(Signac)
library(tidyverse)




  ### 1.0 Load data ------------------------------------------------------------


      # M0 Seurat single-cell data object --------------------------------------

M0 <- qread(
  paste0(
    "../Data/SeuratObjects/", 
    "M0", 
    ".qrds"
  ), 
  nthr = nthr
)




  ### 2.0 Extract sample data --------------------------------------------------

sample_data <-M0@misc$Sample_metrics 
sample_data$Sex <- M0$Sex[match(sample_data$ID, M0$ID)]
setequal(
  unique(paste0(M0@misc$Sample_data$ID, "_", M0@misc$Sample_data$Sex)), 
  paste0(sample_data$ID, "_", sample_data$Sex)
)
sample_data <- sample_data[,c(1,17:19, 2:16)]




  ### 3.0 Generate Random control grouping for all cases -----------------------

sample_data$Rand <- "NA"
sample_data$Case <- factor(sample_data$Case, levels=c("HC", "ALS", "ALS_FTD"))
sample_data$Sex <- factor(sample_data$Sex, levels=c("f", "m"))

Ns <- table(sample_data$Case) %>% 
  setNames(., c("R1", "R2", "R3")) %>% 
  "%*%"(t(prop.table(table(sample_data$Case)))) %>% 
  round(0) %>% 
  t()

i <- which(colSums(Ns)!=table(sample_data$Case))
lapply(i, FUN=function(i){Ns[i,i] <<- Ns[i,i] + table(sample_data$Case)[i]-colSums(Ns)[i]})

Ns %>% 
  addmargins()

rm(i)
Ns

Rand_levels <- c("R1", "R2", "R3")

NsHC <- sample_data %>% 
  filter(Case=="HC") %>% 
  {table(.$Sex)} %>% 
  prop.table  %>% as.numeric() %*%t(as.numeric(Ns[,1])) %>% round() %>% t() 

set.seed(1024)
ind_f <- sample(
  which(sample_data$Case=="HC" & sample_data$Sex=="f"), 
  length(which(sample_data$Case=="HC" & sample_data$Sex=="f"))
)
set.seed(1024)
ind_m <- sample(
  which(sample_data$Case=="HC" & sample_data$Sex=="m"), 
  length(which(sample_data$Case=="HC" & sample_data$Sex=="m"))
)

if(!all(c(colSums(NsHC)==c(length(ind_f), length(ind_m))))) message("Error!") 

sample_data$Rand[ind_f] <- rep(Rand_levels, NsHC[,1])
sample_data$Rand[ind_m] <- rep(Rand_levels, NsHC[,2])

rm(NsHC, ind_f, ind_m)

NsALS <- sample_data %>% 
  filter(Case=="ALS") %>% 
  {table(.$Sex)} %>% 
  prop.table  %>% as.numeric() %*%t(as.numeric(Ns[,2])) %>% round() %>% t() 

set.seed(1024) 
ind_f <- sample(
  which(sample_data$Case=="ALS" & sample_data$Sex=="f"), 
  length(which(sample_data$Case=="ALS" & sample_data$Sex=="f"))
)
set.seed(1024)
ind_m <- sample(
  which(sample_data$Case=="ALS" & sample_data$Sex=="m"), 
  length(which(sample_data$Case=="ALS" & sample_data$Sex=="m"))
)

if(!all(c(colSums(NsALS)==c(length(ind_f), length(ind_m))))) message("Error!")
sum(abs(colSums(NsALS)-c(length(ind_f), length(ind_m))))/2
NsALS[1,1] <- NsALS[1,1]-1 
NsALS[1,2] <- NsALS[1,2]+1 
if(!all(c(colSums(NsALS)==c(length(ind_f), length(ind_m))))) message("Error!")

sample_data$Rand[ind_f] <- rep(Rand_levels, NsALS[,1])
sample_data$Rand[ind_m] <- rep(Rand_levels, NsALS[,2])

rm(NsALS, ind_f, ind_m)

table(sample_data$Rand, sample_data$Case)

NsALSFTD <- sample_data %>% 
  filter(Case=="ALS_FTD") %>% 
  {table(.$Sex)} %>% 
  prop.table  %>% as.numeric() %*%t(as.numeric(Ns[,3])) %>% round() %>% t() 

set.seed(1024) 
ind_f <- sample(
  which(sample_data$Case=="ALS_FTD" & sample_data$Sex=="f"), 
  length(which(sample_data$Case=="ALS_FTD" & sample_data$Sex=="f"))
)
set.seed(1024) 
ind_m <- sample(
  which(sample_data$Case=="ALS_FTD" & sample_data$Sex=="m"), 
  length(which(sample_data$Case=="ALS_FTD" & sample_data$Sex=="m"))
)

if(!all(c(colSums(NsALSFTD)==c(length(ind_f), length(ind_m))))) message("Error!")

sample_data$Rand[ind_f] <- rep(Rand_levels, NsALSFTD[,1])
sample_data$Rand[ind_m] <- rep(Rand_levels, NsALSFTD[,2])

rm(NsALSFTD, ind_f, ind_m)

table(sample_data$Case)
table(sample_data$Rand)

table(sample_data$Case, sample_data$Rand)
table(sample_data$Case, sample_data$Rand) %>% 
  apply(MARGIN=2, FUN=prop.table) %>% 
  "*"(100) %>% 
  round(1) %>% 
  addmargins(margin=1)

table(sample_data$Case, sample_data$Sex)
table(sample_data$Rand, sample_data$Sex)

rm(Rand_levels, Ns)

sample_data$Rand <- factor(sample_data$Rand, levels=c("R1", "R2", "R3"))

M0@misc$Sample_data$Rand1_AllCases <- sample_data$Rand[match(M0@misc$Sample_data$ID, sample_data$ID)]
M0@misc$Sample_data %>% 
  group_by(ID) %>% 
    summarise(Rand=dplyr::first(as.character(Rand1_AllCases))) %>% 
      select(Rand) |> 
      table()

rm(sample_data)




  ### 4.0 Generate Random control grouping for ALS -----------------------------

sample_data <- M0@misc$Sample_metrics 
sample_data$Sex <- M0$Sex[match(sample_data$ID, M0$ID)]
setequal(
  unique(paste0(M0@misc$Sample_data$ID, "_", M0@misc$Sample_data$Sex)), 
  paste0(sample_data$ID, "_", sample_data$Sex)
)
sample_data <- sample_data[,c(1,17:19, 2:16)]

sample_data <- sample_data[sample_data$Case %in% c("ALS", "HC"),]


sample_data$Rand <- "R1" 
table(sample_data$Case)

set.seed(1025)
sample_data$Rand[sample(1:62, 30)] <- "R2"  

table(sample_data$Rand)
table(sample_data$Case, sample_data$Rand) 
table(sample_data$Sex, sample_data$Case)
table(sample_data$Sex, sample_data$Rand) 

sample_data$Rand <- factor(sample_data$Rand, levels=c("R1", "R2"))

M0@misc$Sample_data$Rand1_ALS <- sample_data$Rand[match(M0@misc$Sample_data$ID, sample_data$ID)]

table(M0@misc$Sample_data$Rand1_ALS)
M0@misc$Sample_data %>% 
  group_by(ID) %>% 
    summarize(Rand=dplyr::first(Rand1_ALS)) %>% 
      select(Rand) |> 
        table()

rm(sample_data)




  ### 5.0 Generate Random control grouping for ALSFTD --------------------------

sample_data <- M0@misc$Sample_metrics 
sample_data$Sex <- M0$Sex[match(sample_data$ID, M0$ID)]
setequal(
  unique(paste0(M0@misc$Sample_data$ID, "_", M0@misc$Sample_data$Sex)), 
  paste0(sample_data$ID, "_", sample_data$Sex)
)
sample_data <- sample_data[,c(1,17:19, 2:16)]

sample_data <- sample_data[sample_data$Case %in% c("ALS_FTD", "HC"),]

sample_data$Rand <- "R1" 
table(sample_data$Case)

set.seed(43)
sample_data$Rand[sample(1:49, 17)] <- "R2"  

table(sample_data$Rand)
table(sample_data$Case, sample_data$Rand) 
table(sample_data$Sex, sample_data$Case)
table(sample_data$Sex, sample_data$Rand) 

sample_data$Rand <- factor(sample_data$Rand, levels=c("R1", "R2"))

M0@misc$Sample_data$Rand1_ALSFTD <- sample_data$Rand[match(M0@misc$Sample_data$ID, sample_data$ID)]
table(M0@misc$Sample_data$Rand1_ALSFTD)
M0@misc$Sample_data %>% 
  group_by(ID) %>% 
    summarize(Rand=dplyr::first(Rand1_ALSFTD)) %>% 
      select(Rand) |> 
        table()

rm(sample_data)




  ### 5.0 Generate Random control grouping for ALSFTD_C9 -----------------------

sample_data <- M0@misc$Sample_metrics 
sample_data$Sex <- M0$Sex[match(sample_data$ID, M0$ID)]
setequal(
  unique(paste0(M0@misc$Sample_data$ID, "_", M0@misc$Sample_data$Sex)), 
  paste0(sample_data$ID, "_", sample_data$Sex)
)
sample_data <- sample_data[,c(1,17:19, 2:16)]


sample_data <- sample_data[sample_data$Case_Type %in% c("C9_ALS_FTD", "HC"),]


sample_data$Rand <- "R1" 
table(sample_data$Case)

set.seed(42)
sample_data$Rand[sample(1:32, 7)] <- "R2"  

table(sample_data$Rand)
table(sample_data$Case, sample_data$Rand) 
table(sample_data$Sex, sample_data$Case)
table(sample_data$Sex, sample_data$Rand) 

sample_data$Rand <- factor(sample_data$Rand, levels=c("R1", "R2"))

M0@misc$Sample_data$Rand1_ALSFTD_C9 <- sample_data$Rand[match(M0@misc$Sample_data$ID, sample_data$ID)]
table(M0@misc$Sample_data$Rand1_ALSFTD_C9)

M0@misc$Sample_data %>% 
  group_by(ID) %>% 
    summarize(Rand=dplyr::first(Rand1_ALSFTD_C9)) %>% 
      select(Rand) |> 
        table()

rm(sample_data)




  ### 6.0 Generate Random control grouping for ALSFTD_nonC9 -------------------- 

sample_data <- M0@misc$Sample_metrics 
sample_data$Sex <- M0$Sex[match(sample_data$ID, M0$ID)]
setequal(
  unique(paste0(M0@misc$Sample_data$ID, "_", M0@misc$Sample_data$Sex)), 
  paste0(sample_data$ID, "_", sample_data$Sex)
)
sample_data <- sample_data[,c(1,17:19, 2:16)]

table(sample_data$Case_Type)
sample_data <- sample_data[sample_data$Case_Type %in% c("ALS_FTD", "HC"),]


sample_data$Rand <- "R1" 
table(sample_data$Case)

set.seed(42)
sample_data$Rand[sample(1:42, 10)] <- "R2"  

table(sample_data$Rand)
table(sample_data$Case, sample_data$Rand) 
table(sample_data$Sex, sample_data$Case)
table(sample_data$Sex, sample_data$Rand) 

sample_data$Rand <- factor(sample_data$Rand, levels=c("R1", "R2"))

M0@misc$Sample_data$Rand1_ALSFTD_NonC9 <- sample_data$Rand[match(M0@misc$Sample_data$ID, sample_data$ID)]

table(M0@misc$Sample_data$Rand1_ALSFTD_NonC9)
M0@misc$Sample_data %>% 
  group_by(ID) %>% 
    summarize(Rand=dplyr::first(Rand1_ALSFTD_NonC9)) %>% 
      select(Rand) |> 
        table() 

rm(sample_data)




  ### 7.0 Generate Random control grouping for ALSFTD_C9 vs. ALSFTD_NonC9 ------


sample_data <- M0@misc$Sample_metrics 
sample_data$Sex <- M0$Sex[match(sample_data$ID, M0$ID)]
setequal(
  unique(paste0(M0@misc$Sample_data$ID, "_", M0@misc$Sample_data$Sex)), 
  paste0(sample_data$ID, "_", sample_data$Sex)
)
sample_data <- sample_data[,c(1,17:19, 2:16)] 

sample_data <- sample_data[sample_data$Case_Type %in% c("C9_ALS_FTD", "ALS_FTD"),]

sample_data$Rand <- "R1"
set.seed(31415)
sample_data$Rand[sample(1:17,10)] <- "R2"

table(sample_data$Rand)
table(sample_data$Case_Type, sample_data$Rand) 
table(sample_data$Sex, sample_data$Case_Type)
table(sample_data$Sex, sample_data$Rand) 


sample_data$Rand <- factor(sample_data$Rand, levels=c("R1", "R2"))

M0@misc$Sample_data$Rand1_ALSFTD_C9vsALSFTD_NonC9 <- sample_data$Rand[match(M0@misc$Sample_data$ID, sample_data$ID)]

table(M0@misc$Sample_data$Rand1_ALSFTD_C9vsALSFTD_NonC9)
M0@misc$Sample_data %>% 
  group_by(ID) %>% 
  summarize(Rand=dplyr::first(Rand1_ALSFTD_C9vsALSFTD_NonC9)) %>% 
  select(Rand) |> 
  table() 

rm(sample_data)




  ### 8.0 Export M0 Seurat Object ---------------------------------------------- 

qsave(
  M0, 
  paste0(
    "../Data/SeuratObjects/", 
    "M0", 
    ".qrds"
  ), 
  nthr = nthr
)

                