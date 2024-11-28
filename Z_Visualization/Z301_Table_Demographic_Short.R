

  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 

library(qs)
library(tidyverse)
library(magrittr)


qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)




  ### 1.0 Load data ------------------------------------------------------------

M0_RNA <- qread(
  paste0(
    "../Data/SeuratObjects/", 
    "M0_RNA", 
    ".qrds"
  ), 
  nthr=nthr
)




  ### 2.0 Summarize stats ------------------------------------------------------

M0_RNA@meta.data %>% 
  group_by(ID) %>% 
  summarize(Case=dplyr::first(Case)) %>% 
  select(Case) %>% table()

M0_RNA@meta.data %>% 
  group_by(ID) %>% 
  summarize(Case_Type=dplyr::first(Case_Type)) %>% 
  select(Case_Type) %>% table()

M0_RNA@meta.data %>% 
  group_by(ID) %>% 
  summarize(Case_Type=dplyr::first(Case_Type), Sex=dplyr::first(Sex)) %$% 
  table(Case_Type, Sex) %>% 
  prop.table(margin=1) %>% 
  "*"(100) %>% 
  round(0)

M0_RNA@meta.data %>% 
  group_by(ID) %>% 
  summarize(Case_Type=dplyr::first(Case_Type), Age=dplyr::first(Age)) %>% 
  group_by(Case_Type) %>% 
  summarize(Age_Mean=round(mean(Age), 1), Age_SD=round(sd(Age), 1)) 




