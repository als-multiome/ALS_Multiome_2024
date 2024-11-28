

  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 

library(tidyverse)
library(data.table)


qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)




  ### 1.0 Load data ------------------------------------------------------------

sample_data <- fread(
  paste0(
    "../Data/Input/", 
    "Sample_data",
    ".txt"
  ), 
  data.table=FALSE
)




  ### 2.0 Summarize Well Layout ------------------------------------------------

sample_data$Sort <- paste0(sample_data$Well, "_", sample_data$ID) 

sample_data <- sample_data[str_order(sample_data$Sort, numeric=TRUE), ]
sample_data$Chip <- str_split(sample_data$Well, "Well", simplify=TRUE)[,1]
sample_data$Well <- str_split(sample_data$Well, "Well", simplify=TRUE)[,2]
sample_data$SampleGroup <- str_replace_all(sample_data$SampleGroup, "_", "-")
sample_data$SampleGroup <- str_replace_all(sample_data$SampleGroup, "ALSFTD", "ALS-FTD")
sample_data$SampleGroup[sample_data$SampleGroup=="ALS-FTD" & sample_data$C9ORF72==1] <- "C9 ALS-FTD"

sample_data <- sample_data %>% 
  rename("Disease group" = SampleGroup) %>% 
  select(Chip, Well, "Disease group", ID, Sex)




  ### 3.0 Descriptive statistics -----------------------------------------------

table(table(sample_data$ID) > 1)
  



  ### 4.0 Export Data ----------------------------------------------------------

write.table(
  sample_data, 
  paste0(
    "../Data/Visualization/Tables/", 
    "Well_Layout", 
    ".txt"
  ), quote = FALSE, 
  row.names = FALSE, 
  sep = "\t"
)



