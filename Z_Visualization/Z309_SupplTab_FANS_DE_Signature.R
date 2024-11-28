
  ###  0.0 Load libraries ------------------------------------------------------ 
source("~/ALS_Brain_Multiome.Rcfg") 

library(qs)
library(tidyverse)

  
qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)

 


  ### 1.0 Load data ------------------------------------------------------------

FANS_DE_TDP43_MAST <- qread("../Data/DE/FANS/FANS_DE_TDP43_MAST.qrds", nthr) 




  ### 2.0 Format table ---------------------------------------------------------


Tab <- FANS_DE_TDP43_MAST
Tab <- Tab %>% 
  select(c(gene, p_value, model_log2FC, ci.hi, ci.lo, fdr)) %>% 
  filter(fdr < 0.05)




  ### 3.0 Export table ---------------------------------------------------------

write.csv(
  Tab, 
  "../Data/Visualization/Tables/SupplTab_FANS_TDP43_Signature.csv", 
  quote = FALSE, 
  row.names = FALSE
)
