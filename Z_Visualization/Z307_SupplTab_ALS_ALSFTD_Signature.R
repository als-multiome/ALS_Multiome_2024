
  ###  0.0 Load libraries ------------------------------------------------------ 
source("~/ALS_Brain_Multiome.Rcfg") 

library(qs)
library(tidyverse)

  
qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)

 


  ### 1.0 Load data ------------------------------------------------------------

ALS_ALSFTD_Sign <- qread("../Data/Annotations/Signatures/RNA/Signature_RNA_ALS_ALSFTD_WNN_L25.qrds") 




  ### 2.0 Format table ---------------------------------------------------------

Tab <- ALS_ALSFTD_Sign
Tab <- Tab %>% 
  select(c(Gene, ID_10X, GENETYPE, Signif_In))




  ### 3.0 Export table ---------------------------------------------------------

write.csv(
  Tab, 
  "../Data/Visualization/Tables/SupplTab_ALS_ALSFTD_Signature.csv", 
  quote = FALSE, 
  row.names = FALSE
)
