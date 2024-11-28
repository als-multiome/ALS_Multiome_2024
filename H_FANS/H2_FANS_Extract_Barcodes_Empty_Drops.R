library(DropletUtils) 
library(Matrix) 
library(Seurat)
raw = Read10X(
  "./Input/raw_feature_bc_matrix/"
)

EmptyDrops <- emptyDrops(
  raw
)
EmptyDrops$is.cell=!is.na(EmptyDrops$FDR)
write.table(
  EmptyDrops@rownames[EmptyDrops$is.cell], 
  file="./Input/EmptyDrops_barcodes.tsv", 
  sep="\t", 
  quote=FALSE, 
  row.names=FALSE,
  col.names=FALSE
)
