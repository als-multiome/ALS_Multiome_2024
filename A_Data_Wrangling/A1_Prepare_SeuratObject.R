  ### 0.0 Load libraries -------------------------------------------------------

source("~/ALS_Brain_Multiome.Rcfg") 

library(Seurat) 
library(Signac)
library(tidyverse)
library(qs)




  ### 1.0 Load data ------------------------------------------------------------

Seurat_Objects  = setNames(
  read.csv("../Data/cfg/Files_List.txt", header=FALSE, sep="\t")[,2], 
  nm=read.csv("../Data/cfg/Files_List.txt", header=FALSE, sep="\t")[,1]
)

M0 <- qread(
  Seurat_Objects["Merged_Integrated"], 
  nthr=nthr
)




  ### 2.0 Add a CellId (Unique Barcode) ----------------------------------------

M0@meta.data$CellId <- rownames(M0@meta.data)
any(duplicated(M0$CellId))
all(rownames(M0@meta.data)==M0$CellId)
all(Cells(M0)==M0$CellId)
all(colnames(M0@assays$RNA)==M0$CellId)
all(colnames(M0@assays$ATAC)==M0$CellId)

M0@meta.data <- M0@meta.data[,c(ncol(M0@meta.data), 1:(ncol(M0@meta.data)-1))]




  ### 3.0 Strip of redundant cell-wise metadata --------------------------------



    ## 3.1 Move archival cell-type annotations to Misc data --------------------

M0$Idents151223 <- Idents(M0)
M0@misc$Archived_CellTypeAnnotations <- M0@meta.data[
  ,colnames(M0@meta.data) %in% c(
    "CellId", 
    "clusterATAC", 
    "deepsort_BrainAllen", 
    "deepsort_BrainMerged", 
    "deepsort_BrainNature", 
    "Idents0504", 
    "CellType_Level5", 
    "Idents1704", 
    "CellType_Level5_ATAC", 
    "CellType_Level5_Weighted", 
    "CellType_Level4", 
    "CellType_Level3", 
    "CellType_Level2", 
    "CellType_Level1", 
    "CellType_Level4_ATAC", 
    "CellType_Level3_ATAC", 
    "CellType_Level2_ATAC", 
    "CellType_Level1_ATAC", 
    "CellType_Level4_Weighted", 
    "CellType_Level3_Weighted", 
    "CellType_Level2_Weighted", 
    "CellType_Level1_Weighted",  
    "Celltype_Level5_ATACn", 
    "Celltype_Level5_RNAn", 
    "Celltype_Level5_Weightedn", 
    "Idents151223"
  ) 
] 


all(rownames(M0@misc$Archived_CellTypeAnnotations)==Cells(M0))

M0@meta.data <- M0@meta.data[
  ,-which(colnames(M0@meta.data) %in% c(
    "clusterATAC", 
    "deepsort_BrainAllen", 
    "deepsort_BrainMerged", 
    "deepsort_BrainNature", 
    "Idents0504", 
    "Idents1704", 
    "CellType_Level4", 
    "CellType_Level3", 
    "CellType_Level2", 
    "CellType_Level1", 
    "CellType_Level4_ATAC", 
    "CellType_Level3_ATAC", 
    "CellType_Level2_ATAC", 
    "CellType_Level1_ATAC", 
    "CellType_Level4_Weighted", 
    "CellType_Level3_Weighted", 
    "CellType_Level2_Weighted", 
    "CellType_Level1_Weighted"
    )
  )
]



    ## 3.2 Move GEX QC data to Misc data ---------------------------------------

M0@misc$GEX_QC <- M0@meta.data[
  ,colnames(M0@meta.data) %in% c(
    "CellId", 
    "gex_barcode", 
    "excluded_reason", 
    "gex_raw_reads", 
    "gex_mapped_reads", 
    "gex_conf_intergenic_reads",         
    "gex_conf_exonic_reads", 
    "gex_conf_intronic_reads", 
    "gex_conf_exonic_unique_reads", 
    "gex_conf_exonic_antisense_reads", 
    "gex_conf_exonic_dup_reads", 
    "gex_exonic_umis", 
    "gex_conf_intronic_unique_reads", 
    "gex_conf_intronic_antisense_reads", 
    "gex_conf_intronic_dup_reads", 
    "gex_intronic_umis", 
    "gex_conf_txomic_unique_reads", 
    "gex_umis_count", 
    "gex_genes_count", 
    "percent_mito", 
    "percent_ribo", 
    "percent_hb", 
    "percent_plat", 
    "nCount_RNA", 
    "nFeature_RNA", 
    "nCount_RNACopy", 
    "pct_chrY", 
    "pct_chrX", 
    "FlagRemove", 
    "MADFlag"
  )
]

all(rownames(M0@misc$GEX_QC)==Cells(M0))


M0@meta.data <- M0@meta.data[
  ,-which(colnames(M0@meta.data) %in% c(
    "gex_raw_reads", 
    "gex_mapped_reads", 
    "gex_conf_intergenic_reads",         
    "gex_conf_exonic_reads", 
    "gex_conf_intronic_reads", 
    "gex_conf_exonic_unique_reads", 
    "gex_conf_exonic_antisense_reads", 
    "gex_conf_exonic_dup_reads", 
    "gex_exonic_umis", 
    "gex_conf_intronic_unique_reads", 
    "gex_conf_intronic_antisense_reads", 
    "gex_conf_intronic_dup_reads", 
    "gex_intronic_umis", 
    "gex_conf_txomic_unique_reads", 
    "gex_umis_count", 
    "gex_genes_count", 
    "percent_mito", 
    "percent_ribo", 
    "percent_hb", 
    "percent_plat", 
    "pct_chrY", 
    "pct_chrX"
    )
  ) 
]



    ## 3.3 Move ATAC QC data to Misc data --------------------------------------

M0@misc$ATAC_QC <- M0@meta.data[
  ,colnames(M0@meta.data) %in% c(
    "CellId", 
    "atac_barcode", 
    "atac_raw_reads", 
    "atac_unmapped_reads", 
    "atac_lowmapq", 
    "atac_dup_reads", 
    "atac_chimeric_reads", 
    "atac_mitochondrial_reads", 
    "atac_fragments", 
    "atac_TSS_fragments", 
    "atac_peak_region_fragments", 
    "atac_peak_region_cutsites", 
    "fragments_Count", 
    "mononucleosome_Count", 
    "nucleosome_free_Count", 
    "ATACreads_Count", 
    "nFeature_ATAC", 
    "FRiP", 
    "nucleosome_signal", 
    "nucleosome_percentile", 
    "TSS.enrichment", 
    "TSS.percentile", 
    "high.tss", 
    "nucleosome_group", 
    "FlagRemove", 
    "MADFlag"
  )
]

all(rownames(M0@misc$GEX_QC)==Cells(M0)) 

M0@meta.data <- M0@meta.data[
  ,-which(colnames(M0@meta.data) %in% c(
    "atac_raw_reads", 
    "atac_unmapped_reads", 
    "atac_lowmapq", 
    "atac_dup_reads", 
    "atac_chimeric_reads", 
    "atac_mitochondrial_reads", 
    "atac_fragments", 
    "atac_TSS_fragments", 
    "atac_peak_region_fragments", 
    "atac_peak_region_cutsites", 
    "mononucleosome_Count", 
    "nucleosome_free_Count", 
    "ATACreads_Count", 
    "FRiP", 
    "nucleosome_signal", 
    "nucleosome_percentile", 
    "TSS.enrichment", 
    "TSS.percentile", 
    "high.tss", 
    "nucleosome_group"
  )
  ) 
]



    ## 3.4 Move QC filter data to Misc data ------------------------------------

M0@misc$QC <- M0@meta.data[
  ,colnames(M0@meta.data) %in% c(
    "CellId", 
    "is_cell", 
    "excluded_reason", 
    "multiplets", 
    "multipletsATAC", 
    "DoubletFinderClass", 
    "DoubletFinderWeighted", 
    "DoubletFinderScore",  
    "FlagRemove", 
    "MADFlag" 
  )
]

all(rownames(M0@misc$QC)==Cells(M0)) 

M0@meta.data <- M0@meta.data[
  ,-which(colnames(M0@meta.data) %in% c(
    "is_cell", 
    "excluded_reason", 
    "multiplets", 
    "multipletsATAC", 
    "DoubletFinderClass", 
    "DoubletFinderWeighted", 
    "DoubletFinderScore",  
    "FlagRemove", 
    "MADFlag" 
    )
  ) 
]


View(M0@meta.data) 



    ## 3.5 Remove preliminary IDs ----------------------------------------------

M0@meta.data <- M0@meta.data[,-which(colnames(M0@meta.data)=="ID")]
M0@meta.data <- M0@meta.data[,-which(colnames(M0@meta.data)=="orig.ident")]




  ### 4.0 Save M0 object -------------------------------------------------------

qsave(
  M0, 
  paste0(
    "../Data/SeuratObjects/", 
    "M0", 
    ".qrds"
  ), 
  nthr=nthr
) 



