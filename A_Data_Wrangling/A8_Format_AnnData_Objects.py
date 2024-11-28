
  ### 0.0 Import libraries -----------------------------------------------------


import matplotlib.pyplot as plt 
import scanpy as sc
import mudata as mu
import pertpy as pt 
import pandas as pd 
import arviz as az 




  ### 1.0 Load data ------------------------------------------------------------

M0_RNA = sc.read_h5ad("../Data/AnnDataObjects/M0_RNA_sce.h5ad")
M0_RNA

M0_RNA_Meta = pd.read_csv("../Data/Input/M0_RNA_metadata.csv", index_col=0)




  ### 2.0 Replace AnnData.obs with Pd.DataFrame --------------------------------

M0_RNA.obs.Cohort.value_counts(dropna=False)
M0_RNA.obs.Age.value_counts(dropna=False) 

all(M0_RNA.obs.index == M0_RNA_Meta.index)
all(M0_RNA.obs.CellId == M0_RNA_Meta.CellId)

M0_RNA_Meta["ident"] = M0_RNA_Meta.WNN_L4
all(M0_RNA.obs.columns.values == M0_RNA_Meta.columns.values)

M0_RNA.obs = M0_RNA_Meta
M0_RNA.obs.Cohort.value_counts(dropna=False)
M0_RNA.obs.Age.value_counts(dropna=False) 




  ### 3.0 Save AnnData objects -------------------------------------------------

M0_RNA.write(
  filename="../Data/AnnDataObjects/M0_RNA.h5ad"
)










