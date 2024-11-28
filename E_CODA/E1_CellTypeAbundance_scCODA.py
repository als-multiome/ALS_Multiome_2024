###


  ### 0.0 Import libraries -----------------------------------------------------

import matplotlib.pyplot as plt 
import scanpy as sc
import mudata as mu
import pertpy as pt 
import pandas as pd 
import arviz as az 



  ### 1.0 Load data ------------------------------------------------------------



    ## 1.1 M0_RNA AnnData Object -----------------------------------------------

M0_RNA = sc.read_h5ad("../Data/AnnDataObjects/M0_RNA.h5ad")
M0_RNA




  ### 2.0 scCODA WNN -----------------------------------------------------------



    ## 2.1 WNN_L1 --------------------------------------------------------------

sccoda_WNN_L1_model = pt.tl.Sccoda()
sccoda_WNN_L1_data = sccoda_WNN_L1_model.load(
    M0_RNA, 
    type="cell_level", 
    generate_sample_level=True, 
    cell_type_identifier="WNN_L1", 
    sample_identifier="ID", 
    covariate_obs=["Case", "Sex"],
)

print(sccoda_WNN_L1_data)
print(sccoda_WNN_L1_data["coda"].X)
print(sccoda_WNN_L1_data["coda"].obs)



    ## 2.1 ALSvsHC -------------------------------------------------------------

sccoda_WNN_L1_data.mod["coda_ALS"] = sccoda_WNN_L1_data["coda"][
    sccoda_WNN_L1_data["coda"].obs["Case"].isin(["ALS", "HC"])
].copy()
print(sccoda_WNN_L1_data["coda_ALS"])

sccoda_WNN_L1_model.plot_boxplots(sccoda_WNN_L1_data, modality_key="coda", feature_name="Case", add_dots=False)
plt.show()

sccoda_WNN_L1_model.plot_stacked_barplot(sccoda_WNN_L1_data, modality_key="coda", feature_name="Case")
plt.show()
sccoda_WNN_L1_model.plot_stacked_barplot(sccoda_WNN_L1_data, modality_key="coda", feature_name="Sex") 
plt.show()
sccoda_WNN_L1_model.plot_stacked_barplot(sccoda_WNN_L1_data, modality_key="coda", feature_name="Cohort")
plt.show()

sccoda_WNN_L1_model.plot_boxplots(sccoda_WNN_L1_data, modality_key="coda_ALS", feature_name="Sex", add_dots=False)
plt.show()

sccoda_WNN_L1_model.plot_boxplots(sccoda_WNN_L1_data, modality_key="coda_ALS", feature_name="Cohort", add_dots=False)
plt.show()

sccoda_WNN_L1_data = sccoda_WNN_L1_model.prepare(
    sccoda_WNN_L1_data,
    modality_key="coda",
    formula="C(Case, Treatment('HC'))", 
    reference_cell_type="Non-neuronal",
)
sccoda_WNN_L1_data["coda"]

sccoda_WNN_L1_model.run_nuts(sccoda_WNN_L1_data, modality_key="coda")
sccoda_WNN_L1_data["coda"] 

sccoda_WNN_L1_model.set_fdr(sccoda_WNN_L1_data, modality_key="coda", est_fdr=0.1)
sccoda_WNN_L1_model.summary(sccoda_WNN_L1_data, modality_key="coda") 



###############

CellTypes = sccoda_WNN_L1_data["coda_ALS"].var.index
CellTypes_results = pd.DataFrame(index=CellTypes, columns=["Times_credible"]).fillna(0)

for celltype in CellTypes:
    print(f"Cell type used as reference: {celltype}")

    sccoda_WNN_L1_data = sccoda_WNN_L1_model.prepare(
        sccoda_WNN_L1_data,
        modality_key="coda_ALS",
        formula="Case",
        reference_cell_type=celltype,
    )
    sccoda_WNN_L1_model.run_nuts(sccoda_WNN_L1_data, modality_key="coda_ALS")

    effects = sccoda_WNN_L1_model.credible_effects(sccoda_WNN_L1_data, modality_key="coda_ALS")
    effects.index = effects.index.droplevel(level=0)

    CellTypes_results["Times_credible"] += effects.astype("int")



CellTypes_results["Pct_credible"] = CellTypes_results["Times_credible"] / len(cell_types)
CellTypes_results["Is_credible"] = CellTypes_results["Pct_credible"] > 0.5


  ### 3.0 WNN_L25 --------------------------------------------------------------

sccoda_WNN_L25_model = pt.tl.Sccoda()
sccoda_WNN_L25_data = sccoda_WNN_L25_model.load(
    M0_RNA, 
    type="cell_level", 
    generate_sample_level=True, 
    cell_type_identifier="WNN_L25", 
    sample_identifier="ID", 
    covariate_obs=["Case" , "Case_Type", "Sex"],
)

print(sccoda_WNN_L25_data)
print(sccoda_WNN_L25_data["coda"].X)
print(sccoda_WNN_L25_data["coda"].obs)



    ## 2.1 ALSvsHC -------------------------------------------------------------

sccoda_WNN_L25_data.mod["coda_ALS"] = sccoda_WNN_L25_data["coda"][
    sccoda_WNN_L25_data["coda"].obs["Case"].isin(["ALS", "HC"])
].copy()
print(sccoda_WNN_L25_data["coda_ALS"])

sccoda_WNN_L25_model.plot_boxplots(sccoda_WNN_L25_data, modality_key="coda", feature_name="Case", add_dots=False)
plt.show()

sccoda_WNN_L25_model.plot_stacked_barplot(sccoda_WNN_L25_data, modality_key="coda", feature_name="Case")
plt.show()
sccoda_WNN_L25_model.plot_stacked_barplot(sccoda_WNN_L25_data, modality_key="coda", feature_name="Sex") 
plt.show()
sccoda_WNN_L25_model.plot_stacked_barplot(sccoda_WNN_L25_data, modality_key="coda", feature_name="Cohort")
plt.show()

sccoda_WNN_L25_model.plot_boxplots(sccoda_WNN_L25_data, modality_key="coda_ALS", feature_name="Sex", add_dots=False)
plt.show()

sccoda_WNN_L25_model.plot_boxplots(sccoda_WNN_L25_data, modality_key="coda_ALS", feature_name="Case", add_dots=False)
plt.show()

sccoda_WNN_L25_data = sccoda_WNN_L25_model.prepare(
    sccoda_WNN_L25_data,
    modality_key="coda",
    formula="C(Case, Treatment('HC'))"
)
sccoda_WNN_L25_data["coda"]

sccoda_WNN_L25_model.run_nuts(sccoda_WNN_L25_data, modality_key="coda")
sccoda_WNN_L25_data["coda"] 

sccoda_WNN_L1_model.set_fdr(sccoda_WNN_L1_data, modality_key="coda", est_fdr=0.1)
sccoda_WNN_L1_model.summary(sccoda_WNN_L1_data, modality_key="coda") 



###############

CellTypes = sccoda_WNN_L1_data["coda_ALS"].var.index
CellTypes_results = pd.DataFrame(index=CellTypes, columns=["Times_credible"]).fillna(0)

for celltype in CellTypes:
    print(f"Cell type used as reference: {celltype}")

    sccoda_WNN_L1_data = sccoda_WNN_L1_model.prepare(
        sccoda_WNN_L1_data,
        modality_key="coda_ALS",
        formula="Case",
        reference_cell_type=celltype,
    )
    sccoda_WNN_L1_model.run_nuts(sccoda_WNN_L1_data, modality_key="coda_ALS")

    effects = sccoda_WNN_L1_model.credible_effects(sccoda_WNN_L1_data, modality_key="coda_ALS")
    effects.index = effects.index.droplevel(level=0)

    CellTypes_results["Times_credible"] += effects.astype("int")



CellTypes_results["Pct_credible"] = CellTypes_results["Times_credible"] / len(cell_types)
CellTypes_results["Is_credible"] = CellTypes_results["Pct_credible"] > 0.5


  ### 4.0 sALS-FTD vs C9 ALS-FTD -----------------------------------------------

sccoda_WNN_L25_data.mod["coda_C9vsNonC9"] = sccoda_WNN_L25_data["coda"][
    sccoda_WNN_L25_data["coda"].obs["Case_Type"].isin(["C9_ALS_FTD", "ALS_FTD"])
].copy()
print(sccoda_WNN_L25_data["coda_C9vsNonC9"])

sccoda_WNN_L25_model.plot_boxplots(sccoda_WNN_L25_data, modality_key="coda_C9vsNonC9", feature_name="Case_Type", add_dots=False)
plt.show()

sccoda_WNN_L25_model.plot_boxplots(sccoda_WNN_L25_data, modality_key="coda_C9vsNonC9", feature_name="Sex", add_dots=False)
plt.show()

sccoda_WNN_L25_data = sccoda_WNN_L25_model.prepare(
    sccoda_WNN_L25_data,
    modality_key="coda_C9vsNonC9",
    formula="Case_Type + Sex"
)
sccoda_WNN_L25_data["coda_C9vsNonC9"]

sccoda_WNN_L25_model.run_nuts(sccoda_WNN_L25_data, modality_key="coda_C9vsNonC9")
sccoda_WNN_L25_data["coda_C9vsNonC9"] 

sccoda_WNN_L25_model.set_fdr(sccoda_WNN_L25_data, modality_key="coda_C9vsNonC9", est_fdr=0.1)
sccoda_WNN_L25_model.summary(sccoda_WNN_L25_data, modality_key="coda_C9vsNonC9") 

### All cell types as reference levels 

# Run scCODA with each cell type as the reference
cell_types = sccoda_WNN_L25_data["coda_C9vsNonC9"] .var.index
results_cycle = pd.DataFrame(index=cell_types, columns=["Times_credible_Case_Type", "Times_credible_Sex"]).fillna(0)

for ct in cell_types:
    print(f"Reference: {ct}")
    sccoda_WNN_L25_data = sccoda_WNN_L25_model.prepare(sccoda_WNN_L25_data,modality_key="coda_C9vsNonC9",formula="Case_Type + Sex",reference_cell_type=ct)
      
    sccoda_WNN_L25_model.run_nuts(sccoda_WNN_L25_data, modality_key="coda_C9vsNonC9")
    cred_eff = sccoda_WNN_L25_model.credible_effects(sccoda_WNN_L25_data, modality_key="coda_C9vsNonC9")
    cred_eff.index = cred_eff.index.droplevel(level=0)
    results_cycle["Times_credible_Case_Type"] += cred_eff.astype("int").iloc[0:12]
    results_cycle["Times_credible_Sex"] += cred_eff.astype("int").iloc[12:24]
    results_cycle.to_csv("../Data/Coda/scCODA/C9_ALS_FTD_vs_ALS_FTD/WNN_L25_CaseType_Sex_AllReferences.csv")


  # Calculate percentages
results_cycle["Pct_credible_Case_Type"] = results_cycle["Times_credible_Case_Type"] / len(cell_types)
results_cycle["Pct_credible_Sex"] = results_cycle["Times_credible_Sex"] / len(cell_types)
results_cycle["Is_credible_Case_Type"] = results_cycle["Pct_credible_Case_Type"] > 0.5
results_cycle["Is_credible_Sex"] = results_cycle["Pct_credible_Sex"] > 0.5

results_cycle.to_csv("../Data/Coda/scCODA/C9_ALS_FTD_vs_ALS_FTD/WNN_L25_CaseType_Sex_AllReferences.csv")
print(results_cycle)



  ### 5.0 ALS vs HC ------------------------------------------------------------

sccoda_WNN_L25_data.mod["coda_ALS"] = sccoda_WNN_L25_data["coda"][
    sccoda_WNN_L25_data["coda"].obs["Case"].isin(["ALS", "HC"])
].copy()
print(sccoda_WNN_L25_data["coda_ALS"])

sccoda_WNN_L25_model.plot_boxplots(sccoda_WNN_L25_data, modality_key="coda_ALS", feature_name="Case", add_dots=False)
plt.show()

sccoda_WNN_L25_model.plot_boxplots(sccoda_WNN_L25_data, modality_key="coda_ALS", feature_name="Sex", add_dots=False)
plt.show()

sccoda_WNN_L25_data = sccoda_WNN_L25_model.prepare(
    sccoda_WNN_L25_data,
    modality_key="coda_ALS",
    formula="Case + Sex"
)

sccoda_WNN_L25_data["coda_ALS"]


### All cell types as reference levels 

# Run scCODA with each cell type as the reference
cell_types = sccoda_WNN_L25_data["coda_ALS"] .var.index
results_cycle = pd.DataFrame(index=cell_types, columns=["Times_credible_Case", "Times_credible_Sex"]).fillna(0)

for ct in cell_types:
    print(f"Reference: {ct}")
    sccoda_WNN_L25_data = sccoda_WNN_L25_model.prepare(sccoda_WNN_L25_data,modality_key="coda_ALS",formula="Case + Sex",reference_cell_type=ct)
      
    sccoda_WNN_L25_model.run_nuts(sccoda_WNN_L25_data, modality_key="coda_ALS")
    cred_eff = sccoda_WNN_L25_model.credible_effects(sccoda_WNN_L25_data, modality_key="coda_ALS")
    cred_eff.index = cred_eff.index.droplevel(level=0)
    results_cycle["Times_credible_Case"] += cred_eff.astype("int").iloc[0:12]
    results_cycle["Times_credible_Sex"] += cred_eff.astype("int").iloc[12:24]
    results_cycle.to_csv("../Data/Coda/scCODA/ALS_vs_HC/WNN_L25_CaseType_Sex_AllReferences.csv")


  # Calculate percentages
results_cycle["Pct_credible_Case"] = results_cycle["Times_credible_Case"] / len(cell_types)
results_cycle["Pct_credible_Sex"] = results_cycle["Times_credible_Sex"] / len(cell_types)
results_cycle["Is_credible_Case"] = results_cycle["Pct_credible_Case"] > 0.5
results_cycle["Is_credible_Sex"] = results_cycle["Pct_credible_Sex"] > 0.5

results_cycle.to_csv("../Data/Coda/scCODA/ALS_vs_HC/WNN_L25_CaseType_Sex_AllReferences.csv")
print(results_cycle)

del(results_cycle, ct, cell_types) 

  ### 6.0 ALS_FTD vs HC --------------------------------------------------------
  
sccoda_WNN_L25_data.mod["coda_ALS_FTD"] = sccoda_WNN_L25_data["coda"][
    sccoda_WNN_L25_data["coda"].obs["Case"].isin(["ALS_FTD", "HC"])
].copy()
print(sccoda_WNN_L25_data["coda_ALS_FTD"])

sccoda_WNN_L25_model.plot_boxplots(sccoda_WNN_L25_data, modality_key="coda_ALS_FTD", feature_name="Case", add_dots=False)
plt.show()

sccoda_WNN_L25_model.plot_boxplots(sccoda_WNN_L25_data, modality_key="coda_ALS_FTD", feature_name="Sex", add_dots=False)
plt.show()

sccoda_WNN_L25_data = sccoda_WNN_L25_model.prepare(
    sccoda_WNN_L25_data,
    modality_key="coda_ALS_FTD",
    formula="Case + Sex"
)

sccoda_WNN_L25_data["coda_ALS_FTD"]


### All cell types as reference levels 

# Run scCODA with each cell type as the reference
cell_types = sccoda_WNN_L25_data["coda_ALS_FTD"] .var.index
results_cycle = pd.DataFrame(index=cell_types, columns=["Times_credible_Case", "Times_credible_Sex"]).fillna(0)

for ct in cell_types:
    print(f"Reference: {ct}")
    sccoda_WNN_L25_data = sccoda_WNN_L25_model.prepare(sccoda_WNN_L25_data,modality_key="coda_ALS_FTD",formula="Case + Sex",reference_cell_type=ct)
      
    sccoda_WNN_L25_model.run_nuts(sccoda_WNN_L25_data, modality_key="coda_ALS_FTD")
    cred_eff = sccoda_WNN_L25_model.credible_effects(sccoda_WNN_L25_data, modality_key="coda_ALS_FTD")
    cred_eff.index = cred_eff.index.droplevel(level=0)
    results_cycle["Times_credible_Case"] += cred_eff.astype("int").iloc[0:12]
    results_cycle["Times_credible_Sex"] += cred_eff.astype("int").iloc[12:24]
    results_cycle.to_csv("../Data/Coda/scCODA/ALS_FTD_vs_HC/WNN_L25_Case_Sex_AllReferences.csv")


  # Calculate percentages
results_cycle["Pct_credible_Case"] = results_cycle["Times_credible_Case"] / len(cell_types)
results_cycle["Pct_credible_Sex"] = results_cycle["Times_credible_Sex"] / len(cell_types)
results_cycle["Is_credible_Case"] = results_cycle["Pct_credible_Case"] > 0.5
results_cycle["Is_credible_Sex"] = results_cycle["Pct_credible_Sex"] > 0.5

results_cycle.to_csv("../Data/Coda/scCODA/ALS_FTD_vs_HC/WNN_L25_Case_Sex_AllReferences.csv")
print(results_cycle)

del(results_cycle, ct, cell_types) 


  ### 7.0 ALS_FTD vs ALS ---------------------------------------------------------
  
sccoda_WNN_L25_data.mod["coda_ALS_ALS_FTD"] = sccoda_WNN_L25_data["coda"][
    sccoda_WNN_L25_data["coda"].obs["Case"].isin(["ALS_FTD", "ALS"])
].copy()
print(sccoda_WNN_L25_data["coda_ALS_ALS_FTD"])

sccoda_WNN_L25_model.plot_boxplots(sccoda_WNN_L25_data, modality_key="coda_ALS_ALS_FTD", feature_name="Case", add_dots=False)
plt.show()

sccoda_WNN_L25_model.plot_boxplots(sccoda_WNN_L25_data, modality_key="coda_ALS_ALS_FTD", feature_name="Sex", add_dots=False)
plt.show()

sccoda_WNN_L25_data = sccoda_WNN_L25_model.prepare(
    sccoda_WNN_L25_data,
    modality_key="coda_ALS_ALS_FTD",
    formula="Case + Sex"
)

sccoda_WNN_L25_data["coda_ALS_ALS_FTD"]


### All cell types as reference levels 

# Run scCODA with each cell type as the reference
cell_types = sccoda_WNN_L25_data["coda_ALS_ALS_FTD"] .var.index
results_cycle = pd.DataFrame(index=cell_types, columns=["Times_credible_Case", "Times_credible_Sex"]).fillna(0)

for ct in cell_types:
    print(f"Reference: {ct}")
    sccoda_WNN_L25_data = sccoda_WNN_L25_model.prepare(sccoda_WNN_L25_data,modality_key="coda_ALS_ALS_FTD",formula="Case + Sex",reference_cell_type=ct)
      
    sccoda_WNN_L25_model.run_nuts(sccoda_WNN_L25_data, modality_key="coda_ALS_ALS_FTD")
    cred_eff = sccoda_WNN_L25_model.credible_effects(sccoda_WNN_L25_data, modality_key="coda_ALS_ALS_FTD")
    cred_eff.index = cred_eff.index.droplevel(level=0)
    results_cycle["Times_credible_Case"] += cred_eff.astype("int").iloc[0:12]
    results_cycle["Times_credible_Sex"] += cred_eff.astype("int").iloc[12:24]
    results_cycle.to_csv("../Data/Coda/scCODA/ALS_FTD_vs_ALS/WNN_L25_Case_Sex_AllReferences.csv")


  # Calculate percentages
results_cycle["Pct_credible_Case"] = results_cycle["Times_credible_Case"] / len(cell_types)
results_cycle["Pct_credible_Sex"] = results_cycle["Times_credible_Sex"] / len(cell_types)
results_cycle["Is_credible_Case"] = results_cycle["Pct_credible_Case"] > 0.5
results_cycle["Is_credible_Sex"] = results_cycle["Pct_credible_Sex"] > 0.5

results_cycle.to_csv("../Data/Coda/scCODA/ALS_FTD_vs_ALS/WNN_L25_Case_Sex_AllReferences.csv")
print(results_cycle)

