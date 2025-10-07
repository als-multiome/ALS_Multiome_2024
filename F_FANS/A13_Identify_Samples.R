# A13_Identify_Samples.R 


  ### 0.0 Load libraries ------------------------------------------------------- 

library(qs2)
library(tidyverse)
library(data.table)
library(Seurat)
  
  

        
  ### 1.0 Load data ------------------------------------------------------------
  
Samples <- qs_read(
  "../Data/FANS/Samples/Samples.qs2"
)

Sample_Clusters <- qs_read(
  "../Data/FANS/Samples/Sample_Clusters.qs2"
)

Sample_Clusters$TDP43 <- Samples$TDP43[match(Sample_Clusters$Sample, Samples$ID)]

FANS <- qs_read(
  "../Data/FANS/SeuratObjects/FANS_Merged_QCFiltered_SexAdded.qs2"
)




  ### 3.0 Define auxilliary functions ------------------------------------------

print_genotype_similarity <- function(Subject, Query, QUAL=NULL){
  
  SubjClusts <- as.character(0:(ncol(Subject)-10))
  sapply(
    SubjClusts,
    function(i){
      Subject[[paste0("C", i)]] <<- str_split(Subject[[i]], pattern=":", simplify=TRUE)[,1]
    }
  )
  colnames(Subject) <- str_replace_all(colnames(Subject), "#", "")
  Subject <- Subject %>% 
    select(c("CHROM", "POS", "REF", "ALT", "QUAL"), paste0("C", SubjClusts))
  
  sapply(
    paste0("C", SubjClusts), 
    function(o){
      Subject <<- Subject[Subject[[o]] != "./.",]
    }
  )
  
  
  
  QueryClusts <- as.character(0:(ncol(Query)-10))
  sapply(
    QueryClusts,
    function(i){
      Query[[paste0("C", i)]] <<- str_split(Query[[i]], pattern=":", simplify=TRUE)[,1]
    }
  )
  colnames(Query) <- str_replace_all(colnames(Query), "#", "")
  Query <- Query %>% 
    select(c("CHROM", "POS", "REF", "ALT", "QUAL"), paste0("C", QueryClusts))
  
  sapply(
    paste0("C", QueryClusts), 
    function(o){
      Query <<- Query[Query[[o]] != "./.",]
    }
  )
  
  
  # Identify common SNPs 
  
  Subject$SNP <- paste0(Subject$CHROM, "_", Subject$POS)
  Query$SNP <- paste0(Query$CHROM, "_", Query$POS)
  
  common_SNPs <- intersect(Subject$SNP, Query$SNP)
  Subject <- Subject %>% 
    filter(SNP %in% common_SNPs)
  
  Query <- Query %>% 
    filter(SNP %in% common_SNPs)
  
  all(Subject$SNP == Query$SNP)
  table(Subject$REF == Query$REF)
  
  ind <- which(Subject$REF != Query$REF)
  if(length(ind) > 0){
    Subject <- Subject[-ind,]
    Query <- Query[-ind,]
  }
  
  all(Subject$SNP == Query$SNP)
  table(Subject$ALT == Query$ALT)
  
  ind <- which(Subject$ALT != Query$ALT)
  if(length(ind) > 0){
    Subject <- Subject[-ind,]
    Query <- Query[-ind,]
  }
  
  if(!is.null(QUAL)){
    Subject$QUAL <- as.numeric(Subject$QUAL) 
    Query$QUAL <- as.numeric(Query$QUAL)
    
    ind <- which(Subject$QUAL < QUAL) 
    Subject <- Subject[-ind,]
    Query <- Query[-ind,] 
    
    ind <- which(Query$QUAL < QUAL) 
    Subject <- Subject[-ind,]
    Query <- Query[-ind,]
    
  }
  
  sapply(
    paste0("C", QueryClusts), 
    function(C1){
      tmp <- sapply(
        paste0("C", SubjClusts), 
        function(C0){
          return(
            sum(Subject[[C0]] == Query[[C1]])
          )  
        }
      ) 
      return(
        setNames(
          tmp, 
          0:(length(tmp)-1)
        )
      )
    }
  )    
  
  
}




  ### 3.0 Load VCF Data --------------------------------------------------------


# Read Data 
FC21 <- fread(
  "../Data/FANS/Souporcell_Output/FC21/cluster_genotypes.vcf", 
  skip = "#CHROM"
)

FC22 <- fread(
  "../Data/FANS/Souporcell_Output/FC22/cluster_genotypes.vcf", 
  skip = "#CHROM"
) 

FC23 <- fread(
  "../Data/FANS/Souporcell_Output/FC23/cluster_genotypes.vcf", 
  skip = "#CHROM"
) 
FC24 <- fread(
  "../Data/FANS/Souporcell_Output/FC24/cluster_genotypes.vcf", 
  skip = "#CHROM"
) 

FC25 <- fread(
  "../Data/FANS/Souporcell_Output/FC25/cluster_genotypes.vcf", 
  skip = "#CHROM"
) 

FC26 <- fread(
  "../Data/FANS/Souporcell_Output/FC26/cluster_genotypes.vcf", 
  skip = "#CHROM"
)

FC27 <- fread(
  "../Data/FANS/Souporcell_Output/FC27/cluster_genotypes.vcf", 
  skip = "#CHROM"
) 

FC28 <- fread(
  "../Data/FANS/Souporcell_Output/FC28/cluster_genotypes.vcf", 
  skip = "#CHROM"
)

FC2 <- fread(
  "../Data/FANS/Souporcell_Output/FC2/cluster_genotypes.vcf", 
  skip = "#CHROM"
)

FC3 <- fread(
  "../Data/FANS/Souporcell_Output/FC3/cluster_genotypes.vcf", 
  skip = "#CHROM"
)

FC4 <- fread(
  "../Data/FANS/Souporcell_Output/FC4/cluster_genotypes.vcf", 
  skip = "#CHROM"
)

Chip2Well4 <- fread(
  "../../../../../../Bioinformatic_Data/Cellranger_Output/Human/scMultiom/3_Multiom_M4_FTD/2_Souporcell_GEX_ATAC/GEX/Chip2Well4/cluster_genotypes.vcf", 
  skip = "#CHROM"
)


Chip1Well1 <- fread(
  "../../../../../../Bioinformatic_Data/Cellranger_Output/Human/scMultiom/3_Multiom_M4_FTD/2_Souporcell_GEX_ATAC/GEX/Chip1Well1/cluster_genotypes.vcf", 
  skip = "#CHROM"
) 

Chip1Well2 <- fread(
  "../../../../../../Bioinformatic_Data/Cellranger_Output/Human/scMultiom/3_Multiom_M4_FTD/2_Souporcell_GEX_ATAC/GEX/Chip1Well2/cluster_genotypes.vcf", 
  skip = "#CHROM"
) 

Chip2Well3 <- fread(
  "../../../../../../Bioinformatic_Data/Cellranger_Output/Human/scMultiom/3_Multiom_M4_FTD/2_Souporcell_GEX_ATAC/GEX/Chip2Well3/cluster_genotypes.vcf", 
  skip = "#CHROM"
) 




  ### 4.0 Identify donors and add to Samples data ------------------------------

Sample_Clusters$Donor <- "NA"
Sample_Clusters$Donor[Sample_Clusters$Sample=="FC2"] <- "FANS_1"
Sample_Clusters$Donor[Sample_Clusters$Sample=="FC3"] <- "FANS_2"
Sample_Clusters$Donor[Sample_Clusters$Sample=="FC4"] <- "FANS_3"

Sample_Clusters$Donor[Sample_Clusters$Sample_Cluster=="FC21_0"] <- "FANS_4"
Sample_Clusters$Donor[Sample_Clusters$Sample_Cluster=="FC21_1"] <- "FANS_5"
Sample_Clusters$Donor[Sample_Clusters$Sample_Cluster=="FC21_2"] <- "FANS_6"
Sample_Clusters$Donor[Sample_Clusters$Sample_Cluster=="FC21_3"] <- "FANS_7"
Sample_Clusters$Donor[Sample_Clusters$Sample_Cluster=="FC21_4"] <- "FANS_8"
Sample_Clusters$Donor[Sample_Clusters$Sample_Cluster=="FC21_5"] <- "FANS_9"
Sample_Clusters$Donor[Sample_Clusters$Sample_Cluster=="FC25_0"] <- "FANS_10"
Sample_Clusters$Donor[Sample_Clusters$Sample_Cluster=="FC25_1"] <- "FANS_11"
Sample_Clusters$Donor[Sample_Clusters$Sample_Cluster=="FC27_0"] <- "FANS_12"
Sample_Clusters$Donor[Sample_Clusters$Sample_Cluster=="FC27_1"] <- "FANS_13"

table(
  Sample_Clusters$TDP43, 
  Sample_Clusters$Donor!="NA"
)



    ## 2.1 Identify FC22 Samples -----------------------------------------------

print_genotype_similarity(Subject = FC21, Query = FC22, QUAL = 100)
  # Cluster FC21_0 is Cluster FC22_2, m, donor identified by inspecting known SNPs 
  # Cluster FC22_0 is Cluster FC22_1 ,f, donor identified   
  # Cluster FC21_4 is Cluster FC22_0 and FC22_1  

print_genotype_similarity(
  Query = Chip2Well3, 
  Subject = FC21
)  

print_genotype_similarity(
  Subject = Chip2Well3, 
  Query = FC22
) 

Sample_Clusters$Donor[
  Sample_Clusters$Sample_Cluster=="FC22_2"
] <- Sample_Clusters$Donor[
  Sample_Clusters$Sample_Cluster == "FC21_0"
] 

Sample_Clusters$Donor[
  Sample_Clusters$Sample_Cluster=="FC22_0"
] <- Sample_Clusters$Donor[
  Sample_Clusters$Sample_Cluster == "FC21_4"
]  

Sample_Clusters$Donor[
  Sample_Clusters$Sample_Cluster=="FC22_1"
] <- Sample_Clusters$Donor[
  Sample_Clusters$Sample_Cluster == "FC21_4"
] 

Sample_Clusters$Donor[
  Sample_Clusters$Sample_Cluster=="FC22_0"
] == Sample_Clusters$Donor[
  Sample_Clusters$Sample_Cluster == "FC22_1"
] 



    ## 2.2 Identify FC23 & FC24 Samples ----------------------------------------

print_genotype_similarity(Subject = FC23, Query = FC24, QUAL = 100)
  # Cluster FC23_0 is Cluster FC24_2, m 

print_genotype_similarity(Subject = Chip1Well2, Query = FC24, QUAL = 100)
print_genotype_similarity(Subject = Chip1Well2, Query = FC23, QUAL = 100)

print_genotype_similarity(Subject = FC21, Query = FC23, QUAL = 100)
print_genotype_similarity(Subject = FC21, Query = FC24, QUAL = 100)
  # Cluster FC23_0 is Cluster FC24_2
  # Cluster FC23_0 and FC24_2 are Cluster 21_1
  # Clusters FC23_1, FC23_2, FC24_0 and FC24_1 are Cluster 21_5


print_genotype_similarity(Subject = FC21, Query = Chip1Well2)
  
# Cluster FC23_0 is Cluster FC21_1 
Sample_Clusters$Donor[
  Sample_Clusters$Sample_Cluster=="FC23_0"
] <- Sample_Clusters$Donor[
  Sample_Clusters$Sample_Cluster == "FC21_1"
]

# Cluster FC23_1 is Cluster FC23_2 and Cluster FC21_5  
Sample_Clusters$Donor[
  Sample_Clusters$Sample_Cluster=="FC23_1"
] <- Sample_Clusters$Donor[
  Sample_Clusters$Sample_Cluster == "FC21_5"
]  

Sample_Clusters$Donor[
  Sample_Clusters$Sample_Cluster=="FC23_2"
] <- Sample_Clusters$Donor[
  Sample_Clusters$Sample_Cluster == "FC21_5"
]   


# Cluster FC24_0 is Cluster FC24_1 and Cluster FC21_5 
Sample_Clusters$Donor[
  Sample_Clusters$Sample_Cluster=="FC24_0"
] <- Sample_Clusters$Donor[
  Sample_Clusters$Sample_Cluster == "FC21_5"
]  

Sample_Clusters$Donor[
  Sample_Clusters$Sample_Cluster=="FC24_1"
] <- Sample_Clusters$Donor[
  Sample_Clusters$Sample_Cluster == "FC21_5"
]  

  # Cluster FC24_2 is Cluster FC21_1 
Sample_Clusters$Donor[
  Sample_Clusters$Sample_Cluster=="FC24_2"
] <- Sample_Clusters$Donor[
  Sample_Clusters$Sample_Cluster == "FC21_1"
]



    ## 2.3 Identify FC25 & FC26 Samples  ---------------------------------------

print_genotype_similarity(Subject = FC25, Query = FC26)
print_genotype_similarity(Subject = FC25, Query = FC27)

  # Cluster FC25_0 is Cluster FC26_0 
  # Cluster FC25_1 is Cluster FC26_1 
Sample_Clusters$Donor[
  Sample_Clusters$Sample_Cluster=="FC26_0"
] <- Sample_Clusters$Donor[
  Sample_Clusters$Sample_Cluster == "FC25_0"
] 

Sample_Clusters$Donor[
  Sample_Clusters$Sample_Cluster=="FC26_1"
] <- Sample_Clusters$Donor[
  Sample_Clusters$Sample_Cluster == "FC25_1"
] 

print_genotype_similarity(Subject = Chip2Well4, Query = FC25)
print_genotype_similarity(Subject = Chip1Well1, Query = FC25)



    ## 2.4 Identify FC27 & FC28 ------------------------------------------------

print_genotype_similarity(Subject = FC27, Query = FC28)
print_genotype_similarity(Subject = FC27, Query = FC26)


  # Cluster FC27_0 is Cluster FC28_1 
Sample_Clusters$Donor[
  Sample_Clusters$Sample_Cluster=="FC28_1"
] <- Sample_Clusters$Donor[
  Sample_Clusters$Sample_Cluster == "FC27_0"
] 

  # Cluster FC27_1 is Cluster FC28_0 
Sample_Clusters$Donor[
  Sample_Clusters$Sample_Cluster=="FC28_0"
] <- Sample_Clusters$Donor[
  Sample_Clusters$Sample_Cluster == "FC27_1"
] 




  ### 3.0 Cross-identify samples from M0 Multiome ------------------------------

print_genotype_similarity(
  Subject = Chip1Well2, 
  Query = FC2
)  

print_genotype_similarity(
  Subject = FC27, 
  Query = FC2
)  

print_genotype_similarity(
  Subject = Chip2Well3, 
  Query = FC22
) 

Sample_Clusters$Donor[Sample_Clusters$Donor=="FANS_13"] <- "FANS_1"



print_genotype_similarity(
  Subject = Chip2Well4, 
  Query = FC3
)  

print_genotype_similarity(
  Subject = FC25, 
  Query = FC3
)   

print_genotype_similarity(
  Subject = FC26, 
  Query = FC3
)   

Sample_Clusters$Donor[Sample_Clusters$Donor=="FANS_11"] <- "FANS_2"



print_genotype_similarity(
  Subject = Chip1Well1, 
  Query = FC4
)  

print_genotype_similarity(
  Subject = FC25, 
  Query = FC4
)   

print_genotype_similarity(
  Subject = FC26, 
  Query = FC4
)   

Sample_Clusters$Donor[Sample_Clusters$Donor=="FANS_10"] <- "FANS_3"







  ### 4.0 Double-check ---------------------------------------------------------

unique(Sample_Clusters$Donor) |> length() 

table(
  Sample_Clusters$Donor, 
  Sample_Clusters$Sex
)

table(
  Sample_Clusters$Donor, 
  Sample_Clusters$TDP43
)


  ### 5.0 Add Sample donor data to Seurat object -------------------------------

setequal(FANS$ID_Cluster, Sample_Clusters$Sample_Cluster)

Sample_Clusters$M0_ID <- NA 
Sample_Clusters$M0_ID[Sample_Clusters$Donor=="FANS_1"] <- "FTD3" 
Sample_Clusters$M0_ID[Sample_Clusters$Donor=="FANS_2"] <- "FTD8" 
Sample_Clusters$M0_ID[Sample_Clusters$Donor=="FANS_3"] <- "FTD1" 
Sample_Clusters$M0_ID[Sample_Clusters$Donor=="FANS_4"] <- "FTD5" 
Sample_Clusters$M0_ID[Sample_Clusters$Donor=="FANS_6"] <- "FTD4" 
Sample_Clusters$M0_ID[Sample_Clusters$Donor=="FANS_7"] <- "FTD6" 
Sample_Clusters$M0_ID[Sample_Clusters$Donor=="FANS_8"] <- "FTD2" 



FANS$Sample_donor <- Sample_Clusters$Donor[match(FANS$ID_Cluster, Sample_Clusters$Sample_Cluster)] 

tmp <- Sample_Clusters %>% 
  group_by(Sample) %>% 
  summarize(Donors = paste0(unique(Donor), collapse=","))

setequal(
  Samples$ID, 
  tmp$Sample
)

Samples$Sample_donors <- tmp$Donors[
  match(
    Samples$ID, 
    tmp$Sample
  )  
]

table(
  FANS$Sample_donor, 
  FANS$Sex
)



  ### 6.0 Save data ------------------------------------------------------------

qs_save(
  FANS, 
  "../Data/FANS/SeuratObjects/FANS_Merged_QCFiltered_SexAdded_SampleDonorAdded.qs2"
)

qs_save(
  Sample_Clusters, 
  "../Data/FANS/Samples/Sample_Clusters.qs2"
)

qs_save(
  Samples, 
  "../Data/FANS/Samples/Samples.qs2"
)


