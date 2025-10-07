# A2_CellBender_Background_Removal.sh 

# All analyses run in the CellBender Docker container 
# us.gcr.io/broad-dsde-methods/cellbender   latest    56439f37d58e 
# Cellbender v0.3.0 


# FC1 
cp ~/FC1/outs/raw_feature_bc_matrix.h5 ./

cellbender remove-background  \
 --input ./raw_feature_bc_matrix.h5 \
 --output FC1_corrected_output.h5 \ 
 --expected-cells 92 \ 
 --total-droplets-included 650  
 --cuda

mv ./* ../Data/FANS/CellBender_Output/FC1/ 



# FC2 
cp ~/FC2/outs/raw_feature_bc_matrix.h5 ./

cellbender remove-background  \
 --input ./raw_feature_bc_matrix.h5 \
 --output FC2_corrected_output.h5 \ 
 --expected-cells 919 \ 
 --total-droplets-included 2500  
 --cuda

mv ./* ../Data/FANS/CellBender_Output/FC2/  


# FC3 
cp ~/FC3/outs/raw_feature_bc_matrix.h5 ./

cellbender remove-background  \
 --input ./raw_feature_bc_matrix.h5 \
 --output FC3_corrected_output.h5 \ 
 --expected-cells 2945 \ 
 --total-droplets-included 6000  
 --cuda

mv ./* ../Data/FANS/CellBender_Output/FC3/  



# FC4 
cp ~/FC4/outs/raw_feature_bc_matrix.h5 ./

cellbender remove-background  \
 --input ./raw_feature_bc_matrix.h5 \
 --output FC4_corrected_output.h5 \ 
 --expected-cells 1519 \ 
 --total-droplets-included 3000  
 --cuda

mv ./* ../Data/FANS/CellBender_Output/FC4/  



# FC6 
cp ~/FC6/outs/raw_feature_bc_matrix.h5 ./

cellbender remove-background  \
 --input ./raw_feature_bc_matrix.h5 \
 --output FC6_corrected_output.h5 \ 
 --expected-cells 240 \ 
 --total-droplets-included 1250   
 --cuda

mv ./* ../Data/FANS/CellBender_Output/FC6/  


# FC21 
cp ~/P1_2_TDP43hi/outs/raw_feature_bc_matrix.h5 ./

cellbender remove-background  \
 --input ./raw_feature_bc_matrix.h5 \
 --output FC21_corrected_output.h5 \ 
 --expected-cells 4068 \ 
 --total-droplets-included 10000   
 --cuda

mv ./* ../Data/FANS/CellBender_Output/FC21/  



# FC22 
cp ~/P1_TDP43low/outs/raw_feature_bc_matrix.h5 ./

cellbender remove-background  \
 --input ./raw_feature_bc_matrix.h5 \
 --output FC22_corrected_output.h5 \
 --expected-cells 209 \
 --total-droplets-included 1250 \
 --cuda 

mv ./* ../Data/FANS/CellBender_Output/FC22/  



# FC23 
cp ~/P2_1_TDP43lo/outs/raw_feature_bc_matrix.h5 ./

cellbender remove-background  \
 --input ./raw_feature_bc_matrix.h5 \
 --output FC23_corrected_output.h5 \
 --expected-cells 1192 \
 --total-droplets-included 3200 \
 --cuda 

mv ./* ../Data/FANS/CellBender_Output/FC23/ 



# FC24 
cp ~/P2_2_TDP43lo/outs/raw_feature_bc_matrix.h5 ./

cellbender remove-background  \
 --input ./raw_feature_bc_matrix.h5 \
 --output FC24_corrected_output.h5 \
 --expected-cells 1156 \
 --total-droplets-included 3500 \
 --cuda 

mv ./* ../Data/FANS/CellBender_Output/FC24/   



# FC25 
cp ~/P3_TDP43hi/outs/raw_feature_bc_matrix.h5 ./

cellbender remove-background  \
 --input ./raw_feature_bc_matrix.h5 \
 --output FC25_corrected_output.h5 \
 --expected-cells 1415 \
 --total-droplets-included 4200 \
 --cuda 

mv ./* ../Data/FANS/CellBender_Output/FC25/   



# FC26 
cp ~/P3_TDP43lo/outs/raw_feature_bc_matrix.h5 ./

cellbender remove-background  \
 --input ./raw_feature_bc_matrix.h5 \
 --output FC26_corrected_output.h5 \
 --expected-cells 853 \
 --total-droplets-included 3200 \
 --cuda 

mv ./* ../Data/FANS/CellBender_Output/FC26/   



# FC27 
cp ~/P4_TDP43hi/outs/raw_feature_bc_matrix.h5 ./

cellbender remove-background  \
 --input ./raw_feature_bc_matrix.h5 \
 --output FC27_corrected_output.h5 \
 --expected-cells 1642 \
 --total-droplets-included 12000 \
 --cuda 

mv ./* ../Data/FANS/CellBender_Output/FC27/   



# FC28 
cp ~/P4_TDP43lo/outs/raw_feature_bc_matrix.h5 ./

cellbender remove-background  \
 --input ./raw_feature_bc_matrix.h5 \
 --output FC28_corrected_output.h5 \
 --expected-cells 2620 \
 --total-droplets-included 12000 \
 --cuda 

mv ./* ../Data/FANS/CellBender_Output/FC28/   


# Re-pack H5 files for compatibility with Seurat v5 

ptrepack --complevel 5 ./FC1/FC1_corrected_output.h5:/matrix ./FC1/FC1_corrected_output_Seurat.h5:/matrix
ptrepack --complevel 5 ./FC2/FC2_corrected_output.h5:/matrix ./FC2/FC2_corrected_output_Seurat.h5:/matrix 
ptrepack --complevel 5 ./FC3/FC3_corrected_output.h5:/matrix ./FC3/FC3_corrected_output_Seurat.h5:/matrix
ptrepack --complevel 5 ./FC4/FC4_corrected_output.h5:/matrix ./FC4/FC4_corrected_output_Seurat.h5:/matrix
ptrepack --complevel 5 ./FC6/FC6_corrected_output.h5:/matrix ./FC6/FC6_corrected_output_Seurat.h5:/matrix 
ptrepack --complevel 5 ./FC21/FC21_corrected_output.h5:/matrix ./FC21/FC21_corrected_output_Seurat.h5:/matrix 
ptrepack --complevel 5 ./FC22/FC22_corrected_output.h5:/matrix ./FC22/FC22_corrected_output_Seurat.h5:/matrix 
ptrepack --complevel 5 ./FC23/FC23_corrected_output.h5:/matrix ./FC23/FC23_corrected_output_Seurat.h5:/matrix 
ptrepack --complevel 5 ./FC24/FC24_corrected_output.h5:/matrix ./FC24/FC24_corrected_output_Seurat.h5:/matrix 
ptrepack --complevel 5 ./FC25/FC25_corrected_output.h5:/matrix ./FC25/FC25_corrected_output_Seurat.h5:/matrix 
ptrepack --complevel 5 ./FC26/FC26_corrected_output.h5:/matrix ./FC26/FC26_corrected_output_Seurat.h5:/matrix 
ptrepack --complevel 5 ./FC27/FC27_corrected_output.h5:/matrix ./FC27/FC27_corrected_output_Seurat.h5:/matrix 
ptrepack --complevel 5 ./FC28/FC28_corrected_output.h5:/matrix ./FC28/FC28_corrected_output_Seurat.h5:/matrix 

ptrepack --complevel 5 ./FC1/FC1_corrected_output_filtered.h5:/matrix ./FC1/FC1_corrected_output_filtered_Seurat.h5:/matrix
ptrepack --complevel 5 ./FC2/FC2_corrected_output_filtered.h5:/matrix ./FC2/FC2_corrected_output_filtered_Seurat.h5:/matrix 
ptrepack --complevel 5 ./FC3/FC3_corrected_output_filtered.h5:/matrix ./FC3/FC3_corrected_output_filtered_Seurat.h5:/matrix
ptrepack --complevel 5 ./FC4/FC4_corrected_output_filtered.h5:/matrix ./FC4/FC4_corrected_output_filtered_Seurat.h5:/matrix
ptrepack --complevel 5 ./FC6/FC6_corrected_output_filtered.h5:/matrix ./FC6/FC6_corrected_output_filtered_Seurat.h5:/matrix 
ptrepack --complevel 5 ./FC21/FC21_corrected_output_filtered.h5:/matrix ./FC21/FC21_corrected_output_filtered_Seurat.h5:/matrix 
ptrepack --complevel 5 ./FC22/FC22_corrected_output_filtered.h5:/matrix ./FC22/FC22_corrected_output_filtered_Seurat.h5:/matrix 
ptrepack --complevel 5 ./FC23/FC23_corrected_output_filtered.h5:/matrix ./FC23/FC23_corrected_output_filtered_Seurat.h5:/matrix 
ptrepack --complevel 5 ./FC24/FC24_corrected_output_filtered.h5:/matrix ./FC24/FC24_corrected_output_filtered_Seurat.h5:/matrix 
ptrepack --complevel 5 ./FC25/FC25_corrected_output_filtered.h5:/matrix ./FC25/FC25_corrected_output_filtered_Seurat.h5:/matrix 
ptrepack --complevel 5 ./FC26/FC26_corrected_output_filtered.h5:/matrix ./FC26/FC26_corrected_output_filtered_Seurat.h5:/matrix 
ptrepack --complevel 5 ./FC27/FC27_corrected_output_filtered.h5:/matrix ./FC27/FC27_corrected_output_filtered_Seurat.h5:/matrix 
ptrepack --complevel 5 ./FC28/FC28_corrected_output_filtered.h5:/matrix ./FC28/FC28_corrected_output_filtered_Seurat.h5:/matrix 



