#!/bin/bash 


  ### 0.0 List CellTypes to compare -------------------------------------------- 

ls ../Data/Wig/Case_WNN_L25 | grep bedgraph | grep ALS | cut -d "_" -f2,3,4,5 | cut -d "." -f1 > CellTypes.txt




  ### 1.0 Setup DaPars Environment ---------------------------------------------

dapars_dir=`cat ../Data/cfg/Files_list.txt | grep -w DaPars_Tool1 | cut -d '"' -f2 | cut -d "~" -f2`
cp -rv ~/${dapars_dir}* ./
mv ./src/* ./; rmdir src 




  ### 2.0 APA Analysis with DaPars ---------------------------------------------

while read -r line; do
   echo -n "Processing CellType: "
   echo "$line"
   date
   
   cp -v ../Data/Wig/Case_WNN_L25/ALS_${line}.bedgraph ./
   cp -v ../Data/Wig/Case_WNN_L25/HC_${line}.bedgraph ./

   echo "Genearting DaPars Configuration File... " 
   
   echo "Annotated_3UTR=hg38_3UTR_annotation.bed" > Conf_${line}.txt
   echo Group1_Tophat_aligned_Wig=HC_${line}.bedgraph >> Conf_${line}.txt
   echo Group2_Tophat_aligned_Wig=ALS_${line}.bedgraph >> Conf_${line}.txt
   echo Output_directory=../Data/APA/DaPars1_M0_Case_WNN_L25/ALSvsHC/${line} >> Conf_${line}.txt
   echo Output_result_file=DaPars1_M0_ALSvsHC_${line} >> Conf_${line}.txt

   echo Num_least_in_group1=1 >> Conf_${line}.txt
   echo Num_least_in_group2=1 >> Conf_${line}.txt
   echo Coverage_cutoff=30 >> Conf_${line}.txt
   echo FDR_cutoff=0.05 >> Conf_${line}.txt
   echo PDUI_cutoff=0.5 >> Conf_${line}.txt
   echo Fold_change_cutoff=0.59 >> Conf_${line}.txt

   echo "Analyzing APA with DaPars..." 
   
   python DaPars_main.py Conf_${line}.txt 

   echo "Done!" 
   date 

   rm ./ALS_${line}.bedgraph
   rm ./HC_${line}.bedgraph
   rm Conf_${line}.txt 

done < CellTypes.txt



