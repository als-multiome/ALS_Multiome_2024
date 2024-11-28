#!/bin/bash 


  ### 0.0 Load Input data ------------------------------------------------------ 

cp ../Data/FANS/BedGraph/FANS_TDP43_Neg.bedgraph ./
cp ../Data/FANS/BedGraph/FANS_TDP43_Pos.bedgraph ./
cp ../Data/FANS/BedGraph/FANS_TDP43_Glia.bedgraph ./




  ### 1.0 Setup DaPars Environment ---------------------------------------------

dapars_dir=`cat ../Data/cfg/Files_list.txt | grep -w DaPars_Tool1 | cut -d '"' -f2 | cut -d "~" -f2`
cp -rv ~/${dapars_dir}* ./
mv ./src/* ./; rmdir src 




  ### 2.0 Generate DaPars Run configuration files ------------------------------

echo "Annotated_3UTR=hg38_3UTR_annotation.bed" > FANS_TDP43Pos_vs_TDP43Neg.txt
echo Group1_Tophat_aligned_Wig=FANS_TDP43_Pos.bedgraph >> FANS_TDP43Pos_vs_TDP43Neg.txt
echo Group2_Tophat_aligned_Wig=FANS_TDP43_Neg.bedgraph >> FANS_TDP43Pos_vs_TDP43Neg.txt
echo Output_directory=../Data/APA/DaPars1_TDP43Pos_vs_TDP43Neg/ >> FANS_TDP43Pos_vs_TDP43Neg.txt
echo Output_result_file=DaPars1_TDP43Pos_vs_TDP43Neg >> FANS_TDP43Pos_vs_TDP43Neg.txt 

echo Num_least_in_group1=1 >> FANS_TDP43Pos_vs_TDP43Neg.txt
echo Num_least_in_group2=1 >> FANS_TDP43Pos_vs_TDP43Neg.txt
echo Coverage_cutoff=30 >> FANS_TDP43Pos_vs_TDP43Neg.txt
echo FDR_cutoff=0.05 >> FANS_TDP43Pos_vs_TDP43Neg.txt
echo PDUI_cutoff=0.5 >> FANS_TDP43Pos_vs_TDP43Neg.txt
echo Fold_change_cutoff=0.59 >> FANS_TDP43Pos_vs_TDP43Neg.txt



echo "Annotated_3UTR=hg38_3UTR_annotation.bed" > FANS_TDP43Pos_vs_TDP43Glia.txt
echo Group1_Tophat_aligned_Wig=FANS_TDP43_Pos.bedgraph >> FANS_TDP43Pos_vs_TDP43Glia.txt
echo Group2_Tophat_aligned_Wig=FANS_TDP43_Glia.bedgraph >> FANS_TDP43Pos_vs_TDP43Glia.txt
echo Output_directory=../Data/APA/DaPars1_TDP43Pos_vs_TDP43Glia/ >> FANS_TDP43Pos_vs_TDP43Glia.txt
echo Output_result_file=DaPars1_TDP43Pos_vs_TDP43Glia >> FANS_TDP43Pos_vs_TDP43Glia.txt 

echo Num_least_in_group1=1 >> FANS_TDP43Pos_vs_TDP43Glia.txt
echo Num_least_in_group2=1 >> FANS_TDP43Pos_vs_TDP43Glia.txt
echo Coverage_cutoff=30 >> FANS_TDP43Pos_vs_TDP43Glia.txt
echo FDR_cutoff=0.05 >> FANS_TDP43Pos_vs_TDP43Glia.txt
echo PDUI_cutoff=0.5 >> FANS_TDP43Pos_vs_TDP43Glia.txt
echo Fold_change_cutoff=0.59 >> FANS_TDP43Pos_vs_TDP43Glia.txt




  ### 3.0 FANS ALSFTD DaPars APA Analysis --------------------------------------

python DaPars_main.py FANS_TDP43Pos_vs_TDP43Neg.txt

python DaPars_main.py FANS_TDP43Pos_vs_TDP43Glia.txt 

rm FANS_TDP43* 
rm -r  DaPars_*
rm hg38*
rm GRCh38_*





