#!/bin/bash 

# A23_Merge_Bams.sh


mkdir Input 

  ### 1.0 All Cells ------------------------------------------------------------



    ## 1.1 Load Input data -----------------------------------------------------

cp -rv ../Data/FANS/Bam/PerSample/AllCells/FC* ./Input/



    ## 1.2 Merge files ---------------------------------------------------------

samtools merge -@ $nthr \
  ./Input/FC3.bam \
  ./Input/FC4.bam \
  ./Input/FC21.bam \
  ./Input/FC25.bam \
  ./Input/FC27.bam \
  -o ../Data/FANS/Bam/ByTDP43/AllCells/TDP43_High_AllCells.bam

samtools merge -@ $nthr \
  ./Input/FC2.bam \
  ./Input/FC22.bam \
  ./Input/FC23.bam \
  ./Input/FC24.bam \
  ./Input/FC26.bam \
  ./Input/FC28.bam \
  -o ../Data/FANS/Bam/ByTDP43/AllCells/TDP43_Low_AllCells.bam 
  


    ## 1.3 Sort merged Bam files -----------------------------------------------

samtools sort -@ $nthr \
  ../Data/FANS/Bam/ByTDP43/AllCells/TDP43_High_AllCells.bam \
  -O bam \
  -o tmp1.bam 
  
rm ../Data/FANS/Bam/ByTDP43/AllCells/TDP43_High_AllCells.bam
mv tmp1.bam ../Data/FANS/Bam/ByTDP43/AllCells/TDP43_High_AllCells.bam 
  

samtools sort -@ $nthr \
  ../Data/FANS/Bam/ByTDP43/AllCells/TDP43_Low_AllCells.bam \
  -O bam \
  -o tmp2.bam 
  
rm ../Data/FANS/Bam/ByTDP43/AllCells/TDP43_Low_AllCells.bam
mv tmp2.bam ../Data/FANS/Bam/ByTDP43/AllCells/TDP43_Low_AllCells.bam  
  


    ## 1.4 Index merged Bam files -----------------------------------------------

samtools index -@ $nthr ../Data/FANS/Bam/ByTDP43/AllCells/TDP43_High_AllCells.bam  
samtools index -@ $nthr ../Data/FANS/Bam/ByTDP43/AllCells/TDP43_Low_AllCells.bam  



    ## 1.5 Calculate Total Mapped reads depth ----------------------------------- 

samtools flagstat -@ $nthr ../Data/FANS/Bam/ByTDP43/AllCells/TDP43_High_AllCells.bam > ../Data/FANS/Bam/ByTDP43/AllCells/Flagstat/TDP43_High_AllCells_Flagstat.txt 
samtools flagstat -@ $nthr ../Data/FANS/Bam/ByTDP43/AllCells/TDP43_Low_AllCells.bam > ../Data/FANS/Bam/ByTDP43/AllCells/Flagstat/TDP43_Low_AllCells_Flagstat.txt 

echo -n TDP43_High_AllCells.bedgraph >> ../Data/FANS/Bam/ByTDP43/AllCells/Flagstat/AllCells_Merged_Sample_Depth.txt
echo -en "\t" >> ../Data/FANS/Bam/ByTDP43/AllCells/Flagstat/AllCells_Merged_Sample_Depth.txt
head -n7 ../Data/FANS/Bam/ByTDP43/AllCells/Flagstat/TDP43_High_AllCells_Flagstat.txt | tail -n1 | cut -d "+" -f1 | cut -d " " -f1 >> ../Data/FANS/Bam/ByTDP43/Flagstat/AllCells_Merged_Sample_Depth.txt

echo -n TDP43_Low_AllCells.bedgraph >> ../Data/FANS/Bam/ByTDP43/AllCells/Flagstat/AllCells_Merged_Sample_Depth.txt
echo -en "\t" >> ../Data/FANS/Bam/ByTDP43/AllCells/Flagstat/AllCells_Merged_Sample_Depth.txt
head -n7 ../Data/FANS/Bam/ByTDP43/AllCells/Flagstat/TDP43_Low_AllCells_Flagstat.txt | tail -n1 | cut -d "+" -f1 | cut -d " " -f1 >> ../Data/FANS/Bam/ByTDP43/Flagstat/AllCells_Merged_Sample_Depth.txt

rm -r ./Input/* 




  ### 2.0 WNN_L15 --------------------------------------------------------------



    ## 2.1 Load Input data -----------------------------------------------------

cp -rv ../Data/FANS/Bam/PerSample/WNN_L15/FC* ./Input/



    ## 2.2 Merge files ---------------------------------------------------------

samtools merge -@ $nthr \
  ./Input/FC3_WNN_L15_Glia.bam \
  ./Input/FC4_WNN_L15_Glia.bam \
  ./Input/FC21_WNN_L15_Glia.bam \
  ./Input/FC25_WNN_L15_Glia.bam \
  ./Input/FC27_WNN_L15_Glia.bam \
  -o ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_High_WNN_L15_Glia.bam

samtools merge -@ $nthr \
  ./Input/FC2_WNN_L15_Glia.bam \
  ./Input/FC22_WNN_L15_Glia.bam \
  ./Input/FC23_WNN_L15_Glia.bam \
  ./Input/FC24_WNN_L15_Glia.bam \
  ./Input/FC26_WNN_L15_Glia.bam \
  ./Input/FC28_WNN_L15_Glia.bam \
  -o ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_Low_WNN_L15_Glia.bam 
  

samtools merge -@ $nthr \
  ./Input/FC3_WNN_L15_Exc_Neurons.bam \
  ./Input/FC4_WNN_L15_Exc_Neurons.bam \
  ./Input/FC21_WNN_L15_Exc_Neurons.bam \
  ./Input/FC25_WNN_L15_Exc_Neurons.bam \
  ./Input/FC27_WNN_L15_Exc_Neurons.bam \
  -o ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_High_WNN_L15_Exc_Neurons.bam

samtools merge -@ $nthr \
  ./Input/FC2_WNN_L15_Exc_Neurons.bam \
  ./Input/FC22_WNN_L15_Exc_Neurons.bam \
  ./Input/FC23_WNN_L15_Exc_Neurons.bam \
  ./Input/FC24_WNN_L15_Exc_Neurons.bam \
  ./Input/FC26_WNN_L15_Exc_Neurons.bam \
  ./Input/FC28_WNN_L15_Exc_Neurons.bam \
  -o ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_Low_WNN_L15_Exc_Neurons.bam 


samtools merge -@ $nthr \
  ./Input/FC4_WNN_L15_Inh_Neurons.bam \
  ./Input/FC21_WNN_L15_Inh_Neurons.bam \
  ./Input/FC25_WNN_L15_Inh_Neurons.bam \
  ./Input/FC27_WNN_L15_Inh_Neurons.bam \
  -o ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_High_WNN_L15_Inh_Neurons.bam

samtools merge -@ $nthr \
  ./Input/FC22_WNN_L15_Inh_Neurons.bam \
  ./Input/FC23_WNN_L15_Inh_Neurons.bam \
  ./Input/FC24_WNN_L15_Inh_Neurons.bam \
  ./Input/FC26_WNN_L15_Inh_Neurons.bam \
  ./Input/FC28_WNN_L15_Inh_Neurons.bam \
  -o ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_Low_WNN_L15_Inh_Neurons.bam 



    ## 2.3 Sort merged Bam files -----------------------------------------------

samtools sort -@ $nthr \
  ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_High_WNN_L15_Glia.bam \
  -O bam \
  -o Glia.bam 
  
rm ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_High_WNN_L15_Glia.bam
mv Glia.bam ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_High_WNN_L15_Glia.bam 
  

samtools sort -@ $nthr \
  ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_Low_WNN_L15_Glia.bam \
  -O bam \
  -o Glia2.bam 
  
rm ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_Low_WNN_L15_Glia.bam
mv Glia2.bam ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_Low_WNN_L15_Glia.bam  
  



samtools sort -@ $nthr \
  ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_High_WNN_L15_Exc_Neurons.bam \
  -O bam \
  -o Exc_Neurons.bam 
  
rm ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_High_WNN_L15_Exc_Neurons.bam
mv Exc_Neurons.bam ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_High_WNN_L15_Exc_Neurons.bam 
  

samtools sort -@ $nthr \
  ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_Low_WNN_L15_Exc_Neurons.bam \
  -O bam \
  -o Exc_Neurons2.bam 
  
rm ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_Low_WNN_L15_Exc_Neurons.bam
mv Exc_Neurons2.bam ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_Low_WNN_L15_Exc_Neurons.bam   




samtools sort -@ $nthr \
  ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_High_WNN_L15_Inh_Neurons.bam \
  -O bam \
  -o Inh_Neurons.bam 
  
rm ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_High_WNN_L15_Inh_Neurons.bam
mv Inh_Neurons.bam ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_High_WNN_L15_Inh_Neurons.bam 
  

samtools sort -@ $nthr \
  ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_Low_WNN_L15_Inh_Neurons.bam \
  -O bam \
  -o Inh_Neurons2.bam 
  
rm ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_Low_WNN_L15_Inh_Neurons.bam
mv Inh_Neurons2.bam ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_Low_WNN_L15_Inh_Neurons.bam  
  


    ## 2.4 Index merged Bam files -----------------------------------------------

samtools index -@ $nthr ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_High_WNN_L15_Glia.bam  
samtools index -@ $nthr ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_Low_WNN_L15_Glia.bam  

samtools index -@ $nthr ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_High_WNN_L15_Exc_Neurons.bam  
samtools index -@ $nthr ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_Low_WNN_L15_Exc_Neurons.bam  

samtools index -@ $nthr ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_High_WNN_L15_Inh_Neurons.bam  
samtools index -@ $nthr ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_Low_WNN_L15_Inh_Neurons.bam  



    ## 2.5 Calculate Total Mapped reads depth ----------------------------------- 

samtools flagstat -@ $nthr ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_High_WNN_L15_Glia.bam > ../Data/FANS/Bam/ByTDP43/WNN_L15/Flagstat/TDP43_High_WNN_L15_Glia_Flagstat.txt 
samtools flagstat -@ $nthr ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_Low_WNN_L15_Glia.bam > ../Data/FANS/Bam/ByTDP43/WNN_L15/Flagstat/TDP43_Low_WNN_L15_Glia_Flagstat.txt 

samtools flagstat -@ $nthr ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_High_WNN_L15_Exc_Neurons.bam > ../Data/FANS/Bam/ByTDP43/WNN_L15/Flagstat/TDP43_High_WNN_L15_Exc_Neurons_Flagstat.txt 
samtools flagstat -@ $nthr ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_Low_WNN_L15_Exc_Neurons.bam > ../Data/FANS/Bam/ByTDP43/WNN_L15/Flagstat/TDP43_Low_WNN_L15_Exc_Neurons_Flagstat.txt 

samtools flagstat -@ $nthr ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_High_WNN_L15_Inh_Neurons.bam > ../Data/FANS/Bam/ByTDP43/WNN_L15/Flagstat/TDP43_High_WNN_L15_Inh_Neurons_Flagstat.txt 
samtools flagstat -@ $nthr ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_Low_WNN_L15_Inh_Neurons.bam > ../Data/FANS/Bam/ByTDP43/WNN_L15/Flagstat/TDP43_Low_WNN_L15_Inh_Neurons_Flagstat.txt 



echo -n TDP43_High_WNN_L15_Glia.bedgraph >> ../Data/FANS/Bam/ByTDP43/WNN_L15/Flagstat/WNN_L15_Glia_Merged_Sample_Depth.txt
echo -en "\t" >> ../Data/FANS/Bam/ByTDP43/WNN_L15/Flagstat/WNN_L15_Glia_Merged_Sample_Depth.txt
head -n7 ../Data/FANS/Bam/ByTDP43/WNN_L15/Flagstat/TDP43_High_WNN_L15_Glia_Flagstat.txt | tail -n1 | cut -d "+" -f1 | cut -d " " -f1 >> ../Data/FANS/Bam/ByTDP43/WNN_L15/Flagstat/WNN_L15_Glia_Merged_Sample_Depth.txt

echo -n TDP43_Low_WNN_L15_Glia.bedgraph >> ../Data/FANS/Bam/ByTDP43/WNN_L15/Flagstat/WNN_L15_Glia_Merged_Sample_Depth.txt
echo -en "\t" >> ../Data/FANS/Bam/ByTDP43/WNN_L15/Flagstat/WNN_L15_Glia_Merged_Sample_Depth.txt
head -n7 ../Data/FANS/Bam/ByTDP43/WNN_L15/Flagstat/TDP43_Low_WNN_L15_Glia_Flagstat.txt | tail -n1 | cut -d "+" -f1 | cut -d " " -f1 >> ../Data/FANS/Bam/ByTDP43/WNN_L15/Flagstat/WNN_L15_Glia_Merged_Sample_Depth.txt


echo -n TDP43_High_WNN_L15_Exc_Neurons.bedgraph >> ../Data/FANS/Bam/ByTDP43/WNN_L15/Flagstat/WNN_L15_Exc_Neurons_Merged_Sample_Depth.txt
echo -en "\t" >> ../Data/FANS/Bam/ByTDP43/WNN_L15/Flagstat/WNN_L15_Exc_Neurons_Merged_Sample_Depth.txt
head -n7 ../Data/FANS/Bam/ByTDP43/WNN_L15/Flagstat/TDP43_High_WNN_L15_Exc_Neurons_Flagstat.txt | tail -n1 | cut -d "+" -f1 | cut -d " " -f1 >> ../Data/FANS/Bam/ByTDP43/WNN_L15/Flagstat/WNN_L15_Exc_Neurons_Merged_Sample_Depth.txt

echo -n TDP43_Low_WNN_L15_Exc_Neurons.bedgraph >> ../Data/FANS/Bam/ByTDP43/WNN_L15/Flagstat/WNN_L15_Exc_Neurons_Merged_Sample_Depth.txt
echo -en "\t" >> ../Data/FANS/Bam/ByTDP43/WNN_L15/Flagstat/WNN_L15_Exc_Neurons_Merged_Sample_Depth.txt
head -n7 ../Data/FANS/Bam/ByTDP43/WNN_L15/Flagstat/TDP43_Low_WNN_L15_Exc_Neurons_Flagstat.txt | tail -n1 | cut -d "+" -f1 | cut -d " " -f1 >> ../Data/FANS/Bam/ByTDP43/WNN_L15/Flagstat/WNN_L15_Exc_Neurons_Merged_Sample_Depth.txt


echo -n TDP43_High_WNN_L15_Inh_Neurons.bedgraph >> ../Data/FANS/Bam/ByTDP43/WNN_L15/Flagstat/WNN_L15_Inh_Neurons_Merged_Sample_Depth.txt
echo -en "\t" >> ../Data/FANS/Bam/ByTDP43/WNN_L15/Flagstat/WNN_L15_Inh_Neurons_Merged_Sample_Depth.txt
head -n7 ../Data/FANS/Bam/ByTDP43/WNN_L15/Flagstat/TDP43_High_WNN_L15_Inh_Neurons_Flagstat.txt | tail -n1 | cut -d "+" -f1 | cut -d " " -f1 >> ../Data/FANS/Bam/ByTDP43/WNN_L15/Flagstat/WNN_L15_Inh_Neurons_Merged_Sample_Depth.txt

echo -n TDP43_Low_WNN_L15_Inh_Neurons.bedgraph >> ../Data/FANS/Bam/ByTDP43/WNN_L15/Flagstat/WNN_L15_Inh_Neurons_Merged_Sample_Depth.txt
echo -en "\t" >> ../Data/FANS/Bam/ByTDP43/WNN_L15/Flagstat/WNN_L15_Inh_Neurons_Merged_Sample_Depth.txt
head -n7 ../Data/FANS/Bam/ByTDP43/WNN_L15/Flagstat/TDP43_Low_WNN_L15_Inh_Neurons_Flagstat.txt | tail -n1 | cut -d "+" -f1 | cut -d " " -f1 >> ../Data/FANS/Bam/ByTDP43/WNN_L15/Flagstat/WNN_L15_Inh_Neurons_Merged_Sample_Depth.txt



rm -r ./Input/* 


