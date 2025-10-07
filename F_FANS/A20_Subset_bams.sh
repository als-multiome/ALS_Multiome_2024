#!/bin/bash 

# A20_Subset_Bams.sh 

CR_Output=`cat ../Data/cfg/Files_list.txt | grep FANS_Cellranger_Output_Dir | cut -d $'\t' -f2`
CR_Output="${CR_Output/\$HOME/$HOME}"



for file in ../Data/FANS/Bam/Barcodes/BySample/AllCells/*
do
  
  Sample=`basename $file | cut -d "_" -f1` 
  echo "`date +"%T"`: Copying Bam for Sample $Sample..." 
  
  ID=`cat ../Data/Input/FANS_Seq_TDP43_Samples.txt | grep ${Sample} | cut -d $'\t' -f2`
  cp ${CR_Output}/${ID}/outs/possorted_genome_bam.bam ./${Sample}.bam 
  cp ${CR_Output}/${ID}/outs/possorted_genome_bam.bam.bai ./${Sample}.bam.bai 
  
  echo "`date +"%T"`: Subsetting Bam for Sample $Sample..." 
  
  subset-bam_linux --bam ./${Sample}.bam \
    --cell-barcodes $file \
    --cores 44 \
    --out-bam ./tmp.bam 
  
  
  echo "`date +"%T"`: Indexing Bam for Sample $Sample..." 
  samtools index -@44 -b ./tmp.bam
  
  echo "`date +"%T"`: Moving Bam File and Bam File Index for Sample $Sample..."   
  mv ./tmp.bam ../Data/FANS/Bam/PerSample/${Sample}.bam 
  mv ./tmp.bam.bai ../Data/FANS/Bam/PerSample/${Sample}.bam.bai  
  
   echo "`date +"%T"`: Cleaning Up for Sample $Sample..."     
  rm ./${Sample}.bam
  rm ./${Sample}.bam.bai

  
  unset Sample
  unset ID
  
  echo Done! 
  echo 
  echo 
  echo 
  
done