#!/bin/bash 

ls ../Data/FANS/Bam/ByTDP43/WNN_L4/*bam > BamFilesWNNL4.txt
Total=`wc -l BamFilesWNNL4.txt | cut -d " " -f1` 
while read -r Bamfile; do 
  
  echo 
  echo "####################################################"
  echo 
  
  i=`grep -n ${Bamfile} BamFilesWNNL4.txt | cut -d ":" -f1` 
  echo "Iteration ${i}/${Total}..."
  
  Sample=`echo ${Bamfile} | cut -d "/" -f7,8,9,10,11 | cut -d "." -f1` 
  echo "`date +"%T"`: Bedgraphing ${Bamfile}..." 
  
  bedtools genomecov -ibam ${Bamfile} -bga -split -trackline > ../Data/FANS/Bedgraph/WNN_L4/${Sample}.bedgraph
  bamCoverage -b ${Bamfile} -p 46 -bs 5 --normalizeUsing CPM -o ../Data/FANS/BigWig/WNN_L4/${Sample}.bigWig 

done < BamFilesWNNL4.txt


