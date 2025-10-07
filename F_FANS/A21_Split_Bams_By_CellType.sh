#!/bin/bash 

# A21_Split_Bams_By_CellType.sh

  # Prepare barcodes 
mkdir Barcodes
cp -rv ../Data/FANS/Bam/Barcodes/ByCell/* ./Barcodes
ls ./Barcodes > Barcodes_Files.txt

  # Split files in 200 Barcode chunks 

while read -r barcode_file;
do 
  
  echo "########################################"
  
  Sample=`echo $barcode_file | cut -d "_" -f1` 
  Part=`echo $barcode_file | cut -d "_" -f5 | cut -d "." -f1`
  
  i=`grep -n $barcode_file Barcodes_Files.txt | cut -d ":" -f1`
  echo "Iteration ${i}..."
  
  echo "`date +"%T"`: Splitting Bam for Sample ${Sample}, Part ${Part}..." 

  sinto filterbarcodes -b ../Input/${Sample}.bam \
    -c ./Barcodes/${barcode_file} \
    -p 3 \
    --outdir ../Data/FANS/Bam/PerSample
    
   echo "`date +"%T"`: Done! ..." 
   echo 
   echo $barcode_file >> ./Files_Processed.txt
  
done < Barcodes_Files.txt 


rm ./Barcodes_Files 
rm ./Files_Processed.txt 


ls ../Data/FANS/Bam/PerSample > Files.txt

while read -r bamfile;
do 
  
  echo "########################################"
  
  i=`grep -n $bamfile Files.txt | cut -d ":" -f1`
  echo "Iteration ${i}..."
  
  echo "`date +"%T"`: Sorting Bam ${bamfile}..." 

  samtools sort -@22 \
    ../Data/FANS/Bam/PerSample/${bamfile} \
    -O bam \
    -o ./${bamfile}
    
  rm ../Data/FANS/Bam/PerSample/${bamfile} 
  mv ./${bamfile} ../Data/FANS/Bam/PerSample/${bamfile} 
  
  echo "`date +"%T"`: Indexing Bam ${bamfile}..." 
  samtools index -@22 \
    ../Data/FANS/Bam/PerSample/${bamfile} 
    
  echo "`date +"%T"`: Done! ..." 
  echo 
  echo $bamfile >> ./Files_Processed.txt
  
done < Files.txt

rm ./Files.txt 
rm ./Files_Processed.txt 
