#!/bin/bash 

# A25_Regtools_PerCell.sh

  # Prepare filelist 
../Data/FANS/Bam/PerCell/*.bam | grep -v '\.bai$' > Bamfiles.txt

  

while read -r bamfile;
do 
  
  echo "########################################"
  
  Sample=`echo $bamfile | cut -d "/" -f6 | cut -d "." -f1`

  regtools junctions extract -a 6 \
    -m 30 \
    -M 500000 \
    -s FR \
    -o ./Output/${Sample}_JC.txt \
    ${bamfile}

    
done < Bamfiles.txt 

echo "`date +"%T"`: Done"

