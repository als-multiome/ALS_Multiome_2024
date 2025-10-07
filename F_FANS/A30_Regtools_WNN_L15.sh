#!/bin/bash 

# A30_Regtools_WNN_L15.sh

  # Prepare filelist 
ls ../Data/FANS/Bam/ByTDP43/WNN_L15/*.bam | grep -v '\.bai$' > Bamfiles.txt

while read -r bamfile;
do 
  
  echo "########################################"
  
  Sample=`echo $bamfile | cut -d "/" -f7 | cut -d "." -f1`

  regtools junctions extract -a 6 \
    -m 30 \
    -M 500000 \
    -s FR \
    -o ./Output/${Sample}_JC.txt \
    ${bamfile}

    
done < Bamfiles.txt 

echo "`date +"%T"`: Done" 

rm Bamfiles.txt


  # WNN_L25 
  
ls ../Data/FANS/Bam/ByTDP43/WNN_L25/*.bam | grep -v '\.bai$' > Bamfiles.txt

while read -r bamfile;
do 
  
  echo "########################################"
  
  Sample=`echo $bamfile | cut -d "/" -f7 | cut -d "." -f1`

  regtools junctions extract -a 6 \
    -m 30 \
    -M 500000 \
    -s FR \
    -o ./Output/${Sample}_JC.txt \
    ${bamfile}

    
done < Bamfiles.txt 

echo "`date +"%T"`: Done" 

rm Bamfiles.txt 


  # WNN_L4 
  
ls ../Data/FANS/Bam/ByTDP43/WNN_L4/*.bam | grep -v '\.bai$' > Bamfiles.txt

while read -r bamfile;
do 
  
  echo "########################################"
  
  Sample=`echo $bamfile | cut -d "/" -f7 | cut -d "." -f1`

  regtools junctions extract -a 6 \
    -m 30 \
    -M 500000 \
    -s FR \
    -o ./Output/${Sample}_JC.txt \
    ${bamfile}

    
done < Bamfiles.txt 

echo "`date +"%T"`: Done" 

rm Bamfiles.txt
