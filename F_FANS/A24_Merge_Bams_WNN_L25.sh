#!/bin/bash 

# A23_Merge_Bams_WNN_L25.sh




  ### 1.0 WNN_L25 --------------------------------------------------------------



    ## 2.1 List WNN_L25 Cell Types ---------------------------------------------
    
#ls ../Data/FANS/Bam/PerSample/WNN_L25/ | cut -d "_" -f4,5,6 | cut -d "." -f1 | sort | uniq > CellTypes.txt
echo Oligodendrocytes > CellTypes.txt 
echo OPC >> CellTypes.txt


    ## 2.2. Sample groups ------------------------------------------------------

High=(FC3 FC4 FC21 FC25 FC27)
Low=(FC2 FC22 FC23 FC24 FC26 FC28)


    ## 2.3 Merge files ---------------------------------------------------------

while read -r CellType; do
    
    echo "`date +"%T"`: Processing $CellType..." 
    for list_name in High Low; do
    files=()

    # Get the list of IDs from the named array
    eval "ids=(\"\${${list_name}[@]}\")"

    for id in "${ids[@]}"; do
            file="../Data/FANS/Bam/PerSample/WNN_L25/${id}_WNN_L25_${CellType}.bam"
            echo -n "Looking for file ${file}..."
            if [[ -f "$file" ]]; then
                echo " found."
                files+=("$file") 
            else 
              echo " NOT FOUND!"
            fi
            
        done

        # If any files matched, run merge 
        if [[ ${#files[@]} -gt 0 ]]; then
            
            echo "`date +"%T"`: Merging files for $CellType in $list_name: ${files[*]}"
            samtools merge -@ $nthr \
              "${files[@]}" \
              -o ../Data/FANS/Bam/ByTDP43/WNN_L25/TDP43_${list_name}_WNN_L25_${CellType}.bam
            
            
            echo "###"
            echo "`date +"%T"`: Sorting file for TDP43_${list_name}, ${CellType}..."
            
            samtools sort -@ $nthr \
              ../Data/FANS/Bam/ByTDP43/WNN_L25/TDP43_${list_name}_WNN_L25_${CellType}.bam \
              -O bam \
              -o TDP43_${list_name}_${CellType}.bam 
              
            rm ../Data/FANS/Bam/ByTDP43/WNN_L25/TDP43_${list_name}_WNN_L25_${CellType}.bam
            mv TDP43_${list_name}_${CellType}.bam  ../Data/FANS/Bam/ByTDP43/WNN_L25/TDP43_${list_name}_WNN_L25_${CellType}.bam 
            
            echo "###"
            echo "`date +"%T"`: Indexing file for TDP43_${list_name}, ${CellType}..."
           
            samtools index -@ $nthr ../Data/FANS/Bam/ByTDP43/WNN_L25/TDP43_${list_name}_WNN_L25_${CellType}.bam 
            
            
            echo "###"
            echo "`date +"%T"`: Getting stats for TDP43_${list_name}, ${CellType}..."
           
    
            samtools flagstat -@ $nthr ../Data/FANS/Bam/ByTDP43/WNN_L25/TDP43_${list_name}_WNN_L25_${CellType}.bam  > ../Data/FANS/Bam/ByTDP43/WNN_L25/TDP43_${list_name}_WNN_L25_${CellType}_Flagstat.txt
        
        else
            echo "No files found for $CellType in $list_name. Looked for:"
            for id in "${ids[@]}"; do
                echo "  ../Data/FANS/Bam/PerSample/WNN_L25_${id}_${CellType}.bam"
            done
            echo "Skipping merge."
                
        fi
    done
done < CellTypes.txt

rm CellTypes.txt


mkdir ../Data/FANS/Bam/ByTDP43/WNN_L25/Flagstat 
mv -v ../Data/FANS/Bam/ByTDP43/WNN_L25/TDP43*Flagstat* ../Data/FANS/Bam/ByTDP43/WNN_L25/Flagstat/ 


ls ../Data/FANS/Bam/ByTDP43/WNN_L25/TDP43_High*bam | cut -d "_" -f6,7,8 | \
  cut -d "." -f1 | sort | uniq > CellTypesHigh.txt

while read -r CellType; do 

  echo -n TDP43_High_WNN_L25_${CellType}.bedgraph >> ../Data/FANS/Bam/ByTDP43/WNN_L25/Flagstat/WNN_L25_Merged_Sample_Depth.txt
  echo -en "\t" >> ../Data/FANS/Bam/ByTDP43/WNN_L25/Flagstat/WNN_L25_Merged_Sample_Depth.txt
  head -n7 ../Data/FANS/Bam/ByTDP43/WNN_L25/Flagstat/TDP43_High_WNN_L25_${CellType}_Flagstat.txt | tail -n1 | cut -d "+" -f1 | cut -d " " -f1 >> ../Data/FANS/Bam/ByTDP43/WNN_L25/Flagstat/WNN_L25_Merged_Sample_Depth.txt

done < CellTypesHigh.txt


ls ../Data/FANS/Bam/ByTDP43/WNN_L25/TDP43_Low*bam | cut -d "_" -f6,7,8 | \
  cut -d "." -f1 | sort | uniq > CellTypesLow.txt

while read -r CellType; do 

  echo -n TDP43_Low_WNN_L25_${CellType}.bedgraph >> ../Data/FANS/Bam/ByTDP43/WNN_L25/Flagstat/WNN_L25_Merged_Sample_Depth.txt
  echo -en "\t" >> ../Data/FANS/Bam/ByTDP43/WNN_L25/Flagstat/WNN_L25_Merged_Sample_Depth.txt
  head -n7 ../Data/FANS/Bam/ByTDP43/WNN_L25/Flagstat/TDP43_Low_WNN_L25_${CellType}_Flagstat.txt | tail -n1 | cut -d "+" -f1 | cut -d " " -f1 >> ../Data/FANS/Bam/ByTDP43/WNN_L25/Flagstat/WNN_L25_Merged_Sample_Depth.txt

done < CellTypesLow.txt

cat ../Data/FANS/Bam/ByTDP43/WNN_L25/Flagstat/WNN_L25_Merged_Sample_Depth.txt | cut -d $'\t' \
  -f1 | cut -d "_" -f5,6,7 | cut -d "." -f1 | sort | uniq > CellTypes.txt

while read -r CellType; do
  grep $CellType ../Data/FANS/Bam/ByTDP43/WNN_L25/Flagstat/WNN_L25_Merged_Sample_Depth.txt > ../Data/FANS/Bam/ByTDP43/WNN_L25/Flagstat/WNN_L25_${CellType}_Merged_Sample_Depth.txt
done < CellTypes.txt 

rm ./CellTypes.txt ./CellTypesLow.txt ./CellTypesHigh.txt
