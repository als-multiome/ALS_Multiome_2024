A3_Souporcell_Doublet_Detection.sh 

#!/bin/bash 

singularity pull --arch amd64 library://wheaton5/souporcell/souporcell:release 
singularity --version
#singularity-ce version 4.1.2 

FANS_CR_Output=`grep FANS_Cellranger_Output_Dir Files_list.txt | cut -d $'\t' -f2`
FANS_CR_Output="${FANS_CR_Output/\$HOME/$HOME}"

Souporcell_Output=`grep FANS_Souporcell_Output_Dir Files_list.txt | cut -d $'\t' -f2`
Souporcell_Output="${Souporcell_Output/\$HOME/$HOME}"

Project_Dir=`grep Project_Dir Files_list.txt | cut -d $'\t' -f2`
Project_Dir="${Project_Dir/\$HOME/$HOME}"

mkdir ./Input/
mkdir ./Output/

cp -rv ${Project_Dir}/Data/Input/FANS_Seq_TDP43_Samples.txt ./
dos2unix FANS_Seq_TDP43_Samples.txt

tail -n $((`wc -l FANS_Seq_TDP43_Samples_Tmp.txt | cut -d " " -f1` - 1)) FANS_Seq_TDP43_Samples_Tmp.txt > FANS_Seq_TDP43_Samples.txt

while read -r line;
do
	Name=`echo $line | cut -d " " -f2`
	Samples=`echo $line | cut -d " " -f5`
	ID=`echo $line | cut -d " " -f1`
	
	echo ############################################# 
	echo 
	echo "Working on sample ${ID}" 
	date 
	echo 

	mkdir ./Input/${ID} 
	mkdir ./Output/${ID} 

	cp -rv ${FANS_CR_Output}/${Name}/outs/possorted_genome_bam.bam  ./Input/${ID} 
	cp -rv ${FANS_CR_Output}/${Name}/outs/possorted_genome_bam.bam.bai  ./Input/${ID}  
 	cp -rv ${FANS_CR_Output}/${Name}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz  ./Input/${ID}  
	gzip -d ./Input/${ID}/barcodes.tsv.gz 

	singularity exec ./souporcell_release.sif souporcell_pipeline.py -i ./Input/${ID}/possorted_genome_bam.bam -b ./Input/${ID}/barcodes.tsv -f ./GenomeReference/genome.fa -t 36 -o ./Output/${ID} -k $Samples

	mv stan_consensus.pickle ./Output/{ID} 
	mv ./Output/${ID} ${Souporcell_Output} 

	rm -r ./Input/${ID}


done < FANS_Seq_TDP43_Samples.txt




