


  ### 1.0 Load Input data ------------------------------------------------------ 

mkdir Input
mkdir Input/FANS_S1 
mkdir Input/FANS_S2
mkdir Input/FANS_S3
mkdir Input/FANS_S4
mkdir Input/FANS_S5 
mkdir Input/FANS_S6
mkdir Input/FANS_S7
mkdir Input/FANS_S8



    ## 1.1 FANS CellRanger Bam files ------------------------------------------- 

FANS_S1=`cat ../Data/cfg/Files_list.txt | grep FANS_S1_Outs | cut -d '"' -f2 | cut -d "~" -f2`
FANS_S2=`cat ../Data/cfg/Files_list.txt | grep FANS_S2_Outs | cut -d '"' -f2 | cut -d "~" -f2`
FANS_S3=`cat ../Data/cfg/Files_list.txt | grep FANS_S3_Outs | cut -d '"' -f2 | cut -d "~" -f2`
FANS_S4=`cat ../Data/cfg/Files_list.txt | grep FANS_S4_Outs | cut -d '"' -f2 | cut -d "~" -f2`
FANS_S5=`cat ../Data/cfg/Files_list.txt | grep FANS_S5_Outs | cut -d '"' -f2 | cut -d "~" -f2`
FANS_S6=`cat ../Data/cfg/Files_list.txt | grep FANS_S6_Outs | cut -d '"' -f2 | cut -d "~" -f2`
FANS_S7=`cat ../Data/cfg/Files_list.txt | grep FANS_S7_Outs| cut -d '"' -f2 | cut -d "~" -f2`
FANS_S8=`cat ../Data/cfg/Files_list.txt | grep FANS_S8_Outs| cut -d '"' -f2 | cut -d "~" -f2`

cp -v ~/${FANS_S1}possorted_genome_bam.bam* ./Input/FANS_S1/
cp -v ~/${FANS_S2}possorted_genome_bam.bam* ./Input/FANS_S2/
cp -v ~/${FANS_S3}possorted_genome_bam.bam* ./Input/FANS_S3/
cp -v ~/${FANS_S4}possorted_genome_bam.bam* ./Input/FANS_S4/
cp -v ~/${FANS_S5}possorted_genome_bam.bam* ./Input/FANS_S5/
cp -v ~/${FANS_S6}possorted_genome_bam.bam* ./Input/FANS_S6/
cp -v ~/${FANS_S7}possorted_genome_bam.bam* ./Input/FANS_S7/
cp -v ~/${FANS_S8}possorted_genome_bam.bam* ./Input/FANS_S8/




  ### 2.0 Generate merged FANS Bam files ---------------------------------------

samtools merge -@ $nthr ./Input/FANS_S3/possorted_genome_bam.bam ./Input/FANS_S4/possorted_genome_bam.bam ./Input/FANS_S5/possorted_genome_bam.bam ./Input/FANS_S6/possorted_genome_bam.bam -o ../Data/FANS/Bam/FANS_TDP43_Pos.bam

samtools merge -@ $nthr ./Input/FANS_S1/possorted_genome_bam.bam ./Input/FANS_S2/possorted_genome_bam.bam -o ../Data/FANS/Bam/FANS_TDP43_Neg.bam

samtools merge -@ $nthr ./Input/FANS_S7/possorted_genome_bam.bam ./Input/FANS_S7/possorted_genome_bam.bam -o ../Data/FANS/Bam/FANS_TDP43_Glia.bam 




  ### 3.0 Index merged Bam files -----------------------------------------------

samtools index -@ $nthr ../Data/FANS/Bam/FANS_TDP43_Pos.bam
samtools index -@ $nthr ../Data/FANS/Bam/FANS_TDP43_Neg.bam
samtools index -@ $nthr ../Data/FANS/Bam/FANS_TDP43_Glia.bam 




  ### 4.0 Calculate Total Mapped reads depth ----------------------------------- 

samtools flagstat -@ $nthr ../Data/FANS/Bam/FANS_TDP43_Pos.bam > ../Data/FANS/Flagstat/FANS_TDP43_Pos_Flagstat.txt
samtools flagstat -@ $nthr ../Data/FANS/Bam/FANS_TDP43_Neg.bam > ../Data/FANS/Flagstat/FANS_TDP43_Neg_Flagstat.txt 
samtools flagstat -@ $nthr ../Data/FANS/Bam/FANS_TDP43_Glia.bam > ../Data/FANS/Flagstat/FANS_TDP43_Glia_Flagstat.txt 

echo -e "wig\tdepth" > ../Data/FANS/Flagstat/FANS_Merged_Sample_Depth.txt
echo -n FANS_TDP43_Pos.bedgraph >> ../Data/FANS/Flagstat/FANS_Merged_Sample_Depth.txt
echo -en "\t" >> ../Data/FANS/Flagstat/FANS_Merged_Sample_Depth.txt
head -n7 ../Data/FANS/Flagstat/FANS_TDP43_Pos_Flagstat.txt | tail -n1 | cut -d "+" -f1 | cut -d " " -f1 >> ../Data/FANS/Flagstat/FANS_Merged_Sample_Depth.txt

echo -n FANS_TDP43_Neg.bedgraph >> ../Data/FANS/Flagstat/FANS_Merged_Sample_Depth.txt
echo -en "\t" >> ../Data/FANS/Flagstat/FANS_Merged_Sample_Depth.txt
head -n7 ../Data/FANS/Flagstat/FANS_TDP43_Neg_Flagstat.txt | tail -n1 | cut -d "+" -f1 | cut -d " " -f1 >> ../Data/FANS/Flagstat/FANS_Merged_Sample_Depth.txt

echo -n FANS_TDP43_Glia.bedgraph >> ../Data/FANS/Flagstat/FANS_Merged_Sample_Depth.txt
echo -en "\t" >> ../Data/FANS/Flagstat/FANS_Merged_Sample_Depth.txt
head -n7 ../Data/FANS/Flagstat/FANS_TDP43_Glia_Flagstat.txt | tail -n1 | cut -d "+" -f1 | cut -d " " -f1 >> ../Data/FANS/Flagstat/FANS_Merged_Sample_Depth.txt

rm -r ./Input
