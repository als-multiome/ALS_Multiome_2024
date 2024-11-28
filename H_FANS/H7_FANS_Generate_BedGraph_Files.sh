#


  ### 1.0 Load Input data ------------------------------------------------------ 

cp ../Data/FANS/Bam/FANS_TDP43_Neg.bam ./ 
cp ../Data/FANS/Bam/FANS_TDP43_Pos.bam ./ 
cp ../Data/FANS/Bam/FANS_TDP43_Glia.bam ./ 




  ### 2.0 Generate FANS Data Bedgraph files ------------------------------------ 

bedtools genomecov -ibam ./FANS_TDP43_Neg.bam -bga -split -trackline > ../Data/FANS/BedGraph/FANS_TDP43_Neg.bedgraph 
bedtools genomecov -ibam ./FANS_TDP43_Pos.bam -bga -split -trackline > ../Data/FANS/BedGraph/FANS_TDP43_Pos.bedgraph 
bedtools genomecov -ibam ./FANS_TDP43_Glia.bam -bga -split -trackline > ../Data/FANS/BedGraph/FANS_TDP43_Glia.bedgraph 




  ### 3.0 Generate Normalized BigWig files ------------------------------------- 

bamCoverage -b ./FANS_TDP43_Neg.bam -p 46 -bs 5 --normalizeUsing CPM -o ../Data/FANS/BigWig/FANS_TDP43_Neg_CPM.bigWig
bamCoverage -b ./FANS_TDP43_Pos.bam -p 46 -bs 5 --normalizeUsing CPM -o ../Data/FANS/BigWig/FANS_TDP43_Pos_CPM.bigWig
bamCoverage -b ./FANS_TDP43_Glia.bam -p 46 -bs 5 --normalizeUsing CPM -o ../Data/FANS/BigWig/FANS_TDP43_Glia_CPM.bigWig

rm FANS_TDP43*bam 
