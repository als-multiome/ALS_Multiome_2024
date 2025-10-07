# 


bedtools genomecov -ibam ../Data/FANS/Bam/ByTDP43/TDP43_High_AllCells.bam -bga -split -trackline > ./TDP43_High_AllCells.bedgraph
bedtools genomecov -ibam ../Data/FANS/Bam/ByTDP43/TDP43_Low_AllCells.bam -bga -split -trackline > ./TDP43_Low_AllCells.bedgraph
bamCoverage -b ../Data/FANS/Bam/ByTDP43/TDP43_High_AllCells.bam -p 46 -bs 5 --normalizeUsing CPM -o ./TDP43_High_AllCells.bigWig
bamCoverage -b ../Data/FANS/Bam/ByTDP43/TDP43_Low_AllCells.bam -p 46 -bs 5 --normalizeUsing CPM -o ../Data/FANS/BigWig/TDP43_Low_AllCells.bigWig  


bedtools genomecov -ibam ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_High_WNN_L15_Glia.bam -bga -split -trackline > ../Data/FANS/Bedgraph/WNN_L15/TDP43_High_WNN_L15_Glia.bedgraph
bedtools genomecov -ibam ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_Low_WNN_L15_Glia.bam -bga -split -trackline > ../Data/FANS/Bedgraph/WNN_L15/TDP43_Low_WNN_L15_Glia.bedgraph

bedtools genomecov -ibam ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_High_WNN_L15_Exc_Neurons.bam -bga -split -trackline > ../Data/FANS/Bedgraph/WNN_L15/TDP43_High_WNN_L15_Exc_Neurons.bedgraph
bedtools genomecov -ibam ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_Low_WNN_L15_Exc_Neurons.bam -bga -split -trackline > ../Data/FANS/Bedgraph/WNN_L15/TDP43_Low_WNN_L15_Exc_Neurons.bedgraph

bedtools genomecov -ibam ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_High_WNN_L15_Inh_Neurons.bam -bga -split -trackline > ../Data/FANS/Bedgraph/WNN_L15/TDP43_High_WNN_L15_Inh_Neurons.bedgraph
bedtools genomecov -ibam ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_Low_WNN_L15_Inh_Neurons.bam -bga -split -trackline > ../Data/FANS/Bedgraph/WNN_L15/TDP43_Low_WNN_L15_Inh_Neurons.bedgraph


bamCoverage -b ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_High_WNN_L15_Glia.bam -p 46 -bs 5 --normalizeUsing CPM -o ../Data/FANS/BigWig/WNN_L15/TDP43_High_WNN_L15_Glia.bigWig
bamCoverage -b ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_Low_WNN_L15_Glia.bam -p 46 -bs 5 --normalizeUsing CPM -o ../Data/FANS/BigWig/WNN_L15/TDP43_Low_WNN_L15_Glia.bigWig 

bamCoverage -b ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_High_WNN_L15_Exc_Neurons.bam -p 46 -bs 5 --normalizeUsing CPM -o ../Data/FANS/BigWig/WNN_L15/TDP43_High_WNN_L15_Exc_Neurons.bigWig
bamCoverage -b ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_Low_WNN_L15_Exc_Neurons.bam -p 46 -bs 5 --normalizeUsing CPM -o ../Data/FANS/BigWig/WNN_L15/TDP43_Low_WNN_L15_Exc_Neurons.bigWig 

bamCoverage -b ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_High_WNN_L15_Inh_Neurons.bam -p 46 -bs 5 --normalizeUsing CPM -o ../Data/FANS/BigWig/WNN_L15/TDP43_High_WNN_L15_Inh_Neurons.bigWig
bamCoverage -b ../Data/FANS/Bam/ByTDP43/WNN_L15/TDP43_Low_WNN_L15_Inh_Neurons.bam -p 46 -bs 5 --normalizeUsing CPM -o ../Data/FANS/BigWig/WNN_L15/TDP43_Low_WNN_L15_Inh_Neurons.bigWig 

