

  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 


library(qs) 
library(tidyverse) 




  
  ### 1.0 Load data ------------------------------------------------------------ 

Chromosomes <- qread(
  paste0(
    "../Data/Annotations/Genomic_Regions/",
    "Chromosomes", 
    ".qrds"
  ), 
  nthr=nthr
)

Peaks <- qread(
  paste0(
    "../Data/Annotations/Genomic_Regions/", 
    "Peaks", 
    "qrds"
  ), 
  nthr=nthr
)




  ### 2.0 Plot ATAC Peak GC content --------------------------------------------

ggplot(Chromosomes) + 
  aes(GC_Content*100, PeakGC, label=Chromosome) + 
  geom_abline(slope=1, intercept = 0, color="#FF0000") + 
  geom_smooth(method="lm") + 
  geom_point(aes(fill=Chromosome), size=4, pch=21) +  
  geom_text_repel(max.overlaps = 200) + 
  scale_x_continuous(limits=c(35,50)) + 
  scale_y_continuous(limits=c(35,55)) + 
  xlab("\n Chromosome GC content [ % ]") + 
  ylab("Peaks GC content [ % ] \n") + 
  theme_classic()  + 
  theme(
    axis.title = element_text(size=12, color="black", face="bold"), 
    axis.text.x = element_text(color="black", face="bold.italic", angle=45, hjust = 1),  
    axis.text.y = element_text(color="black", face="bold.italic"), 
    legend.position="Null"
  )

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/SupplFig_ATAC_Peak_QC/", 
    "ATAC_GRs_GC_Content", 
    ".pdf"
  ), 
  width = 128, 
  height = 128, 
  units="mm"
)




  ### 3.0 Plot ATAC Peak numbers pro chromosome --------------------------------

ggplot(Chromosomes) + 
  aes(Length, N_Peaks, label=Chromosome) + 
  geom_smooth(method="lm") + 
  geom_point(aes(fill=Chromosome), size=4, pch=21) +  
  geom_text_repel() + 
  xlab("\n Chromosome length [Mbp]") + 
  ylab("Number of peaks \n") + 
  theme_classic()  + 
  theme(
    axis.title = element_text(size=12, color="black", face="bold"), 
    axis.text.x = element_text(color="black", face="bold.italic", angle=45, hjust = 1),  
    axis.text.y = element_text(color="black", face="bold.italic"), 
    legend.position="Null"
  )

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/SupplFig_ATAC_Peak_QC/", 
    "ATAC_GRs_Chromosomes", 
    ".pdf"
  ), 
  width = 128, 
  height = 128, 
  units="mm"
)




  ### 4.0 Plot ATAC Peak width -------------------------------------------------

ggplot(as.data.frame(mcols(Peaks))) + 
  geom_violin(aes(x="All", y=sequence.length), fill="olivedrab3") + 
  geom_boxplot(aes(x="All", y=sequence.length), outlier.shape = NA, width=0.1) + 
  stat_summary(aes(x="All", y=sequence.length),
               geom="point",
               fun="mean",
               shape=21,
               fill="orangered2",
               size=5
  ) + 
  ylab("Peak width [bp]\n ") + 
  theme_classic() + 
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(), 
    axis.line = element_line(color="black"), 
    axis.text = element_text(size=10, color="black", face="bold.italic"), 
    axis.title = element_text(size=14, color="black", face="bold.italic")
  )

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/SupplFig_ATAC_Peak_QC/", 
    "ATAC_GRs_Length_Violin", 
    ".pdf"
  ), 
  width = 64, 
  height = 128, 
  units="mm"
)




  ### 5.0 Plot ATAC Peak width by Chromosome -----------------------------------

ggplot(as.data.frame(mcols(Peaks))) + 
  geom_boxplot(aes(x=Chromosome, y=sequence.length, fill=Chromosome), outlier.shape = NA) + 
  stat_summary(aes(x=Chromosome, y=sequence.length),
               geom="point",
               fun="mean",
               shape=21,
               fill="orangered2",
               size=5
  ) + 
  ylab("Peak width [bp]\n ") + 
  theme_classic() + 
  theme(
    axis.title.x = element_blank(), 
    axis.line = element_line(color="black"), 
    axis.text.x = element_text(size=10, angle=45, hjust = 1, color="black", face="bold.italic"),
    axis.text.y = element_text(size=10, color="black", face="bold.italic"),  
    
    axis.title = element_text(size=14, color="black", face="bold.italic"), 
    legend.position = "Null"
  )

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/SupplFig_ATAC_Peak_QC/", 
    "ATAC_GRs_Length_Per_Chromosome_Boxplot", 
    ".pdf"
  ), 
  width = 192, 
  height = 128, 
  units="mm"
)


  ### 6.0 Plot ATAC Peaks TFIDF counts Violin Plot -----------------------------

ggplot(as.data.frame(mcols(Peaks))) + 
  geom_violin(aes(x="All", y=TFIDFcounts), fill="olivedrab3") + 
  geom_boxplot(aes(x="All", y=TFIDFcounts), outlier.shape = NA, width=0.1) + 
  stat_summary(aes(x="All", y=TFIDFcounts),
               geom="point",
               fun="mean",
               shape=21,
               fill="orangered2",
               size=5
  ) + 
  ylab("Peak counts [normalized]\n ") + 
  theme_classic() + 
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(), 
    axis.line = element_line(color="black"), 
    axis.text = element_text(size=10, color="black", face="bold.italic"), 
    axis.title = element_text(size=14, color="black", face="bold.italic")
  )

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/SupplFig_ATAC_Peak_QC/", 
    "ATAC_GRs_TFIDF_Counts_Violin", 
    ".pdf"
  ), 
  width = 64, 
  height = 128, 
  units="mm"
)




  ### 7.0 Plot ATAC Peaks TFIDF counts per Chromosome Boxplot ------------------

ggplot(as.data.frame(mcols(Peaks))) + 
  geom_boxplot(aes(x=Chromosome, y=TFIDFcounts, fill=Chromosome), outlier.shape = NA) + 
  stat_summary(aes(x=Chromosome, y=TFIDFcounts),
               geom="point",
               fun="mean",
               shape=21,
               fill="orangered2",
               size=5
  ) + 
  ylab("Peak counts [ TFIDF norm. ]\n ") + 
  theme_classic() + 
  theme(
    axis.title.x = element_blank(), 
    axis.line = element_line(color="black"), 
    axis.text.x = element_text(size=10, angle=45, hjust = 1, color="black", face="bold.italic"),
    axis.text.y = element_text(size=10, color="black", face="bold.italic"),  
    
    axis.title = element_text(size=14, color="black", face="bold.italic"), 
    legend.position = "Null"
  )

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/SupplFig_ATAC_Peak_QC/", 
    "ATAC_GRs_Length_Per_Chromosome_Boxplot", 
    ".pdf"
  ), 
  width = 192, 
  height = 128, 
  units="mm"
)




  ### 8.0 Plot ATAC Peaks TFIDF counts by peak length --------------------------

ggplot(as.data.frame(mcols(Peaks))) + 
  geom_bin2d(
    aes(
      sequence.length,
      TFIDFcounts
    ),
    bins=400
  ) + 
  geom_smooth(
    aes(
      sequence.length,
      TFIDFcounts
    ), 
    method="lm", 
    col="red"
  ) + 
  scale_fill_viridis_c() + 
  theme_classic() 

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/SupplFig_ATAC_Peak_QC/", 
    "ATAC_GRs_TFIDF_Counts_by_GR_Length", 
    ".pdf"
  ), 
  width = 128, 
  height = 128, 
  units="mm"
)




  ### 9.0 Plot ATAC Peaks TFIDF counts by peak GC content ----------------------

ggplot(as.data.frame(mcols(Peaks))) + 
  geom_bin2d(
    aes(
      GC.percent,
      TFIDFcounts
    ),
    bins=400
  ) + 
  geom_smooth(
    aes(
      GC.percent,
      TFIDFcounts
    ), 
    method="lm", 
    col="red"
  ) + 
  scale_fill_viridis_c() + 
  theme_classic() 

ggsave(
  filename=paste0(
    "../Data/Visualization/Figures/SupplFig_ATAC_Peak_QC/", 
    "ATAC_GRs_TFIDF_Counts_by_GR_GC_Content", 
    ".pdf"
  ), 
  width = 128, 
  height = 128, 
  units="mm"
)

ggplot(data.frame(counts=M0_TFIDF@assays$ATAC@data@x)) + 
  geom_density(aes(counts)) + 
  scale_x_continuous(expand=c(0,0.0,0,0.1)) + 
  scale_y_continuous(expand=c(0,0.0,0,0.1)) + 
  theme_classic() 


