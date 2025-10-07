# A5_ScDblFinder_Detect_Doublets.R 

  ### 0.0 Load libraries ------------------------------------------------------- 

library(data.table)
library(Seurat)
library(qs2)
library(tidyverse)
library(scDblFinder)

  
    
  
  ### 1.0 Load data ------------------------------------------------------------
  
Samples <- qs_read(
  "../Data/FANS/Samples/Samples.qs2"
)

FANS_Seurats <- qs_read(
  "../Data/FANS/SeuratObjects/FANS_Single_SeuratObjects_Filtered.qs2"
)




  ### 2.0 Detect doublets with ScDblFinder -------------------------------------

n_ArtDoubl <- c(
  100, 200, 500, 1000, 
  2000, 5000, 10000, 20000, 
  40000, 50000, 80000, 100000
)  

FANS_ScDblFinder_SCEs <- lapply(
  FANS_Seurats, 
  function(Seurat){
    print(paste0("Working on sample ", Seurat$orig.ident[1], "... ")) 
    return(
      lapply(
        setNames(n_ArtDoubl, nm = n_ArtDoubl), 
        function(n_ArtDoubl){
          return(
            scDblFinder(
              sce = Seurat@assays$RNA$counts, 
              clusters = TRUE, 
              artificialDoublets = n_ArtDoubl
              )
          )
        }
      )
    )
  }
)




  ### 3.0 Plot ArtDoub/Doublet rates and choose threshold ---------------------- 

Doublet_Rates <- data.frame(
  n_ArtDoubl = n_ArtDoubl  
)
lapply(
  names(FANS_ScDblFinder_SCEs), 
  function(Sample){
    Doublet_Rates[[Sample]] <<- sapply(
      FANS_ScDblFinder_SCEs[[Sample]], 
      function(x){
        tab <- table(x$scDblFinder.class, useNA = 'always')
        return(
          tab['doublet']/sum(tab)
        )
      }
    )
  }
)

Doublet_Rates <- Doublet_Rates %>% 
  pivot_longer(
    cols = starts_with("FC"), 
    names_to = "Sample", 
    values_to = "Pct_Doublets", 
    values_drop_na = FALSE
  )

Doublet_Rates %>% 
  mutate(
    Sample = factor(
      Sample, 
      levels = str_sort(unique(Doublet_Rates$Sample), numeric = TRUE)
    )
  ) %>% 
  ggplot() + 
  aes(n_ArtDoubl, Pct_Doublets*100, fill = Sample) + 
  geom_vline(xintercept = 5000, col="#FF0000AA") + 
  geom_line(aes(col=Sample)) + 
  geom_point(pch=21, size = 2) + 
  scale_x_continuous(trans="log10") + 
  scale_y_continuous(limits = c(0,40)) + 
  theme_classic() 




  ### 3.0 Detect doublets with 5000 artificial doublets ------------------------

lapply(
  names(FANS_Seurats), 
  function(Seurat){
    SCE = scDblFinder(
      sce = FANS_Seurats[[Seurat]]@assays$RNA$counts, 
      clusters = TRUE, 
      artificialDoublets = 5000
    ) 
    message(
      paste0(
        "Sorting concordant: ", 
        all(
          rownames(SCE@colData) == rownames(FANS_Seurats[[Seurat]]@meta.data)
        )
      )
    ) 
    FANS_Seurats[[Seurat]]$scDblFinder.class <<- SCE$scDblFinder.class[
      match(
        rownames(FANS_Seurats[[Seurat]]@meta.data), 
        rownames(SCE@colData)
      )
    ]
  }
)

lapply(
  FANS_Seurats, 
  function(x){
    return(
      x@meta.data %>% 
        select(scDblFinder.class, nFeature_RNA, nCount_RNA)
    )
  }
) %>%  
  bind_rows(.id = "Sample") %>% 
  mutate(scDblFinder.class = str_replace_all(scDblFinder.class, "singlet", "Singlet")) %>% 
  mutate(scDblFinder.class = str_replace_all(scDblFinder.class, "doublet", "Doublet"))-> tmp 
  
ggplot() + 
  geom_point(
    aes(nCount_RNA, nFeature_RNA, fill = scDblFinder.class), 
    data = (tmp %>% filter(scDblFinder.class == "Singlet")), 
    pch = 21, 
    alpha = 0.6
  ) + 
  geom_point(
    aes(nCount_RNA, nFeature_RNA, fill = scDblFinder.class), 
    data = (tmp %>% filter(scDblFinder.class == "Doublet")), 
    pch = 21, 
    size = 2, 
    alpha = 1
  ) + 
  scale_x_continuous(trans="log10") + 
  scale_y_continuous(trans="log10") + 
  scale_fill_manual(
    values = c(
     "Singlet"= "olivedrab2", 
     "Doublet" = "salmon"
    )
  ) + 
  xlab("Counts") + 
  ylab("Features") + 
  theme_classic() + 
  theme(
    legend.title = element_blank()
  )

ggsave(
  "../Data/FANS/Visualization/scDblFinder_Doublets_Counts_Features_Plot_Log10.pdf", 
  dpi = 300, 
  width = 7, 
  height = 6.34
) 

ggplot() + 
  geom_point(
    aes(nCount_RNA, nFeature_RNA, fill = scDblFinder.class), 
    data = (tmp %>% filter(scDblFinder.class == "Singlet")), 
    pch = 21, 
    alpha = 0.6
  ) + 
  geom_point(
    aes(nCount_RNA, nFeature_RNA, fill = scDblFinder.class), 
    data = (tmp %>% filter(scDblFinder.class == "Doublet")), 
    pch = 21, 
    size = 2, 
    alpha = 1
  ) + 
  scale_fill_manual(
    values = c(
      "Singlet"= "olivedrab2", 
      "Doublet" = "salmon"
    )
  ) + 
  xlab("Counts") + 
  ylab("Features") + 
  theme_classic() + 
  theme(
    legend.title = element_blank()
  )

ggsave(
  "../Data/FANS/Visualization/scDblFinder_Doublets_Counts_Features_Plot_Linear.pdf", 
  dpi = 300, 
  width = 7, 
  height = 6.34
)


tmp <- tmp %>% 
  group_by(Sample) %>% 
  summarize(
    Pct_Doublets = sum(scDblFinder.class == "Doublet")*100/length(Sample)
  )



Samples$PctDoublets_scDblFinder <- tmp$Pct_Doublets[
  match(
    Samples$ID, 
    tmp$Sample
  )
]

 


  ### 5.0 Save data ------------------------------------------------------------

qs_save(
  FANS_Seurats, 
  "../Data/FANS/SeuratObjects/FANS_Single_SeuratObjects_Filtered_scDblFinder_Marked.qs2"
)


qs_save(
  Samples, 
  "../Data/FANS/Samples/Samples.qs2"
)



