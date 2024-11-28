

  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 


library(qs) 
library(tidyverse) 
library(ggrepel)
library(patchwork)
library(ggpubr)

  


  ### 1.0 Load data ------------------------------------------------------------ 



    ## 1.1 Load ATAC DA Results without Co-Variates ----------------------------

Results_NoCov <- qread(
  paste0(
    "../Data/DA/WNN/AllCase/", 
    "DESeq_Results_Index", 
    ".qrds"
  ), 
  nthr=nthr
)



    ## 1.2 Load ATAC DA SVA 1-15 Results ---------------------------------------

for (i in 1:15){
  assign(
    x=paste0(
      "Results_", 
      i, 
      "_SVs"
      ),
    value = qread(
      paste0(
        "../Data/DA/WNN/AllCase/SVA/1_15_SVs/", 
        "DESeq_Results_", 
        i, 
        "SVs_Index", 
        ".qrds"
      ), 
      nthr=nthr
    )
  )
}


    ## 1.3 Load ColDicts -------------------------------------------------------

readxl::excel_sheets("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx")

ColDict_WNN_L1 <- setNames(
  object = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "WNN_L1")$Color, 
  nm = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "WNN_L1")$WNN_L1
)

ColDict_WNN_L2 <- setNames(
  object = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "WNN_L2")$Color, 
  nm = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "WNN_L2")$WNN_L2
)


ColDict_WNN_L3 <- setNames(
  object = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "WNN_L3")$Color, 
  nm = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "WNN_L3")$WNN_L3
)

ColDict_Case <- setNames(
  object = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "Case")$Color, 
  nm = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "Case")$Case
)

ColDict_CaseType <- setNames(
  object = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "Case_Type")$Color, 
  nm = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "Case_Type")$Case_Type
)




  ### 2.0 Define plotting functions --------------------------------------------

plot_covariate_results <- function(
    results_list, 
    comparisons, 
    X_Labels, 
    point_labels, 
    CellTypeLevel, 
    CellType, 
    stat, 
    returnData=FALSE, 
    exlude_labels_from_line=NULL, 
    labels_to_plot=NULL, 
    fill_scale=NULL, 
    col_scale=NULL, 
    ylab=NULL, 
    xlab=NULL
){
  
  df <- expand.grid(
    Comparison=comparisons, 
    Labels=X_Labels
  ) 
  
  df$CellTypeLevel <- CellTypeLevel 
  df$CellType <- CellType 
  
  df <- df[,c(3,4,2,1)] 
  
  df$PointLabels <- unlist(lapply(results_list, FUN=function(x){return(x[[point_labels]][x$CellTypeLevel==CellTypeLevel & x$CellType==CellType & x$Comparison %in% comparisons])})) 
  df[[stat]] <- unlist(lapply(results_list, FUN=function(x){return(x[[stat]][x$CellTypeLevel==CellTypeLevel & x$CellType==CellType & x$Comparison %in% comparisons])})) 
  
  if(returnData) return(df) 
  
  p1 <- ggplot(df) + 
    aes(Labels, .data[[stat]], fill=Comparison, label=paste0(PointLabels, " SVs")) + 
    geom_line(aes(group=Comparison, col=Comparison), data=if(is.null(exlude_labels_from_line)) {df} else {df[!df$Labels %in% exlude_labels_from_line,]}) + 
    geom_point(pch=21, size=3) + 
#   geom_text_repel(data = if(is.null(labels_to_plot)) {df} else {df[df$Labels %in% labels_to_plot,]}) + 
    {if(!is.null(fill_scale)) scale_fill_manual(values=fill_scale)} +  
    {if(!is.null(col_scale)) scale_color_manual(values=col_scale)} + 
    {if(!is.null(xlab)) xlab(xlab)} + 
    {if(!is.null(ylab)) ylab(ylab)} + 
    theme_classic() + 
    theme(
      axis.title.x = element_blank(), 
      axis.title.y = element_text(size=12, face="bold", color="#000000"),
      axis.text.x = element_text(size=10, face="bold", color="#000000", angle=45, hjust=1), 
      axis.text.y = element_text(size=10, face="bold", color="#000000", hjust=1) 
      
    ) 
  
  return(p1)
}




  ### 3.0 Summarize data -------------------------------------------------------

Results <- list() 

Results[["NoCov"]] <- Results_NoCov
Results[["NoCov"]]$SVs <- 0  

for (i in 1:15){
  Results[[paste0(i, "SVs")]] <- get(
    paste0(
      "Results_", 
      i, 
      "_SVs"
    )
  )

 Results[[paste0(i, "SVs")]]$SVs <- Results[[paste0(i, "SVs")]]$Covariates
}
rm(i)


p_ALS_All <- plot_covariate_results(
  results_list = Results, 
  comparisons = c("All_Cases", "Rand"), 
  X_Labels = c("NoCov", "1SV", "2SV", "3SV", "4SV", "5SV", "6SV", "7SV", "8SV", "9SV", "10SV", "11SV", "12SV", "13SV", "14SV", "15SV"), 
  CellTypeLevel = "WNN_L25", 
  CellType = "Inh_PVALB", 
  point_labels="SVs",
  stat="All_q_0.05_ALS", 
  returnData=FALSE, 
  fill_scale=c(as.character(ColDict_Case["ALS"]), "#BBBBBB"), 
  col_scale=c(as.character(ColDict_Case["ALS"]), "#BBBBBB"), 
  ylab="# of DARs (q<0.05)\n"
) + 
  theme(legend.position = "bottom", legend.title = element_blank()) + 
  scale_y_continuous(limits=c(0,1000)) + 
  scale_fill_manual(name="Comparison", labels=c("ALSvsHC", "Random control"), values=c(as.character(ColDict_Case["ALS"]), "#BBBBBB")) + 
  scale_color_manual(name="Comparison", labels=c("ALSvsHC", "Random control"), values=c(as.character(ColDict_Case["ALS"]), "#BBBBBB"))

p_ALSFTD_All <- plot_covariate_results(
    results_list = Results, 
    comparisons = c("All_Cases", "Rand"), 
    X_Labels = c("NoCov", "1SV", "2SV", "3SV", "4SV", "5SV", "6SV", "7SV", "8SV", "9SV", "10SV", "11SV", "12SV", "13SV", "14SV", "15SV"), 
    CellTypeLevel = "WNN_L25", 
    CellType = "Inh_PVALB", 
    point_labels="SVs",
    stat="All_q_0.05_ALSFTD", 
    returnData=FALSE, 
    fill_scale=c(as.character(ColDict_Case["ALS_FTD"]), "#BBBBBB"), 
    col_scale=c(as.character(ColDict_Case["ALS_FTD"]), "#BBBBBB"), 
    ylab="# of DARs (q<0.05)\n"
  ) + 
    theme(legend.position = "bottom", legend.title = element_blank()) + 
    scale_y_continuous(limits=c(0,1000))  + 
    scale_fill_manual(name="Comparison", labels=c("ALS_FTDvsHC", "Random control"), values=c(as.character(ColDict_Case["ALS_FTD"]), "#BBBBBB")) + 
    scale_color_manual(name="Comparison", labels=c("ALS_FTDvsHC", "Random control"), values=c(as.character(ColDict_Case["ALS_FTD"]), "#BBBBBB"))

    
p_ALS_All + p_ALSFTD_All
