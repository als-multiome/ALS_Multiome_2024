

  ### 0.0 Load libraries ------------------------------------------------------- 

source("~/ALS_Brain_Multiome.Rcfg") 

library(qs)  
library(Seurat)
library(Signac)
library(tidyverse)
library(data.table)
library(ggraph)
library(pbapply)


qs::set_trust_promises(TRUE)
qs::set_trust_promises(TRUE)




  ### 1.0 Load data ------------------------------------------------------------

M0 <- qread(
  paste0(
    "../Data/SeuratObjects/",  
    "M0",
    ".qrds"
  ),
  nthr=nthr 
)

ColDict_WNN_L4 <- setNames(
  object = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "WNN_L4")$Color, 
  nm = readxl::read_xlsx("../Data/Visualization/ALS_Brain_Multiome_ColDicts.xlsx", sheet = "WNN_L4")$WNN_L4
)




  ### 2.0 Define auxiliary functions -------------------------------------------


BuildClusterTree2 <- function (object, assay = NULL, features = NULL, dims = NULL, 
                               reduction = "pca", graph = NULL, slot = "data", reorder = FALSE, 
                               reorder.numeric = FALSE, verbose = TRUE, .dist="euclidean", .method="complete") 
{
  if (!PackageCheck("ape", error = FALSE)) {
    stop(cluster.ape, call. = FALSE)
  }
  assay <- assay %||% DefaultAssay(object = object)
  if (!is.null(x = graph)) {
    idents <- levels(x = object)
    nclusters <- length(x = idents)
    data.dist <- matrix(data = numeric(length = 1L), nrow = nclusters, 
                        ncol = nclusters, dimnames = list(idents, idents))
    graph <- object[[graph]]
    cxi <- CellsByIdentities(object = object)
    cpairs <- na.omit(object = unique(x = t(x = apply(X = expand.grid(1:nclusters, 
                                                                      1:nclusters)[, c(2, 1)], MARGIN = 1, FUN = function(x) {
                                                                        if (length(x = x) == length(x = unique(x = x))) {
                                                                          return(sort(x = x))
                                                                        }
                                                                        return(c(NA, NA))
                                                                      }))))
    if (verbose) {
      pb <- txtProgressBar(style = 3, file = stderr())
    }
    for (i in 1:nrow(x = cpairs)) {
      i1 <- cpairs[i, ][1]
      i2 <- cpairs[i, ][2]
      graph.sub <- graph[cxi[[idents[i1]]], cxi[[idents[i2]]]]
      d <- mean(x = graph.sub)
      if (is.na(x = d)) {
        d <- 0
      }
      data.dist[i1, i2] <- d
      if (verbose) {
        setTxtProgressBar(pb = pb, value = i/nrow(x = cpairs))
      }
    }
    if (verbose) {
      close(con = pb)
    }
    diag(x = data.dist) <- 1
    if (.dist=="cosine") {data.dist <- as.dist(lsa::cosine(data.dist))} else{
      data.dist <- dist(x = data.dist, method=.dist)
    }
  }
  else if (!is.null(x = dims)) {
    my.lapply <- ifelse(test = verbose, yes = pblapply, no = lapply)
    embeddings <- Embeddings(object = object, reduction = reduction)[, 
                                                                     dims]
    data.dims <- my.lapply(X = levels(x = object), FUN = function(x) {
      cells <- WhichCells(object = object, idents = x)
      if (length(x = cells) == 1) {
        cells <- c(cells, cells)
      }
      temp <- colMeans(x = embeddings[cells, ])
    })
    data.dims <- do.call(what = "cbind", args = data.dims)
    colnames(x = data.dims) <- levels(x = object)
    if(.dist=="cosine"){data.dist <- as.dist(lsa::cosine(x=data.dims))}else{
      data.dist <- dist(x = t(x = data.dims), method=.dist)
    }
  }
  else {
    features <- features %||% VariableFeatures(object = object)
    features <- intersect(x = features, y = rownames(x = object))
    data.avg <- AverageExpression(object = object, assays = assay, 
                                  features = features, slot = slot, verbose = verbose)[[1]]
    if(.dist=="cosine"){data.dist <- as.dist(lsa::cosine(x=data.avg[features,]))}else{
      data.dist <- dist(x = t(x = data.avg[features, ]), method=.dist)
    }
  }
  data.tree <- ape::as.phylo(x = hclust(d = data.dist, method=.method))
  Tool(object = object) <- data.tree
  if (reorder) {
    if (verbose) {
      message("Reordering identity classes and rebuilding tree")
    }
    old.ident.order <- levels(x = object)
    data.tree <- Tool(object = object, slot = "BuildClusterTree")
    all.desc <- GetDescendants(tree = data.tree, node = (data.tree$Nnode + 
                                                           2))
    all.desc <- old.ident.order[all.desc[all.desc <= (data.tree$Nnode + 
                                                        1)]]
    Idents(object = object) <- factor(x = Idents(object = object), 
                                      levels = all.desc, ordered = TRUE)
    if (reorder.numeric) {
      new.levels <- sort(x = unique(x = as.integer(x = Idents(object = object))))
      Idents(object = object) <- factor(x = as.integer(x = Idents(object = object)), 
                                        levels = new.levels)
      object[["tree.ident"]] <- as.integer(x = Idents(object = object))
    }
    object <- BuildClusterTree(object = object, assay = assay, 
                               features = features, dims = dims, reduction = reduction, 
                               graph = graph, slot = slot, reorder = FALSE, verbose = verbose)
  }
  return(object)
}


PlotClusterTree2 <- function (object, direction = "downwards", nodelabels=FALSE, ...) 
{
  if (!rlang::is_installed("ape")) {
    stop(cluster.ape, call. = FALSE)
  }
  if (is.null(x = Tool(object = object, slot = "BuildClusterTree2"))) {
    stop("Phylogenetic tree does not exist, build using BuildClusterTree")
  }
  data.tree <- Tool(object = object, slot = "BuildClusterTree2")
  ape::plot.phylo(x = data.tree, direction = direction, show.node.label=FALSE, ...)
  if(nodelabels){ape::nodelabels()} 
}




  ### 3.0 Hierarchical cell-type clustering based on RNA expression ------------

M0$WNN_L4 <- droplevels(M0$WNN_L4)
Idents(M0) <- M0$WNN_L4


RNA_Tree <- BuildClusterTree2(
  M0, 
  assay="RNA",
  slot="data", 
  features=Features(M0, assay = "RNA", layer="data"), 
  .dist = "manhattan", 
  .method = "average",
  verbose = TRUE
)

data.tree <- Tool(object = RNA_Tree, slot = "BuildClusterTree2")

ColDict_All <- c(
  setNames(ColDict_WNN_L4, nm=str_replace_all(names(ColDict_WNN_L4), "_", "-")),
  setNames(rep("#FFFFFF", 69), nm=paste0("Node", 1:69))
)

(p1 <- ggraph(data.tree, "dendrogram") + 
  geom_edge_elbow(flipped=FALSE) + 
  geom_node_label(aes(label = name, fill=name), hjust = 0, repel = FALSE, fontface="bold.italic") + 
  scale_fill_manual(values=ColDict_All) + 
  coord_flip() + 
  scale_y_reverse(expand=c(0.5,0,1,0)) + 
  theme_void() + NoLegend() 
)

ggsave(
  filename = paste0(
    "../Data/Visualization/Figures/SupplFig_Hierarchical_CellTypes_Clustering/", 
    "Hierarchical_CellType_Clustering_Dendrogram_RNA", 
    ".pdf"
  ),
  p1,
  dpi=300, width = 3600, height = 4800, units = "px"
)

rm(data.tree, p1, RNA_Tree, ColDict_All)

sum(rowSums(M0@assays$RNA)>0)
dim(M0)
length(unique(M0@misc$Sample_data$ID))



  ### 4.0 Hierarchical cell-type clustering based on ATAC accessibility --------

M0$WNN_L4 <- droplevels(M0$WNN_L4)
Idents(M0) <- M0$WNN_L4
DefaultAssay(M0) <- "ATAC"

ATAC_Tree <- BuildClusterTree2(
  M0, 
  assay="ATAC",
  slot="data", 
  features=Features(M0, assay = "ATAC", layer="data"), 
  .dist = "manhattan", 
  .method = "average", 
  verbose = TRUE
)

data.tree <- Tool(object = ATAC_Tree, slot = "BuildClusterTree2")


ColDict_All <- c(
  setNames(ColDict_WNN_L4, nm=str_replace_all(names(ColDict_WNN_L4), "_", "-")),
  setNames(rep("#FFFFFF", 69), nm=paste0("Node", 1:69))
)

(p1 <- ggraph(data.tree, "dendrogram") + 
    geom_edge_elbow(flipped=FALSE) + 
    geom_node_label(aes(label = name, fill=name), hjust = 0, repel = FALSE, fontface="bold.italic") + 
    scale_fill_manual(values=ColDict_All) + 
    coord_flip() + 
    scale_y_reverse(expand=c(0.5,0,1,0)) + 
    theme_void() + NoLegend() 
)

ggsave(
  filename = paste0(
    "../Data/Visualization/Figures/SupplFig_Hierarchical_CellTypes_Clustering/", 
    "Hierarchical_CellType_Clustering_Dendrogram_ATAC", 
    ".pdf"
  ),
  p1,
  dpi=300, width = 3600, height = 4800, units = "px"
)



rm(data.tree, p1, ATAC_Tree, ColDict_All)

sum(rowSums(M0@assays$ATAC)>0)
dim(M0)
length(unique(M0@misc$Sample_data$ID))


  ### 5.0 Hierarchical cell-type clustering based on RNA & ATAC Latent space ---- 



    ## 5.1 Generate common reductions space with RNA and ATAC scVI latents ----- 

rownames(M0@reductions$scvirna) 
rownames(M0@reductions$scviatac)
all(rownames(M0@reductions$scvirna)==rownames(M0@reductions$scviatac))

common_reduction = cbind(M0@reductions$scvirna@cell.embeddings, M0@reductions$scviatac@cell.embeddings)

all(Cells(M0)==rownames(common_reduction))

M0@reductions$Latents <- M0@reductions$scvirna
M0@reductions$Latents@cell.embeddings <- common_reduction

rm(common_reduction)



    ## 5.2 Generate and plot hierarchical clustering tree ----------------------

Common_Tree <- BuildClusterTree2(
  M0, 
  dims=TRUE, 
  reduction = "Latents", 
  .dist = "manhattan", 
  .method = "average", 
  verbose = TRUE
)



Common_Tree <- Common_Tree@tools$BuildClusterTree


ColDict_All <- c(
  setNames(ColDict_WNN_L4, nm=names(ColDict_WNN_L4)),
  setNames(rep("#FFFFFF", 69), nm=paste0("Node", 1:69))
)

(p1 <- ggraph(Common_Tree, "dendrogram") + 
    geom_edge_elbow(flipped=FALSE) + 
    geom_node_label(aes(label = name, fill=name), hjust = 0, repel = FALSE, fontface="bold.italic") + 
    scale_fill_manual(values=ColDict_All) + 
    coord_flip() + 
    scale_y_reverse(expand=c(0.5,0,1,0)) + 
    theme_void() + NoLegend() 
)

ggsave(
  filename = paste0(
    "../Data/Visualization/Figures/SupplFig_Hierarchical_CellTypes_Clustering/", 
    "Hierarchical_CellType_Clustering_Dendrogram_Latents_Common_Space", 
    ".pdf"
  ),
  p1,
  dpi=300, width = 3600, height = 4800, units = "px"
)


rm(Common_Tree, p1, ColDict_All)





