# /usr/bin/env R

library(anndata)
library(tidyverse)
library(optparse)
library(rliger)
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)

option_list <- list(
  make_option(c("-i", "--input_h5Seurat"),
    type = "character", default = NULL,
    help = "Path to input preprocessed h5Seurat file"
  ),
  make_option(c("-o", "--out_h5Seurat"),
    type = "character", default = NULL,
    help = "Output fastMNN from Seurat wrappers integrated h5Seurat file"
  ),
  make_option(c("-p", "--out_UMAP"),
    type = "character", default = NULL,
    help = "Output UMAP after fastMNN integration"
  ),
  make_option(c("-b", "--batch_key"),
    type = "character", default = NULL,
    help = "Batch key identifier to integrate"
  ),
  make_option(c("-s", "--species_key"),
    type = "character", default = NULL,
    help = "Species key identifier"
  ),
  make_option(c("-c", "--cluster_key"),
    type = "character", default = NULL,
    help = "Cluster key for UMAP plotting"
  )
)

# parse input
opt <- parse_args(OptionParser(option_list = option_list))

input_h5Seurat <- opt$input_h5Seurat
out_h5Seurat <- opt$out_h5Seurat
out_UMAP <- opt$out_UMAP
batch_key <- opt$batch_key
species_key <- opt$species_key
cluster_key <- opt$cluster_key

## create Seurat object via h5Seurat

# Convert(input_h5ad, dest = "h5seurat", overwrite = TRUE)
# input_h5Seurat <- gsub("h5ad", "h5Seurat", input_h5ad)
obj <- LoadH5Seurat(input_h5Seurat)

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj, split.by = batch_key, do.center = FALSE)

# LIGER
obj <- RunOptimizeALS(obj, k = 30, lambda = 5, split.by = batch_key)
obj <- RunQuantileNorm(obj, split.by = batch_key)

obj <- FindNeighbors(obj, reduction = "iNMF", dims = 1:30)
obj <- FindClusters(obj, resolution = 0.4)
# Dimensional reduction and plotting
obj <- RunUMAP(obj, dims = 1:ncol(obj[["iNMF"]]), reduction = "iNMF", n_neighbors = 15L,  min_dist = 0.3)


# have to convert all factor to character, or when later converting to h5ad, the factors will be numbers
i <- sapply(obj@meta.data, is.factor)
obj@meta.data[i] <- lapply(obj@meta.data[i], as.character)


SaveH5Seurat(obj,
  filename = out_h5Seurat,
  overwrite = TRUE
)

Convert(out_h5Seurat, "h5ad", overwrite = TRUE)
# iNMF embedding will be in .obsm['iNMF']

pdf(out_UMAP, height = 6, width = 10)
DimPlot(obj, reduction = "umap", group.by = species_key, shuffle = TRUE, label = TRUE)
DimPlot(obj, reduction = "umap", group.by = batch_key, shuffle = TRUE, label = TRUE)
DimPlot(obj, reduction = "umap", group.by = cluster_key, shuffle = TRUE, label = TRUE)

dev.off()
