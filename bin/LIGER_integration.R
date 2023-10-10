# /usr/bin/env R

# © EMBL-European Bioinformatics Institute, 2023
# Yuyao Song <ysong@ebi.ac.uk>

library(optparse)
library(rliger)
library(Seurat)
library(SeuratWrappers)

option_list <- list(
  make_option(c("-i", "--input_rds"),
    type = "character", default = NULL,
    help = "Path to input preprocessed rds file"
  ),
  make_option(c("-o", "--out_rds"),
    type = "character", default = NULL,
    help = "Output fastMNN from Seurat wrappers integrated rds file"
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

input_rds <- opt$input_rds
out_rds <- opt$out_rds
out_UMAP <- opt$out_UMAP
batch_key <- opt$batch_key
species_key <- opt$species_key
cluster_key <- opt$cluster_key

## create Seurat object via rds

# Convert(input_h5ad, dest = "rds", overwrite = TRUE)
# input_rds <- gsub("h5ad", "rds", input_h5ad)
obj <- readRDS(input_rds)

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


saveRDS(obj,
  file = out_rds
)

# iNMF embedding will be in .obsm['iNMF']

pdf(out_UMAP, height = 6, width = 10)
DimPlot(obj, reduction = "umap", group.by = species_key, shuffle = TRUE, label = TRUE)
DimPlot(obj, reduction = "umap", group.by = batch_key, shuffle = TRUE, label = TRUE)
DimPlot(obj, reduction = "umap", group.by = cluster_key, shuffle = TRUE, label = TRUE)

dev.off()
