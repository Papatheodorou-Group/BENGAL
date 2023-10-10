# /usr/bin/env R

library(Seurat)
library(optparse)


option_list <- list(
  make_option(c("-i", "--input_rds"),
    type = "character", default = NULL,
    help = "Path to input preprocessed rds file"
  ),
  make_option(c("-o", "--out_rds"),
    type = "character", default = NULL,
    help = "Output Seurat CCA integrated rds file"
  ),
  make_option(c("-p", "--out_UMAP"),
    type = "character", default = NULL,
    help = "Output UMAP after Seurat CCA integration"
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
input <- readRDS(input_rds)

object.list <- SplitObject(input, split.by = batch_key)

# normalize and identify variable features for each dataset independently
object.list <- lapply(X = object.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = object.list)

input.anchors <- FindIntegrationAnchors(object.list = object.list, anchor.features = features)

input.combined <- IntegrateData(anchorset = input.anchors)

DefaultAssay(input.combined) <- "integrated"

input.combined <- ScaleData(input.combined, verbose = FALSE)
input.combined <- RunPCA(input.combined, npcs = 50, verbose = FALSE)
input.combined <- RunUMAP(input.combined, reduction = "pca", dims = 1:30, n_neighbors = 15L,  min_dist = 0.3)
input.combined <- FindNeighbors(input.combined, reduction = "pca", dims = 1:30)
input.combined <- FindClusters(input.combined, resolution = 0.4)

saveRDS(input.combined,
  file = out_rds
)

pdf(out_UMAP, height = 6, width = 10)
DimPlot(input.combined, reduction = "umap", group.by = species_key, shuffle = TRUE, label = TRUE)
DimPlot(input.combined, reduction = "umap", group.by = batch_key, shuffle = TRUE, label = TRUE)
DimPlot(input.combined, reduction = "umap", group.by = cluster_key, shuffle = TRUE, label = TRUE)

dev.off()
