#!/usr/bin/env R

# © EMBL-European Bioinformatics Institute, 2023
# Yuyao Song <ysong@ebi.ac.uk>
library(anndata)
library(mclust)
library(cluster)
library(optparse)
library(kBET)
library(tidyverse)
library(FNN)


option_list <- list(
  make_option(c("-i", "--input_h5ad"),
    type = "character", default = NULL,
    help = "Path to input integrated h5ad file"
  ),
  make_option(c("-o", "--out_csv"),
    type = "character", default = NULL,
    help = "Output csv file with various batch effect measurements"
  ),
  make_option(c("-m", "--method"),
    type = "character", default = NULL,
    help = "Integration method used, have an effect on use embedding or use count"
  ),
  make_option(c("-b", "--batch_key"),
    type = "character", default = NULL,
    help = "Batch key identifier to integrate"
  ),
  make_option(c("-c", "--cluster_key"),
    type = "character", default = NULL,
    help = "Cell type cluster key"
  ),
  make_option(c("-s", "--species_key"),
    type = "character", default = NULL,
    help = "Species key identifier"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

input_h5ad <- opt$input_h5ad
out_csv <- opt$out_csv
method <- opt$method
batch_key <- opt$batch_key
cluster_key <- opt$cluster_key
species_key <- opt$species_key

message("read in h5ad")
ad <- read_h5ad(input_h5ad)

## identify shared cell types for measuring batch effect
species_all = levels(factor(ad$obs[[species_key]]))
cell_types = list()
for(species in species_all){
    cell_types[[species]] = levels(factor(ad$obs[ad$obs[[species_key]] == species, cluster_key]))
}

ct_shared=Reduce(intersect, cell_types)
message(paste0("shared cell types include ", paste(ct_shared, collapse = ", ")))

ad <- ad[ad$obs[[cluster_key]] %in% ct_shared, ]

# get the matrix to compute batch effect
# these methods return embedding
if (method == "harmony") {
  data <- data.frame(ad$obsm[["X_pca_harmony"]], row.names = ad$obs_names)
  do_pca <- FALSE
}
if (method == "scanorama") {
  data <- data.frame(ad$obsm[["X_scanorama"]], row.names = ad$obs_names)
  do_pca <- FALSE
}
if (method == "scVI") {
  data <- data.frame(ad$obsm[["X_scVI"]], row.names = ad$obs_names)
  do_pca <- FALSE
}
if (method == "LIGER") {
  data <- data.frame(ad$obsm[["X_iNMF"]], row.names = ad$obs_names)
  do_pca <- FALSE
}
if (method == "rligerUINMF") {
  data <- data.frame(ad$obsm[["X_inmf"]], row.names = ad$obs_names)
  do_pca <- FALSE
}
if (method == "fastMNN") {
  data <- data.frame(ad$obsm[["X_mnn"]], row.names = ad$obs_names)
  do_pca <- FALSE
}

if (method == "SAMap") {
  data <- data.frame(ad$obsm[["wPCA"]], row.names = ad$obs_names)
  do_pca <- FALSE
}

# seurat return a pseudo-count matrix that is after the integration
if (method %in% c("seuratCCA", "seuratRPCA", "unintegrated")) {
  data <- data.frame(as.matrix((ad$X)), row.names = ad$obs_names)
  do_pca <- TRUE
} # sparse matrix to dense

batch <- ad$obs[[batch_key]]

# if only 2 batches use pval in linear model, if multiple batches use pval in ANOVA F test in PC regression
if (length(levels(factor(as.character(batch)))) > 2) {
  use_pval <- "p.value.F.test"
} else {
  use_pval <- "p.value.lm"
}

kbet_per_ct <- data.frame()
for (ct in levels(factor(ad$obs[[cluster_key]]))) {
  message(ct)
  ad_ct <- ad[ad$obs[[cluster_key]] == ct, ]
  data_ct <- data[ad_ct$obs_names, ]
  batch <- ad_ct$obs[[batch_key]]
  k0 <- floor(mean(table(batch))) # neighbourhood size: mean batch size
  knn <- get.knn(data_ct, k = k0, algorithm = "cover_tree")
  # now run kBET with pre-defined nearest neighbours.
  batch.estimate <- suppressWarnings(kBET(data_ct, batch, k = k0, knn = knn, plot = FALSE, do.pca = do_pca))

  if (any(is.na(batch.estimate))) {
    message(paste("cell type ", ct, " have less than 10 cells, skip cell-type kBET"), sep = "")
    next
  }

  add <- batch.estimate$summary["mean", ]
  add$cell_type <- ct

  if(nrow(kbet_per_ct) == 0){

      kbet_per_ct <- add %>% pivot_longer(cols = 1:3, names_to = "measure", values_to = "value")

  } else {
  kbet_per_ct <- rbind(kbet_per_ct, add %>% pivot_longer(cols = 1:3, names_to = "measure", values_to = "value"))
  }
}


message("write summary")

summary <- kbet_per_ct %>% pivot_wider(id_cols = cell_type, names_from = measure, values_from = value)
summary$method <- method
summary$input_h5ad <- input_h5ad
write_csv(summary, out_csv)
