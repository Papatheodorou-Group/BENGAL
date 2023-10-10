# /usr/bin/env R

# © EMBL-European Bioinformatics Institute, 2023
# Yuyao Song <ysong@ebi.ac.uk>

library(optparse)
library(rliger)
library(scCustomize) # liger to seurat function keep_metadata
library(SeuratDisk)


option_list <- list(
  make_option(c("--metadata"),
    type = "character", default = NULL,
    help = "Path to a file indicate species-rds mapping, tab-seperated"
  ),
  make_option(c("--basename"),
      type = "character", default = NULL,
      help = "Basename of file to save"
  ),
  make_option(c("--out_dir"),
    type = "character", default = NULL,
    help = "output dir to write rds files"
  ),
  make_option(c("--cluster_key"),
    type = "character", default = NULL,
    help = "paths to rliger UINMF ready metadata"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

metadata_path <- opt$metadata
basename <- opt$basename
out_dir <- opt$out_dir
cluster_key <- opt$cluster_key


metadata = read.table(metadata_path, sep = '\t', header=TRUE)
metadata = as.data.frame(metadata)

for(type in colnames(metadata)[-1]){

message(type)

obj_list = list()
for(i in seq(1, nrow(metadata))){

    obj_list[[i]] = readRDS(metadata[i, type])
}

liger <- seuratToLiger(obj_list, remove.missing = FALSE, renormalize = FALSE, names = metadata$species)

meta_all = data.frame()
keep_cols = Reduce(intersect, lapply(obj_list, FUN = function(x) colnames(x@meta.data)))

for(i in seq(1, nrow(metadata))){
    message(i)
    obj_list[[i]]@meta.data$cell_id = rownames(obj_list[[i]]@meta.data)
    use = obj_list[[i]]@meta.data[, c("cell_id", keep_cols)]
    meta_all = rbind(meta_all, use)

}

meta_new = meta_all[match(rownames(liger@cell.data), rownames(meta_all)), ]
meta_new = cbind(meta_new, liger@cell.data)
meta_new = meta_new[match(rownames(liger@cell.data), rownames(meta_new)), ]

liger@cell.data <- meta_new

species.liger <- normalize(liger)
species.liger <- rliger::selectGenes(species.liger, var.thres = 0.3, unshared = TRUE, unshared.datasets = list(1, 2), unshared.thresh = 0.3)
species.liger <- scaleNotCenter(species.liger)
species.liger <- optimizeALS(species.liger, lambda = 5, use.unshared = TRUE, thresh = 1e-10, k = 30)
species.liger <- quantile_norm(species.liger, ref_dataset = metadata$species[1])
species.liger <- louvainCluster(species.liger)
species.liger <- runUMAP(species.liger)

seurat_obj <- scCustomize::Liger_to_Seurat(
  species.liger,
  nms = names(species.liger@H),
  renormalize = TRUE,
  use.liger.genes = TRUE,
  by.dataset = FALSE,
  keep_meta = TRUE,
  reduction_label = "umap", # in line with the X_umap in scanpy
  seurat_assay = "RNA"
)

k <- sapply(seurat_obj@meta.data, is.factor)
seurat_obj@meta.data[k] <- lapply(seurat_obj@meta.data[k], as.character) # known bug in seurat to h5ad for factors

message("save pdf")
pdf(paste0(out_dir,"/", basename, "_", type, "_rligerUINMF_integrated_UMAP.pdf"), height = 10, width = 12)
DimPlot(object = seurat_obj, reduction = 'umap', group.by = c('species', 'cell_ontology_mapped'), ncol = 1)
dev.off()


message("save seurat object")

saveRDS(seurat_obj,
  file = paste0(out_dir,"/", basename, "_", type, "_rligerUINMF_integrated.rds")
)

message(paste0("finish", type))

}

