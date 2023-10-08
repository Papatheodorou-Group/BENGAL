# /usr/bin/env R

library(sceasy)
library(anndata)
library(optparse)
library(SeuratDisk)


option_list <- list(
  make_option(c("-i", "--input_file"),
    type = "character", default = NULL,
    help = "Path to input file for convrting"
  ),
  make_option(c("-o", "--output_file"),
    type = "character", default = NULL,
    help = "Output file after conversion"
  ),
  make_option(c("-t", "--type"),
    type = "character", default = NULL,
    help = "Conversion type, choose between anndata_to_seurat or seurat_to_anndata"
  )
)

# parse input
opt <- parse_args(OptionParser(option_list = option_list))

input_file <- opt$input_file
output_file <- opt$output_file
type <- opt$type


if(type == 'anndata_to_seurat'){

    message(paste0("from anndata to seurat, input: ", input_file))

    obj <- anndata::read_h5ad(input_file)

    sceasy::convertFormat(obj, from="anndata", to="seurat",
                       outFile=output_file)
} else if (type == 'seurat_to_anndata'){

    message(paste0("from seurat to anndata, input: ", input_file))
    
    obj <- readRDS(input_file)

    sceasy::convertFormat(obj, from="seurat", to="anndata",
                       outFile=output_file)

} else {

    warning("type must be either anndata_to_seurat or seurat_to_anndata")
}
