# /usr/bin/env R

library(optparse)

# © EMBL-European Bioinformatics Institute, 2023
# Yuyao Song <ysong@ebi.ac.uk>

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
  ),
    make_option(c("--conda_path"),
    type = "character", default = NULL,
    help = "Conda for python executable to use for reticulate, important to match the prepared conda env!"
  )
)

# parse input
opt <- parse_args(OptionParser(option_list = option_list))

input_file <- opt$input_file
output_file <- opt$output_file
type <- opt$type
conda_path <- opt$conda_path

# set sys env before loading reticulate
Sys.setenv(RETICULATE_PYTHON=paste0(conda_path, "/bin/python3"))
Sys.setenv(RETICULATE_PYTHON_ENV=conda_path)

library(reticulate)
library(sceasy)
library(anndata)


if(type == 'anndata_to_seurat'){

    message(paste0("from anndata to seurat, input: ", input_file))

    sceasy::convertFormat(input_file, from="anndata", to="seurat",
                       outFile=output_file)
} else if (type == 'seurat_to_anndata'){

    message(paste0("from seurat to anndata, input: ", input_file))
    obj <- readRDS(input_file)
    dir_name = dirname(input_file)
    ## it is bizzar that from seurat to anndata needs loading object, but the other way works on-disk
    sceasy::convertFormat(obj, from="seurat", to="anndata",
                       outFile=output_file)

} else {

    warning("type must be either anndata_to_seurat or seurat_to_anndata")
}
