#!/usr/bin/env R

##########
## concatenate anndata object by homology for LIGER UINMF pipeline
## ysong@ebi.ac.uk for CrossSpeciesIntegration pipeline
#########

# © EMBL-European Bioinformatics Institute, 2023
# Yuyao Song <ysong@ebi.ac.uk>

library(optparse)
library(anndata)
library(dplyr)
library(purrr)
library(readr)
library(magrittr)
library(tibble)
library(biomaRt)

option_list <- list(
  make_option(c("--metadata"),
    type = "character", default = NULL,
    help = "Path to a file indicate species-h5ad mapping, tab-seperated"
  ),
  make_option(c("--homology_tbl"),
      type = "character", default = NULL,
      help = "Ensembl homology tbl output"
  ),
  make_option(c("--out_dir"),
    type = "character", default = NULL,
    help = "output dir to write h5ad files"
  ),
  make_option(c("--metadata_output"),
    type = "character", default = NULL,
    help = "paths to rliger UINMF ready metadata"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

input_metadata <- opt$metadata
homology_tbl_csv <- opt$homology_tbl
out_dir <- opt$out_dir
metadata_output <- opt$metadata_output


message("create seurat objects for rliger")


genes_main_chr <- list(
  hsapiens = as.character(c(lapply(seq(1, 22, by = 1), as.character), "X", "Y")),
  sscrofa = as.character(c(lapply(seq(1, 18, by = 1), as.character), "X", "Y")),
  mmusculus = as.character(c(lapply(seq(1, 19, by = 1), as.character), "X", "Y")),
  mmulatta = as.character(c(lapply(seq(1, 20, by = 1), as.character), "X", "Y")),
  mfascicularis = as.character(c(lapply(seq(1, 20, by = 1), as.character), "X")),
  dmelanogaster = c("2L", "2R", "3L", "3R", "4", "X", "Y"),
  xtropicalis = as.character(c(lapply(seq(1, 10, by = 1), as.character))),
  drerio = as.character(c(lapply(seq(1, 25, by = 1), as.character)))
)

metadata <- read_tsv(input_metadata, col_names = FALSE)
species_list=metadata[['X1']]
message(paste0("start concatenating anndata object from ", paste(species_list, collapse = ', ')))


adatas = list()
for( species in metadata$X1){
    adatas[[species]] = read_h5ad(as.character(metadata[which(metadata$X1 == species), 'X2']))
}

for(species_now in species_list){
    adatas[[species_now]]$var[[paste0(species_now, "_homolog_ensembl_gene")]] = adatas[[species_now]]$var_names
}


count_all <- data.frame()
species_1 = species_list[1]


## version = '105'
mart <- useEnsembl("ensembl", version = "105", dataset = as.character(paste(species_1, "_gene_ensembl", sep = "")))

# get genes in the main chrs of the first species
genes_species_1 <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"), mart = mart) %>%
  dplyr::filter(chromosome_name %in% (genes_main_chr[species_1] %>% purrr::flatten_chr()))
message(paste("start querying ensembl biomaRt for gene homology"))

avail_attributes <- listAttributes(mart) %>% filter(grepl((paste(species_list[-1],collapse="|")), name)) %>% filter(grepl("homo", name)) %>%
filter(!grepl("Query protein or transcript ID", description))
biomartCacheClear()

homology_tbl = getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position",
    avail_attributes$name), mart = mart, filters = "ensembl_gene_id",
  values = genes_species_1[["ensembl_gene_id"]])


#names(homology_tbl)[colnames(homology_tbl) == 'Gene stable ID'] = paste0(species_1, "_homolog_ensembl_gene")
#names(homology_tbl)[colnames(homology_tbl) == 'Gene name'] = paste0(species_1, "_homolog_associated_gene_name")
#names(homology_tbl)[colnames(homology_tbl) == 'Chromosome/scaffold name'] = paste0(species_1, "_homolog_chromosome")
#names(homology_tbl)[colnames(homology_tbl) == 'Gene start (bp)'] = paste0(species_1, "_homolog_chrom_start")
#names(homology_tbl)[colnames(homology_tbl) == 'Gene end (bp)'] = paste0(species_1, "_homolog_chrom_end")

#names(homology_tbl)[match(avail_attributes[,"description"], names(homology_tbl))] = avail_attributes[,"name"]

write_csv(homology_tbl, file = homology_tbl_csv)


message("start building one2one")
homology_tbl[paste0(species_1, "_homolog_associated_gene_name")] = homology_tbl$external_gene_name
homology_tbl[paste0(species_1, "_homolog_ensembl_gene")] = homology_tbl$ensembl_gene_id
homology_tbl[paste0(species_1, "_homolog_chromosome")] = homology_tbl$chromosome_name
homology_tbl[paste0(species_1, "_homolog_chrom_start")] = homology_tbl$start_position
homology_tbl[paste0(species_1, "_homolog_chrom_end")] = homology_tbl$end_position
one2one = homology_tbl %>% filter_at(vars(ends_with("homolog_orthology_type")), all_vars(. == 'ortholog_one2one')) %>%
distinct(get(paste0(species_1, "_homolog_ensembl_gene")), `.keep_all` = TRUE)

adatas_one2one = list()

for(species_now in names(adatas)){
    adatas_one2one[[species_now]] = adatas[[species_now]][, tolower(adatas[[species_now]]$var_names) %in% tolower(one2one[[paste0(species_now, "_homolog_ensembl_gene")]])]
}


for(species_now in species_list[-1]){
    message(species_now)
     one2one_now = one2one %>% filter(tolower(get(paste0(species_now, "_homolog_ensembl_gene"))) %in% tolower(adatas_one2one[[species_now]]$var_names)) %>%
distinct(get(paste0(species_1, "_homolog_ensembl_gene")), `.keep_all` = TRUE)
     adatas_one2one[[species_now]] = adatas_one2one[[species_now]][, tolower(adatas_one2one[[species_now]]$var_names) %in% tolower(one2one_now[[paste0(species_now, "_homolog_ensembl_gene")]])]
     one2one_now = one2one_now[match(tolower(adatas_one2one[[species_now]]$var_names), tolower(one2one_now[[paste0(species_now, '_homolog_ensembl_gene')]])), ]
     adatas_one2one[[species_now]]$var[[paste0(species_now, '_homolog_ensembl_gene')]] = adatas_one2one[[species_now]]$var_names
     adatas_one2one[[species_now]]$var[[paste0(species_1, '_homolog_ensembl_gene')]] = one2one_now[[paste0(species_1, '_homolog_ensembl_gene')]]
     adatas_one2one[[species_now]]$var_names = one2one_now[[paste0(species_1, '_homolog_ensembl_gene')]]
     adatas_one2one[[species_now]] = adatas_one2one[[species_now]][, !duplicated(adatas_one2one[[species_now]]$var_names)]
}

adatas_one2one_uinmf = list()

for(species_now in names(adatas)){

    adata_unshared = adatas[[species_now]][, !(adatas[[species_now]]$var[[paste0(species_now, '_homolog_ensembl_gene')]] %in% adatas_one2one[[species_now]]$var[[paste0(species_now, '_homolog_ensembl_gene')]])]
    adatas_one2one_uinmf[[species_now]] = concat(list(adatas_one2one[[species_now]], adata_unshared), axis = 1L, join = 'inner', merge = 'first', index_unique = '-')

    adatas_one2one_uinmf[[species_now]] = adatas_one2one_uinmf[[species_now]][, !duplicated(adatas_one2one_uinmf[[species_now]]$var_names)]

    i <- sapply(adatas_one2one_uinmf[[species_now]]$obs, is.factor)
    adatas_one2one_uinmf[[species_now]]$obs[i] <- lapply(adatas_one2one_uinmf[[species_now]]$obs[i], as.character)
    j <- sapply(adatas_one2one_uinmf[[species_now]]$var, is.factor)
    adatas_one2one_uinmf[[species_now]]$var[j] <- lapply(adatas_one2one_uinmf[[species_now]]$var[j], as.character)

    write_h5ad(anndata = adatas_one2one_uinmf[[species_now]],
               filename = paste0(out_dir, "/", species_now, "_one2one_only_ligerUINMF.h5ad"), compression = "gzip")
    message(species_now)

    }

message("finish building one2one")


message("start building many2many")
message("start higher expression level")


many2many = homology_tbl %>% filter_at(vars(ends_with("homolog_orthology_type")), all_vars(. != 'ortholog_one2one'))


for(species_now in species_list){
    print(species_now)
    many2many = many2many %>% filter(!is.na(get(paste0(species_now, "_homolog_ensembl_gene"))) & get(paste0(species_now, "_homolog_ensembl_gene")) != "")
    print(dim(many2many))
    many2many = many2many %>%
    filter(tolower(get(paste0(species_now, "_homolog_ensembl_gene"))) %in% tolower(adatas[[species_now]]$var[[paste0(species_now, "_homolog_ensembl_gene")]]))
    print(dim(many2many))
}

many2many_copy <- many2many %>% rowid_to_column("index")
adata_many2many = AnnData(shape = list(0, 0))


adatas_many2many_all = list()
while (nrow(many2many_copy) > 0) {

    if(nrow(many2many_copy) < 6000 & nrow(many2many_copy) > 5990) message(nrow(many2many_copy))
    
    if(nrow(many2many_copy) < 4000 & nrow(many2many_copy) > 3990) message(nrow(many2many_copy))

    if(nrow(many2many_copy) < 1000 & nrow(many2many_copy) > 990) message(nrow(many2many_copy))


    dd <- many2many_copy %>%
        filter(get(paste0(species_1, "_homolog_ensembl_gene"))  == levels(factor(many2many_copy[[paste0(species_1, "_homolog_ensembl_gene")]]))[1])

    genes_now = dd %>% dplyr::select(ends_with("_homolog_ensembl_gene")) %>% flatten() %>% unique() %>% as.character()

    gene_group <- many2many_copy %>%
        dplyr::filter_at(vars(ends_with("_homolog_ensembl_gene")), any_vars(. %in% genes_now))

    many2many_copy <- many2many_copy %>%
    filter(!index %in% gene_group$index)

    adatas_many2many = list()

    for(species_now in species_list){
            adatas_many2many[[species_now]] = adatas[[species_now]][, tolower(adatas[[species_now]]$var_names) %in% tolower(gene_group[[paste0(species_now, "_homolog_ensembl_gene")]])]
            keep_row = adatas_many2many[[species_now]]$var %>%
            arrange(desc(mean_counts)) %>%
            slice(1)
            adatas_many2many[[species_now]] = adatas_many2many[[species_now]][, which(tolower(adatas_many2many[[species_now]]$var_names) == tolower(rownames(keep_row)))]
        }

    new_name = adatas_many2many[[species_1]]$var_names
    #message(new_name)

    for(species_now in species_list[-1]){
    adatas_many2many[[species_now]]$var[[paste0(species_1, "_homolog_ensembl_gene")]] = new_name
    rownames(adatas_many2many[[species_now]]$var) = new_name
}

    for(species_now in species_list){
    if (is.null(adatas_many2many_all[[species_now]])) {
    adatas_many2many_all[[species_now]] <- adatas_many2many[[species_now]]
  } else {
    adatas_many2many_all[[species_now]] <- concat(list(adatas_many2many_all[[species_now]], adatas_many2many[[species_now]]), axis = 1L, join = "outer", merge = "first", index_unique = '-')
  }

}

  rm(adatas_many2many)
}


adatas_one2one_and_many_expr = list()

for(species_now in names(adatas)){

    adatas_one2one_and_many_expr[[species_now]] = concat(list(adatas_one2one[[species_now]], adatas_many2many_all[[species_now]]), axis = 1L, join = 'inner', merge = 'first', index_unique = '-')
    message(species_now)

    }

rm(adatas_many2many_all)

adatas_one2one_higher_expr_uinmf = list()

for(species_now in names(adatas)){

    adata_unshared = adatas[[species_now]][, !(adatas[[species_now]]$var[[paste0(species_now, '_homolog_ensembl_gene')]] %in% adatas_one2one_and_many_expr[[species_now]]$var[[paste0(species_now, '_homolog_ensembl_gene')]])]
    adatas_one2one_higher_expr_uinmf[[species_now]] = concat(list(adatas_one2one_and_many_expr[[species_now]], adata_unshared), axis = 1L, join = 'inner', merge = 'first', index_unique = '-')

    adatas_one2one_higher_expr_uinmf[[species_now]]= adatas_one2one_higher_expr_uinmf[[species_now]][, !duplicated(adatas_one2one_higher_expr_uinmf[[species_now]]$var_names)]

    i <- sapply(adatas_one2one_higher_expr_uinmf[[species_now]]$obs, is.factor)
    adatas_one2one_higher_expr_uinmf[[species_now]]$obs[i] <- lapply(adatas_one2one_higher_expr_uinmf[[species_now]]$obs[i], as.character)
    j <- sapply(adatas_one2one_higher_expr_uinmf[[species_now]]$var, is.factor)
    adatas_one2one_higher_expr_uinmf[[species_now]]$var[j] <- lapply(adatas_one2one_higher_expr_uinmf[[species_now]]$var[j], as.character)

    write_h5ad(anndata = adatas_one2one_higher_expr_uinmf[[species_now]],
               filename = paste0(out_dir, "/", species_now, "_many_higher_expr_ligerUINMF.h5ad"), compression = "gzip")
    message(species_now)

    }
message("finish building higher expression level")

rm(adatas_one2one_and_many_expr)
rm(adatas_one2one_higher_expr_uinmf)


message("start higher homology confidence")

## available attributes to indicate confidence of homology
adatas_many2many_homo_all = list()

## mind the order
order = metadata$X1
avail_ordered = c()
for (attr in c("orthology_confidence", "homolog_goc_score", "homolog_wga_coverage")){

    avail_homo = c(avail_attributes$name[grepl(attr, avail_attributes$name)])

    for (i in seq(1, length(order))){
   avail_ordered = c(avail_ordered, avail_homo[grepl(order[i], avail_homo)])
  }
}

avail_homo = avail_ordered

many2many_copy_homo <- many2many %>% rowid_to_column("index")
adata_many2many_homo = AnnData(shape = list(0, 0))

while (nrow(many2many_copy_homo) > 0) {
    

    if(nrow(many2many_copy_homo) < 6000 & nrow(many2many_copy_homo) > 5990) message(nrow(many2many_copy_homo))

    if(nrow(many2many_copy_homo) < 4000 & nrow(many2many_copy_homo) > 3990) message(nrow(many2many_copy_homo))

    if(nrow(many2many_copy_homo) < 1000 & nrow(many2many_copy_homo) > 990) message(nrow(many2many_copy_homo))


    dd <- many2many_copy_homo %>%
        filter(get(paste0(species_1, "_homolog_ensembl_gene"))  == levels(factor(many2many_copy_homo[[paste0(species_1, "_homolog_ensembl_gene")]]))[1])

    genes_now = dd %>% dplyr::select(ends_with("_homolog_ensembl_gene")) %>% flatten() %>% unique() %>% as.character()

    ## find a gene group with many-to-many relationships
    gene_group <- many2many_copy_homo %>%
        dplyr::filter_at(vars(ends_with("_homolog_ensembl_gene")), any_vars(. %in% genes_now))

    if(nrow(gene_group) == 1) {
        gene_keep = gene_group
    } else {
    gene_keep = gene_group %>%
        dplyr::arrange(
        sapply(avail_homo, FUN = function(x) get(x)) ## keep the member of group with highest overall confidence
        ) %>%
            slice(n())
    }

    many2many_copy_homo <- many2many_copy_homo %>%
            filter(!index %in% gene_group$index)

    adatas_many2many_homo = list()
    for(species_now in species_list){

    adatas_many2many_homo[[species_now]] = adatas[[species_now]][, adatas[[species_now]]$var_names %in% gene_keep[[paste0(species_now, "_homolog_ensembl_gene")]]]

}

    new_name = adatas_many2many_homo[[species_1]]$var_names
    #message(new_name)
    for(species_now in species_list[-1]){
    adatas_many2many_homo[[species_now]]$var[[paste0(species_1, "_homolog_ensembl_gene")]] = new_name
    rownames(adatas_many2many_homo[[species_now]]$var) = new_name
}

    for(species_now in species_list){
    if (is.null(adatas_many2many_homo_all[[species_now]])) {
    adatas_many2many_homo_all[[species_now]] <- adatas_many2many_homo[[species_now]]
  } else {
    adatas_many2many_homo_all[[species_now]] <- concat(list(adatas_many2many_homo_all[[species_now]], adatas_many2many_homo[[species_now]]), axis = 1L, join = "outer", merge = "first", index_unique = '-')
  }
}
  rm(adatas_many2many_homo)

}
adatas_one2one_and_many_homo = list()
for(species_now in names(adatas)){

    adatas_one2one_and_many_homo[[species_now]] = concat(list(adatas_one2one[[species_now]], adatas_many2many_homo_all[[species_now]]), axis = 1L, join = 'inner', merge = 'first', index_unique = '-')
    message(species_now)

    }

rm(adatas_many2many_homo_all)


adatas_one2one_higher_homo_uinmf = list()

for(species_now in names(adatas)){

    adata_unshared = adatas[[species_now]][, !(adatas[[species_now]]$var[[paste0(species_now, '_homolog_ensembl_gene')]] %in% adatas_one2one_and_many_homo[[species_now]]$var[[paste0(species_now, '_homolog_ensembl_gene')]])]
    adatas_one2one_higher_homo_uinmf[[species_now]] = concat(list(adatas_one2one_and_many_homo[[species_now]], adata_unshared), axis = 1L, join = 'inner', merge = 'first', index_unique = '-')

    adatas_one2one_higher_homo_uinmf[[species_now]]= adatas_one2one_higher_homo_uinmf[[species_now]][, !duplicated(adatas_one2one_higher_homo_uinmf[[species_now]]$var_names)]

    i <- sapply(adatas_one2one_higher_homo_uinmf[[species_now]]$obs, is.factor)
    adatas_one2one_higher_homo_uinmf[[species_now]]$obs[i] <- lapply(adatas_one2one_higher_homo_uinmf[[species_now]]$obs[i], as.character)
    j <- sapply(adatas_one2one_higher_homo_uinmf[[species_now]]$var, is.factor)
    adatas_one2one_higher_homo_uinmf[[species_now]]$var[j] <- lapply(adatas_one2one_higher_homo_uinmf[[species_now]]$var[j], as.character)

    write_h5ad(anndata = adatas_one2one_higher_homo_uinmf[[species_now]],
               filename = paste0(out_dir, "/", species_now, "_many_higher_homology_conf_ligerUINMF.h5ad"), compression = "gzip")
    message(species_now)

    }
message("finish higher homology confidence")

rm(adatas_one2one_and_many_homo)
rm(adatas_one2one_higher_homo_uinmf)



metadata = data.frame('species' = names(adatas))
metadata$one2one_only = paste0(out_dir, "/", metadata$species, "_one2one_only_ligerUINMF.rds")
metadata$many_higher_expr = paste0(out_dir, "/", metadata$species, "_many_higher_expr_ligerUINMF.rds")
metadata$many_higher_homology_conf = paste0(out_dir, "/", metadata$species, "_many_higher_homology_conf_ligerUINMF.rds")
write_tsv(metadata, file = metadata_output, col_names = TRUE)
