#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Updated upon containerization of BENGAL envs
// Yuyao Song <ysong@ebi.ac.uk>
// Oct 2023

// data requirements: batch_key, species_key, sample_key, mapped in config file, present in adata.obs
// data requirements: mean_count in adata.var from scanpy QC
// raw h5ad file naming: <species>
// adata.var_names are ensembl gene ids

log.info """
         ===========================================================


         Cross-species integration and assessment - nextflow pipeline
         Use singularity containers for cluster execution
         - check inout format
         - concatenate input anndata
         Author: ysong@ebi.ac.uk
         Initial date: Mar 2022
         Latest date: Oct 2023

         ===========================================================

         """
         .stripIndent()



process validate_adata_input {

    // labels for cluster options and containers
    // no need to change here, adjust as per dataset in config file

    label 'validate'
    label 'regular_resource'

    input:
    tuple val(basename), path(metadata)

    script:
    """
    python ${projectDir}/bin/validate_input.py ${metadata} --batch_key ${params.batch_key} --species_key ${params.species_key} --cluster_key ${params.cluster_key}

    """

}


process concat_by_homology {

    label 'concat'
    label 'regular_resource'

    publishDir "${params.results}/results/h5ad_homology_concat", mode: 'copy'

    input:

    tuple val(basename), path(metadata)

    output:

    path "${basename}_one2one_only.h5ad"
    path "${basename}_many_higher_expr.h5ad"
    path "${basename}_many_higher_homology_conf.h5ad"
    path "${basename}_homology_tbl.csv"


    script:
    """
    Rscript ${projectDir}/bin/concat_by_homology_multiple_species_by_gene_id.R --metadata ${metadata} \
    --one2one_h5ad ${basename}_one2one_only.h5ad --many_higher_expr_h5ad ${basename}_many_higher_expr.h5ad \
    --many_higher_homology_conf_h5ad ${basename}_many_higher_homology_conf.h5ad --homology_tbl ${basename}_homology_tbl.csv
    """
}

process concat_by_homology_rliger_uinmf {

    label 'concat'
    label 'regular_resource'

    publishDir "${params.results}/results/rligerUINMF/h5ad_homology_concat", mode: 'copy'

    input:

    tuple val(basename), path(metadata)

    output:
    path "homology_tbl.csv"
    path "${basename}_liger.tsv"
    //path "*_ligerUINMF.h5ad"

    script:
    """
    mkdir -p ${params.results}/results/rligerUINMF/h5ad_homology_concat && \
    Rscript ${projectDir}/bin/concat_by_homology_rligerUINMF_multiple_species.R --metadata ${metadata} \
    --out_dir ${params.results}/results/rligerUINMF/h5ad_homology_concat  --homology_tbl homology_tbl.csv \
    --metadata_output ${basename}_liger.tsv

    """
}



workflow {

    metadata_ch = Channel.fromPath(params.input_metadata)
                                         .map { file -> tuple(file.baseName, file) }
    metadata_ch.view()
    validate_adata_input(metadata_ch)
    concat_by_homology(metadata_ch)
    concat_by_homology_rliger_uinmf(metadata_ch)


}
