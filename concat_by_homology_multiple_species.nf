#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// data requirements: batch_key, species_key, sample_key, mapped in config file, present in adata.obs
// data requirements: mean_count in adata.var from scanpy QC
// raw h5ad file naming: <species>
// adata.var_names are ensembl gene ids

log.info """
         ===========================================================


         Cross-species integration and assessment - nextflow pipeline
         - Integrate scRNA-seq data from multiple species - this workflow
         - SCCAF projection on integrated data to assess integration quality
         - harmony, scanorama, scVI - python based
         - Seurat CCA, Seurat RPCA, fastMNN, LIGER, LIGER-UINMF - r based
         Author: ysong@ebi.ac.uk
         Mar 2022

         ===========================================================

         """
         .stripIndent()


process concat_by_homology {

    publishDir "${params.results}/results/h5ad_homology_concat", mode: 'copy'
    cpus 12
    queue 'research'
    memory '100GB'
    conda "${projectDir}/envs/r_based_concat_integration.yml"
    cache 'lenient'

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

    publishDir "${params.results}/results/rligerUINMF/h5ad_homology_concat", mode: 'copy'
    cpus 12
    queue 'research'
    memory '100GB'
    conda "${projectDir}/envs/r_based_concat_integration.yml"
    cache 'lenient'

    input:

    tuple val(basename), path(metadata)

    output:
    path "homology_tbl.csv"
    path "*_ligerUINMF.h5ad"

    script:
    """
    mkdir -p ${params.results}/results/rligerUINMF/h5ad_homology_concat && \
    Rscript ${projectDir}/bin/concat_by_homology_rligerUINMF_multiple_species.R --metadata ${metadata} \
    --out_dir ${params.results}/results/rligerUINMF/h5ad_homology_concat  --homology_tbl homology_tbl.csv \
    --metadata_output ${baseName}.tsv

    """
}







workflow {

    metadata_ch = Channel.fromPath(params.input_metadata)
                                         .map { file -> tuple(file.baseName, file) }
    concat_by_homology(metadata_ch)
    concat_by_homology_rliger_uinmf(metadata_ch)


}
