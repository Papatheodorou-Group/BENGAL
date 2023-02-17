#!/usr/bin/env nextflow


nextflow.enable.dsl=2


// data requirements: batch_key, species_key, sample_key in adata.obs
// data requirements: mean_count in adata.var from scanpy QC
// name the raw h5ad file using the same string as in homology table
// adata.var have column 'gene_name' as gene_name, be unique!
// adata.obs and adata.var cannot contain NA and the column names cannot contain "-"
// adata.var.cluster_key factor is ordered to be the same across two datasets with the same cell type // todo

log.info """
         ===========================================================

         Cross-species integration and assessment - nextflow pipeline
         - Integrate scRNA-seq data from multiple species
         - SCCAF projection on integrated data to assess integration quality - this workflow
         - harmony, scanorama, scVI - python based
         - Seurat CCA, Seurat RPCA, fastMNN, LIGER, LIGER-UINMF - r based
         - SAMap
         Author: ysong@ebi.ac.uk
         Mar 2022

         ===========================================================

         """
         .stripIndent()


metadata = channel.fromPath(params.input_metadata)

all_integrated_h5ad_mapped_ch = channel
                .fromPath(params.integrated_h5ad)
                .map { file -> tuple(file.baseName.split("_")[-2], file.baseName, file) }

all_orig_h5ad_mapped_ch = channel
                .fromPath(params.integrated_h5ad)
                .map { file -> (params.homology_concat_h5ad - "*.h5ad" + file.baseName - file.baseName.split("_")[-2] - '__integrated'  + ".h5ad") }
                .unique()

all_orig_h5ad_with_base_mapped_ch = channel
                .fromPath(params.homology_concat_h5ad)
                .map { file -> tuple(params.task_name, file) }
                .unique()


all_integrated_and_orig_h5ad_mapped_ch = channel
                .fromPath(params.integrated_h5ad)
                .map { file -> tuple(file.baseName.split("_")[-2], file.baseName, (params.homology_concat_h5ad - "*.h5ad" + file.baseName - file.baseName.split("_")[-2] - '__integrated'  + ".h5ad"), file) }



// method name, file basename, file

process copy_for_rliger {


    publishDir "${params.results}/results/h5ad_homology_concat", mode: 'copy'
    cpus 1
    queue 'research'
    memory '10GB'
    conda "${projectDir}/envs/py_based_integration.yml"
    cache 'lenient'

    input:
    tuple val(basename), path(unintegrated_h5ad)

    output:
    path "*"


    shell:
    '''
    rliger_file=$(echo !{unintegrated_h5ad}  | sed "s/!{basename}/rliger_uinmf_metadata/g") && echo ${rliger_file} && cp !{unintegrated_h5ad} ${rliger_file}
    '''



}



process sccaf_assessment {
    publishDir "${params.results}/results/per_species", mode: 'copy'
    cpus 8
    queue 'production'
    memory '20GB'
    conda "${projectDir}/envs/py_based_integration.yml"


    input:
    path(metadata)

   output:
    path '*'

    script:
    """
    python ${projectDir}/bin/sccaf_assessment_metadata.py ${metadata} ${metadata}_SCCAF_AUC ${metadata}_SCCAF_accuracy_summary \
    --use_embedding True --embedding_key X_pca \
    --integration_method unintegrated \
    --cluster_key ${params.cluster_key}
    """

}


process sccaf_projection {

    publishDir "${params.results}/results/${method}/cross_species/SCCAF_projection/", mode: 'copy'

    cpus 12
    queue 'production'
    memory '50GB'
    conda "${projectDir}/envs/py_based_integration.yml"
    cache 'lenient'

    input:
    tuple val(method), val(basename), path(cross_species_integrated_h5ad)

    output:
    path '*'

    script:
    """
    python ${projectDir}/bin/sccaf_projection_multiple_species.py \
    --species_key ${params.species_key} --cluster_key ${params.cluster_key_sccaf} --projection_key ${params.projection_key_sccaf} \
    --integration_method ${method} ${cross_species_integrated_h5ad} \
    ${basename}_SCCAF_projection_result.h5ad \
    ${basename}_SCCAF_projection_figures.pdf \
    ${basename}_SCCAF_accuracy_summary.csv
    """
}



process batch_metrics{

    publishDir "${params.results}/batch_metrics/cross_species", mode: 'copy'
    cpus 4
    queue 'research'
    memory '100GB'
    conda "${projectDir}/envs/scib-pipeline-R4.0.yml"
    cache 'lenient'

    input:
    tuple val(method), val(basename), path(unintegrated_h5ad), path(cross_species_integrated_h5ad)

    output:
    path "*"

    
    script:
    """
    python ${projectDir}/bin/scIB_metrics_individual.py \
    ${cross_species_integrated_h5ad} \
    ${unintegrated_h5ad} \
    ${basename}_scIB_metrics.csv ${basename}_cell_type_basw.csv \
    ${basename}_orig_scIB_metrics.csv ${basename}_orig_cell_type_basw.csv \
    ${basename}_scIB.h5ad \
    --integration_method ${method} --batch_key ${params.batch_key} \
    --species_key ${params.species_key} --cluster_key ${params.cluster_key} \
    --num_cores 4

    """



}



process trajectory_metrics{

    publishDir "${params.results}/batch_metrics/cross_species", mode: 'copy'
    cpus 4
    queue 'research'
    conda "${projectDir}/envs/scib-pipeline-R4.0.yml"
    memory '100GB'
    cache 'lenient'

    input:
    tuple val(method), val(basename), path(unintegrated_h5ad), path(cross_species_integrated_h5ad)

    output:
    path "*"

    script:
    """
    python ${projectDir}/bin/scIB_trajectory.py \
    ${cross_species_integrated_h5ad} ${unintegrated_h5ad} \
    ${basename}_trajectory_metrics.csv \
    ${basename}_trajectory_scIB.h5ad \
    --integration_method ${method} --batch_key ${params.batch_key} --species_key ${params.species_key} --cluster_key ${params.cluster_key} --root_cell ${params.root_cell}

    """

}


workflow {
    integrated_h5ad = all_integrated_h5ad_mapped_ch
    integrated_and_orig_h5ad = all_integrated_and_orig_h5ad_mapped_ch
    orig_h5ad = all_orig_h5ad_mapped_ch
    metadata_ch = metadata
    all_orig_h5ad_with_base_mapped_ch.view()
    copy_for_rliger(all_orig_h5ad_with_base_mapped_ch)
    batch_metrics(integrated_and_orig_h5ad)
    trajectory_metrics(integrated_and_orig_h5ad)
}
