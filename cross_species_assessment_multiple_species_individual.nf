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


metadata_ch = channel.fromPath(params.input_metadata)

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


all_integrated_rds_mapped_ch = channel
                .fromPath(params.integrated_rds)
                .map { file -> tuple(file.baseName, file, file.getParent()) }
// method name, file basename, file

process copy_for_rliger {

    label 'regular_resource'

    publishDir "${params.results}/results/h5ad_homology_concat", mode: 'copy'

    input:
    tuple val(basename), path(unintegrated_h5ad)

    output:
    val true, emit: signal

    shell:
    '''
    if ! [ -n "\$(find !{params.results}/results/h5ad_homology_concat -type f -regex '.*liger.*')" ]
    then 
        rliger_file=$(echo !{unintegrated_h5ad}  | sed "s/!{basename}/!{basename}_liger/g") && echo ${rliger_file} && cp !{unintegrated_h5ad} ${rliger_file}
    fi
    '''

}

process convert_format_rds {

    label 'convert'
    label 'regular_resource'

    //publishDir "${out_dir}/", mode: 'copy'

    input: 
    tuple val(basename), path(integrated_rds), path(out_dir)

    //output:
    //path "*.h5ad"
    //directly write to results out_dir, not elegant, but works

    script:

    """
    Rscript ${projectDir}/bin/convert_format.R \
    -i ${integrated_rds} -o ${out_dir}/${basename}.h5ad -t seurat_to_anndata \
    --conda_path ${params.sceasy_conda}

    """


}



process sccaf_assessment {

    label 'regular_resource'
    label 'sccaf_based'

    publishDir "${params.results}/results/per_species", mode: 'copy'

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

    label 'regular_resource'
    label 'sccaf_based'


    publishDir "${params.results}/results/${method}/cross_species/SCCAF_projection/", mode: 'copy'


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

    label 'regular_intg_resource'
    label 'scIB_based'

    publishDir "${params.results}/batch_metrics/cross_species", mode: 'copy'

    input:
    
    tuple val(ready), val(method), val(basename), path(unintegrated_h5ad), path(cross_species_integrated_h5ad)

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
    --num_cores 4 --conda_path ${params.scib_conda}

    """



}



process trajectory_metrics{

    label 'regular_resource'
    label 'scIB_based'

    publishDir "${params.results}/batch_metrics/cross_species", mode: 'copy'

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

    //all_integrated_rds_mapped_ch.view()

    convert_format_rds(all_integrated_rds_mapped_ch)

    copy_for_rliger(all_orig_h5ad_with_base_mapped_ch)
    sccaf_assessment(metadata_ch)
    sccaf_projection(all_integrated_h5ad_mapped_ch)
    ch_test = copy_for_rliger.out.signal.combine(all_integrated_and_orig_h5ad_mapped_ch).unique()
    batch_metrics(ch_test)
    //trajectory_metrics(iall_integrated_and_orig_h5ad_mapped_ch)
}
