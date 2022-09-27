#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """
         ===========================================================

         Cross-species integration and assessment - nextflow pipeline
         - Integrate scRNA-seq data from multiple species
         - Batch metrics and biology conservation metrics to assess integration quality - this workflow
         - harmony, scanorama, scVI - python based
         - Seurat CCA, Seurat RPCA, fastMNN, LIGER, LIGER-UINMF - r based
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

all_integrated_and_orig_h5ad_mapped_ch = channel
                .fromPath(params.integrated_h5ad)
                .map { file -> tuple(file.baseName.split("_")[-2], file.baseName, (params.homology_concat_h5ad - "*.h5ad" + file.baseName - file.baseName.split("_")[-2] - '__integrated'  + ".h5ad"), file) }



// method name, file basename, file




process sccaf_assessment {
    publishDir "${params.results}/results/per_species", mode: 'copy'
    cpus 8
    queue 'production'
    memory '20GB'
    // run in singularity container for stability

    input:
    path(metadata)

    output:
    path '*_SCCAF_accuracy_summary_*.csv'
    path '*_SCCAF_AUC_*.png'

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
    cache 'lenient'
    // run in singularity container for stability

    input:
    tuple val(method), val(basename), path(cross_species_integrated_h5ad)

    output:
    path "${basename}_SCCAF_projection_result.h5ad"
    path "${basename}_SCCAF_projection_figures.pdf"
    path "${basename}_SCCAF_accuracy_summary.csv"

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


process kBET{

    publishDir "${params.results}/batch_metrics/cross_species", mode: 'copy'
    cpus 8
    queue 'bigmem'
    memory '300GB'
    conda "${projectDir}/envs/r_based_concat_integration.yml"
    cache 'lenient'

    input:
    tuple val(method), val(basename), path(cross_species_integrated_h5ad)

    output:
    path "${basename}_kBET.csv"

    script:
    """
    Rscript ${projectDir}/bin/kbet_multiple_species.R \
    --input_h5ad ${cross_species_integrated_h5ad} --out_csv ${basename}_kBET.csv \
    --method ${method} --batch_key ${params.batch_key} --species_key ${params.species_key} --cluster_key ${params.cluster_key}

    """

}

process batch_metrics{

    publishDir "${params.results}/batch_metrics/cross_species", mode: 'copy'
    cpus 8
    queue 'bigmem'
    memory '300GB'
    conda "${projectDir}/envs/scib-pipeline-R4.0.yml"
    cache 'lenient'

    input:
    tuple val(method), val(basename), path(unintegrated_h5ad), path(cross_species_integrated_h5ad)

    output:
    path "${basename}_batch_metrics.csv"
    path "${basename}_batch_metrics.csv_nmi_opt_cluster.csv"
    path "${basename}_scIB.h5ad"
    path "${basename}_cell_type_silhouette.csv"


    script:
    """
    python ${projectDir}/bin/scIB_metrics.py \
    ${cross_species_integrated_h5ad} ${unintegrated_h5ad} \
    ${basename}_batch_metrics.csv ${basename}_cell_type_silhouette.csv \
    ${basename}_scIB.h5ad \
    --integration_method ${method} --batch_key ${params.batch_key} --species_key ${params.species_key} --cluster_key ${params.cluster_key}

    """

}


process trajectory_metrics{

    publishDir "${params.results}/batch_metrics/cross_species", mode: 'copy'
    cpus 4
    queue 'research'
    memory '100GB'
    conda "${projectDir}/envs/scib-pipeline-R4.0.yml"
    cache 'lenient'

    input:
    tuple val(method), val(basename), path(unintegrated_h5ad), path(cross_species_integrated_h5ad)

    output:
    path "${basename}_trajectory_metrics.csv"
    path "${basename}_trajectory_scIB.h5ad"

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
    integrated_and_orig_h5ad.view()

    return
    sccaf_assessment(metadata_ch)
    sccaf_projection(integrated_h5ad)
    kBET(integrated_h5ad)
    batch_metrics(integrated_and_orig_h5ad)
    trajectory_metrics(integrated_and_orig_h5ad)
}
