#!/usr/bin/env nextflow
nextflow.enable.dsl=2


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


process harmony_integration {

    publishDir "${params.results}/results/harmony/cross_species/integrated_h5ad", mode: 'copy'
    cpus 12
    queue 'research'
    memory '150GB'
    conda "${projectDir}/envs/py_based_integration.yml"
    cache 'lenient'

    input:
    tuple val(baseName), path(file)

    output:
    path "${baseName}_harmony_integrated.h5ad", emit: harmony_cross_species_h5ad_ch
    path "${baseName}_harmony_integrated_UMAP.png", emit: harmony_cross_species_umap_ch

    script:
    """
    python ${projectDir}/bin/harmony_integration.py \
    --batch_key ${params.batch_key} --species_key ${params.species_key} --cluster_key ${params.cluster_key} \
    ${file} ${baseName}_harmony_integrated.h5ad ${baseName}_harmony_integrated_UMAP.png
    """
}


process scanorama_integration {

    publishDir "${params.results}/results/scanorama/cross_species/integrated_h5ad", mode: 'copy'
    cpus 12
    queue 'research'
    memory '150GB'
    conda "${projectDir}/envs/py_based_integration.yml"
    cache 'lenient'

    input:
    tuple val(baseName), path(file)

    output:
    path "${baseName}_scanorama_integrated.h5ad", emit: scanorama_cross_species_h5ad_ch
    path "${baseName}_scanorama_integrated_UMAP.png", emit: scanorama_cross_species_umap_ch

    script:
    """
    python ${projectDir}/bin/scanorama_integration.py \
    --batch_key ${params.batch_key} --species_key ${params.species_key} --cluster_key ${params.cluster_key} \
    ${file} ${baseName}_scanorama_integrated.h5ad ${baseName}_scanorama_integrated_UMAP.png
    """
}


// scVI use gpu version
process scvi_integration {

    publishDir "${params.results}/results/scVI/cross_species/integrated_h5ad", mode: 'copy'
    queue 'gpu'
    clusterOptions ' -gpu "num=2:j_exclusive=no" -P gpu -n 4 '
    memory '50GB'
    conda "${projectDir}/envs/scVI_integration.yml"
    cache 'lenient'

    input:
    tuple val(baseName), path(file)

    output:
    path "${baseName}_scVI_integrated.h5ad", emit: scvi_cross_species_h5ad_ch
    path "${baseName}_scVI_integrated_UMAP.png", emit: scvi_cross_species_umap_ch

    script:
    """
    python ${projectDir}/bin/scvi_integration.py \
    --batch_key ${params.batch_key} --species_key ${params.species_key} --cluster_key ${params.cluster_key} \
    ${file} ${baseName}_scVI_integrated.h5ad ${baseName}_scVI_integrated_UMAP.png
    """
}


process seurat_CCA_integration{

    publishDir "${params.results}/results/seuratCCA/cross_species/integrated_h5ad/", mode: 'copy'
    cpus 12
    queue 'bigmem'
    memory '350GB'
    conda "${projectDir}/envs/r_based_concat_integration.yml"
    cache 'lenient'

    input:
    tuple val(baseName), path(file)

    output:
    path "${baseName}_seuratCCA_integrated.h5seurat", emit: seuratCCA_cross_species_h5seurat_ch
    path "${baseName}_seuratCCA_integrated_UMAP.pdf" , emit: seuratCCA_cross_species_umap_ch

    script:
    """
    Rscript ${projectDir}/bin/seurat_CCA_integration.R -i ${file} -o ${baseName}_seuratCCA_integrated.h5seurat -p ${baseName}_seuratCCA_integrated_UMAP.pdf -b ${params.batch_key} -s ${params.species_key} -c ${params.cluster_key}
    """
}


process seurat_RPCA_integration{

    publishDir "${params.results}/results/seuratRPCA/cross_species/integrated_h5ad", mode: 'copy'
    cpus 12
    queue 'research'
    memory '200GB'
    conda "${projectDir}/envs/r_based_concat_integration.yml"
    cache 'lenient'

    input:
    tuple val(baseName), path(file)

    output:
    path "${baseName}_seuratRPCA_integrated.h5seurat" , emit: seuratRPCA_cross_species_h5seurat_ch
    path "${baseName}_seuratRPCA_integrated_UMAP.pdf" , emit: seuratRPCA_cross_species_umap_ch

    script:
    """
    Rscript ${projectDir}/bin/seurat_RPCA_integration.R -i ${file} -o ${baseName}_seuratRPCA_integrated.h5seurat -p ${baseName}_seuratRPCA_integrated_UMAP.pdf -b ${params.batch_key} -s ${params.species_key} -c ${params.cluster_key}
    """
}

process fastMNN_integration{

    publishDir "${params.results}/results/fastMNN/cross_species/integrated_h5ad", mode: 'copy'
    cpus 8
    queue 'research'
    memory '200GB'
    conda "${projectDir}/envs/r_based_concat_integration.yml"
    cache 'lenient'

    input:
    tuple val(baseName), path(file)

    output:
    path "${baseName}_fastMNN_integrated.h5seurat" , emit: fastMNN_cross_species_h5seurat_ch
    path "${baseName}_fastMNN_integrated_UMAP.pdf" , emit: fastMNN_cross_species_UMAP_ch

    script:
    """
    Rscript ${projectDir}/bin/fastMNN_integration.R -i ${file} -o ${baseName}_fastMNN_integrated.h5seurat -p ${baseName}_fastMNN_integrated_UMAP.pdf -b ${params.batch_key} -s ${params.species_key} -c ${params.cluster_key}
    """
}


process liger_integration{

    publishDir "${params.results}/results/LIGER/cross_species/integrated_h5ad", mode: 'copy'
    cpus 8
    queue 'research'
    memory '200GB'
    conda "${projectDir}/envs/r_based_concat_integration.yml"
    cache 'lenient'

    input:
    tuple val(baseName), path(file)

    output:
    path "${baseName}_LIGER_integrated.h5seurat", emit: liger_cross_species_h5seurat_ch
    path "${baseName}_LIGER_integrated_UMAP.pdf", emit: liger_cross_species_UMAP_ch

    script:
    """
    Rscript ${projectDir}/bin/LIGER_integration.R -i ${file} -o ${baseName}_LIGER_integrated.h5seurat -p ${baseName}_LIGER_integrated_UMAP.pdf -b ${params.batch_key} -s ${params.species_key} -c ${params.cluster_key}
    """
}


process rligerUINMF_integration{

    publishDir "${params.results}/results/rligerUINMF/cross_species/integrated_h5ad", mode: 'copy'
    cpus 8
    queue 'research'
    memory '200GB'
    conda "${projectDir}/envs/r_based_concat_integration.yml"
    cache 'lenient'

    input:
    tuple val(baseName), path(metadata)

    output:
    path "rliger_uinmf_metadata_*_rligerUINMF_integrated.h5seurat"
    path "rliger_uinmf_metadata_*_rligerUINMF_integrated_UMAP.pdf"

    script:
    """
    mkdir -p ${params.results}/results/rligerUINMF/cross_species/integrated_h5ad && \
    Rscript ${projectDir}/bin/rliger_integration_UINMF_multiple_species.R \
    --basename ${baseName} \
    --metadata ${params.liger_metadata} \
    --out_dir ${params.results}/results/rligerUINMF/cross_species/integrated_h5ad  \
    --cluster_key ${params.cluster_key}
    """
}


workflow {

    all_homology_h5ad_mapped_ch = Channel.fromPath(params.homology_concat_h5ad)
                                         .map { file -> tuple(file.baseName, file) }

    all_homology_h5seurat_mapped_ch = Channel.fromPath(params.homology_concat_h5seurat)
                                         .map { file -> tuple(file.baseName, file) }

    concatenated_h5ad = all_homology_h5ad_mapped_ch
    concatenated_h5ad.view()

    harmony_integration(concatenated_h5ad)
    scanorama_integration(concatenated_h5ad)
    scvi_integration(concatenated_h5ad)

    concatenated_h5seurat = all_homology_h5seurat_mapped_ch
    seurat_CCA_integration(concatenated_h5seurat)
    seurat_RPCA_integration(concatenated_h5seurat)
    fastMNN_integration(concatenated_h5seurat)
    liger_integration(concatenated_h5seurat)

    liger_metadata = Channel.fromPath(params.liger_metadata)
                                         .map { file -> tuple(file.baseName, file) }
    liger_metadata.view()
    rligerUINMF_integration(liger_metadata)

}
