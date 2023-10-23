#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// © EMBL-European Bioinformatics Institute, 2023

// Updated upon containerization of BENGAL envs
// Yuyao Song <ysong@ebi.ac.uk>
// Oct 2023


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
    
    label 'scanpy_based'
    label 'regular_intg_resource'

    publishDir "${params.results}/results/harmony/cross_species/integrated_h5ad", mode: 'copy'
   
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

    label 'scanpy_based'
    label 'regular_intg_resource'

    publishDir "${params.results}/results/scanorama/cross_species/integrated_h5ad", mode: 'copy'
    
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


// scVI uses GPU
process scvi_integration {

    label 'scvi_based'
    label 'GPU_intg_resource'

    publishDir "${params.results}/results/scVI/cross_species/integrated_h5ad", mode: 'copy'

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
// scANVI uses GPU
process scanvi_integration {


    label 'scvi_based'
    label 'GPU_intg_resource'

    publishDir "${params.results}/results/scANVI/cross_species/integrated_h5ad", mode: 'copy'

    input:
    tuple val(baseName), path(file)

    output:
    path "${baseName}_scANVI_integrated.h5ad", emit: scanvi_cross_species_h5ad_ch
    path "${baseName}_scANVI_integrated_UMAP.png", emit: scanvi_cross_species_umap_ch

    script:
    """
    python ${projectDir}/bin/scANVI_integration.py \
    --batch_key ${params.batch_key} --species_key ${params.species_key} --cluster_key ${params.cluster_key} \
    ${file} ${baseName}_scANVI_integrated.h5ad ${baseName}_scANVI_integrated_UMAP.png ${baseName}_scANVI_integrated_mde.png
    """
}



process seurat_CCA_integration{


    label 'seurat_based'
    label 'bigmem_intg_resource'

    publishDir "${params.results}/results/seuratCCA/cross_species/integrated_h5ad/", mode: 'copy'

    input:
    tuple val(baseName), path(file)

    output:
    path "${baseName}_seuratCCA_integrated.rds", emit: seuratCCA_cross_species_rds_ch
    path "${baseName}_seuratCCA_integrated_UMAP.pdf" , emit: seuratCCA_cross_species_umap_ch

    script:
    """
    Rscript ${projectDir}/bin/seurat_CCA_integration.R -i ${file} -o ${baseName}_seuratCCA_integrated.rds -p ${baseName}_seuratCCA_integrated_UMAP.pdf -b ${params.batch_key} -s ${params.species_key} -c ${params.cluster_key}
    """
}


process seurat_RPCA_integration{

    label 'seurat_based'
    label 'bigmem_intg_resource'

    publishDir "${params.results}/results/seuratRPCA/cross_species/integrated_h5ad", mode: 'copy'
   
    input:
    tuple val(baseName), path(file)

    output:
    path "${baseName}_seuratRPCA_integrated.rds" , emit: seuratRPCA_cross_species_rds_ch
    path "${baseName}_seuratRPCA_integrated_UMAP.pdf" , emit: seuratRPCA_cross_species_umap_ch

    script:
    """
    Rscript ${projectDir}/bin/seurat_RPCA_integration.R -i ${file} -o ${baseName}_seuratRPCA_integrated.rds -p ${baseName}_seuratRPCA_integrated_UMAP.pdf -b ${params.batch_key} -s ${params.species_key} -c ${params.cluster_key}
    """
}

process fastMNN_integration{

    label 'R_based'
    label 'regular_intg_resource'

    publishDir "${params.results}/results/fastMNN/cross_species/integrated_h5ad", mode: 'copy'

    input:
    tuple val(baseName), path(file)

    output:
    path "${baseName}_fastMNN_integrated.rds" , emit: fastMNN_cross_species_rds_ch
    path "${baseName}_fastMNN_integrated_UMAP.pdf" , emit: fastMNN_cross_species_UMAP_ch

    script:
    """
    Rscript ${projectDir}/bin/fastMNN_integration.R -i ${file} -o ${baseName}_fastMNN_integrated.rds -p ${baseName}_fastMNN_integrated_UMAP.pdf -b ${params.batch_key} -s ${params.species_key} -c ${params.cluster_key}
    """
}


process liger_integration{

    label 'R_based'
    label 'regular_intg_resource'

    publishDir "${params.results}/results/LIGER/cross_species/integrated_h5ad", mode: 'copy'

    input:
    tuple val(baseName), path(file)

    output:
    path "${baseName}_LIGER_integrated.rds", emit: liger_cross_species_rds_ch
    path "${baseName}_LIGER_integrated_UMAP.pdf", emit: liger_cross_species_UMAP_ch

    script:
    """
    Rscript ${projectDir}/bin/LIGER_integration.R -i ${file} -o ${baseName}_LIGER_integrated.rds -p ${baseName}_LIGER_integrated_UMAP.pdf -b ${params.batch_key} -s ${params.species_key} -c ${params.cluster_key}
    """
}


process rligerUINMF_integration{

    label 'R_based'
    label 'regular_intg_resource'

    publishDir "${params.results}/results/rligerUINMF/cross_species/integrated_h5ad", mode: 'copy'

    input:
    tuple val(baseName), path(metadata)

    output:
    path "*_rligerUINMF_integrated.rds"
    path "*_rligerUINMF_integrated_UMAP.pdf"

    script:
    """
    mkdir -p ${params.results}/results/rligerUINMF/cross_species/integrated_h5ad && \
    Rscript ${projectDir}/bin/rliger_integration_UINMF_multiple_species.R \
    --basename ${baseName} \
    --metadata ${params.liger_metadata} \
    --out_dir .  \
    --cluster_key ${params.cluster_key}
    """
}

// nextflow will search for the output file in the cwd 
// so set the out_dir in the script to cwd, then nextflow will copy the results to the publishDir anyway

workflow {

    all_homology_h5ad_mapped_ch = Channel.fromPath(params.homology_concat_h5ad)
                                         .map { file -> tuple(file.baseName, file) }

    all_homology_rds_mapped_ch = Channel.fromPath(params.homology_concat_rds)
                                         .map { file -> tuple(file.baseName, file) }

    concatenated_h5ad = all_homology_h5ad_mapped_ch
    concatenated_h5ad.view()

    harmony_integration(concatenated_h5ad)
    scanorama_integration(concatenated_h5ad)
    scvi_integration(concatenated_h5ad)
    scanvi_integration(concatenated_h5ad)

    concatenated_rds = all_homology_rds_mapped_ch
    seurat_CCA_integration(concatenated_rds)
    seurat_RPCA_integration(concatenated_rds)
    fastMNN_integration(concatenated_rds)
    liger_integration(concatenated_rds)

    liger_metadata = Channel.fromPath(params.liger_metadata)
                                         .map { file -> tuple(file.baseName, file) }
    liger_metadata.view()
    rligerUINMF_integration(liger_metadata)

}
