## BENGAL: BENchmarking strateGies for cross-species integrAtion of singLe-cell RNA sequencing data ##

Author&maintainer: Yuyao Song <ysong@ebi.ac.uk>

A Nextflow DSL2 pipeline to perform cross-species single-cell RNA-seq data integration and assessment of integration results.

**On Oct 2023, Yuyao updated the pipeline for containerization, improvements in anndata/seurat conversion, updates in scrips and updates in Nextflow.**

## This repo includes:

1) the Nextflow pipeline for cross-species integration of scRNA-seq data using various strategies
2) containers paths, dockerfiles used and conda environments
3) an example config file for running the Nextflow pipeline

## System requirements
#### Hardware:
This workflow is written to be executed on HPC clusters with LSF job scheduler. It could be easily adapted to other schedulers by changing job resource syntax in the nextflow config file. If the GPU inplementation of scVI/scANVI is to be used (beneficial for speeding up the integration on large datasets), GPU computing nodes are required, please refer to [scVI-tools site](https://scvi-tools.org/) for respective setups.

#### OS:
Development of this workflow was done on Rocky Linux 8.5 (RHEL), while in theory this can be run on any linux distribution. To run the GPU inplementation of scVI/scANVI we used Nvidia Tesla V100 GPUs. 

## Installation:

#### Clone the source code of BENGAL:
`git clone -b main git@github.com:Functional-Genomics/BENGAL.git`

**If nextflow or singularity is not installed in your cluster, install them. This can take some efforts and it might worth discussing with cluster IT managers. Please refer to [nextflow documentation](https://www.nextflow.io/docs/latest/getstarted.html) or [singularity documentation](https://singularity-tutorial.github.io/01-installation/).** 


## Inputs
The nextflow script takes one input file: a tab-seperated metadata file mapping species to the paths of raw count AnnData objects, in the form of .h5ad files. See example: `example_metadata_nf.tsv`. 

The config file defines project directories and parameters. See example: `config/example.config`

**Please change the metadata file, and the directories and parameters in the config file for your own application as appropriate**
**Read the instrictions in the config file and change relevant entries is crucial for running the pipeline**

#### Input Requirements:
*These requirements will be checked in the first process of the pipeline.*

The raw count AnnData objects need to have the following row or column annotations. Note that the exact column name of each key is specified in the config file.

1) a `species_key` in adata.obs to store species identity. Naming should be in line with the short name in ENSEMBL, such as hsapiens; mmusculus; drerio etc.
2) a `cluster_key` in adata.obs to store cell types. If assessment is performed, this column will be used to match homologous cell types across species. Preferably, use [Cell Ontology annotation](https://obofoundry.org/ontology/cl.html). 
3) `mean_counts` in adata.var computed by `sc.pp.calculate_qc_metrics` from [scanpy](https://github.com/scverse/scanpy).

The .var_names of the raw count AnnData file should be ENSEMBL gene ids.
The .X of the raw count AnnData file should be stored in dense matrix format, if SeuratDisk is used for .h5ad/.h5seurat conversion.


## Run instructions

#### Perpare the conda environment for anndata/seurat conversion. 
In principle, you can use any program to perform the conversion. Since Oct 2023 we now use [sceasy](https://github.com/cellgeni/sceasy). We also no longer use h5seurat format due to challenges in converting to/from anndata. 

It didn't seem so necessary to containerize this process so we provide a light conda environment that is compatible with other parts of the pipeline. [Mamba](https://github.com/mamba-org/mamba) is recommended as a faster substitute for conda. 

I personally perfer creating a conda env independent of nextflow and then point nextflow to the absolute path of the conda env. This way saves the running time of the pipeline and make reuse of the same env and debug easier.

First create a conda environment for the conversion:

`conda env create -f envs/sceasy.yml`

Then put the path of your sceasy conda environment into the config file in the indicated place.

#### Perpare the conda environment for running scVI/scANVI and scIB

These two parts are also not containerized since the conda env is relatively easy to set up while the respective container will be very heavy. In the future we might consider containerizing it if necessary.

`conda env create -f envs/scvi.yml`

`conda env create -f envs/scib.yml`

Then put the path of your scvi and scib conda environments into the config file in the indicated place. These env files are just created as I followed the installation instruction from [scvi](https://docs.scvi-tools.org/en/stable/installation.html) and [scib](https://scib.readthedocs.io/en/stable/installation.html) under Python 3.10.10, so if you encounter any issues, feel free to create your own evns based on their instructions.

#### Pull the containers used in BENGAL. 

We now provide a few containers to help execute the pipeline (well deserved yay due to the complexity of building them...). Please pull these containers into a local dir and specify in the config file. Here we assume you use [singularity](https://sylabs.io/) to run these containers on a HPC cluster.

1. Concatenate anndata files cross-species: `singularity pull bengal_concat.sif docker://yysong123/bengal_concat:4.2.0`
2. Python based integration: `singularity pull bengal_py.sif docker://yysong123/bengal_py:1.9.2`
3. Seurat/R based integration: `singularity pull bengal_seurat.sif docker://yysong123/bengal_seurat:4.3.0`
4. SCCAF assessment for ALCS: `singularity pull bengal_sccaf.sif docker://yysong123/bengal_sccaf:0.0.11` 

### To run BENGAL:
In a bash shell, check your metadata/config files are set and run:

1) `conda activate nextflow && nextflow -C config/example.config run concat_by_homology_multiple_species.nf`. Add flag `-with-trace -with-report report.html` if you want nextflow run stats.
2) `nextflow -C config/example.config run cross_species_integration_multiple_species.nf`
3) `nextflow -C config/example.config run cross_species_assessment_multiple_species_individual.nf`

Note: add resume flag `-resume` as appropriate to avoid re-calculation of the same data during multiple runs.

## Outputs

1) Concatenated raw count AnnData objects containing cells from all species, in the form of .h5ad files. Objects are concatenated by matching genes between species using gene homology annotation from ENSEMBL.  
2) Integration result from different algorithms including: [fastMNN](https://bioconductor.org/packages/release/bioc/html/batchelor.html), [harmony](https://github.com/slowkow/harmonypy), [LIGER](https://github.com/welch-lab/liger), [LIGER-UINMF](https://github.com/welch-lab/liger), [scanorama](https://github.com/brianhie/scanorama), [scVI](https://scvi-tools.org/), [SeuratV4CCA](https://satijalab.org/seurat/) and [SeuratV4RPCA](https://satijalab.org/seurat/), in the form of AnnData (.h5ad) or Seurat (.h5seurat) objects.
3) Respective UMAP visualizations with species; batches or cell types color code.
4) Assessment metrics for each integrated results. There are 4 batch correction metrics and 6 biology conservation metrics. Plots associated with the metrics are also generated for visual inspection. 
5) Cross-species cell type annotation transfer results with [SCCAF](https://github.com/SCCAF/sccaf).

Estimated execution time: ~6h for integrated dataset with 100,000 cells using resources specified in the .nf scripts.

## Citation

The publication in which we described and applied BENGAL is [here](https://www.nature.com/articles/s41467-023-41855-w). Please cite it if you use BENGAL.

Song, Y., Miao, Z., Brazma, A. et al. Benchmarking strategies for cross-species integration of single-cell RNA sequencing data. Nat Commun 14, 6495 (2023). https://doi.org/10.1038/s41467-023-41855-w

The BENGAL pipeline used upon publication of the paper is archived in zenodo:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8268784.svg)](https://doi.org/10.5281/zenodo.8268784)

LICENSE: GPLv3 license

NOTE: we moved this git repo from Functional-Genomics/BENGAL to Papatheodorou-Group/BENGAL on 23 Oct 2023, but redirection should happen automatically.
