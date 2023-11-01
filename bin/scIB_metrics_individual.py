#!/usr/bin/env python3

# © EMBL-European Bioinformatics Institute, 2023
# Yuyao Song <ysong@ebi.ac.uk>


import os

import click
import pandas as pd
import scanpy as sc
import random
import scib
import numpy

# set R for kBET
import os

## set seed 
random.seed(123)
numpy.random.seed(456)

@click.command()
@click.argument("input_h5ad", type=click.Path(exists=True))
@click.argument("unintegrated_h5ad", type=click.Path(exists=False), default=None)
@click.argument("out_integrated_metrics", type=click.Path(exists=False), default=None)
@click.argument("out_integrated_basw", type=click.Path(exists=False), default=None)
@click.argument("out_orig_metrics", type=click.Path(exists=False), default=None)
@click.argument("out_orig_basw", type=click.Path(exists=False), default=None)
@click.argument("out_h5ad", type=click.Path(exists=False), default=None)
@click.option(
    "--species_key", type=str, default=None, help="Species key to distinguish species"
)
@click.option(
    "--batch_key",
    type=str,
    default=None,
    help="Batch key on which integration is performed",
)
@click.option("--integration_method", type=str, default=None, help="Integration method")
@click.option(
    "--cluster_key",
    type=str,
    default=None,
    help="Cluster key in species one to use as labels to transfer to species two",
)
@click.option(
    "--num_cores",
    type=int,
    default=1,
    help="Number of cores for graph LISI scores",
)
@click.option(
    "--conda_path",
    type=str,
    default=None,
    help="scIB conda path",
)


def run_scIB_metrics(
    input_h5ad,
    unintegrated_h5ad,
    out_integrated_metrics,
    out_integrated_basw,
    out_orig_metrics,
    out_orig_basw,
    out_h5ad,
    species_key,
    batch_key,
    cluster_key,
    integration_method,
    num_cores,
    conda_path
):
    

    os.environ['R_HOME'] = f"{conda_path}/lib/R"
    os.environ['PATH'] = f"{conda_path}/bin:" + os.environ['PATH']
    os.environ['R_LIBS_USER'] = f"{conda_path}/lib/R/library"
    os.environ['R_LIBS'] = f"{conda_path}lib/R/library"
    os.environ['LD_LIBRARY_PATH'] = f"{conda_path}/lib"

    # dictionary for method properties
    embedding_keys = {
        "harmony": "X_pca_harmony",
        "scanorama": "X_scanorama",
        "scVI": "X_scVI",
        "scANVI": "X_scANVI",
        "LIGER": "X_inmf",
        "rligerUINMF": "X_inmf",
        "fastMNN": "X_mnn",
    }
    use_embeddings = {
        "harmony": True,
        "scanorama": True,
        "scVI": True,
        "scANVI": True,
        "LIGER": True,
        "rligerUINMF": True,
        "fastMNN": True,
        "seuratCCA": False,
        "seuratRPCA": False,
        "unintegrated": False,
    }
    from_h5seurat = {
        "harmony": False,
        "scanorama": False,
        "scVI": False,
        "scANVI": False,
        "LIGER": True,
        "rligerUINMF": True,
        "fastMNN": True,
        "seuratCCA": True,
        "seuratRPCA": True,
        "unintegrated": False,
    }
    sc.set_figure_params(dpi_save=200, frameon=False, figsize=(10, 5))

    click.echo("Read anndata")
    input_ad = sc.read_h5ad(input_h5ad)
    orig_ad = sc.read_h5ad(unintegrated_h5ad)
    species_all = input_ad.obs[species_key].astype("category").cat.categories.values

    ## for files from h5seurat sometimes these are not stored as category

    input_ad.obs[species_key] = input_ad.obs[species_key].astype("category")
    input_ad.obs[cluster_key] = input_ad.obs[cluster_key].astype("category")  
    input_ad.obs[batch_key] = input_ad.obs[batch_key].astype("category")
    
    # known bug - fix when convert h5Seurat to h5ad the index name error
    #if from_h5seurat[integration_method] is True:
    #    input_ad.__dict__["_raw"].__dict__["_var"] = (
    #        input_ad.__dict__["_raw"]
    #        .__dict__["_var"]
    #        .rename(columns={"_index": "features"})
    #    )

    use_embedding = use_embeddings[integration_method]
    if use_embedding is True:
        embedding_key = embedding_keys[integration_method]

    # re-calculate on integrated and unintegrated data
    # due to scIB hard-coding, make sure input_ad.obsp['connectivities'], input_ad.uns['neighbours'] are from the embedding
    # for lisi type_='knn'
    # LIGER embedding only have 20 dims


    if use_embedding is True:
        click.echo("calculate KNN graph from embedding " + embedding_key)
        num_pcs = min(input_ad.obsm[embedding_key].shape[1], 20)
        if num_pcs < 20:
            click.echo("using less PCs: " + str(num_pcs))
        sc.pp.neighbors(input_ad, n_neighbors=20, n_pcs=num_pcs, use_rep=embedding_key)
        # compute knn if use embedding
    else:
        click.echo("use PCA to compute KNN graph")
        sc.tl.pca(input_ad, svd_solver="arpack")
        sc.pp.neighbors(input_ad, n_neighbors=20, n_pcs=20, use_rep="X_pca")
        embedding_key = "X_pca"

    # while no embedding, compute PCA and compute knn

    # get neighbour graph from unintegrated data
    sc.pp.normalize_total(orig_ad, target_sum=1e4)
    sc.pp.log1p(orig_ad)
    sc.pp.highly_variable_genes(orig_ad, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pp.scale(orig_ad, max_value=10)
    sc.tl.pca(orig_ad, svd_solver="arpack")
    sc.pp.neighbors(orig_ad, n_neighbors=20, n_pcs=40)
    sc.tl.umap(orig_ad, min_dist=0.3)
    sc.pl.umap(orig_ad, color=[batch_key, species_key, cluster_key])

    click.echo(
        "Start computing various batch metrics using scIB, the integrated file is "
        + input_h5ad
    )

    output_integrated = pd.DataFrame()
    output_integrated_basw = pd.DataFrame()
    output_orig = pd.DataFrame()
    output_orig_basw = pd.DataFrame()

    click.echo("PC regression")

    output_integrated.loc["PCR", "value"] = scib.metrics.pc_regression(
        input_ad.obsm[embedding_key], covariate=input_ad.obs[species_key], n_comps=50
    )
    output_orig.loc["PCR", "value"] = scib.metrics.pc_regression(
        orig_ad.obsm["X_pca"],
        covariate=orig_ad.obs[species_key],
        pca_var=orig_ad.uns["pca"]["variance"],
    )

    click.echo("Silhouette batch")

    integrated_basw = scib.metrics.silhouette_batch(
        input_ad,
        batch_key=batch_key,
        group_key=cluster_key,
        embed=embedding_key,
        metric="euclidean",
        return_all=True,
        scale=True,
        verbose=False,
    )

    output_integrated.loc["bASW", "value"] = integrated_basw[0]

    orig_basw = scib.metrics.silhouette_batch(
        orig_ad,
        batch_key=batch_key,
        group_key=cluster_key,
        embed="X_pca",
        metric="euclidean",
        return_all=True,
        scale=True,
        verbose=False,
    )

    output_orig.loc["bASW", "value"] = orig_basw[0]

    click.echo("Graph connectivity")

    output_integrated.loc["GC", "value"] = scib.metrics.graph_connectivity(
        input_ad, label_key=cluster_key
    )
    output_orig.loc["GC", "value"] = scib.metrics.graph_connectivity(
        orig_ad, label_key=cluster_key
    )

    # click.echo("graph iLISI")

    # output_integrated.loc["iLISI", "value"] = scib.metrics.ilisi_graph(
    #     input_ad,
    #     batch_key,
    #     type_="embed",
    #     use_rep=embedding_key,
    #     k0=50,
    #     subsample=None,
    #     n_cores=num_cores,
    #     scale=True,
    #     verbose=True,
    # )

    # output_orig.loc["iLISI", "value"] = scib.metrics.ilisi_graph(
    #     orig_ad,
    #     batch_key,
    #     type_="full",
    #     k0=50,
    #     subsample=None,
    #     n_cores=num_cores,
    #     scale=True,
    #     verbose=True,
    # )

    click.echo("kBET")

    output_integrated.loc["kBET", "value"] = scib.metrics.kBET(
        input_ad,
        batch_key=batch_key,
        label_key=cluster_key,
        type_="knn", ## for equal treatment
        scaled=True,
        return_df=False,
        verbose=True,
    )

    output_orig.loc["kBET", "value"] = scib.metrics.kBET(
        orig_ad,
        batch_key=batch_key,
        label_key=cluster_key,
        type_="full",
        scaled=True,
        return_df=False,
        verbose=True,
    )

    # click.echo("cLISI")

    # output_integrated.loc["cLISI", "value"] = scib.metrics.clisi_graph(
    #     input_ad,
    #     label_key=cluster_key,
    #     type_="embed",
    #     use_rep=embedding_key,
    #     k0=50,
    #     subsample=None,
    #     scale=True,
    #     n_cores=num_cores,
    #     verbose=True,
    # )

    # output_orig.loc["cLISI", "value"] = scib.metrics.clisi_graph(
    #     orig_ad,
    #     label_key=cluster_key,
    #     type_="full",
    #     k0=50,
    #     subsample=None,
    #     scale=True,
    #     n_cores=num_cores,
    #     verbose=True,
    # )

    click.echo("NMI")
    click.echo("clustering optimization with leiden")

    scib.me.cluster_optimal_resolution(
        input_ad,
        cluster_key="cluster",
        label_key=cluster_key,
        cluster_function=sc.tl.leiden,
    )
    scib.me.cluster_optimal_resolution(
        orig_ad,
        cluster_key="cluster",
        label_key=cluster_key,
        cluster_function=sc.tl.leiden,
    )

    output_integrated.loc["NMI", "value"] = scib.me.nmi(
        input_ad, cluster_key="cluster", label_key=cluster_key
    )
    output_integrated.loc["ARI", "value"] = scib.me.ari(
        input_ad, cluster_key="cluster", label_key=cluster_key
    )

    output_orig.loc["NMI", "value"] = scib.me.nmi(
        orig_ad, cluster_key="cluster", label_key=cluster_key
    )
    output_orig.loc["ARI", "value"] = scib.me.ari(
        orig_ad, cluster_key="cluster", label_key=cluster_key
    )

    click.echo("Silhouette cell type")

    output_integrated.loc["cASW", "value"] = scib.me.silhouette(
        input_ad, label_key=cluster_key, embed=embedding_key
    )

    output_orig.loc["cASW", "value"] = scib.me.silhouette(
        orig_ad, label_key=cluster_key, embed="X_pca"
    )
    
    click.echo("Isolated label F1")

    output_integrated.loc["iso_F1", "value"] = scib.me.isolated_labels_f1(
        input_ad, embed=None, batch_key=batch_key, label_key=cluster_key
    )

    output_orig.loc["iso_F1", "value"] = scib.me.isolated_labels_f1(
        orig_ad, embed=None, batch_key=batch_key, label_key=cluster_key
    )

    click.echo("write results")
    click.echo("metrics of integrated data")

    output_integrated.loc["input_h5ad", "value"] = input_h5ad
    output_integrated.loc["unintegrated_h5ad", "value"] = unintegrated_h5ad
    output_integrated.loc["species_key", "value"] = species_key
    output_integrated.loc["batch_key", "value"] = batch_key
    output_integrated.loc["cluster_key", "value"] = cluster_key
    output_integrated.loc["integration_method", "value"] = integration_method

    output_integrated.T.to_csv(out_integrated_metrics)

    click.echo("writing clustering optimized integrated data")
    input_ad.write(out_h5ad, compression="gzip")

    click.echo("metric of unintegrated data")
    output_orig.loc["input_h5ad", "value"] = unintegrated_h5ad
    output_orig.loc["unintegrated_h5ad", "value"] = unintegrated_h5ad
    output_orig.loc["species_key", "value"] = species_key
    output_orig.loc["batch_key", "value"] = batch_key
    output_orig.loc["cluster_key", "value"] = cluster_key
    output_orig.loc["integration_method", "value"] = integration_method

    output_orig.T.to_csv(out_orig_metrics)

    click.echo("cell type bASW of integrated data")

    integrated_basw[1]["input_h5ad"] = input_h5ad
    integrated_basw[1]["unintegrated_h5ad"] = unintegrated_h5ad
    integrated_basw[1]["species_key"] = species_key
    integrated_basw[1]["batch_key"] = batch_key
    integrated_basw[1]["cluster_key"] = cluster_key
    integrated_basw[1]["integration_method"] = integration_method
    integrated_basw[1].to_csv(out_integrated_basw)

    orig_basw[1]["input_h5ad"] = unintegrated_h5ad
    orig_basw[1]["unintegrated_h5ad"] = unintegrated_h5ad
    orig_basw[1]["species_key"] = species_key
    orig_basw[1]["batch_key"] = batch_key
    orig_basw[1]["cluster_key"] = cluster_key
    orig_basw[1]["integration_method"] = integration_method
    orig_basw[1].to_csv(out_orig_basw)

    click.echo("finish scIB metrics calculation")


if __name__ == "__main__":
    run_scIB_metrics()
