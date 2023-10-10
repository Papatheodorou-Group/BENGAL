# © EMBL-European Bioinformatics Institute, 2023
# Yuyao Song <ysong@ebi.ac.uk>

K_default=20

import pandas as pd
import numpy as np
import scanpy as sc
import random
import numpy
import click

random.seed(123)
numpy.random.seed(456)


@click.command()
@click.argument("input_h5ad", type=click.Path(exists=True))
@click.argument("unintegrated_h5ad", type=click.Path(exists=False), default=None)
@click.argument("out_integrated_metrics", type=click.Path(exists=False), default=None)
@click.argument("out_orig_metrics", type=click.Path(exists=False), default=None)
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
def run_alignment_score(
    input_h5ad,
    unintegrated_h5ad,
    out_integrated_metrics,
    out_orig_metrics,
    species_key,
    batch_key,
    cluster_key,
    integration_method,
):
    # dictionary for method properties
    embedding_keys = {
        "harmony": "X_pca_harmony",
        "scanorama": "X_scanorama",
        "scVI": "X_scVI",
        "scANVI": "X_scANVI",
        "LIGER": "X_iNMF",
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
        "SAMap": False,
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
        "SAMap": False,
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
    if from_h5seurat[integration_method] is True:
        input_ad.__dict__["_raw"].__dict__["_var"] = (
            input_ad.__dict__["_raw"]
            .__dict__["_var"]
            .rename(columns={"_index": "features"})
        )

    use_embedding = use_embeddings[integration_method]
    if use_embedding is True:
        embedding_key = embedding_keys[integration_method]

    # re-calculate on integrated and unintegrated data
    # due to scIB hard-coding, make sure input_ad.obsp['connectivities'], input_ad.uns['neighbours'] are from the embedding
    # for lisi type_='knn'
    # LIGER embedding only have 20 dims

    if integration_method == "SAMap":
        click.echo("use SAMap KNN graph")
        # do nothing
    elif use_embedding is True:
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
    output_orig = pd.DataFrame()

    click.echo("alignment score")
    
    def q(x):
        return np.array(list(x))

    def avg_as(ad):
        x = q(ad.obs[batch_key])
        xu = np.unique(x)
        a = np.zeros((xu.size,xu.size))
        for i in range(xu.size):
            for j in range(xu.size):
                if i!=j:
                    a[i,j] = ad.obsp['connectivities'][x==xu[i],:][:,x==xu[j]].sum(1).A.flatten().mean() / K_default
        return pd.DataFrame(data=a,index=xu,columns=xu)


    output_integrated.loc["Alignment_score", "value"] = avg_as(input_ad)[avg_as(input_ad) != 0].mean().mean()
    
    output_orig.loc["Alignment_score", "value"] = avg_as(orig_ad)[avg_as(orig_ad) != 0].mean().mean()

    output_integrated.loc["input_h5ad", "value"] = input_h5ad
    output_integrated.loc["unintegrated_h5ad", "value"] = unintegrated_h5ad
    output_integrated.loc["species_key", "value"] = species_key
    output_integrated.loc["batch_key", "value"] = batch_key
    output_integrated.loc["cluster_key", "value"] = cluster_key
    output_integrated.loc["integration_method", "value"] = integration_method

    output_integrated.T.to_csv(out_integrated_metrics)

    click.echo("metric of unintegrated data")
    output_orig.loc["input_h5ad", "value"] = unintegrated_h5ad
    output_orig.loc["unintegrated_h5ad", "value"] = unintegrated_h5ad
    output_orig.loc["species_key", "value"] = species_key
    output_orig.loc["batch_key", "value"] = batch_key
    output_orig.loc["cluster_key", "value"] = cluster_key
    output_orig.loc["integration_method", "value"] = integration_method

    output_orig.T.to_csv(out_orig_metrics)

    
    

if __name__ == "__main__":
    run_alignment_score()