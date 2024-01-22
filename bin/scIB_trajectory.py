#!/usr/bin/env python3

# © EMBL-European Bioinformatics Institute, 2023
# Yuyao Song <ysong@ebi.ac.uk>


import click
import pandas as pd
import scanpy as sc
import numpy as np
import scib


@click.command()

@click.argument("input_h5ad", type=click.Path(exists=True))
@click.argument("unintegrated_h5ad", type=click.Path(exists=False), default=None)
@click.argument("out_csv", type=click.Path(exists=False), default=None)
@click.argument("out_h5ad", type=click.Path(exists=False), default=None)
@click.option('--species_key', type=str, default=None, help="Species key to distinguish species")
@click.option('--batch_key', type=str, default=None, help="Batch key on which integration is performed")
@click.option('--integration_method', type=str, default=None, help="Integration method")
@click.option('--cluster_key', type=str, default=None, help="Cluster key in species one to use as labels to transfer to species two")
@click.option('--root_cell', type=str, default=None, help="Root cell in trajectory, should be one category in adata.obs[cluster_key]")


def run_scIB_trajectory(input_h5ad, unintegrated_h5ad, out_csv,  out_h5ad, species_key, batch_key, cluster_key, integration_method, root_cell):
    # dictionary for method properties
    embedding_keys={"scANVI": "X_scANVI", "harmony": "X_pca_harmony", "scanorama": "X_scanorama", "scVI": "X_scVI", "LIGER": "X_iNMF", "rligerUINMF":"X_inmf", "fastMNN": "X_mnn", "SAMap": "wPCA" }
    use_embeddings={"scANVI": True, "harmony": True, "scanorama": True, "scVI": True, "LIGER": True, "rligerUINMF":True, "fastMNN": True, "SAMap": False , "seuratCCA": False, "seuratRPCA": False, "unintegrated": False}
    from_h5seurat={"scANVI": False, "harmony": False, "scanorama": False, "scVI": False, "LIGER": True, "rligerUINMF":True, "fastMNN": True, "SAMap": False , "seuratCCA": True, "seuratRPCA": True, "unintegrated": False}
    sc.set_figure_params(dpi_save=200, frameon=False, figsize=(10, 5))
    click.echo("Read anndata")
    input_ad = sc.read_h5ad(input_h5ad)
    orig_ad = sc.read_h5ad(unintegrated_h5ad)

    species_all=input_ad.obs[species_key].astype("category").cat.categories.values
    # known bug - fix when convert h5Seurat to h5ad the index name error
    # if from_h5seurat[integration_method] is True:
    #    input_ad.__dict__['_raw'].__dict__['_var'] = input_ad.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})

    use_embedding=use_embeddings[integration_method]
    if use_embedding is True:
        embedding_key=embedding_keys[integration_method]

    # register color in .uns
    sc.pl.umap(input_ad, color = cluster_key)
    color = dict(zip(input_ad.obs[cluster_key].cat.categories.to_list(), input_ad.uns[cluster_key+'_colors']))

    # get neighbours on integrated and unintegrated data
    # for lisi type_='knn'
    ## LIGER embedding only have 20 dims

    if integration_method == 'SAMap':
        click.echo("use SAMap KNN graph")

        ## do nothing
    elif use_embedding is True:
        click.echo("calculate KNN graph from embedding " + embedding_key)
        num_pcs = min(input_ad.obsm[embedding_key].shape[1], 20)
        if num_pcs < 20:
            click.echo("using less PCs: " + str(num_pcs))
        sc.pp.neighbors(input_ad, n_neighbors=20, n_pcs=num_pcs, use_rep=embedding_key, key_added=integration_method)        
        ## compute knn if use embedding
    else:
        click.echo("use PCA to compute KNN graph")
        sc.tl.pca(input_ad, svd_solver='arpack')
        sc.pp.neighbors(input_ad, n_neighbors=20, n_pcs=20, use_rep='X_pca', key_added=integration_method)
        embedding_key='X_pca'
        ## while no embedding, compute PCA and compute knn

    ## get neighbour graph from unintegrated data
    ## use PCA to get initial neighbours for diffusion map
    ## use diffusion map to compute neighbours - denoise
    if from_h5seurat[integration_method] is True:
        orig_ad.obs_names = orig_ad.obs_names.str.replace("-", "_")
        ## seurat does not convert - to _ which does not match integrated data

    sc.pp.normalize_total(orig_ad, target_sum=1e4)
    sc.pp.log1p(orig_ad)
    sc.pp.highly_variable_genes(orig_ad, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pp.scale(orig_ad, max_value=10)
    sc.tl.pca(orig_ad, svd_solver='arpack')
    sc.pp.neighbors(orig_ad, n_neighbors=20, n_pcs=40, use_rep='X_pca')
    sc.tl.diffmap(orig_ad)
    sc.pp.neighbors(orig_ad, n_neighbors=10, use_rep='X_diffmap')

    ## set root cell for diffusion map pseudo-time
    orig_ad.uns['iroot'] = np.flatnonzero(orig_ad.obs[cluster_key]  == root_cell)[0]
    sc.tl.dpt(orig_ad)

    ## integrated data
    ## use integrated embedding to get initial neighbours for diffusion map
    ## use diffusion map neighbours
    if integration_method == 'SAMap':
        click.echo("use SAMap KNN graph for diffmap")
        sc.tl.diffmap(input_ad, neighbors_key=None)

    else:
        click.echo("calculate diffmap for " + integration_method)
        sc.tl.diffmap(input_ad, neighbors_key=integration_method)
    ## neighbours with diffmap need not have a key
    sc.pp.neighbors(input_ad, n_neighbors=10, use_rep='X_diffmap')
    input_ad.uns['iroot'] = np.flatnonzero(input_ad.obs[cluster_key]  == root_cell)[0]
    sc.tl.dpt(input_ad)

    ## scIB trajectory conservation score
    ## per-species conservation: set batch_key='species'
    score_batch = scib.metrics.trajectory_conservation(adata_pre=orig_ad, adata_post=input_ad, label_key=cluster_key, pseudotime_key='dpt_pseudotime', batch_key=species_key)
    score = scib.metrics.trajectory_conservation(adata_pre=orig_ad, adata_post=input_ad, label_key=cluster_key, pseudotime_key='dpt_pseudotime', batch_key=None)

    output=pd.DataFrame()

    output.loc['trajectory_conservation_score_batch', 'value'] = score_batch
    output.loc['trajectory_conservation_score_none', 'value'] = score
    output.loc['input_h5ad', 'value'] = input_h5ad
    output.loc['unintegrated_h5ad', 'value'] = unintegrated_h5ad
    output.loc['species_key', 'value'] = species_key
    output.loc['batch_key', 'value'] = batch_key
    output.loc['cluster_key', 'value'] = cluster_key
    output.loc['root_cell', 'value'] = root_cell
    output.loc['integration_method', 'value'] = integration_method
    output.T.to_csv(out_csv)
    input_ad.write(out_h5ad, compression = 'gzip')
    click.echo("finish trajectory conservation metrics")


if __name__ == '__main__':
    run_scIB_trajectory()
