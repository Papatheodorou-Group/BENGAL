#!/usr/bin/env python3

# © EMBL-European Bioinformatics Institute, 2023
# Yuyao Song <ysong@ebi.ac.uk>


import click
import pandas as pd
import scanpy as sc
import scib


@click.command()

@click.argument("input_h5ad", type=click.Path(exists=True))
@click.argument("unintegrated_h5ad", type=click.Path(exists=False), default=None)
@click.argument("out_csv", type=click.Path(exists=False), default=None)
@click.argument("out_silhouette", type=click.Path(exists=False), default=None)
@click.argument("out_h5ad", type=click.Path(exists=False), default=None)
@click.option('--species_key', type=str, default=None, help="Species key to distinguish species")
@click.option('--batch_key', type=str, default=None, help="Batch key on which integration is performed")
@click.option('--integration_method', type=str, default=None, help="Integration method")
@click.option('--cluster_key', type=str, default=None, help="Cluster key in species one to use as labels to transfer to species two")


def run_scIB_metrics(input_h5ad, unintegrated_h5ad, out_csv, out_silhouette, out_h5ad, species_key, batch_key, cluster_key, integration_method):
    # dictionary for method properties
    embedding_keys={"harmony": "X_pca_harmony", "scanorama": "X_scanorama", "scVI": "X_scVI", "LIGER": "X_iNMF", "rligerUINMF":"X_inmf", "fastMNN": "X_mnn", "SAMap": "wPCA" }
    use_embeddings={"harmony": True, "scanorama": True, "scVI": True, "LIGER": True, "rligerUINMF":True, "fastMNN": True, "SAMap": True , "seuratCCA": False, "seuratRPCA": False, "unintegrated": False}
    from_h5seurat={"harmony": False, "scanorama": False, "scVI": False, "LIGER": True, "rligerUINMF":True, "fastMNN": True, "SAMap": False , "seuratCCA": True, "seuratRPCA": True, "unintegrated": False}
    sc.set_figure_params(dpi_save=200, frameon=False, figsize=(10, 5))
    click.echo("Read anndata")
    input_ad = sc.read_h5ad(input_h5ad)
    orig_ad = sc.read_h5ad(unintegrated_h5ad)
    species_all=input_ad.obs[species_key].astype("category").cat.categories.values
    # known bug - fix when convert h5Seurat to h5ad the index name error
    if from_h5seurat[integration_method] is True:
        input_ad.__dict__['_raw'].__dict__['_var'] = input_ad.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})

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
        sc.pp.neighbors(input_ad, n_neighbors=20, n_pcs=20, use_rep=embedding_key)
        ## compute knn if use embedding
    else:
        click.echo("use PCA to compute KNN graph")
        sc.tl.pca(input_ad, svd_solver='arpack')
        sc.pp.neighbors(input_ad, n_neighbors=20, n_pcs=20, use_rep='X_pca')
        embedding_key='X_pca'
    ## while no embedding, compute PCA and compute knn

    ## get neighbour graph from unintegrated data
    sc.pp.normalize_total(orig_ad, target_sum=1e4)
    sc.pp.log1p(orig_ad)
    sc.pp.highly_variable_genes(orig_ad, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pp.scale(orig_ad, max_value=10)
    sc.tl.pca(orig_ad, svd_solver='arpack')
    sc.pp.neighbors(orig_ad, n_neighbors=20, n_pcs=40)

    click.echo("Start computing various batch metrics using scIB "+input_h5ad)
    click.echo("silhouette_batch")

    ## silhouette_batch
    silb = scib.metrics.silhouette_batch(input_ad, batch_key = batch_key, group_key = cluster_key, embed = embedding_key, metric='euclidean',
                                         return_all=True, scale=True, verbose=True)
    click.echo("PC regression")

    ## pcr

    pcr = scib.metrics.pcr_comparison(adata_pre = orig_ad, adata_post = input_ad, covariate = batch_key,
                                      embed = embedding_key, n_comps=50, scale=True, verbose=True)
    ## click.echo("HVG overlap")
    ## I don't understand how methods that generate an embedding change HVG
    ## hvg_overlap

    ## hvg = scib.metrics.hvg_overlap(adata_pre = orig_ad, adata_post = input_ad, batch = batch_key, n_hvg=500, verbose=True)

    ## sil on embedding
    click.echo("silhouette on embedding")

    sil = scib.metrics.silhouette(adata = input_ad, group_key = cluster_key, embed = embedding_key, metric='euclidean', scale=True)

    ## nmi and ari cluster vs label, graph_conn, clisi and ilisi
    click.echo("NMI, ARI, grapth connectivity, cLISI and iLISI")

    res = scib.metrics.metrics(adata =  orig_ad, adata_int = input_ad, batch_key=batch_key, label_key=cluster_key, embed=embedding_key,
                               cluster_key='cluster_scIB', cluster_nmi=out_csv+"_nmi_opt_cluster.csv", type_='knn',
                               isolated_labels_asw_=False,
                               silhouette_=False,
                               hvg_score_=False,
                               graph_conn_=True,
                               pcr_=False,
                               isolated_labels_f1_=True,
                               nmi_=True,
                               ari_=True,
                               kBET_=True,
                               ilisi_=True,
                               clisi_=True,
                               cell_cycle_=False, trajectory_=False)

    ## collect results and write output
    output=pd.DataFrame()
    output.loc['NMI_cluster/label', 'value'] = res.loc['NMI_cluster/label', 0]
    output.loc['ARI_cluster/label', 'value'] = res.loc['ARI_cluster/label', 0]
    output.loc['iLISI', 'value'] = res.loc['iLISI', 0]
    output.loc['cLISI', 'value'] = res.loc['cLISI', 0]
    output.loc['graph_conn', 'value'] = res.loc['graph_conn', 0]
    #output.loc['HVG_overlap', 'value'] = hvg
    output.loc['pcr', 'value'] = pcr
    output.loc['silhouette', 'value'] = sil
    output.loc['silhouette_batch', 'value'] = silb[0]
    output.loc['input_h5ad', 'value'] = input_h5ad
    output.loc['unintegrated_h5ad', 'value'] = unintegrated_h5ad
    output.loc['species_key', 'value'] = species_key
    output.loc['batch_key', 'value'] = batch_key
    output.loc['cluster_key', 'value'] = cluster_key
    output.loc['integration_method', 'value'] = integration_method

    output.T.to_csv( out_csv)

    silb[1]['input_h5ad'] = input_h5ad
    silb[1]['unintegrated_h5ad'] = unintegrated_h5ad
    silb[1]['species_key'] = species_key
    silb[1]['batch_key'] = batch_key
    silb[1]['cluster_key'] = cluster_key
    silb[1]['integration_method'] = integration_method
    silb[1].to_csv( out_silhouette)

    input_ad.write(out_h5ad, compression = 'gzip')
    click.echo("finish batch metrics")


if __name__ == '__main__':
    run_scIB_metrics()
