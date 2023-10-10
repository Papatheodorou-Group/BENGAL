#/usr/bin/env python3

# © EMBL-European Bioinformatics Institute, 2023
# Yuyao Song <ysong@ebi.ac.uk>

import click
import matplotlib.pyplot as plt
import scanpy as sc


@click.command()
@click.argument("input_h5ad", type=click.Path(exists=True))
@click.argument("out_h5ad", type=click.Path(exists=False), default=None)
@click.argument("out_umap", type=click.Path(exists=False), default=None)
@click.option('--batch_key', type=str, default=None, help="Batch key in identifying HVG and harmony integration")
@click.option('--species_key', type=str, default=None, help="Species key to distinguish species")
@click.option('--cluster_key', type=str, default=None, help="Cluster key in species one to use as labels to transfer to species two")


def run_scanorama(input_h5ad, out_h5ad, out_umap, batch_key, species_key, cluster_key):
    click.echo('Start scanorama integration')
    sc.set_figure_params(dpi_save=300, frameon=False, figsize=(10, 6))
    adata = sc.read_h5ad(input_h5ad)
    adata.var_names_make_unique()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    click.echo("HVG")
    sc.pp.highly_variable_genes(adata, batch_key=batch_key)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata, use_rep='X_pca', n_neighbors=15, n_pcs=40)
    sc.tl.umap(adata, min_dist=0.3) ## to match min_dist in seurat
    adata.obsm['X_umapraw'] = adata.obsm['X_umap']
    click.echo("Scanorama")
    sc.external.pp.scanorama_integrate(adata, key=species_key, basis='X_pca', adjusted_basis='X_scanorama')
    sc.pp.neighbors(adata, use_rep='X_scanorama', key_added = 'scanorama', n_neighbors=15, n_pcs=40)
    sc.tl.umap(adata, neighbors_key = 'scanorama', min_dist=0.3) ## to match min_dist in seurat
    sc.pl.umap(adata, color=[batch_key, species_key, cluster_key], ncols=1)
    plt.savefig(out_umap, dpi=300,  bbox_inches='tight')
    adata.obsm['X_umapscanorama'] = adata.obsm['X_umap']
    click.echo("Save output")
    adata.write(out_h5ad)
    click.echo("Done scanorama")


if __name__ == '__main__':
    run_scanorama()

