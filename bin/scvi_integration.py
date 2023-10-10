#/usr/bin/env python3

# © EMBL-European Bioinformatics Institute, 2023
# Yuyao Song <ysong@ebi.ac.uk>

import click
import matplotlib.pyplot as plt
import scanpy as sc
import scvi



@click.command()
@click.argument("input_h5ad", type=click.Path(exists=True))
@click.argument("out_h5ad", type=click.Path(exists=False), default=None)
@click.argument("out_umap", type=click.Path(exists=False), default=None)
@click.option('--batch_key', type=str, default=None, help="Batch key in identifying HVG and harmony integration")
@click.option('--species_key', type=str, default=None, help="Species key to distinguish species")
@click.option('--cluster_key', type=str, default=None, help="Cluster key in species one to use as labels to transfer to species two")


def run_scVI(input_h5ad, out_h5ad, out_umap, batch_key, species_key, cluster_key):
    click.echo('Start scVI integration - use cpu mode')
    sc.set_figure_params(dpi_save=300, frameon=False, figsize=(10, 6))
    adata = sc.read_h5ad(input_h5ad)
    adata.var_names_make_unique()
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=2000,
        ##layer="counts",
        batch_key=batch_key,
        subset=True
    )

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    #sc.pp.scale(adata, max_value=10)
    #sc.tl.pca(adata)
    #sc.pp.neighbors(adata)
    #sc.tl.umap(adata)
    #adata.obsm['X_umapraw'] = adata.obsm['X_umap']

    click.echo("setup scVI model")
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=batch_key)
    vae = scvi.model.SCVI(adata, n_layers=2, n_latent=40, gene_likelihood="nb")
    vae.train()
    adata.obsm["X_scVI"] = vae.get_latent_representation()
    sc.pp.neighbors(adata, use_rep="X_scVI", key_added='scVI', n_neighbors=15, n_pcs=40)
    sc.tl.umap(adata, neighbors_key='scVI', min_dist=0.3) ## to match min_dist in seurat
    sc.pl.umap(adata, neighbors_key='scVI', color=[batch_key, species_key, cluster_key], ncols=1)
    plt.savefig(out_umap, dpi=300,  bbox_inches='tight')
    adata.obsm['X_umapscVI'] = adata.obsm['X_umap']
    click.echo("Save output")
    adata.write(out_h5ad)
    click.echo("Done scVI")


if __name__ == '__main__':
    run_scVI()

