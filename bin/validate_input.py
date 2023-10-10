#/usr/bin/env python3

# © EMBL-European Bioinformatics Institute, 2023
# Yuyao Song <ysong@ebi.ac.uk>

# ysong update pydantic to v2 for containerization Oct 2023

import click
from pydantic import StringConstraints, Field, ConfigDict, BaseModel, ValidationError
import pandas as pd
import numpy as np
import scanpy as sc
from typing import List
from typing_extensions import Annotated

@click.command()
@click.argument("input_metadata", type=click.Path(exists=True))
@click.option('--batch_key', type=str, default=None, help="Column in adata.obs for batch label")
@click.option('--species_key', type=str, default=None, help="Column in adata.obs for species label")
@click.option('--cluster_key', type=str, default=None, help="Column in adata.obs for cell type label")

def validate_adata_input(input_metadata, batch_key, cluster_key, species_key):

    ## define classes to validate
    meta = pd.read_csv(input_metadata, sep = '\t', header=None)

    class adata_obs_for_csi(BaseModel):
        species_key: pd.Series
        cluster_key: pd.Series
        batch_key: pd.Series
        model_config = ConfigDict(arbitrary_types_allowed=True)

    class adata_X_for_csi(BaseModel):
        X: np.ndarray
        model_config = ConfigDict(arbitrary_types_allowed=True)

    class adata_var_for_csi(BaseModel):
        mean_counts: pd.Series
        var_names: List[Annotated[str, Field(pattern="^ENS[A-Z]{3}[GP][0-9]{11}$|^ENSG[0-9]{11}(\.[0-9]+_[A-Z]+_[A-Z]+)?$")]] # require ensembl gene id as var_names 
        model_config = ConfigDict(arbitrary_types_allowed=True)

    ## validate

    for i in range(0, meta.shape[0]):
        print('validating input for species '+ meta.iloc[i, 0])
        print('input anndata path is '+ meta.iloc[i, 1])

        species_now = meta.iloc[i, 0]

        print('read in adata')
        ad_now = sc.read_h5ad(meta.iloc[i, 1])

        ## validate required fields in adata.obs
        for key in {species_key, batch_key, cluster_key}:
            try:
                ad_now.obs[key]
            except KeyError as ke:
                print(str(ke) + ' does not exist in ' + species_now + ' adata.obs')
                raise
            else:
                try:
                    adata_obs_for_csi(species_key = ad_now.obs[species_key], batch_key = ad_now.obs[batch_key], cluster_key = ad_now.obs[cluster_key] )
                    print('obs test pass for ' + species_now)
                except ValidationError as ve:
                    print(ve.json())
                    raise
        ## validate adata.X is a dense matrix in np.array and all positive values (i.e. non-scaled data)
        try:
            adata_X_for_csi(X = ad_now.X)
            print('count matrix test pass for ' + species_now)
        except ValidationError as e:
            print(e.json())
        else:
            if not all(ad_now.X.flatten() >= 0):
                raise ValueError('values in adata.X are not all positive, please make sure raw count matrix is in adata.X')

        ## validate required fields in adata.var, and adata.var_names are ensembl gene ids
        try:
            ad_now.var['mean_counts']
        except KeyError as ke:
            print(str(ke) + ' does not exist in ' + species_now + ' adata.var')
            raise
        else:
            try:
                adata_var_for_csi(mean_counts = ad_now.var['mean_counts'], var_names = ad_now.var_names.values.tolist())
                print('var test pass for ' + species_now)
            except ValidationError as e:
                print(e.json())
                raise

    print('input validation complete')

if __name__ == '__main__':
    validate_adata_input()
