#/usr/bin/env python3

# © EMBL-European Bioinformatics Institute, 2023
# Yuyao Song <ysong@ebi.ac.uk>


import click
import matplotlib.pyplot as plt
from sklearn import metrics
import pandas as pd
import scanpy as sc
from SCCAF import *

@click.command()
@click.argument("input_metadata", type=click.Path(exists=True))
@click.argument("out_auc", type=click.Path(exists=False), default=None)
@click.argument("out_acc_csv", type=click.Path(exists=False), default=None)
@click.option('--integration_method', type=str, default=None, help="Integration method")
@click.option('--use_embedding', type=bool, default=False, help="Whether use embedding for SCCAF assessment, default False to use count matrix")
@click.option('--embedding_key', type=str, default=None, help="If use embedding, the embedding key in input_h5ad.obsm to calculate SCCAF assessment")
@click.option('--cluster_key', type=str, default=None, help="Cluster key in input_h5ad.obs to use as the label to calculate SCCAF assessment")

def run_sccaf_assessment(input_metadata, out_auc, use_embedding, embedding_key, cluster_key, out_acc_csv, integration_method):

    def get_cell_type_auc(clf, y_test, y_prob):
        rc_aucs = [] #AUC
        rp_aucs = [] # AUC from recall precision
        fprs = [] #FPR
        tprs = [] #TPR
        prss = [] #Precision
        recs = [] #Recall
        for i, cell_type in enumerate(clf.classes_):
            fpr, tpr, _ = metrics.roc_curve(y_test == cell_type, y_prob[:, i])
            prs, rec, _ = metrics.precision_recall_curve(y_test == cell_type, y_prob[:, i])
            fprs.append(fpr)
            tprs.append(tpr)
            prss.append(prs)
            recs.append(rec)
            rc_aucs.append(metrics.auc(fpr, tpr))
            rp_aucs.append(metrics.auc(rec, prs))
        tbl = pd.DataFrame(data=list(zip(clf.classes_, rp_aucs, rc_aucs)), columns=['cell_type', "ROC_AUC", "PR_AUC"])
        return tbl

    meta = pd.read_csv(input_metadata, delimiter='\t', header=None)
    for i in range(0, meta.shape[0]):
        species=meta.loc[i, 0]
        input_h5ad=meta.loc[i, 1]
        print(species)
        input_ad = sc.read_h5ad(input_h5ad)

        sc.pp.normalize_total(input_ad, target_sum=1e4)
        sc.pp.log1p(input_ad)
        sc.pp.highly_variable_genes(input_ad, min_mean=0.0125, max_mean=3, min_disp=0.5)
        input_ad.raw = input_ad
        sc.pp.scale(input_ad, max_value=10)
        sc.tl.pca(input_ad, svd_solver='arpack')
        sc.pp.neighbors(input_ad, n_neighbors=10, n_pcs=40)
        sc.tl.umap(input_ad, min_dist=0.3)
        sc.pl.umap(input_ad, color = [cluster_key])
        sc.set_figure_params(dpi_save=300, frameon=False, figsize=(10, 5))
        ##input_ad = input_ad[~(input_ad.obs.cell_type=='T_cell'), ]

        if use_embedding is True and not embedding_key in input_ad.obsm.keys():
            raise ValueError("`input_ad.obsm['%s']` doesn't exist. Please assign the embedding key if use embedding."%(embedding_key))
        if use_embedding is True and embedding_key is not None:
            input_matrix = input_ad.obsm[embedding_key]
        else:
            input_matrix = input_ad.X

        colors = input_ad.uns[cluster_key+"_colors"]
        y_prob, y_pred, y_test, clf, cvsm, acc = SCCAF_assessment(input_matrix, input_ad.obs[cluster_key], n=200)
        aucs = plot_roc(y_prob, y_test, clf, cvsm=cvsm, acc=acc, colors=colors)

        plt.savefig(out_auc+"_"+species+".png", dpi=200, bbox_inches='tight')

        tbl1=get_cell_type_auc(clf, y_test, y_prob)

        tbl1[['test_acc']] = acc
        tbl1[['CV_acc']] = cvsm
        tbl1[['type_label']] = 'original'
        tbl1['from_species'] = species
        tbl1['to_species'] = species
        tbl1['integration_method'] = integration_method
        tbl1['input_file'] = input_h5ad
        tbl1['key_use'] = cluster_key
        tbl1['adj_rand_score'] = 'NaN'
        tbl1['pct_cell_type_kept'] = 'NaN'

        tbl1.to_csv(out_acc_csv+"_"+species+".csv", index=False, header=True)

if __name__ == '__main__':
    run_sccaf_assessment()
