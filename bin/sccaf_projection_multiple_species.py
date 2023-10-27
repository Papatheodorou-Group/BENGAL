#/usr/bin/env python3

# © EMBL-European Bioinformatics Institute, 2023
# Yuyao Song <ysong@ebi.ac.uk>


import click
from typing import List

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import scanpy as sc
import anndata as ad
from anndata import AnnData
import itertools
from SCCAF import *
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn import metrics


@click.command()

@click.argument("input_h5ad", type=click.Path(exists=True))
@click.argument("out_projection_h5ads", type=click.Path(exists=False), default=None)
@click.argument("out_figures", type=click.Path(exists=False), default=None)
@click.argument("out_acc_csv", type=click.Path(exists=False), default=None)
@click.option('--species_key', type=str, default=None, help="Species key to distinguish species")
@click.option('--integration_method', type=str, default=None, help="Integration method")
@click.option('--cluster_key', type=str, default=None, help="Cluster key in species to use to assess label preservation")
@click.option('--projection_key', type=str, default=None, help="Projection key in species one to use as labels to transfer to species two")


def run_sccaf_projection(input_h5ad, species_key,  cluster_key, projection_key, integration_method, out_projection_h5ads, out_figures, out_acc_csv):
# dictionary for method properties
    embedding_keys={"harmony": "X_pca_harmony", "scanorama": "X_scanorama", "scVI": "X_scVI", "LIGER": "X_inmf", "rligerUINMF":"X_inmf", "fastMNN": "X_mnn", "SAMap": "wPCA", "scANVI": "X_scANVI"}
    use_embeddings={"harmony": True, "scanorama": True, "scVI": True, "LIGER": True, "rligerUINMF":True, "fastMNN": True, "SAMap": True , "seuratCCA": False, "seuratRPCA": False, "scANVI": True, "unintegrated": False}
    from_h5seurat={"harmony": False, "scanorama": False, "scVI": False, "LIGER": True, "rligerUINMF":True, "fastMNN": True, "SAMap": False , "seuratCCA": True, "seuratRPCA": True, "scANVI": False, "unintegrated": False}
    sc.set_figure_params(dpi_save=200, frameon=False, figsize=(11, 6))
    input_ad = sc.read_h5ad(input_h5ad)
    species_all=input_ad.obs[species_key].astype("category").cat.categories.values
    acc_summary=pd.DataFrame()
     
    if integration_method == 'unintegrated':

    ## raw concatinated data needs preprocessing
        sc.pp.normalize_total(input_ad, target_sum=1e4)
        sc.pp.log1p(input_ad)
        sc.pp.highly_variable_genes(input_ad, min_mean=0.0125, max_mean=3, min_disp=0.5)
        input_ad.raw = input_ad
        sc.pp.scale(input_ad, max_value=10)
        sc.tl.pca(input_ad, svd_solver='arpack')
        sc.pp.neighbors(input_ad, n_neighbors=10, n_pcs=40)
        sc.tl.umap(input_ad, min_dist=0.3)

    # register color in .uns
    sc.pl.umap(input_ad, color = cluster_key)
    sc.pl.umap(input_ad, color = projection_key)

    click.echo("Start SCCAF projection workflow with input: "+input_h5ad)
    use_embedding=use_embeddings[integration_method]
    if use_embedding is True:
        embedding_key=embedding_keys[integration_method]


    # known bug - fix when convert h5Seurat to h5ad the index name error
    #if from_h5seurat[integration_method] is True:
    #    input_ad.__dict__['_raw'].__dict__['_var'] = input_ad.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
    #    click.echo("From h5seurat")

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

    ##def unique_combinations(elements: list[str]) -> list[tuple[str, str]]:
    def unique_combinations(elements) -> list:
    ## Precondition: `elements` does not contain duplicates.
    ## Postcondition: Returns unique combinations of length 2 from `elements`.

        return list(itertools.combinations(elements, 2))

    species_pairs=unique_combinations(species_all)

    with PdfPages(out_figures) as pdf:
        for species_pair in species_pairs:
            species_one=species_pair[0]
            species_two=species_pair[1]

            adx1 = input_ad[input_ad.obs[species_key]==species_one, :]
            adx2 = input_ad[input_ad.obs[species_key]==species_two, :]

            if use_embedding is True and not embedding_key in input_ad.obsm.keys():
                raise ValueError("`adata.obsm['%s']` doesn't exist. Please assign the embedding key if use embedding."%(embedding_key))
            if use_embedding is True and embedding_key is not None:
                input_matrix_1 = adx1.obsm[embedding_key]
                input_matrix_2 = adx2.obsm[embedding_key]
            else:
                input_matrix_1 = adx1.X
                input_matrix_2 = adx2.X
            click.echo("running projection between "+species_one+" and "+species_two)

            colors1 = adx1.uns[cluster_key+"_colors"]

            y_prob1, y_pred1, y_test1, clf1, cvsm1, acc1 = SCCAF_assessment(input_matrix_1, adx1.obs[cluster_key], n=200)
            aucs1 = plot_roc(y_prob1, y_test1, clf1, cvsm=cvsm1, acc=acc1, colors=colors1, title=species_one+"_"+integration_method+"_original_label")
            pdf.attach_note("AUC of original label in "+species_one+" by "+integration_method)
            pdf.savefig(dpi=200, bbox_inches='tight')
            plt.close()

            tbl1=get_cell_type_auc(clf1, y_test1, y_prob1)

            tbl1[['test_acc']] = acc1
            tbl1[['CV_acc']] = cvsm1
            tbl1[['type_label']] = 'original'
            tbl1['from_species'] = species_one
            tbl1['to_species'] = species_one
            tbl1['integration_method'] = integration_method
            tbl1['input_file'] = input_h5ad
            tbl1['key_use'] = cluster_key
            tbl1['adj_rand_score'] = 'NaN'
            tbl1['pct_cell_type_kept'] = 'NaN'


            colors2 = adx2.uns[cluster_key+"_colors"]

            y_prob2, y_pred2, y_test2, clf2, cvsm2, acc2 = SCCAF_assessment(input_matrix_2, adx2.obs[cluster_key], n=200)
            aucs2 = plot_roc(y_prob2, y_test2, clf2, cvsm=cvsm2, acc=acc2, colors=colors2, title=species_two+"_"+integration_method+"_original_label")
            pdf.attach_note("AUC of original label in "+species_two+" by "+integration_method)
            pdf.savefig(dpi=200, bbox_inches='tight')
            plt.close()

            tbl2=get_cell_type_auc(clf2, y_test2, y_prob2)

            tbl2['test_acc'] = acc2
            tbl2[['CV_acc']] = cvsm2
            tbl2[['type_label']] = 'original'
            tbl2['from_species'] = species_two
            tbl2['to_species'] = species_two
            tbl2['integration_method'] = integration_method
            tbl2['input_file'] = input_h5ad
            tbl2['key_use'] = cluster_key
            tbl2['adj_rand_score'] = 'NaN'
            tbl2['pct_cell_type_kept'] = 'NaN'


            # run projection only on shared cell types
            adx1 = input_ad[input_ad.obs[species_key]==species_one, :]
            adx2 = input_ad[input_ad.obs[species_key]==species_two, :]

            a = adx1.obs[projection_key].cat.categories.tolist()
            b = adx2.obs[projection_key].cat.categories.tolist()

            shared_ct =  list(set(a) & set(b))
            ##shared_ct.remove('T_cell')
            adx1 = adx1[adx1.obs[projection_key].isin(shared_ct), :]
            adx2 = adx2[adx2.obs[projection_key].isin(shared_ct), :]

            color = dict(zip(input_ad.obs[projection_key].cat.categories.to_list(), input_ad.uns[projection_key+'_colors']))

            color_palette = {key: value for key, value in color.items() if key in shared_ct}

            if use_embedding is True and not embedding_key in input_ad.obsm.keys():
                raise ValueError("`adata.obsm['%s']` doesn't exist. Please assign the embedding key if use embedding."%(embedding_key))
            if use_embedding is True and embedding_key is not None:
                input_matrix_1 = adx1.obsm[embedding_key]
                input_matrix_2 = adx2.obsm[embedding_key]
            else:
                input_matrix_1 = adx1.X
                input_matrix_2 = adx2.X


            ## species one label inferred from species two
            y_prob2, y_pred2, y_test2, clf2, cvsm2, acc2 = SCCAF_assessment(input_matrix_2, adx2.obs[projection_key], n=200)
            adx1.obs['logit_inferred'] = clf2.predict(input_matrix_1)
            pct_1 = len(set(adx1.obs['logit_inferred'].astype("category").cat.categories.tolist()) & set(adx1.obs[projection_key].astype("category").cat.categories.tolist())) / len(adx1.obs[projection_key].astype("category").cat.categories.tolist())

            umap1 = sc.pl.umap(adx1, color=[projection_key,"logit_inferred"], palette=color_palette, title = [species_one+' original label', species_one+' inferred label by '+species_two], ncols=2, frameon=False, show=False)
            pdf.attach_note("UMAP of species " + species_one + " original and transferred label")
            pdf.savefig(dpi=200, bbox_inches='tight')
            plt.close()

            ars_1  = adjusted_rand_score(adx1.obs['logit_inferred'], adx1.obs[projection_key])

            # sometimes the projection doesn't work at all and there will be only one cell type left after the projection
            if len(set(adx1.obs['logit_inferred'].astype("category").cat.categories.tolist())) <= 2:
                click.echo("projection doesn't work, return NA")
                tbl3=pd.DataFrame(data=list(zip(['logit_not_working'], ["NaN"], ["NaN"])), columns=['cell_type', "ROC_AUC", "PR_AUC"])

                tbl3['test_acc'] = 'NaN'
                tbl3[['CV_acc']] = 'NaN'
                tbl3[['type_label']] = 'logit_inferred'
                tbl3['from_species'] = species_two
                tbl3['to_species'] = species_one
                tbl3['integration_method'] = integration_method
                tbl3['input_file'] = input_h5ad
                tbl3['key_use'] = projection_key
                tbl3['adj_rand_score'] = ars_1
                tbl3['pct_cell_type_kept'] = pct_1

            else:
                click.echo("projection worked")
                colors3 = adx1.uns['logit_inferred_colors']
                y_prob3, y_pred3, y_test3, clf3, cvsm3, acc3 = SCCAF_assessment(input_matrix_1, adx1.obs['logit_inferred'],n=200)
                aucs3 = plot_roc(y_prob3, y_test3, clf3, cvsm=cvsm3, acc=acc3, colors=colors3)
                pdf.attach_note("AUC of logit inferred label in "+species_one)
                pdf.savefig(dpi=200, bbox_inches='tight')
                plt.close()

                tbl3=get_cell_type_auc(clf3, y_test3, y_prob3)

                tbl3['test_acc'] = acc3
                tbl3[['CV_acc']] = cvsm3
                tbl3[['type_label']] = 'logit_inferred'
                tbl3['from_species'] = species_two
                tbl3['to_species'] = species_one
                tbl3['integration_method'] = integration_method
                tbl3['input_file'] = input_h5ad
                tbl3['key_use'] = projection_key
                tbl3['adj_rand_score'] = ars_1
                tbl3['pct_cell_type_kept'] = pct_1



            y_prob1, y_pred1, y_test1, clf1, cvsm1, acc1 = SCCAF_assessment(input_matrix_1, adx1.obs[projection_key], n=200)

            adx2.obs['logit_inferred'] = clf1.predict(input_matrix_2)
            pct_2 = len(set(adx2.obs['logit_inferred'].astype("category").cat.categories.tolist()) & set(adx2.obs[projection_key].astype("category").cat.categories.tolist())) / len(adx2.obs[projection_key].astype("category").cat.categories.tolist())
            umap2 = sc.pl.umap(adx2, color=[projection_key,"logit_inferred"], palette=color_palette, title = [species_two+' original label', species_two+' inferred label by '+species_one], ncols=2, frameon=False, show=False)
            pdf.attach_note("UMAP of species " + species_two + " original and transferred label")
            pdf.savefig(dpi=200, bbox_inches='tight')
            plt.close()

            ars_2  = adjusted_rand_score(adx2.obs['logit_inferred'], adx2.obs[projection_key])

            if len(set(adx2.obs['logit_inferred'].astype("category").cat.categories.tolist())) <= 2:
                click.echo("projection doesn't work, return NA")
                tbl4=pd.DataFrame(data=list(zip(['logit_not_working'], ["NaN"], ["NaN"])), columns=['cell_type', "ROC_AUC", "PR_AUC"])

                tbl4['test_acc'] = 'NaN'
                tbl4[['CV_acc']] = 'NaN'
                tbl4[['type_label']] = 'logit_inferred'
                tbl4['from_species'] = species_one
                tbl4['to_species'] = species_two
                tbl4['integration_method'] = integration_method
                tbl4['input_file'] = input_h5ad
                tbl4['key_use'] = projection_key
                tbl4['adj_rand_score'] = ars_2
                tbl4['pct_cell_type_kept'] = pct_2

            else:
                click.echo("projection worked")
                colors4 = adx2.uns['logit_inferred_colors']
                y_prob4, y_pred4, y_test4, clf4, cvsm4, acc4 = SCCAF_assessment(input_matrix_2, adx2.obs['logit_inferred'],n=200)
                aucs4 = plot_roc(y_prob4, y_test4, clf4, cvsm=cvsm4, acc=acc4, colors=colors4)
                pdf.attach_note("AUC of logit inferred label in " + species_two)
                pdf.savefig(dpi=200, bbox_inches='tight')
                plt.close()

                tbl4=get_cell_type_auc(clf4, y_test4, y_prob4)

                tbl4['test_acc'] = acc4
                tbl4[['CV_acc']] = cvsm4
                tbl4[['type_label']] = 'logit_inferred'
                tbl4['from_species'] = species_one
                tbl4['to_species'] = species_two
                tbl4['integration_method'] = integration_method
                tbl4['input_file'] = input_h5ad
                tbl4['key_use'] = projection_key
                tbl4['adj_rand_score'] = ars_2
                tbl4['pct_cell_type_kept'] = pct_2

            acc_summary = pd.concat([acc_summary, tbl1, tbl2, tbl3, tbl4], axis=0, ignore_index=True)
            adata_out = ad.concat([adx1, adx2], join='inner', merge='same', label='sccaf_projection', keys=[species_one+"_from_"+species_one+"_"+species_two+"_"+integration_method, species_two+"_from_"+species_one+"_"+species_two+"_"+integration_method])
                
            adata_out.write(out_projection_h5ads+"_"+species_one+"_"+species_two+".h5ad", compression = 'gzip')
            click.echo("finish projection between "+species_one+" and "+species_two)

    click.echo("write acc summary for all")
    acc_summary.to_csv(out_acc_csv, index=False)

if __name__ == '__main__':
    run_sccaf_projection()


