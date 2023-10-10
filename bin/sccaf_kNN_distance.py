#/usr/bin/env python3

# © EMBL-European Bioinformatics Institute, 2023
# Yuyao Song <ysong@ebi.ac.uk>

import pandas as pd
import numpy as np
import click
import scanpy as sc
from scipy.sparse import issparse

# for reading/saving clf model

from sklearn.preprocessing import MinMaxScaler, MaxAbsScaler
from sklearn.model_selection import cross_val_score, train_test_split
from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC

from sklearn.neighbors import KNeighborsClassifier

np.random.seed(seed=123)

@click.command()
@click.argument("input_h5ad", type=click.Path(exists=True))
@click.argument("out_acc_csv", type=click.Path(exists=False), default=None)
@click.option('--integration_method', type=str, default=None, help="Integration method")
@click.option('--cluster_key', type=str, default=None, help="Cluster key in input_h5ad.obs to use as the label to calculate SCCAF assessment")
@click.option('--n_neighbors', type=int, default=15, help="Number of neighbours for kNN calculation and assessment")


    
    
def anndata_knn_acc(input_h5ad, cluster_key, out_acc_csv, integration_method, n_neighbors=15):
    
    click.echo("Read input h5ad")
    input_ad = sc.read_h5ad(input_h5ad)
    
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
        "per_species": False,
    }

    
    use_embedding = use_embeddings[integration_method]
    
    if use_embedding is True:
        embedding_key = embedding_keys[integration_method]
        
    ## prepare the connectivity matrix
    
    if integration_method == "SAMap":
        click.echo("use SAMap KNN graph")
        # do nothing
    elif use_embedding is True:
        click.echo(f"Calculate KNN graph from embedding {embedding_key}")
        num_pcs = min(input_ad.obsm[embedding_key].shape[1], 40)
        if num_pcs < 40:
            click.echo(f"using {str(num_pcs)} PCs")
        sc.pp.neighbors(input_ad, n_neighbors=n_neighbors, n_pcs=num_pcs, use_rep=embedding_key, key_added=integration_method)
        
        # compute knn if use embedding
    elif integration_method == "per_species":
        click.echo('calculate NN graph from raw count for per species data')
        sc.pp.normalize_total(input_ad, target_sum=1e4)
        sc.pp.log1p(input_ad)
        sc.pp.highly_variable_genes(input_ad, min_mean=0.0125, max_mean=3, min_disp=0.5)
        input_ad.raw = input_ad
        sc.pp.scale(input_ad, max_value=10)
        sc.tl.pca(input_ad, svd_solver='arpack')
        sc.pp.neighbors(input_ad, n_neighbors=n_neighbors, n_pcs=40, use_rep="X_pca")
        embedding_key = "X_pca"
    elif integration_method == 'unintegrated':
        click.echo("use PCA to compute KNN graph for unintegrated data")
        
        sc.pp.normalize_total(input_ad, target_sum=1e4)
        sc.pp.log1p(input_ad)
        sc.pp.highly_variable_genes(input_ad, min_mean=0.0125, max_mean=3, min_disp=0.5)
        sc.pp.scale(input_ad, max_value=10)
        sc.tl.pca(input_ad, svd_solver="arpack")
        sc.pp.neighbors(input_ad, n_neighbors=n_neighbors, n_pcs=40, use_rep="X_pca")
        embedding_key = "X_pca"
        
    else:
        click.echo("use PCA to compute KNN graph for corrected counts output")
        sc.pp.scale(input_ad, max_value=10)
        sc.tl.pca(input_ad, svd_solver="arpack")
        sc.pp.neighbors(input_ad, n_neighbors=n_neighbors, n_pcs=40, use_rep="X_pca")
        embedding_key = "X_pca"
    
    # take the distance matrix as input for kNN classifier
    # in theory the "distance" slot is correct, while if we use "connectivities" the output is exactly the same
    # this is because a non-weighted majority vote is performed, essentially, the most frequent label of the neighbours
    # so the weight or distance does not matter numerically. I could also use weight='distance'
    if use_embedding and f"{integration_method}_distances" in input_ad.obsp.keys():
        X_use = input_ad.obsp[f"{integration_method}_distances"]
        distances_use = f"{integration_method}_distances"
        
    elif not use_embedding and "distances" in input_ad.obsp.keys():
        X_use = input_ad.obsp["distances"]
        distances_use = "distances"
    elif integration_method == 'SAMap':
        X_use = input_ad.obsp['connectivities']
        distances_use = "connectivities_SAMap"
        # in SAMap, connectivities is just 1 if neighbor, 0 if not
    else:
        raise ValueError("Error in neighbours calculation")
     
        # minmaxscaler requires dense input
        
    ## use half data for training and half for testing
    
    click.echo(f"testing on k={n_neighbors}")
        
    y_prob, y_pred, y_test, clf, cvsm, acc = self_projection(X_use, input_ad.obs[cluster_key], n=0, fraction=0.5, classifier="KNN", K_knn=n_neighbors)
    
    pd.DataFrame(data={"test_accuracy":acc, "CV_accuracy":cvsm, "integration_method": integration_method, "input_file": input_h5ad, "connectivities": distances_use}, index=[0]).to_csv(out_acc_csv, index=False)
    
    #return X_use, y_prob, y_pred, y_test, clf, cvsm, acc


def train_test_split_per_type_square(X, y, frac=0.5):
    """
    This function is identical to train_test_split, but can split the data either based on number of cells or by fraction.

    Input
    -----
    X: `numpy.array` or sparse matrix
        the feature matrix
    y: `list of string/int`
        the class assignments
    n: `int` optional (default: 100)
        maximum number sampled in each label
    fraction: `float` optional (default: 0.8)
        Fraction of data included in the training set. 0.5 means use half of the data for training,
        if half of the data is fewer than maximum number of cells (n).

    return
    -----
    X_train, X_test, y_train, y_test
    """
    df = pd.DataFrame(y)
    df.index = np.arange(len(y))
    df.columns = ['class']
    c_idx = df.groupby('class').apply(lambda x: x.sample(frac=frac)).index.get_level_values(None)
    d_idx = ~np.isin(np.arange(len(y)), c_idx)
    # the test matrix is n_test by n_indexed(trained)
    return X[c_idx, :][:, c_idx], X[d_idx, :][:, c_idx], y[c_idx], y[d_idx]


def inverse_weight(arr):
    result = arr
    for i in range(len(arr)):
        for j in range(len(arr[i])):
                result[i][j] = 1 / arr[i][j]
    return result

## make distance huge when not neighbours
def replace_zeros(arr):
    for i in range(len(arr)):
        for j in range(len(arr[i])):
            if arr[i][j] == 0:
                arr[i][j] = 1e5
    return arr


def self_projection(X,
                    cell_types,
                    classifier="LR",
                    K_knn=15,
                    metric_knn='precomputed',
                    penalty='l1',
                    sparsity=0.5,
                    fraction=0.5,
                    random_state=1,
                    solver='liblinear',
                    n=0,
                    cv=5,
                    whole=False,
                    n_jobs=None):
    # n = 100 should be good.
    """
    This is the core function for running self-projection.

    Input
    -----
    X: `numpy.array` or sparse matrix
        the expression matrix, e.g. ad.raw.X.
    cell_types: `list of String/int`
        the cell clustering assignment
    classifier: `String` optional (defatul: 'LR')
        a machine learning model in "LR" (logistic regression), "KNN" (k-nearest neighbour),\
        "RF" (Random Forest), "GNB"(Gaussion Naive Bayes), "SVM" (Support Vector Machine) and "DT"(Decision Tree).
    K_knn:  `int` optional (default: 15)
        the "k" in a KNN classifier if used
    metric_knn: `String` optional (default: 'precomputed')
        the distance metric for KNN classifier, default is to use a precomputed distance matrix in sklearn.neighbors.KNeighborsClassifier
    penalty: `String` optional (default: 'l2')
        the standardization mode of logistic regression. Use 'l1' or 'l2'.
    sparsity: `fload` optional (default: 0.5)
        The sparsity parameter (C in sklearn.linear_model.LogisticRegression) for the logistic regression model.
    fraction: `float` optional (default: 0.5)
        Fraction of data included in the training set. 0.5 means use half of the data for training,
        if half of the data is fewer than maximum number of cells (n).
    random_state: `int` optional (default: 1)
        random_state parameter for logistic regression.
    n: `int` optional (default: 100)
        Maximum number of cell included in the training set for each cluster of cells.
        only fraction is used to split the dataset if n is 0.
    cv: `int` optional (default: 5)
        fold for cross-validation on the training set.
        0 means no cross-validation.
    whole: `bool` optional (default: False)
        if measure the performance on the whole dataset (include training and test).
    n_jobs: `int` optional, number of threads to use with the different classifiers (default: None - unlimited).

    return
    -----
    y_prob, y_pred, y_test, clf
    y_prob: `matrix of float`
        prediction probability
    y_pred: `list of string/int`
        predicted clustering of the test set
    y_test: `list of string/int`
        real clustering of the test set
    clf: the classifier model.
    """
    # split the data into training and testing
    
    if classifier == 'KNN':
        scaler = MaxAbsScaler()
        
        # for kNN classifier, we need to normalize the data to unify the distance measurement between different algorithms
        X = scaler.fit_transform(X)
        assert X.shape[0] == X.shape[1], "Matrix is not square, is a connectivity matrix used as input for kNN classifier?"

    if issparse(X):
        X = X.todense()
        click.echo("to dense after maxabs scaling")
        
    X_train, X_test, y_train, y_test = train_test_split_per_type_square(X, cell_types, frac=fraction)
    X_train = replace_zeros(np.array(X_train))
    X_test = replace_zeros(np.array(X_test))

        # fraction means test size
    # set the classifier

    if classifier == 'LR':
        clf = LogisticRegression(random_state=1, penalty=penalty, C=sparsity, multi_class="ovr", solver=solver)
    elif classifier == 'KNN':
        clf = clf = KNeighborsClassifier(n_neighbors=K_knn, metric=metric_knn, n_jobs=n_jobs, weights=inverse_weight)


    # mean cross validation score
    cvsm = 0
    if cv > 0:
        cvs = cross_val_score(clf, X_train, np.array(y_train), cv=cv, scoring='accuracy', n_jobs=n_jobs)
        cvsm = cvs.mean()
        print("Mean CV accuracy: %.4f" % cvsm)

    # accuracy on cross validation and on test set
    clf.fit(X_train, y_train)
    accuracy = clf.score(X_train, y_train)
    print("Accuracy on the training set: %.4f" % accuracy)
    accuracy_test = clf.score(X_test, y_test)
    print("Accuracy on the hold-out set: %.4f" % accuracy_test)

    # accuracy of the whole dataset
    if whole:
        accuracy = clf.score(X, cell_types)
        print("Accuracy on the whole set: %.4f" % accuracy)

    # get predicted probability on the test set
    y_prob = None
    #if not classifier in ['SH', 'PCP']:
    y_prob = clf.predict_proba(X_test)
    y_pred = clf.predict(X_test)

    return y_prob, y_pred, y_test, clf, cvsm, accuracy_test

if __name__ == '__main__':
    anndata_knn_acc()
