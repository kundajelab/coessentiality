import numpy as np, pandas as pd, time, os, subprocess, scipy as sp, diffmap as dm
import importlib, matplotlib.pyplot as plt, sklearn.covariance as skcov
import scipy.io, sklearn.metrics
from scipy.sparse import csr_matrix

# Use two single-cell libraries to call approx nearest neighbors and UMAP.
import scanpy as sc, anndata as adata

itime = time.time()

# Read in GLS -log(p) values.
pairwise_distances = np.load('GLS_p.npy')
gene_names = np.ravel(pd.read_csv('genes.txt', header=None))

# Read in clusterONE analysis results.
cone_clusts = pd.read_csv('clusterOne_clusters.tsv', sep="\t", index_col=0, header=0)
clustdata = cone_clusts.iloc[:,10:].fillna(value='')
nclusts = clustdata.shape[0]
cluster_ndces = {}
for gname in gene_names:
    cluster_ndces[gname] = np.where(clustdata == gname)[0]
    if len(cluster_ndces)%500 == 0:
        print("{} genes finished. \tTime: {}".format(len(cluster_ndces), time.time() - itime))
clust_mmships = np.zeros((len(gene_names), nclusts))
for i in range(len(gene_names)):
    clust_mmships[i, cluster_ndces[gene_names[i]]] = 1
clust_mmships = sp.sparse.csr_matrix(clust_mmships)
sp.sparse.save_npz('clusterone_memberships.npz', clust_mmships)
# Build and save pairwise Jaccard similarities between genes, according to the clusterings given.
gg_jaccsims = dm.pairwise_jaccard_graph(clust_mmships)
gm_inc = dm.build_knn(gg_jaccsims, k=10, symmetrize_type='inclusive')

# Use GLS -log(p) values between each pair of genes (the (genes x genes) matrix GLS_log_p) as the adjacency matrix of the GLS graph.
a = -np.log(GLS_log_p)
a[np.isinf(a)] = 0
GLS_pvals_100 = dm.build_knn(a, k=100, symmetrize_type='inclusive')
GLS_pvals_10 = dm.build_knn(GLS_pvals_100, k=10, symmetrize_type='inclusive')
sp.sparse.save_npz('GLS_pvals_10NN.npz', GLS_pvals_10)

# Construct the combined graph
frac_CO_graph = 0.99
frac_GLS_graph = 1-frac_CO_graph
CO_graph = gm_inc
GLS_graph = GLS_pvals_10

adj_mat = sp.sparse.csr_matrix(
    (frac_CO_graph * CO_graph) + 
    (frac_GLS_graph * GLS_graph)
)
n_cmps = 100
reduced_dim = 50
sffix = "_GLS01_CO99"
vizdf_filename = "vizdf{}.csv".format(sffix)

# emap_naive, eigvals = dm.diffmap_proj(adj_mat, t=0, n_comps=reduced_dim, embed_type='naive', return_eigvals=True)
# print("Laplacian eigenmap computed. Time: {}".format(time.time() - itime))
emap_heat = dm.diffmap_proj(adj_mat, n_comps=n_cmps, n_dims=40, min_energy_frac=0.9, embed_type='diffmap', return_eigvals=False)
print("Diffusion components computed. Time: {}".format(time.time() - itime))
ann_heat = adata.AnnData(X=emap_heat[:, :40])

sc.pp.neighbors(ann_heat)
print(time.time() - itime)
sc.tl.umap(ann_heat)
print(time.time() - itime)
heat_umap = np.array(ann_heat.obsm['X_umap'])
vizdf = pd.DataFrame(data=heat_umap, columns=['hUMAP_x', 'hUMAP_y'])
vizdf['gene_names'] = gene_names
vizdf.to_csv(vizdf_filename, sep="\t", index=False)