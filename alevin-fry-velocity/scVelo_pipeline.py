import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import scvelo as scv
import matplotlib
import scipy
import json
import os

# Processing alevin-fry count matrix
frydir = "pancreas_quant_res"
e2n_path = "data/geneid_to_name.txt"
meta_info = json.load(open(os.path.sep.join([frydir, "meta_info.json"])))
ng = meta_info['num_genes']
usa_mode = meta_info['usa_mode']

if usa_mode:
    print("processing input in USA mode, will return A+S as the spliced count, and U as the unspliced count")
else:
    print("please follow previous steps to generate the ount matrix in the USA mode")
    assert(False)

af_raw = sc.read_mtx(os.path.sep.join([frydir, "alevin", "quants_mat.mtx"]))
ng = int(ng/3)
e2n = dict([ l.rstrip().split() for l in open(e2n_path).readlines()])
var_names = [ l.rstrip() for l in open(os.path.sep.join([frydir, "alevin", "quants_mat_cols.txt"])).readlines()][:ng]
var_names = [e2n[e] for e in var_names]

obs_names = [ l.rstrip() for l in open(os.path.sep.join([frydir, "alevin", "quants_mat_rows.txt"])).readlines() ]

x = af_raw.X
spliced = x[:,range(0,ng)] + x[:,range(2*ng,3*ng)]
unspliced = x[:,range(ng, 2*ng)]

# creating AnnData using spliced and unspliced count matrix
adata = anndata.AnnData(X = spliced, 
                        layers = dict(spliced = spliced, 
                                    unspliced = unspliced))
adata.obs_names = obs_names
adata.var_names = var_names
adata.var_names_make_unique()

# get embeddings
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.tsne(adata)
sc.tl.umap(adata, n_components = 2)

# housekeeping
matplotlib.use('AGG')
scv.settings.set_figure_params('scvelo')

# get the proportion of spliced and unspliced count
scv.utils.show_proportions(adata)

# filter cells and genes, then normalize expression values
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000,enforce=True)

# scVelo pipeline
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(adata, n_jobs = 11)
scv.tl.velocity(adata, mode = 'dynamical')
scv.tl.velocity_graph(adata)
adata.write('pancreas_full_dim_scvelo.h5ad', compression='gzip')
scv.pl.velocity_embedding_stream(adata, basis='X_umap', save="pancreas_full_dim")


# overwriting AnnData
adata = anndata.AnnData(X = spliced, 
                        layers = dict(spliced = spliced, 
                                    unspliced = unspliced))
# assign obs_names(cell IDs) and var_names(gene names)
adata.obs_names = obs_names
adata.var_names = var_names
adata.var_names_make_unique()

example_adata = scv.datasets.pancreas()
subset_adata = adata[example_adata.obs_names, np.intersect1d(example_adata.var_names, adata.var_names)]
example_adata = example_adata[subset_adata.obs_names, subset_adata.var_names]
subset_adata.obs = example_adata.obs
subset_adata.obsm["X_umap"] = example_adata.obsm["X_umap"]

# housekeeping
matplotlib.use('AGG')
scv.settings.set_figure_params('scvelo')

# get the proportion of spliced and unspliced count
scv.utils.show_proportions(subset_adata)

# filter cells and genes, then normalize expression values
scv.pp.filter_and_normalize(subset_adata, min_shared_counts=20, n_top_genes=2000,enforce=True)

# scVelo pipeline
scv.pp.moments(subset_adata, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(subset_adata, n_jobs = 11)
scv.tl.velocity(subset_adata, mode = 'dynamical')
scv.tl.velocity_graph(subset_adata)
subset_adata.write('pancreas_trimmed_scvelo.h5ad', compression='gzip')
scv.pl.velocity_embedding_stream(subset_adata, basis='X_umap', save="pancreas_trimmed")
