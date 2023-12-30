
#!/usr/bin/env python
# coding: utf-8

import ncem as nc
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq
import pandas as pd 

sc.settings.set_figure_params(dpi=80)
import warnings
warnings.filterwarnings("ignore")

sc.logging.print_header()

print(f"ncem=={nc.__version__}")

from ncem.interpretation import InterpreterInteraction
from ncem.data import get_data_custom, customLoader


ad = sq.read.vizgen(path='./',
    counts_file="./cell_by_gene.csv",
    meta_file="./cell_metadata.csv")

annot = pd.read_csv('./ARKFrozen62PFC_annotation.csv',
                   sep=',', index_col=0)

annot = annot[['type']]

ad = ad[annot.index,]

ad.obs = ad.obs.join(annot)

ad = ad[ad.obs.type.isin(ad.obs.type.value_counts()[ad.obs.type.value_counts() > 100].index),:]

ad.layers["counts"] = ad.X.copy()
sc.pp.normalize_total(ad, inplace=True)
sc.pp.log1p(ad)
sc.pp.pca(ad)
sc.pp.neighbors(ad)

ad.obs['Cluster'] = ad.obs['type']

ncem = InterpreterInteraction()

ad.obs['donor'] = 'unknown1'
ad.obs['library_id'] = 'unknown2'

ad.obs['donor'] = ad.obs['donor'].astype('category')
ad.obs['library_id'] = ad.obs['library_id'].astype('category')
ad.obs['library_id'].value_counts()

ncem.data = customLoader(adata=ad, 
                         cluster='Cluster', 
                         patient='donor', 
                         library_id='library_id', 
                         radius=40)
get_data_custom(interpreter=ncem)


var_decomp = ncem.data.compute_variance_decomposition()

ncem.data.variance_decomposition(var_decomp, figsize=(5,3), save='./ARKFrozen62PFC_r40')

var_decomp.mean(axis=0)[['intra cell type variance', 'inter cell type variance', 'gene variance']]

print('Completed variance decomposition')

# Type coupling analysis

ncem.get_sender_receiver_effects()

## Magnitude
type_coupling = ncem.type_coupling_analysis_circular(edge_attr='magnitude', 
                                                     figsize=(15,12), 
                                                     de_genes_threshold=25,
                                                     text_space=1.28,
                                                    save='./ARKFrozen62PFC_r40_magnitude')

type_coupling.to_csv("./ARKFrozen62PFC_r40_magnitude_all.csv", sep=',', index=False)
type_coupling[(type_coupling['de_genes_abs'] > 25) & (type_coupling['magnitude'] > 1)].to_csv('./ARKFrozen62PFC_r40_magnitude.csv', sep=',', index=None)

## DEGs
type_coupling = ncem.type_coupling_analysis_circular(edge_attr='de_genes', 
                                                     figsize=(15,12), 
                                                     de_genes_threshold=25,
                                                     text_space=1.28,
                                                     save='./ARKFrozen62PFC_r40_DEG')


