#!/usr/bin/env python3

import os
import tempfile
import argparse
import pandas as pd
import numpy as np
import scanpy as sc
import scvi
import seaborn as sns
import matplotlib.pyplot as plt
import torch
from scvi.external import CellAssign

parser = argparse.ArgumentParser()

parser.add_argument('--adata',
                    type=str,
                    help="Path to adata.h5ad file")

parser.add_argument('--markers',
                    type=str,
                    help='Path to marker genes .csv')

args = parser.parse_args()

scvi.settings.seed = 123
print("Last run with scvi-tools version:", scvi.__version__)

adata = sc.read_h5ad(args.adata)

lib_size = adata.X.sum(axis=1)
adata = adata[lib_size > 0]

lib_size = np.asarray(adata.X.sum(axis=1))
adata.obs["size_factor"] = lib_size / np.mean(lib_size)

sc.pl.umap(adata, color=["louvain"], save="_teste1.png")

celltype_markers = pd.read_csv(args.markers,
                               index_col=0,)

celltype_markers['Unknown'] = 0

celltype_markers = celltype_markers.loc[celltype_markers.index.isin(adata.var.index), :]

bdata = adata[:, celltype_markers.index].copy()

sc.pl.umap(bdata, color=["louvain"], save="_teste2.png")

scvi.external.CellAssign.setup_anndata(bdata,
                                       layer="counts",
                                       size_factor_key="size_factor")

model = CellAssign(bdata, 
                   celltype_markers)

model.train()

predictions = model.predict()

bdata.obs["cellassign_predictions"] = predictions.idxmax(axis=1).values

sc.pl.umap(bdata,
           color=["cellassign_predictions"],
           save = "_scverse.png")