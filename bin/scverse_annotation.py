#!/usr/bin/env python3

import os
import tempfile
import argparse
import pandas as pd
import numpy as np
import scanpy as sc
import scvi
import seaborn as sns
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

scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)

sc.set_figure_params(figsize=(6, 6), frameon=False)
sns.set_theme()
torch.set_float32_matmul_precision("high")
save_dir = tempfile.TemporaryDirectory()

adata = sc.read(args.adata)

celltype_markers = pd.read_csv(args.markers,
                               index_col=0,)

adata.var_names_make_unique()
adata.obs_names_make_unique()

celltype_markers = celltype_markers.loc[:, celltype_markers.columns.isin(adata.var.index)]

bdata = adata[:, adata.var_names.isin(celltype_markers.index)]

lib_size = np.asarray(bdata.X.sum(1))
bdata.obs["size_factor"] = lib_size / np.mean(lib_size)

print(bdata)

scvi.external.CellAssign.setup_anndata(bdata,
                                       size_factor_key="size_factor")

model = CellAssign(bdata, 
                   celltype_markers)

model.train()

predictions = model.predict()

bdata.obs["cellassign_predictions"] = predictions.idxmax(axis=1).values

# celltype is the original CellAssign prediction
sc.pl.umap(bdata,
           color=["cellassign_predictions"],
           save = "_scverse.png")
