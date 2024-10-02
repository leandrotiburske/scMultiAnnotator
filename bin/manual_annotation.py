#!/usr/bin/env python3

#
# Load required libraries
#

import pandas as pd
import numpy as np
import anndata as an
import scanpy as sc
import urllib.request
import argparse
import logging
import urllib.request

parser = argparse.ArgumentParser()

parser.add_argument('--dataset',
                    type=str,
                    help="Path to the .h5ad dataset file.")

parser.add_argument('--markers',
                    type=str,
                    help="Path to .csv file with markers.")

parser.add_argument('--clustering',
                    type=str,
                    help="Clustering method used")

args = parser.parse_args()

adata = sc.read(args.dataset)

markers_csv = pd.read_csv(args.markers,                              
                          index_col=0,)

#markers_csv = markers_csv.loc[:, markers_csv.columns.isin(adata.var.index)]


markers_dict = dict()

for cell in markers_csv.columns:
    cell_markers = pd.DataFrame(markers_csv.loc[:,cell])
    cell_markers = cell_markers.loc[cell_markers[cell] == 1]
    markers_dict[cell] = list(cell_markers.index)

marker_genes_in_data = dict()
for ct, markers in markers_dict.items():
    markers_found = list()
    for marker in markers:
        if marker in adata.var.index:
            markers_found.append(marker)
    marker_genes_in_data[ct] = markers_found

cells = list(marker_genes_in_data.keys())

for ct in cells:
    print(f"{ct.upper()}:")  # print cell subtype name
    print(marker_genes_in_data[ct])
    sc.pl.umap(
        adata,
        save = "_{}_markers.png".format(ct),
        color=marker_genes_in_data[ct],
        vmin=0,
        vmax="p99",  # set vmax to the 99th percentile of the gene count instead of the maximum, to prevent outliers from making expression in other cells invisible. Note that this can cause problems for extremely lowly expressed genes.
        sort_order=False,  # do not plot highest expression on top, to not get a biased view of the mean expression among cells
        frameon=False,
        cmap="Reds",  # or choose another color map e.g. from here: https://matplotlib.org/stable/tutorials/colors/colormaps.html

    )


celltype_markers = {
    ct: [m for m in ct_markers if m in adata.var.index]
    for ct, ct_markers in markers_dict.items()
    if ct in cells
}

sc.pl.dotplot(
    adata,
    groupby = args.clustering,
    var_names = celltype_markers,
    standard_scale = "var",  # standard scale: normalize each gene to range from 0 to 1
    save = "markers.png"
)