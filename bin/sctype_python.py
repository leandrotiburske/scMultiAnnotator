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

parser.add_argument('--clustering',
                    type=str,
                    help="Clustering method to be used: 'leiden' or 'louvain'. Default: 'leiden'")

parser.add_argument('--tissue',
                    type=str,
                    help="One of ScType's tissues.",
                    nargs="*")

args = parser.parse_args()

tissue = ' '.join(args.tissue)

# Fetch the script from the URL
url = "https://raw.githubusercontent.com/kris-nader/sc-type-py/main/sctype_py.py"
response = urllib.request.urlopen(url)
script = response.read().decode()

# Execute the script
exec(script)

adata = sc.read_h5ad(args.dataset)

adata.raw = adata

# Scale dataset
sc.pp.scale(adata,
            max_value=10)

scaled_data = pd.DataFrame(adata.X)

# Change column indexes
scaled_data.columns =adata.var_names

# Change the row indexes
scaled_data.index = adata.obs_names
scaled_data=scaled_data.T

scRNAseqData = scaled_data
gs_list = gene_sets_prepare(path_to_db_file="https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx",
                            cell_type = tissue)

es_max = sctype_score(scRNAseqData = scRNAseqData,
                      scaled = True,
                      gs = gs_list['gs_positive'],
                      gs2 = gs_list['gs_negative'])

unique_clusters = adata.obs[args.clustering].unique()

# Apply the function to each unique cluster and combine the results into a DataFrame
cL_results = pd.concat([process_cluster(cluster,
                                        adata,
                                        es_max,
                                        args.clustering) for cluster in unique_clusters])

# Group by cluster and select the top row based on scores
sctype_scores = cL_results.groupby('cluster').apply(lambda x: x.nlargest(1, 'scores')).reset_index(drop=True)

# Set low-confidence clusters to "Unknown"
sctype_scores.loc[sctype_scores['scores'] < sctype_scores['ncells'] / 4, 'type'] = 'Unknown'

# Iterate over unique clusters
adata.obs['sctype_classification'] = ""
for cluster in sctype_scores['cluster'].unique():
    # Filter sctype_scores for the current cluster
    cl_type = sctype_scores[sctype_scores['cluster'] == cluster]
    # Get the type for the current cluster
    cl_type_value = cl_type['type'].iloc[0]
    # Update 'sctype_classification' in pbmc.obs for cells belonging to the current cluster
    adata.obs.loc[adata.obs[args.clustering] == cluster, 'sctype_classification'] = cl_type_value

# Plot the UMAP with sctype_classification as labels
sc.pl.umap(adata, color='sctype_classification',
           title='UMAP with sctype_classification',
           save = "sctype.png")

sctype_scores.to_csv(path_or_buf="sctype_scores.tsv",
                     sep="\t")