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

#
# Parse arguments
#

parser = argparse.ArgumentParser()

parser.add_argument('--counts',
                    type=str,
                    help="Path to CellRanger's output.")

parser.add_argument('--resolution',
                    type=float,
                    help="Resolution for clustering cells. Default: 0.8")

parser.add_argument('--min_genes',
                    type=int,
                    help="Minimum number of genes each cell must have. Default: 200")

parser.add_argument('--min_cells',
                    type=int,
                    help="Minimum number of cells that each gene must be present. Default: 3")

parser.add_argument('--organism',
                    type=str,
                    help="Organism that is being analysed: human | mouse. Default: human")

parser.add_argument('--clustering',
                    type=str,
                    help="Clustering method to be used: 'leiden' or 'louvain'. Default: 'leiden'")

parser.add_argument('--tissue',
                    type=str,
                    help="One of ScType's tissues.",
                    nargs="*")

args = parser.parse_args()

tissue = ' '.join(args.tissue)

#
# Quality control of single-cell data
#

# Read dataset
adata=sc.read_10x_mtx(args.counts)

# Calculate percentages of mitochondrial genes: "MT-" for human, "Mt-" for mouse
if args.organism == "human":
    # Mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # Ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    # Hemoglobin genes
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

elif args.organism == "mouse":
    # Mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("Mt-")
    # Ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
    # Hemoglobin genes
    adata.var["hb"] = adata.var_names.str.contains("^Hb[^(p)]")

else:
    logging.ERROR("Organism name is not valid. Provide 'human' or 'mouse'.")
    exit(1)

sc.pp.calculate_qc_metrics(adata,
                           qc_vars=["mt", "ribo", "hb"],
                           inplace=True,
                           log1p=True)

# Plot quality control metrics
sc.pl.violin(adata,
             ["n_genes_by_counts",
               "total_counts", 
               "pct_counts_mt"],
               jitter=0.4,
               multi_panel=True,
               save="_qc.png")

# Filter cells according to provided thresholds (minimim number of genes per cell and cell per genes)
sc.pp.filter_cells(adata, 
                   min_genes = args.min_genes)

sc.pp.filter_genes(adata, 
                   min_cells = args.min_cells)

adata.layers["counts"] = adata.X.copy()

#
# Process dataset
#

# Normalise
sc.pp.normalize_total(adata,
                      target_sum=1e4)

# Apply X = log(X + 1), where X is each value in the expression matrix
sc.pp.log1p(adata)

# Find highly variable genes
sc.pp.highly_variable_genes(adata, 
                            n_top_genes=2000, 
                            flavor="seurat_v3",
                            layer="counts")
adata.raw = adata

# Scale and run PCA
sc.pp.scale(adata,
            max_value=10)

scaled_data = pd.DataFrame(adata.X)
# Change column indexes
scaled_data.columns =adata.var_names
# Change the row indexes
scaled_data.index = adata.obs_names
scaled_data=scaled_data.T

# Perform dimensionality reduction by Principal Component Analysis (PCA)
sc.tl.pca(adata,
          zero_center=False)

# Find the 20 closest cells to each cell in the PCA space
sc.pp.neighbors(adata, 
                n_neighbors=20,
                n_pcs=10,
                use_rep="X_pca")

# Find clusters based on neighbors either by Louvain or Leiden (default)
if args.clustering == "leiden":
    sc.tl.leiden(adata, 
                 resolution = args.resolution)
    
elif args.clustering == "louvain":
    sc.tl.louvain(adata, 
                 resolution = args.resolution)

# Perform Uniform Manifold Approximation and Projection (UMAP) to better visualise scRNAseq data
sc.tl.umap(adata,min_dist=0.3)

# Save UMAP
sc.pl.umap(adata, color=[args.clustering], save="_clusters.png")

# Export dataset
adata.write("adata.h5ad")

#
# Find cluster markers
#
sc.tl.rank_genes_groups(adata,
                        groupby = args.clustering,
                        method = "wilcoxon",
                        key_added = "Markers")

# Filter relevant genes
sc.tl.filter_rank_genes_groups(adata,
                               min_in_group_fraction=0.2,
                               max_out_group_fraction=0.2,
                               key="Markers",
                               key_added="Markers_filtered",)

# Create dotplot with marker genes
sc.pl.rank_genes_groups_dotplot(adata,
                                groupby = args.clustering,
                                standard_scale = "var",
                                n_genes = 5,
                                key = "Markers_filtered",
                                save = "markers.png")

# Get dataframe with all marker genes
markers_df = sc.get.rank_genes_groups_df(adata,
                                         key = "Markers_filtered",
                                         group = None,)

# Write dataframe
markers_df.to_csv(path_or_buf="clusters_markers.tsv",
                  sep="\t")


import urllib.request

# Fetch the script from the URL
url = "https://raw.githubusercontent.com/kris-nader/sc-type-py/main/sctype_py.py"
response = urllib.request.urlopen(url)
script = response.read().decode()

# Execute the script
exec(script)


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