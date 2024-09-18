import pandas as pd
import numpy as np
import anndata as an
import scanpy as sc
import urllib.request
import argparse
import logging

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

args = parser.parse_args()

# Fetch the script from the URL
url = "https://raw.githubusercontent.com/kris-nader/sc-type-py/main/sctype_py.py"
response = urllib.request.urlopen(url)
script = response.read().decode()

# Execute the script
exec(script)

adata=sc.read_10x_mtx(args.counts)

# mitochondrial genes, "MT-" for human, "Mt-" for mouse
if args.organism == "human":
    # mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

elif args.organism == "mouse":
    # mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("Mt-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
    # hemoglobin genes
    adata.var["hb"] = adata.var_names.str.contains("^Hb[^(p)]")

else:
    logging.ERROR("Organism name is not valid. Provide 'human' or 'mouse'.")
    exit(1)

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
)

sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
    save="_qc.png"
)

sc.pp.filter_cells(adata, 
                   min_genes = args.min_genes)

sc.pp.filter_genes(adata, 
                   min_cells = args.min_cells)


adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_total(adata,
                      target_sum=1e4)

sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, 
                            n_top_genes=2000, 
                            flavor="seurat_v3",
                            layer="counts")
adata.raw = adata

# Scale and run PCA
sc.pp.scale(adata,
            max_value=10)

scaled_data = pd.DataFrame(adata.X)
# change column indexes
scaled_data.columns =adata.var_names
# Change the row indexes
scaled_data.index = adata.obs_names
scaled_data=scaled_data.T

sc.tl.pca(adata,
          zero_center=False)

sc.pp.neighbors(adata, 
                n_neighbors=20,
                n_pcs=10,
                use_rep="X_pca")

if args.clustering == "leiden":
    sc.tl.leiden(adata, 
                 resolution = args.resolution)
    
elif args.clustering == "louvain":
    sc.tl.louvain(adata, 
                 resolution = args.resolution)

sc.tl.umap(adata,min_dist=0.3)

sc.pl.umap(adata, color=[args.clustering], save="_clusters.png")