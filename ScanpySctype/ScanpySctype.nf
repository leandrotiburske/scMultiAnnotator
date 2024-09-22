process ScanpySctype {

    // Set verbosity

    if (params.verbose == true) {
        debug true
    }

    else {
        debug false
    }

    // Name output directory
    publishDir "scanpySctype", mode: 'symlink'

    // Input: Directory with FASTQ files and a reference transcriptome
    input:
    path counts

    // Consider all files as output
    output:
    path "figures"
    path "clusters_markers.tsv"
    path "sctype_scores.tsv"
    path "adata.h5ad", emit: dataset

    script:
    """
    python3.9 $projectDir/ScanpySctype/scanpy_processing.py \
        --counts /home/leandro/Downloads/filtered_gene_bc_matrices/hg19 \
        --resolution $params.resolution \
        --min_genes $params.min_genes \
        --min_cells $params.min_cells \
        --organism $params.organism \
        --clustering $params.clustering \
        --tissue $params.tissue
    """

}

