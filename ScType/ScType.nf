process ScType {

    // Set verbosity

    if (params.verbose == true) {
        debug true
    }

    else {
        debug false
    }

    // Name output directory
    publishDir "scType", mode: 'symlink'

    // Input: Directory with FASTQ files and a reference transcriptome
    input:
    path adata

    // Consider all files as output
    output:
    path "figures"
    path "sctype_scores.tsv"

    script:
    """
    sctype_python.py \
        --clustering $params.clustering \
        --tissue $params.tissue \
        --dataset $adata
    """

}

