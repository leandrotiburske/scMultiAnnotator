process ScVerse {

    // Set verbosity

    if (params.verbose == true) {
        debug true
    }

    else {
        debug false
    }

    //conda workflow.projectDir/'ScVerse/scverse.yaml'


    // Name output directory
    publishDir "scVerse", mode: 'symlink'

    // Input: Directory with FASTQ files and a reference transcriptome
    input:
    path adata

    // Consider all files as output
    output:
    path "*"

    script:
    """
    scverse_annotation.py \
        --adata $adata \
        --markers $params.markers
    """

}