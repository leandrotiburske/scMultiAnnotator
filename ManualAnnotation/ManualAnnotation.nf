process ManualAnnotation {

    // Set verbosity

    if (params.verbose == true) {
        debug true
    }

    else {
        debug false
    }

    // Name output directory
    publishDir "ManualAnnotation", mode: 'symlink'

    // Input: .h5ad file from ProcessCounts's output
    input:
    path dataset

    // Use all produced files as output
    output:
    path "*"

    script:
    """
    manual_annotation.py \
        --dataset $dataset \
        --markers $params.markers \
        --clustering $params.clustering
    """

}