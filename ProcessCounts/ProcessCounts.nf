process ProcessCounts {

    // Set verbosity

    if (params.verbose == true) {
        debug true
    }

    else {
        debug false
    }

    // Name output directory
    publishDir "processCounts", mode: 'symlink'

    // Input: Directory with FASTQ files and a reference transcriptome
    input:
    path counts

    // Consider all files as output
    output:
    path "*"

    script:
    """
    python3.10 $projectDir/ProcessCounts/scanpy,py
    """

}

