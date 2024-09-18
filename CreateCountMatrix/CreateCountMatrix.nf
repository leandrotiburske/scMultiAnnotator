process CreateCountMatrix {

    // Set verbosity

    if (params.verbose == true) {
        debug true
    }

    else {
        debug false
    }

    // Name output directory
    publishDir "createCountMatrix", mode: 'symlink'

    // Input: Directory with FASTQ files and a reference transcriptome
    input:
    path fastq
    val reference_transcriptome

    // Consider all files as output
    output:
    path "*"

    script:

    // If no reference transcriptome is provided, use the human transcriptome
    if (reference_transcriptome == "none") {
        """
        echo "Downloading human reference transcriptome..."
        
        wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
        tar -zxvf refdata-gex-GRCh38-2020-A.tar.gz

        cellranger count --id=$params.experiment_name \
            --create-bam false \
            --fastqs=$fastq \
            --transcriptome=refdata-gex-GRCh38-2020-A
    """

    // Run CellRanger when reference transcriptome is provided
    } else {
        """
        cellranger count --id=$params.experiment_name \
            --create-bam false \
            --fastqs=$fastq \
            --transcriptome=$reference_transcriptome
        """
    }

}