#!/usr/bin/env nextflow

// Define workflow parameters

params.counts = null
params.fastq = null
params.reference_transcriptome = "none"
params.organism = "human"
params.tissue = null
params.markers = null
params.experiment_name = "cellranger_output"
params.resolution = 0.8
params.min_genes = 200
params.min_cells = 3
params.clustering = "leiden"
params.verbose = false
params.help = null

// Help message
if (params.help) {
    log.info """\

    scType Nextflow: 
    single-cell processing and annotation workflow

    ----------------------------------------------

    Usage:

        nextflow run

    Parameters:

        --counts: Path to the directory with the barcodes, genes and matrix. Default: null
        --organism: Either `human` or `mouse`. Default: human
        --tissue: One of ScType's available tissues. Default: null
        --markers: Path to a csv file with markers according to ScVerse's specifications. Default: null
        --resolution: Louvain or Leiden's clustering resolution. Default: 0.8
        --min_genes: Minimum number of genes for cell filtering. Default: 200
        --min_cells: Minimum number of cells for gene filtering. Default: 3
        --clustering: Clustering method. Either `leiden` or `louvain`. Default: leiden
        --verbose: true or false. Default: false
    """
    .stripIdent()

    exit(0)

}

// Print parameters
log.info """\

scType Nextflow: 
single-cell processing and annotation workflow

----------------------------------------------

"""


// Include processes

//include { CreateCountMatrix } from './CreateCountMatrix/CreateCountMatrix.nf'
include { ProcessCounts } from './ProcessCounts/ProcessCounts.nf'
include { ManualAnnotation } from './ManualAnnotation/ManualAnnotation.nf'
include { ScType } from './ScType/ScType.nf'
include { ScVerse } from './ScVerse/ScVerse.nf'

// Setup workflow

workflow {

    ProcessCounts(params.counts)
    ManualAnnotation(ProcessCounts.out.dataset)
    ScType(ProcessCounts.out.dataset)
    scverse_ch = ScVerse(ProcessCounts.out.dataset)
    //countMatrix_ch = CreateCountMatrix(params.fastq, params.reference_transcriptome)

}

// Print message when workflow is completed
workflow.onComplete {
    log.info ( workflow.success ? "\nDone!\n" : "Something went wrong" )
}