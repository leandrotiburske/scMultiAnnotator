#!/usr/bin/env nextflow

// Define workflow parameters

params.teste_pbmc = "/home/leandro/Downloads/filtered_gene_bc_matrices/hg19/"
params.fastq = null
params.reference_transcriptome = "none"
params.organism = "human"
params.experiment_name = "cellranger_output"
params.resolution = 0.8
params.min_genes = 200
params.min_cells = 3
params.clustering = "leiden"
params.verbose = true
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

    Options:
    """
    .stripIdent()

    exit(0)

}

// Print parameters
log.info """\

scType Nextflow: 
single-cell processing and annotation workflow

----------------------------------------------

fastq: ${params.fastq}
verbose: ${params.verbose}

"""


// Include processes

include { CreateCountMatrix } from './CreateCountMatrix/CreateCountMatrix.nf'
include { ProcessCounts } from './ProcessCounts/ProcessCounts.nf'


// Setup workflow

workflow {

    processedCounts = ProcessCounts(params.teste_pbmc)
    countMatrix = CreateCountMatrix(params.fastq, params.reference_transcriptome)

}

// Print message when workflow is completed
workflow.onComplete {
    log.info ( workflow.success ? "\nDone!\n" : "Something went wrong" )
}