# scMultiAnnotator

## Description

scMultiAnnotator helps performing single-cell annotation according to its [gold standards](). This Nextflow-orchestrated pipeline performs three annotation algorithms with different methodologies between them: 

1. Manual annotation: Scanpy plots to help the user manually annotate cell clusters using custom cell markers;

2. Scoring-based method: Python implementation of ScType, which uses its own cell markers database;

3. Probabilistic model: Annotate cell-types without the need of reference single-cell data by using custom cell markers.

The user can provide either `.fastq` files (in this case, `cellranger count` is performed) or a count matrix as input. 

## Table of contents

## Intallation

1. Nextflow

    First, install Nextflow according to its [documentation](https://www.nextflow.io/docs/latest/install.html):

    ```
    # If `java -version` < 11, reinstall Java

    curl -s https://get.sdkman.io | bash
    sdk install java 17.0.10-tem
    java -version
    ```

    ```
    # Install Nextflow

    curl -s https://get.nextflow.io | bash
    chmod +x nextflow
    sudo mv nextflow /usr/local/bin
    nextflow info
    ```

2. Docker

    If Docker is not installed yet, follow the steps on [this tutorial]().

## Usage

```
nextflow run leandrotiburske/scMultiAnnotator \
                                                --fastq ${PWD}/samples.fastq \
                                                --reference_transcriptome ${PWD}/reference \
                                                --organism mouse \
                                                --tissue "Immune system" \
                                                --markers ${PWD}/markers.csv
```

## Example

## Future improvements