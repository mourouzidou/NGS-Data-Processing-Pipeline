## NGS Data Processing Pipeline for Sequence Retrieval, Trimming, and Analysis

This project includes a pipeline for handling Next-Generation Sequencing (NGS) data for :

* retrieving sequence data from public repositories
* trimming low-quality reads
* performing downstream analysis
  
The pipeline is automated through a series of shell scripts.

## Features

- **Sequence Retrieval**: Automatically download raw sequencing data from public databases like SRA (Sequence Read Archive).
- **Quality Trimming**: Perform quality control and trimming of reads to remove low-quality bases and adapters.
- **NGS Analysis**: Execute downstream analysis steps, such as alignment, variant calling, or assembly, depending on the project requirements.

### 1. Retrieve Sequence Data

Use the `getSRR_1.sh` script to download raw sequence data from SRA. 
The script will download the corresponding sequencing data in FASTQ format.

```bash
bash getSRR_1.sh <SRR_ID>
```
### 2. Quality Trimming

This script processes RNA-Seq data by performing read trimming, quality control, read mapping, and gene expression quantification.

``` bash

bash trimming_2.sh <input_fastq> <output_fastq>

``` 

### 3. NGS Analysis

This script processes genomic data from the 1000 Genomes Project by extracting, downloading, aligning, and analyzing sequence data.

``` bash

bash ngs.sh <list_of_identifiers> <reference_genome.fasta>
```
<list_of_identifiers>: A file containing a list of sample identifiers to process.
<reference_genome.fasta>: The reference genome file to use for alignment.








