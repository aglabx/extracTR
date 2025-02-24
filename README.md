# extracTR

## Introduction

extracTR is a tool for identifying and analyzing tandem repeats in genomic sequences. It works with raw sequencing data (FASTQ) or assembled genomes (FASTA), using k-mer based approaches to detect repetitive patterns efficiently.

## Features

- Efficient tandem repeat detection from raw sequencing data
- Support for single-end and paired-end FASTQ files
- Support for genome assemblies in FASTA format
- Customizable parameters for fine-tuning repeat detection
- Output in easy-to-analyze CSV format
- Multi-threaded processing for improved performance

## Requirements

- Python 3.7 or later
- Jellyfish 2.3.0 or later
- Conda (for easy environment management)

## Installation

We recommend installing extracTR in a separate Conda environment to manage dependencies effectively.

1. Create a new Conda environment:

```bash
conda create -n extractr_env python=3.9
```

2. Activate the environment:

```bash
conda activate extractr_env
```

3. Install Jellyfish:

```bash
conda install -c bioconda jellyfish
```

4. Install extracTR using pip:

```bash
pip install extracTR
```

To deactivate the environment when you're done:

```bash
conda deactivate
```

## Usage

Before running extracTR, ensure that you have removed adapters from your sequencing reads and activated the Conda environment:

```bash
conda activate extractr_env
```

Basic usage:

For paired-end FASTQ files:
```bash
extracTR -1 reads_1.fastq -2 reads_2.fastq -o output_prefix -c 30
```

For single-end FASTQ file:
```bash
extracTR -1 reads.fastq -o output_prefix -c 30
```

For genome assembly in FASTA format:
```bash
extracTR -f genome.fasta -o output_prefix -c 1
```

Advanced usage with custom parameters:

```bash
extracTR -1 reads_1.fastq -2 reads_2.fastq -o output_prefix -t 64 -c 30 -k 25
```

Options:
- `-1, --fastq1`: Input file with forward DNA sequences in FASTQ format
- `-2, --fastq2`: Input file with reverse DNA sequences in FASTQ format (optional for paired-end data)
- `-f, --fasta`: Input genome assembly in FASTA format
- `-o, --output`: Prefix for output files
- `-t, --threads`: Number of threads to use (default: 32)
- `-c, --coverage`: Coverage to use for indexing (required)
- `-k, --k`: K-mer size to use for indexing (default: 23)
- `--lu`: Coverage cutoff for k-mers (default: 100 * coverage)

Note: You must provide either FASTQ file(s) or a FASTA file as input.

## Output

extracTR generates the following output files:
- `{output_prefix}.csv`: Main output file containing detected tandem repeats
- `{output_prefix}.sdat`: Intermediate file with k-mer frequency data
- Additional files for detailed analysis and debugging
