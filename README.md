# extracTR

## Introduction

extracTR is a tool for identifying and analyzing tandem repeats (satellite DNA) in genomic sequences. It works with raw sequencing data (FASTQ) or assembled genomes (FASTA), using k-mer based approaches to detect repetitive patterns efficiently. extracTR can also design FISH probes for detected satellites and enrich monomer variant sequences directly from the de Bruijn graph.

## Features

- Efficient tandem repeat detection from raw sequencing data
- Support for single-end and paired-end FASTQ files
- Support for genome assemblies in FASTA format
- Support for precomputed aindex (`--aindex`)
- FISH probe design for detected satellite monomers
- Monomer variant enrichment via de Bruijn graph cycle search
- IUPAC degenerate consensus generation from variants
- Customizable parameters for fine-tuning repeat detection
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

### Basic usage

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

For precomputed aindex:
```bash
extracTR --aindex /path/to/index_prefix -o output_prefix -c 30
```

### Advanced usage

Custom k-mer size and probe design parameters:
```bash
extracTR -1 reads_1.fastq -2 reads_2.fastq -o output_prefix -t 64 -c 30 -k 25 \
    --probe-length 45 --top-probes 5 --min-gc 0.40 --max-gc 0.60
```

Skip variant enrichment and/or probe design:
```bash
extracTR -1 reads.fastq -o output_prefix -c 30 --skip-variants --skip-probes
```

### Options

#### Input / output
| Option | Description |
|--------|-------------|
| `-1, --fastq1` | Input file with forward reads in FASTQ format |
| `-2, --fastq2` | Input file with reverse reads in FASTQ format (optional) |
| `-f, --fasta` | Input genome assembly in FASTA format |
| `--aindex` | Prefix for a precomputed aindex (skips index computation) |
| `-o, --output` | Prefix for output files (required) |

#### Indexing parameters
| Option | Default | Description |
|--------|---------|-------------|
| `-t, --threads` | 32 | Number of threads for index computation |
| `-c, --coverage` | — | Data coverage (required; set 1 for genome assembly) |
| `-k, --k` | 23 | K-mer size for aindex |
| `--lu` | 100 * coverage | Minimum k-mer frequency cutoff |

#### FISH probe design
| Option | Default | Description |
|--------|---------|-------------|
| `--probe-length` | 40 | Probe length in bp |
| `--top-probes` | 3 | Number of top probes to report per monomer |
| `--min-gc` | 0.35 | Minimum GC content for probes |
| `--max-gc` | 0.65 | Maximum GC content for probes |
| `--skip-probes` | false | Skip the FISH probe design step |

#### Variant enrichment
| Option | Default | Description |
|--------|---------|-------------|
| `--skip-variants` | false | Skip the variant enrichment step |

**Note:** You must provide either FASTQ file(s), a FASTA file, or a precomputed aindex as input.

## Pipeline

extracTR runs the following steps:

1. **Index computation** — Build or load a k-mer frequency index (aindex)
2. **Tandem repeat detection** — Bidirectional greedy walk in the de Bruijn graph to find circular paths (tandem repeats) and linear elements
3. **Save results** — Write detected monomers and dispersed elements to FASTA
4. **Analyze repeat borders** — (placeholder for future development)
5. **Variant enrichment** — For each monomer, search alternative cycles in the de Bruijn graph to find sequence variants. Generate IUPAC degenerate consensus from same-length variants
6. **FISH probe design** — Slide a window across each circular monomer, score candidates by k-mer frequency strength and specificity (CV-based), filter by GC% and Tm, remove overlapping probes

## Output

extracTR generates the following output files:

| File | Description |
|------|-------------|
| `{prefix}.fa` | Predicted tandem repeat monomers (FASTA) |
| `{prefix}_te.fa` | Predicted dispersed / non-circular elements (FASTA) |
| `{prefix}_variants.fa` | Monomer sequence variants from de Bruijn graph |
| `{prefix}_consensus.fa` | IUPAC degenerate consensus sequences |
| `{prefix}_probes.fa` | FISH probe candidates (FASTA with metrics in headers) |
| `{prefix}_probes.tsv` | FISH probe candidates with full metrics (TSV) |

### Probe TSV columns

| Column | Description |
|--------|-------------|
| `probe_id` | Unique probe identifier |
| `source_monomer` | Monomer the probe was designed from |
| `position` | Start position within the monomer |
| `length` | Probe length in bp |
| `sequence` | Probe nucleotide sequence |
| `gc_content` | GC fraction (0-1) |
| `melting_temp` | Estimated melting temperature (°C) |
| `frequency_score` | Normalized mean k-mer frequency (signal strength) |
| `specificity_score` | CV-based specificity (1 = highly specific) |
| `composite_score` | frequency_score * specificity_score |

## Benchmarking

extracTR includes a benchmark module for evaluating repeat detection against TRF ground truth on reference genomes:

```bash
python -m extractr.benchmark \
    --fastq1 reads_1.fastq --fastq2 reads_2.fastq \
    -o benchmark_results -c 30 \
    --ref-index /path/to/genome.23 \
    --ref-header /path/to/genome.header \
    --ref-sdat /path/to/genome.23.sdat \
    --trf /path/to/genome.trf \
    --hl /path/to/srf.fa
```

## License

BSD License
