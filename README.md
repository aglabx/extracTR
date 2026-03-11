# extracTR

## Introduction

extracTR is a tool for identifying and analyzing tandem repeats (satellite DNA) in genomic sequences. It works with raw sequencing data (FASTQ) or assembled genomes (FASTA), using k-mer based approaches to detect repetitive patterns efficiently. extracTR can also design FISH probes for detected satellites and enrich monomer variant sequences directly from the de Bruijn graph.

For genome assemblies, extracTR integrates FasTAN to first identify tandem repeat arrays, then builds a focused k-mer index from those arrays for more accurate satellite detection.

## Features

- Efficient tandem repeat detection from raw sequencing data
- Support for single-end and paired-end FASTQ files (including gzipped)
- Support for genome assemblies in FASTA format with FasTAN pre-processing
- Support for precomputed aindex (`--aindex`)
- FISH probe design for detected satellite monomers
- Monomer variant enrichment via de Bruijn graph cycle search
- IUPAC degenerate consensus generation from variants
- Optional Rust backend for accelerated k-mer graph traversal
- Multi-threaded processing for improved performance

## Requirements

- Python 3.8 or later
- Jellyfish 2.x
- pigz (for gzipped FASTQ input)
- C compiler (gcc/clang) — for building FasTAN and tanbed from source

## Installation

We recommend installing extracTR in a separate Conda environment.

1. Create and activate a Conda environment:

```bash
conda create -n extractr python=3.9
conda activate extractr
```

2. Install Jellyfish and pigz:

```bash
conda install -c bioconda jellyfish
conda install pigz
```

3. Install extracTR:

```bash
pip install extracTR
```

4. Install external tools (FasTAN, tanbed) — required for genome assembly mode:

```bash
extracTR --install-tools
```

This clones and compiles [FasTAN](https://github.com/ad3002/FASTAN) and [tanbed](https://github.com/richarddurbin/alntools) from source. Requires `git`, `make`, and a C compiler.

## Usage

```bash
conda activate extractr
```

### From paired-end reads

```bash
extracTR -1 reads_1.fastq.gz -2 reads_2.fastq.gz -o output_prefix -c 30
```

### From single-end reads

```bash
extracTR -1 reads.fastq.gz -o output_prefix -c 30
```

### From genome assembly

```bash
extracTR -f genome.fasta -o output_prefix -c 1
```

This runs the FasTAN pipeline first (find TR arrays → build index from arrays only → detect satellites). To skip FasTAN and use the direct graph approach on the whole genome:

```bash
extracTR -f genome.fasta -o output_prefix -c 1 --no-fastan
```

### From precomputed aindex

```bash
extracTR --aindex /path/to/index_prefix -o output_prefix -c 30
```

### Advanced usage

Custom k-mer size, threads, and probe design parameters:
```bash
extracTR -1 reads_1.fastq.gz -2 reads_2.fastq.gz -o output_prefix -t 64 -c 30 -k 25 \
    --probe-length 45 --top-probes 5 --min-gc 0.40 --max-gc 0.60
```

Two-tier threshold for detecting satellites with variable k-mer frequencies (e.g. alpha satellite):
```bash
extracTR -1 reads_1.fastq.gz -2 reads_2.fastq.gz -o output_prefix -c 30 --ext-lu 50
```

Skip variant enrichment and/or probe design:
```bash
extracTR -1 reads.fastq -o output_prefix -c 30 --skip-variants --skip-probes
```

### Options

#### Input / output
| Option | Description |
|--------|-------------|
| `-1, --fastq1` | Input file with forward reads in FASTQ format (plain or .gz) |
| `-2, --fastq2` | Input file with reverse reads in FASTQ format (optional, plain or .gz) |
| `-f, --fasta` | Input genome assembly in FASTA format |
| `--aindex` | Prefix for a precomputed aindex (skips index computation) |
| `-o, --output` | Prefix for output files (required) |

#### Indexing parameters
| Option | Default | Description |
|--------|---------|-------------|
| `-t, --threads` | 32 | Number of threads for index computation |
| `-c, --coverage` | — | Data coverage (required; set 1 for genome assembly) |
| `-k, --k` | 23 | K-mer size for aindex |
| `--lu` | 100 * coverage | Minimum k-mer frequency cutoff for seed selection |
| `--ext-lu` | same as lu | Extension threshold (lower than lu to detect satellites with variable k-mer freq) |

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
| `--max-backtracks` | 1000 | Max DFS backtracks per seed (Rust backend only) |

#### Genome assembly mode
| Option | Default | Description |
|--------|---------|-------------|
| `--no-fastan` | false | Skip FasTAN pre-step, use direct graph approach on whole genome |

#### Tools management
| Option | Description |
|--------|-------------|
| `--install-tools` | Install external tools (FasTAN, tanbed) and exit |

## Pipeline

extracTR runs the following steps:

1. **FasTAN pre-processing** (genome mode only) — Run FasTAN to find tandem repeat arrays, extract array sequences, build k-mer index from arrays only
2. **Index computation** — Build or load a k-mer frequency index (jellyfish + aindex). For gzipped FASTQ, uses `pigz -dc` streaming via jellyfish generator
3. **Tandem repeat detection** — Bidirectional greedy walk in the de Bruijn graph to find circular paths (tandem repeats) and linear elements
4. **Save results** — Write detected monomers and dispersed elements to FASTA
5. **Variant enrichment** — For each monomer, search alternative cycles in the de Bruijn graph to find sequence variants. Generate IUPAC degenerate consensus from same-length variants
6. **FISH probe design** — Slide a window across each circular monomer, score candidates by k-mer frequency strength and specificity (CV-based), filter by GC% and Tm, remove overlapping probes

## Output

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
extracTR-benchmark \
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
