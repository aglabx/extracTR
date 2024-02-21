# extracTR

## Introduction

This is a tool

## Installation

```bash
pip install extracTR
```

## Usage

Remember to remove adapters from the reads before running extracTR. Otherwise, some of predicted repeats may be contaminated by adapter sequences due to the fact that the adapter sequences are highly repetitive.

```bash
extracTR -1 reads.trimmed_1.fastq -2 reads.trimmed_2.fastq -o output_prefix
```