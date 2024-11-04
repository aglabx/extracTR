#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 20.02.2024
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import argparse
from .core_functions.index_tools import compute_and_get_index, compute_and_get_index_for_fasta
from .core_functions.tr_finder import tr_greedy_finder_bidirectional

def run_it():

    parser = argparse.ArgumentParser(description="Extract and analyze tandem repeats from raw DNA sequences.")
    parser.add_argument("-1", "--fastq1", help="Input file with DNA sequences in FASTQ format.", default=None)
    parser.add_argument("-2", "--fastq2", help="Input file with DNA sequences in FASTQ format.", default=None)
    parser.add_argument("-f", "--fasta", help="Input genome fasta file", required=False, default=None)
    parser.add_argument("-o", "--output", help="Output file with tandem repeats in CSV format.", required=True)
    parser.add_argument("-t", "--threads", help="Number of threads to use.", default=32, type=int, required=False)
    parser.add_argument("-c", "--coverage", help="Data coverage, set 1 for genome assembly", type=int, required=True)
    parser.add_argument("--lu", help="Minimal repeat kmers coverage [100 * coverage].", default=None, type=int, required=False)
    parser.add_argument("-k", "--k", help="K-mer size to use for aindex.", default=23, type=int, required=False)
    args = parser.parse_args()
    
    settings = {
        "fastq1": args.fastq1,
        "fastq2": args.fastq2,
        "fasta": args.fasta,
        "output": args.output,
        "threads": 32,
        "coverage": 1,
        "lu": args.lu,
        "k": 23,
        "min_fraction_to_continue": 30,
    }

    fastq1 = settings.get("fastq1", None)
    fastq2 = settings.get("fastq2", None)
    fasta = settings.get("fasta", None)
    threads = settings.get("threads", 32)
    coverage = settings.get("coverage", 1)
    lu = settings.get("lu", 100 * coverage)
    prefix = settings.get("output", "test")
    min_fraction_to_continue = settings.get("min_fraction_to_continue", 30)
    k = settings.get("k", 23)

    ### step 1. Compute aindex for reads
    if fastq1 and fastq2:
        kmer2tf, sdat = compute_and_get_index(fastq1, fastq2, prefix, threads, lu=lu)
    elif fasta:
        kmer2tf, sdat = compute_and_get_index_for_fasta(fasta, prefix, threads, lu=lu)
    else:
        raise Exception("No input data")

    ### step 2. Find tandem repeats using circular path in de bruijn graph

    repeats = tr_greedy_finder_bidirectional(sdat, kmer2tf, max_depth=30_000, coverage=coverage, min_fraction_to_continue=min_fraction_to_continue, k=k, lu=lu)

    all_predicted_trs = []
    for i, (status, second_status, next_rid, next_i, seq) in enumerate(repeats):
        if status == "tr":
            seq = seq[:-k]
            # print(status, second_status, next_rid, next_i, len(seq), seq)
            all_predicted_trs.append(seq)
        elif status == "frag":
            pass
        elif status == "zero":
            pass
        elif status == "long":
            pass
        else:
            # print(status, second_status, next_rid, next_i, len(seq), seq)
            raise Exception("Unknown status")
        
    print(f"Predicted {len(all_predicted_trs)} tandem repeats.")

    ### step 3. Save results to CSV

    output_file = f"{prefix}.csv"

    with open(output_file, "w") as fh:
        for i, seq in enumerate(all_predicted_trs):
            fh.write(f">{i}_{len(seq)}bp\n{seq}\n")

    ### step 4. Analyze repeat borders

    ### step 5. Enrich repeats variants

if __name__ == "__main__":
    run_it()