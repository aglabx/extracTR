#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 20.02.2024
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import argparse
import os
import aindex

from .core_functions.index_tools import compute_and_get_index

from .core_functions import tr_greedy_finder

def run_it(settings):

    fastq1 = settings.get("fastq1", None)
    fastq2 = settings.get("fastq2", None)
    threads = settings.get("threads", 32)
    coverage = settings.get("coverage", 1)
    lu = settings.get("lu", 100 * coverage)
    prefix = settings.get("output", "test")

    ### step 1. Compute aindex for reads
    kmer2tf, sdat = compute_and_get_index(fastq1, fastq2, prefix, threads, lu=lu)

    ### step 2. Find tandem repeats using circular path in de bruijn graph

    repeats = tr_greedy_finder(sdat, kmer2tf, max_depth=30_000, coverage=30, min_fraction_to_continue=30, k=23)

    all_predicted_trs = []
    for i, (status, second_status, next_rid, next_i, seq) in enumerate(repeats):
        if status == "tr":
            seq = seq[:-k]
            print(status, second_status, next_rid, next_i, len(seq), seq)
            all_predicted_trs.append(seq)
    print(len(all_predicted_trs))

    ### step 3. Save results to CSV

    output_file = f"{prefix}.csv"

    with open(output_file, "w") as fh:
        for i, seq in enumerate(all_predicted_trs):
            fh.write(f">{i}_{len(seq)}bp\n{seq}\n")

    ### step 4. Analyze repeat borders

    ### step 5. Enrich repeats variants
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract and analyze tandem repeats from raw DNA sequences.")
    parser.add_argument("-1", "--fastq1", help="Input file with DNA sequences in FASTQ format.")
    parser.add_argument("-2", "--fastq2", help="Input file with DNA sequences in FASTQ format.")
    parser.add_argument("-o", "output", help="Output file with tandem repeats in CSV format.")
    args = parser.parse_args()
    
    settings = {
        "fastq1": args.fastq1,
        "fastq2": args.fastq2,
        "output": args.output,
    }

    run_it(settings)