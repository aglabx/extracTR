#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 20.02.2024
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import argparse

from .core_functions.sdat_tools import load_sdat_as_dict
from .core_functions.index_tools import compute_and_get_index
from .extractr import check_dependencies
import aindex
from .core_functions.sdat_tools import compute_abundace_anomaly
from .core_functions.tr_finder import naive_tr_finder, tr_greedy_finder
from collections import Counter
from .core_functions.helpers import get_revcomp
from .core_functions.evaluation import compute_score
from .core_functions.helpers import sc_iter_fasta_brute
from tqdm import tqdm


def get_ref_aindex(ref_index_prefix, header_file, sdat_file, trf_file, lu=1):
    """Load reference aindex and ground truth TRF data.

    Args:
        ref_index_prefix: Prefix for reference aindex files (e.g. /path/to/genome.23)
        header_file: Path to .header file for reference
        sdat_file: Path to reference .sdat file
        trf_file: Path to TRF ground truth file
        lu: Minimum k-mer frequency cutoff
    """
    settings = {
        "index_prefix": ref_index_prefix,
        "aindex_prefix": ref_index_prefix,
        "reads_file": ref_index_prefix.rsplit('.', 1)[0] + ".reads" if '.' in ref_index_prefix else ref_index_prefix + ".reads",
        "max_tf": 10000000,
    }

    ref2tf = aindex.load_aindex(settings, skip_reads=False, skip_aindex=False)
    ref2tf.load_header(header_file)

    chrm2start = {}
    for start, name in ref2tf.headers.items():
        chrm2start[name.split()[0]] = start

    sdat_ref = load_sdat_as_dict(sdat_file, minimal_tf=lu)

    trf_data = []
    trf_our_format = []
    with open(trf_file) as fh:
        for line in tqdm(fh, desc="Loading TRF data"):
            d = line.strip().split("\t")
            trf_data.append(d)
            array = d[14].upper()
            if len(array) < 100000:
                continue
            trf_our_format.append((d[18].split()[0], int(d[6]), int(d[7]), d[13].upper()))

    return ref2tf, chrm2start, sdat_ref, trf_our_format


def run_it(settings):

    ### step 1. Compute aindex for reads

    fastq1 = settings["fastq1"]
    fastq2 = settings["fastq2"]
    output = settings["output"]
    threads = settings.get("threads", 32)
    coverage = settings.get("coverage", 1)
    lu = settings.get("lu", 100 * coverage)
    prefix = settings.get("index", "raw")
    ref_index_prefix = settings["ref_index"]
    header_file = settings["ref_header"]
    sdat_file = settings["ref_sdat"]
    trf_file = settings["trf"]
    debug = settings.get("debug", False)
    k = 23

    check_dependencies(need_index_tools=True)

    kmer2tf, sdat = compute_and_get_index(fastq1, fastq2, prefix, threads, lu=lu, debug=debug)
    ref2tf, chrm2start, sdat_ref, trf_our_format = get_ref_aindex(
        ref_index_prefix, header_file, sdat_file, trf_file, lu=1
    )

    ### Step 4a. Find kmers underrepresented in the assembly
    ### Step 4b. Find kmers overrepresented in the assembly

    all_rep, kmer2abandacy_diff = compute_abundace_anomaly(sdat, sdat_ref, coverage, ref_coverage=1)

    absent_in_ref = [x for x in all_rep if x[3] == 0]
    absent_in_raw = [x for x in all_rep if x[2] == 0]

    result_unerrepresented_file = f"{output}_underrepresented.csv"
    with open(result_unerrepresented_file, "w") as fh:
        for rep in absent_in_ref:
            fh.write(f"{'\t'.join(map(str,rep))}\n")

    result_overrepresented_file = f"{output}_overrepresented.csv"
    with open(result_overrepresented_file, "w") as fh:
        for rep in absent_in_raw:
            fh.write(f"{'\t'.join(map(str,rep))}\n")

    ### step 2. Find tandem repeats using circular path in de bruijn graph

    naive_tr_repeats, _, _, _ = naive_tr_finder(sdat, kmer2tf, min_tf_extension=3000, min_fraction_to_continue=30, k=k)
    repeats = tr_greedy_finder(sdat, kmer2tf, max_depth=30_000, coverage=30, min_fraction_to_continue=30, k=k)
    if debug:
        print(Counter([x[0] for x in repeats]))

    all_predicted_trs_v2 = []
    for i, (status, second_status, next_rid, next_i, seq) in enumerate(repeats):
        if status == "tr":
            seq = seq[:-k]
            if debug:
                print(status, second_status, next_rid, next_i, len(seq), seq)
            all_predicted_trs_v2.append(seq)
    if debug:
        print("Len of all_predicted_trs_v2", len(all_predicted_trs_v2))

    ### 6. Get dataset for evaluation

    MIN_PRED_SIZE = 30000

    trii = 0
    all_predicted_all = []
    for ii, x in enumerate(repeats):
        trii += 1
        all_predicted_all.append(x[3])
        all_predicted_all.append(get_revcomp(x[3]))

    trii = 0
    all_predicted_trs = []
    for ii, x in enumerate(repeats):
        if x[1] == "TR":
            trii += 1
            monomer = x[3]
            if len(monomer) < 5:
                continue
            kmer = (x[3] * k)[:k]
            size = kmer2tf[kmer]//30 * len(x[3])
            if size < MIN_PRED_SIZE:
                continue
            if debug:
                print(trii, f"Size: {size}", x)
            all_predicted_trs.append(x[3])
    if debug:
        print(len(all_predicted_trs))

    evaluation1, missed_repeats_fp1, missed_repeats_fn1 = compute_score(all_predicted_trs_v2, trf_our_format, chrm2start, ref2tf, delta=30_000, min_array_length=100, min_fish_strength=100, locus_length_cutoff=10_000, k=k, debug=debug)
    evaluation2, missed_repeats_fp2, missed_repeats_fn2 = compute_score(all_predicted_trs, trf_our_format, chrm2start, ref2tf, delta=30_000, min_array_length=100, min_fish_strength=100, locus_length_cutoff=10_000, k=k, debug=debug)

    srf_file = settings.get("hl")
    if srf_file:
        srf_reps = []
        for header, seq in sc_iter_fasta_brute(srf_file):
            if isinstance(seq, list):
                seq = seq[0]
            if len(seq) > 15000:
                continue
            srf_reps.append(seq)
        if debug:
            print(len(srf_reps))

        evaluation_srf, missed_repeats_fp_srf, missed_repeats_fn_srf = compute_score(srf_reps, trf_our_format, chrm2start, ref2tf, delta=30_000, min_array_length=100, min_fish_strength=100, locus_length_cutoff=10_000, k=k, debug=debug)
    else:
        evaluation_srf = None

    ### step 3. Save evaluation results to CSV

    output_file = f"{output}.csv"

    evaluation_fields = ["FP", "FN", "TP", "TN", "Accuracy", "P", "R", "F1"]

    with open(output_file, "w") as fh:
        for field in evaluation_fields:
            fh.write(f"{field},")
        fh.write("\n")
        for field in evaluation_fields:
            fh.write(f"{evaluation1[field]},")
        fh.write("\n")
        for field in evaluation_fields:
            fh.write(f"{evaluation2[field]},")
        fh.write("\n")
        if evaluation_srf:
            for field in evaluation_fields:
                fh.write(f"{evaluation_srf[field]},")
            fh.write("\n")

    found_repeats_fasta = f"{output}.fasta"
    with open(found_repeats_fasta, "w") as fh:
        for i, seq in enumerate(repeats):
            fh.write(f">{i}_{len(seq)}bp\n{seq}\n")


def main():
    parser = argparse.ArgumentParser(description="Benchmark on human T2T data.")
    parser.add_argument("--fastq1", help="Input left reads in FASTQ format.", required=True)
    parser.add_argument("--fastq2", help="Input right reads in FASTQ format.")
    parser.add_argument("-o", "--output", help="Output prefix for benchmark results.", required=True)
    parser.add_argument("--rindex", help="(deprecated) Use --ref-index instead.")
    parser.add_argument("--index", help="Index prefix for raw reads.")
    parser.add_argument("--fasta", help="Reference data in fasta format.")
    parser.add_argument("-t", "--threads", help="Number of threads to use.", default=32, type=int)
    parser.add_argument("-c", "--coverage", help="Coverage for aindex.", default=30, type=float)
    parser.add_argument("--lu", help="LU for aindex", default=None, type=int)
    parser.add_argument("--hl", help="Heng Li / SRF fasta for comparison")
    parser.add_argument("--ref-index", help="Prefix for reference aindex (e.g. /path/to/genome.23)", required=True)
    parser.add_argument("--ref-header", help="Path to reference .header file", required=True)
    parser.add_argument("--ref-sdat", help="Path to reference .sdat file", required=True)
    parser.add_argument("--trf", help="Path to TRF ground truth file", required=True)
    parser.add_argument("--debug", help="Show verbose diagnostic output.", action="store_true", default=False)

    args = parser.parse_args()

    settings = {
        "fastq1": args.fastq1,
        "fastq2": args.fastq2,
        "output": args.output,
        "threads": args.threads,
        "coverage": args.coverage,
        "lu": args.lu,
        "rindex": args.rindex,
        "index": args.index,
        "fasta": args.fasta,
        "hl": args.hl,
        "ref_index": args.ref_index,
        "ref_header": args.ref_header,
        "ref_sdat": args.ref_sdat,
        "trf": args.trf,
        "debug": args.debug,
    }

    run_it(settings)


if __name__ == "__main__":
    main()
