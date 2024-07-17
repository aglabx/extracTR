#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 20.02.2024
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import argparse

from extractr.core_functions.sdat_tools import load_sdat_as_dict
from .core_functions.index_tools import compute_and_get_index
import aindex
from .core_functions.sdat_tools import compute_abundace_anomaly
from .core_functions.tr_finder import naive_tr_finder, tr_greedy_finder
from collections import Counter
from trseeker.tools.sequence_tools import get_revcomp
from .core_functions.evaluation import compute_score
from trseeker.tools.sequence_tools import sc_iter_fasta_brute


def get_ref_aindex(fasta_file, ref_prefix, threads, lu=1):
    
    settings = {
        "index_prefix": "/media/eternus1/data/human/raw_reads/GCF_009914755.1_T2T-CHM13v2.0_genomic.L2.23",
        "aindex_prefix": "/media/eternus1/data/human/raw_reads/GCF_009914755.1_T2T-CHM13v2.0_genomic.L2.23",
        "reads_file": "/media/eternus1/data/human/raw_reads/GCF_009914755.1_T2T-CHM13v2.0_genomic.L2.reads",
        "max_tf": 10000000,
    }

    ref2tf = aindex.load_aindex(settings, skip_reads=False, skip_aindex=False)
    ref2tf.load_header("/media/eternus1/data/human/raw_reads/GCF_009914755.1_T2T-CHM13v2.0_genomic.L2.header")

    chrm2start = {}
    for start, name in ref2tf.headers.items():
        chrm2start[name.split()[0]] = start

    sdat_file = f"{ref_prefix}.sdat"
    sdat_file = "/media/eternus1/data/human/raw_reads/GCF_009914755.1_T2T-CHM13v2.0_genomic.L2.23"

    sdat_ref = load_sdat_as_dict(sdat_file, minimal_tf=lu)

    trf_file = "/media/eternus1/nfs/projects/databases/t2t/trf/GCF_009914755.1_T2T-CHM13v2.0_genomic.10kb.trf"
    ground_truth = {}
    trf_data = []
    trf_our_format = []
    with open(trf_file) as fh:
        for line in tqdm(fh):
            d = line.strip().split("\t")
            trf_data.append(d)
            array = d[14].upper()
            if len(array) < 100000:
                continue
            # for i in range(len(array)-k+1):
            #   kmer = array[i:i+k]
            #   ground_truth.setdefault(kmer, {})
            #   trid = len(trf_data) - 1
            #   ground_truth[kmer].setdefault(trid, 0)
            #   ground_truth[kmer][trid] += 1
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
    fasta_file = settings.get("fasta", None)
    ref_prefix = settings.get("rindex", "ref")

    kmer2tf, sdat = compute_and_get_index(fastq1, fastq2, prefix, threads, lu=lu)
    ref2tf, chrm2start, sdat_ref, trf_our_format = get_ref_aindex(fasta_file, ref_prefix, threads, lu=1)

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

    naive_tr_repeats, _, _, _ = naive_tr_finder(sdat, kmer2tf, min_tf_extension=3000, min_fraction_to_continue=30, k=23)
    repeats = tr_greedy_finder(sdat, kmer2tf, max_depth=30_000, coverage=30, min_fraction_to_continue=30, k=23)
    print(Counter([x[0] for x in repeats]))

    all_predicted_trs_v2 = []
    for i, (status, second_status, next_rid, next_i, seq) in enumerate(repeats):
        if status == "tr":
            seq = seq[:-k]
            print(status, second_status, next_rid, next_i, len(seq), seq)
            all_predicted_trs_v2.append(seq)
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
            print(trii, f"Size: {size}", x)
            # all_predicted_trs.append((size, x[3]))
            all_predicted_trs.append(x[3])
            # all_predicted_trs.append(get_revcomp(x[3]))
    print(len(all_predicted_trs))

    evaluation1, missed_repeats_fp1, missed_repeats_fn1 = compute_score(all_predicted_trs_v2, trf_our_format, chrm2start, ref2tf, delta=30_000, min_array_length=100, min_fish_strength=100, locus_length_cutoff=10_000, k=23)
    evaluation2, missed_repeats_fp2, missed_repeats_fn2 = compute_score(all_predicted_trs, trf_our_format, chrm2start, ref2tf, delta=30_000, min_array_length=100, min_fish_strength=100, locus_length_cutoff=10_000, k=23)

    srf_file = "/mnt/data/human/raw_reads/srf.fa"

    srf_reps = []
    for header, seq in sc_iter_fasta_brute(srf_file):
        if isinstance(seq, list):
            seq = seq[0]
        if len(seq) > 15000:
            continue
        srf_reps.append(seq)
    print(len(srf_reps))

    evaluation_srf, missed_repeats_fp_srf, missed_repeats_fn_srf = compute_score(srf_reps, trf_our_format, chrm2start, ref2tf, delta=30_000, min_array_length=100, min_fish_strength=100, locus_length_cutoff=10_000, k=23)

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
        for field in evaluation_fields:
            fh.write(f"{evaluation_srf[field]},")
        fh.write("\n")


    found_repeats_fasta = f"{output}.fasta"
    with open(found_repeats_fasta, "w") as fh:
        for i, seq in enumerate(repeats):
            fh.write(f">{i}_{len(seq)}bp\n{seq}\n")
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Benchmark on human T2T data.")
    parser.add_argument("--fastq1", help="Input left reads in FASTQ format.")
    parser.add_argument("--fastq2", help="Input right reads in  in FASTQ format.")
    parser.add_argument("-o", "--output", help="Output prefix for tandem repeats search results.")
    parser.add_argument("--rindex", help="Index prefix for reference.")
    parser.add_argument("--index", help="Index prefix for raw reads.")
    parser.add_argument("--fasta", help="Reference data in fasta format.")
    parser.add_argument("-t", "--threads", help="Number of threads to use.")
    parser.add_argument("-c", "--coverage", help="Coverage for aindex.")
    parser.add_argument("--lu", help="LU for aindex")
    parser.add_argument("--hl", help="Heng Li fasta")
    parser.add_argument("--tarean", help="TAREAM fasta")

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
        "tarean": args.tarean,
    }

    run_it(settings)