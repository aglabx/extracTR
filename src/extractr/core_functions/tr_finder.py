#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 20.02.2024
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

from collections import defaultdict
from trseeker.tools.sequence_tools import get_revcomp
from tqdm import tqdm

def naive_tr_finder(sdat, kmer2tf, min_tf_extension=3000, min_fraction_to_continue=30, k=23):
    ### TODO: save forks decisions

    repeats = []
    used_kmers = set()

    alphabet = ["A", "C", "T", "G"]
    rid = 0
    kmer2rid = {}
    kmer2repeat = defaultdict(list)
    frag2type = {}
    for (start_kmer, tf) in sdat:
        step = 0
        if start_kmer in used_kmers:
            continue
        used_kmers.add(start_kmer)
        used_kmers.add(get_revcomp(start_kmer))
        frag2type[start_kmer] = "SS"
        frag2type[get_revcomp(start_kmer)] = "SS"
        kmer2rid[start_kmer] = rid
        seq = start_kmer
        profile = [tf]
    while True:
        kmer = seq[-k:]
        ori_tf = kmer2tf[kmer]
        solutions = []
        prefix = kmer[-k+1:]
        for i, nucleotide in enumerate(alphabet):
            ctf = kmer2tf[prefix+nucleotide]
            if ctf < min_tf_extension:
                continue
            fr = 100.*ctf//(ori_tf+1)
            solutions.append((fr, nucleotide))

        solutions.sort()
    
        if solutions and solutions[-1][0] >= min_fraction_to_continue:
            seq += "".join(solutions[-1][1])
            profile.append(solutions)
            step += 1

            if seq[-k:] in used_kmers and seq[-k:] != start_kmer:
                repeated_kmer = seq[-k:]
                prev_type = frag2type[repeated_kmer]

                if len(prev_type) == 2:
                    prev_type = "F"+prev_type

                repeats.append((rid, "FRAG", len(seq), start_kmer, seq))
                for iii in range(len(seq)-k+1):
                    kmer_ = seq[iii:iii+k]
                    kmer2repeat[kmer_].append((rid, iii, prev_type))
                    kmer2repeat[get_revcomp(kmer_)].append((rid, iii, prev_type))
                    if not kmer_ in frag2type or frag2type[kmer_] == "FPP":
                        frag2type[kmer_] = prev_type
                        frag2type[get_revcomp(kmer_)] = prev_type
                rid += 1
                break
            used_kmers.add(seq[-k:])
            used_kmers.add(get_revcomp(seq[-k:]))
            frag2type[seq[-k:]] = "PP"
            frag2type[get_revcomp(seq[-k:])] = "PP"
            kmer2rid[seq[-k:]] = rid
            ### TODO: case ABAC => AB
            if seq[-k:] == start_kmer:
                repeats.append((rid, "TR", len(seq[:step]), seq[:step], profile))
                for iii in range(len(seq)-k+1):
                    kmer2repeat[seq[iii:iii+k]].append((rid, iii, "TR"))
                    kmer2repeat[get_revcomp(seq[iii:iii+k])].append((rid, iii, "TR"))
                    frag2type[seq[iii:iii+k]] = "TR"
                    frag2type[get_revcomp(seq[iii:iii+k])] = "TR"
                rid += 1
                break

            continue
        repeats.append((rid, "TE", len(seq), seq, solutions))
        for iii in range(len(seq)-k+1):
            kmer2repeat[seq[iii:iii+k]].append((rid, iii, "TE"))
            kmer2repeat[get_revcomp(seq[iii:iii+k])].append((rid, iii, "TE"))
            frag2type[seq[iii:iii+k]] = "TE"
            frag2type[get_revcomp(seq[iii:iii+k])] = "TE"
        rid += 1
        break

    return repeats, kmer2rid, kmer2repeat, frag2type

def tr_greedy_finder(sdat, kmer2tf, max_depth=30_000, coverage=30, min_fraction_to_continue=30, k=23):
    # ### Step 5b. Greedy find the most possible circles in the graph not working

    MIN_TF = coverage * 100

    repeats = []
    rid = 0

    alphabet = ["A", "C", "T", "G"]
    rid = 0
    cache = {}
    print("Expected iterations:", len([x for x in sdat if x[1] > MIN_TF]))
    for (start_kmer, tf) in tqdm(sdat):
        if tf < MIN_TF:
            break
        if start_kmer in cache:
            continue
        second_status = None
        status = None
        next_rid = None
        next_i = None
        length = k
        seq = [start_kmer]
        prefix = start_kmer[1:]

        cache[start_kmer] = (rid, 0, length)
        
        while True:
            solutions = []
            for i, nucleotide in enumerate(alphabet):
                ctf = kmer2tf[prefix+nucleotide]
                if ctf < MIN_TF:
                    continue
                if ctf:
                    solutions.append((ctf, nucleotide))  
            if not solutions:
                status = "zero"
                break
            
            solutions.sort()
            ctf, nucleotide = solutions[-1]
            kmer = prefix + nucleotide
            seq.append(nucleotide)
            if start_kmer == kmer:
                status = "tr"
                break

            prefix = kmer[1:]
            length += 1

            if kmer in cache:
                rid, strand, i = cache[kmer]
                if not rid in repeats:
                    second_status = "self"
                else:
                    second_status = repeats[rid][0]
                next_rid = rid
                next_i = i
                status = "frag"
                # seq.append(repeats[rid][-1][i:])
                break
            cache[kmer] = (rid, 0, length)
            if length == max_depth:
                status == "long"
                break
        repeats.append((status, second_status, next_rid, next_i, "".join(seq)))
        rid += 1

    return repeats