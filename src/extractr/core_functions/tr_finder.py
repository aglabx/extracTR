#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 20.02.2024
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

from collections import defaultdict
from trseeker.tools.sequence_tools import get_revcomp
from tqdm import tqdm
from collections import deque


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


def tr_greedy_finder_bidirectional(sdat, kmer2tf, max_depth=30_000, coverage=30, min_fraction_to_continue=30, k=23):
    """
    Extends sequences in both left and right directions except for TRs.
    
    Parameters:
        sdat: List of tuples (kmer, tf) sorted by tf descending.
        kmer2tf: Dictionary mapping k-mer to its tf value.
        max_depth: Maximum length to extend.
        coverage: Coverage threshold.
        min_fraction_to_continue: Minimum fraction to continue extending.
        k: Length of k-mers.
    
    Returns:
        List of repeats with their status and sequences.
    """
    
    MIN_TF = coverage * 100
    
    repeats = []
    rid = 0
    
    alphabet = ["A", "C", "T", "G"]
    cache = {}
    print("Expected iterations:", len([x for x in sdat if x[1] > MIN_TF]))
    
    for (start_kmer, tf) in tqdm(sdat):
        if tf < MIN_TF:
            break
        if start_kmer in cache:
            continue
        
        # Initialize for bidirectional extension
        seq_deque = deque([start_kmer])
        status = None
        second_status = None
        next_rid = None
        next_i = None
        total_length = k
        
        # Mark the start k-mer in cache
        cache[start_kmer] = (rid, 0, total_length)
        
        # Right extension
        right_prefix = start_kmer[1:]
        right_seq = []
        right_length = 0
        right_status = None
        
        while right_length < max_depth:
            solutions = []
            for nucleotide in alphabet:
                new_kmer = right_prefix + nucleotide
                ctf = kmer2tf.get(new_kmer, 0)
                if ctf < MIN_TF:
                    continue
                if ctf:
                    solutions.append((ctf, nucleotide))
            if not solutions:
                right_status = "zero"
                break
            
            # Select the nucleotide with the highest coverage
            solutions.sort(reverse=True)
            ctf, nucleotide = solutions[0]
            new_kmer = right_prefix + nucleotide
            right_seq.append(nucleotide)
            total_length += 1
            
            if new_kmer == start_kmer:
                right_status = "tr"
                break
            
            if new_kmer in cache:
                existing_rid, strand, i = cache[new_kmer]
                if existing_rid not in repeats:
                    second_status = "self"
                else:
                    second_status = repeats[existing_rid][0]
                next_rid = existing_rid
                next_i = i
                right_status = "frag"
                break
            
            # Update for next iteration
            cache[new_kmer] = (rid, 0, total_length)
            right_prefix = new_kmer[1:]
            right_length += 1
        
        # If TR detected in right extension, skip left extension
        if right_status == "tr":
            final_status = "tr"
            full_seq = ''.join(right_seq)
        else:
            # Left extension
            left_suffix = start_kmer[:-1]
            left_seq = []
            left_length = 0
            left_status = None
            
            while left_length < max_depth:
                solutions = []
                for nucleotide in alphabet:
                    new_kmer = nucleotide + left_suffix
                    ctf = kmer2tf.get(new_kmer, 0)
                    if ctf < MIN_TF:
                        continue
                    if ctf:
                        solutions.append((ctf, nucleotide))
                if not solutions:
                    left_status = "zero"
                    break
                
                # Select the nucleotide with the highest coverage
                solutions.sort(reverse=True)
                ctf, nucleotide = solutions[0]
                new_kmer = nucleotide + left_suffix
                left_seq.append(nucleotide)
                total_length += 1
                
                if new_kmer == start_kmer:
                    left_status = "tr"
                    break
                
                if new_kmer in cache:
                    existing_rid, strand, i = cache[new_kmer]
                    if existing_rid not in repeats:
                        second_status = "self"
                    else:
                        second_status = repeats[existing_rid][0]
                    next_rid = existing_rid
                    next_i = i
                    left_status = "frag"
                    break
                
                # Update for next iteration
                cache[new_kmer] = (rid, 0, total_length)
                left_suffix = new_kmer[:-1]
                left_length += 1
            
            # Determine overall status
            if left_status == "tr" or right_status == "tr":
                final_status = "tr"
            elif left_status == "frag" or right_status == "frag":
                final_status = "frag"
            elif left_status == "zero" and right_status == "zero":
                final_status = "zero"
            else:
                final_status = "extended"
            
            # Combine left and right sequences
            full_seq = ''.join(reversed(left_seq)) + start_kmer + ''.join(right_seq)
        
        repeats.append((final_status, second_status, next_rid, next_i, full_seq))
        rid += 1
    
    return repeats
