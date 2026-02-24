#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 20.02.2024
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

from intervaltree import IntervalTree
from tqdm import tqdm

def get_bed_coordinates(x, ref2tf, chrm2start):
  try:
    start_header = ref2tf.get_header(x.begin).split()[0]
    end_header = ref2tf.get_header(x.end).split()[0]
    start_pos = x.begin - chrm2start[start_header]
    end_pos = x.end - chrm2start[end_header]
  except:
    print("ERROR2:", x)
    return None, None, None
  if start_header != end_header:
    print("ERROR1:", start_header, start_pos, end_header, end_pos, rep)
    return None, None, None
  if end_pos < start_pos:
      start_pos, end_pos = end_pos, start_pos
  return start_header, start_pos, end_pos


def get_bed_data(all_loci, ref2tf, chrm2start, locus_length_cutoff=10000):

  data = []
  for rep in all_loci:
    for x in all_loci[rep]:
      start_header, start_pos, end_pos = get_bed_coordinates(x, ref2tf, chrm2start)
      if start_header is None:
        continue
      if abs(end_pos - start_pos) < locus_length_cutoff:
        continue
      data.append((start_header, start_pos, end_pos, rep))

  return data


def compute_stats(data, trf_our_format, name="our", debug=False):
  ref_IT = {}
  for chrmA, startA, endA, repA in trf_our_format:
    ref_IT.setdefault(chrmA, IntervalTree())
    ref_IT[chrmA].addi(startA, endA, repA)

  dat_IT = {}
  for chrmA, startA, endA, repA in data:
    dat_IT.setdefault(chrmA, IntervalTree())
    dat_IT[chrmA].addi(startA, endA, repA)

  false_positive = 0
  false_negative =0
  true_positive = 0
  true_negative = 0
  true_positive_our = 0

  missed_repeats_fp = []
  missed_repeats_fn = []

  for chrmA, startA, endA, repA in trf_our_format:
    if chrmA in dat_IT and dat_IT[chrmA][startA:endA]:
      true_positive += 1
    else:
      false_negative += 1
      missed_repeats_fn.append((chrmA, startA, endA, repA))


  for chrmA, startA, endA, repA in data:
    if chrmA in ref_IT and not ref_IT[chrmA][startA:endA]:
      false_positive += 1
      missed_repeats_fp.append((chrmA, startA, endA, repA))
    else:
      true_positive_our += 1

  p = true_positive/(true_positive+false_positive)
  r = true_positive/(true_positive+false_negative)

  if debug:
    print("FP", false_positive)
    print("FN", false_negative)
    print("TP", true_positive)
    print("TN", true_negative)
    print(f"Accuracy {name} to trf:", true_positive/len(trf_our_format))
    print(f"Accuracy trf to {name}:", true_positive_our/len(data))
    print("P", p)
    print("R", r)
    print("F1", 2 * p * r / (p+r))

  evaluation = {
    "FP": false_positive,
    "FN": false_negative,
    "TP": true_positive,
    "TN": true_negative,
    "Accuracy": true_positive/len(trf_our_format),
    "P": p,
    "R": r,
    "F1": 2 * p * r / (p+r)
  }

  return evaluation, missed_repeats_fp, missed_repeats_fn


### version 2 all kmers
### get all ref loci for repeats

def compute_loci(predicted_trs, ref2tf, delta=1000, k=23, debug=False):
  ### TODO: replace interval tree with chains
  all_loci = {}
  N = len(predicted_trs)
  for rep_id, repeat in enumerate(predicted_trs):
    computed_kmers = set()
    if len(repeat) < k:
      repeat = (repeat * k)[:2*k]
    pos_tree = IntervalTree()
    if debug:
      print(f"{rep_id}/{N}", repeat)
    for i in tqdm(range(len(repeat)-k+1), desc=f"Loci {rep_id}/{N}", disable=not debug):
      key_kmer = repeat[i:i+k]
      if key_kmer in computed_kmers:
        continue
      computed_kmers.add(key_kmer)
      ref_poses = ref2tf.pos(key_kmer)
      for pos in ref_poses:
        pos_tree.addi(pos-delta, pos+delta)
    if debug:
      print(len(pos_tree))
    pos_tree.merge_overlaps()
    if debug:
      print(len(pos_tree))
    all_loci[repeat] = pos_tree
  return all_loci




def compute_loci_chains(predicted_trs, ref2tf, delta=1000, min_array_length=10000, k=23, min_fish_strength=100, debug=False):
  ### TODO: replace interval tree with chains
  all_loci = {}
  N = len(predicted_trs)
  for rep_id, repeat in enumerate(predicted_trs):
    computed_kmers = set()
    if len(repeat) < k:
      repeat = (repeat * k)[:2*k]
    all_positions = []
    if debug:
      print(f"{rep_id}/{N}", repeat)
    for i in tqdm(range(len(repeat)-k+1), desc=f"Chains {rep_id}/{N}", disable=not debug):
      key_kmer = repeat[i:i+k]
      if key_kmer in computed_kmers:
        continue
      computed_kmers.add(key_kmer)
      all_positions += ref2tf.pos(key_kmer)
    all_positions.sort()

    chains = []
    i, j = 0, 1
    fish_strength = 0
    while i < len(all_positions) and j < len(all_positions):
      dist = all_positions[j] - all_positions[i]
      if dist < delta:
        j += 1
        fish_strength += 1
      else:
        if fish_strength > min_fish_strength:
          d = all_positions[j-1] - all_positions[i]
          if d >= min_array_length:
            chains.append((all_positions[i], all_positions[j-1]))
        i = j
        j = i + 1
        fish_strength = 0
    
    if fish_strength > min_fish_strength and i < len(all_positions):
      d = all_positions[j-1] - all_positions[i]
      if d >= min_array_length:
        chains.append((all_positions[i], all_positions[j-1]))
    ## todo: handle last case
    pos_tree = IntervalTree()
    for i, j in chains:
      pos_tree.addi(i, j)
    if debug:
      print(len(pos_tree))
    all_loci[repeat] = pos_tree
  return all_loci


def compute_score(monomers_dataset, trf_our_format, chrm2start, ref2tf, delta=30_000, min_array_length=100,min_fish_strength=100, locus_length_cutoff=10_000, k=23, debug=False):
    all_loci_chains = compute_loci_chains(monomers_dataset, ref2tf, delta=delta, min_array_length=min_array_length, k=k, min_fish_strength=min_fish_strength, debug=debug)
    data = get_bed_data(all_loci_chains, ref2tf, chrm2start, locus_length_cutoff=locus_length_cutoff)
    evaluation, missed_repeats_fp, missed_repeats_fn = compute_stats(data, trf_our_format, debug=debug)
    return evaluation, missed_repeats_fp, missed_repeats_fn