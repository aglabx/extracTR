#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 20.02.2024
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

from trseeker.tools.sequence_tools import get_revcomp

def load_sdat_as_dict(file_name, minimal_tf=100):
    ''' Load kmer tf data from sorted tsv file.
    '''
    data = {}
    with open(file_name) as fh:
        for line in fh:
            kmer, tf = line.strip().split()
            data[kmer] = int(tf)
            if int(tf) < minimal_tf:
                break
    return data

def load_sdat_as_list(file_name, minimal_tf=100):
    ''' Load kmer tf data from sorted tsv file.
    '''
    data = []
    with open(file_name) as fh:
        for line in fh:
            kmer, tf = line.strip().split()
            data.append((kmer, int(tf)))
            if int(tf) < minimal_tf:
                break
    return data

def compute_abundace_anomaly(sdat, ref_sdat, coverage, ref_coverage=1):
    ''' Compute abundance anomaly between two sdat files.
    '''
    all_rep = []
    kmer2abandacy_diff = {}
    for kmer, tf in sdat:
        a = tf//coverage
        b = 0
        if kmer in ref_sdat:
            b = ref_sdat[kmer]//ref_coverage
        v = a/(b+0.000001)
        all_rep.append((v, kmer, a, b))
        kmer2abandacy_diff[kmer] = (v, kmer, a, b)
        kmer2abandacy_diff[get_revcomp(kmer)] = (v, kmer, a, b)
    all_rep.sort(reverse=True)
    return all_rep, kmer2abandacy_diff
    