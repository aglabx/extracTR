#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 13.02.2026
# @author: Iaroslav Chelombitko
# @contact: les4sixstm@gmail.com

"""FISH probe design for satellite DNA monomers.

Designs optimal FISH probes from predicted tandem repeat monomers using
k-mer frequency profiles from the aindex. No reference genome required.
"""

from dataclasses import dataclass, field
from typing import List, Tuple, Dict, Optional
import math


@dataclass
class ProbeCandidate:
    """A candidate FISH probe with quality metrics."""
    sequence: str
    gc_content: float
    melting_temp: float
    frequency_score: float
    specificity_score: float
    composite_score: float
    source_monomer_id: str
    position_in_monomer: int


# IUPAC ambiguity codes for degenerate consensus
_IUPAC = {
    frozenset('A'): 'A', frozenset('C'): 'C',
    frozenset('G'): 'G', frozenset('T'): 'T',
    frozenset('AG'): 'R', frozenset('CT'): 'Y',
    frozenset('GC'): 'S', frozenset('AT'): 'W',
    frozenset('GT'): 'K', frozenset('AC'): 'M',
    frozenset('CGT'): 'B', frozenset('AGT'): 'D',
    frozenset('ACT'): 'H', frozenset('ACG'): 'V',
    frozenset('ACGT'): 'N',
}


def compute_gc_content(seq: str) -> float:
    """Compute GC fraction of a DNA sequence."""
    if not seq:
        return 0.0
    seq = seq.upper()
    gc = sum(1 for c in seq if c in ('G', 'C'))
    return gc / len(seq)


def estimate_melting_temp(seq: str) -> float:
    """Estimate melting temperature using the Wallace rule for short oligos
    and the basic salt-adjusted formula for longer ones.

    For len <= 14: Tm = 2*(A+T) + 4*(G+C)  (Wallace rule)
    For len > 14:  Tm = 81.5 + 41*(G+C)/N - 675/N  (salt-adjusted, ~50mM Na+)
    """
    if not seq:
        return 0.0
    seq = seq.upper()
    n = len(seq)
    gc = sum(1 for c in seq if c in ('G', 'C'))
    at = n - gc

    if n <= 14:
        return 2 * at + 4 * gc
    else:
        return 81.5 + 41.0 * gc / n - 675.0 / n


def compute_kmer_frequency_profile(
    monomer: str, kmer2tf, k: int
) -> List[Tuple[str, int, int]]:
    """Return (kmer, position, frequency) for each k-mer in the monomer.

    The monomer is treated as circular: positions wrap around.
    """
    n = len(monomer)
    circular = monomer + monomer[:k - 1]  # wrap for circularity
    profile = []
    for i in range(n):
        kmer = circular[i:i + k]
        freq = kmer2tf[kmer]
        profile.append((kmer, i, freq))
    return profile


def compute_specificity_score(
    probe_kmer_freqs: List[int], monomer_median_freq: float
) -> float:
    """Compute specificity score for a probe based on coefficient of variation
    of its k-mer frequencies.

    Low CV = all k-mers have similar frequency = probe is specific to one repeat.
    High CV = some k-mers are shared with other elements = less specific.

    Returns a value in [0, 1] where 1 is maximally specific (CV=0).
    """
    if not probe_kmer_freqs or monomer_median_freq <= 0:
        return 0.0

    n = len(probe_kmer_freqs)
    mean_freq = sum(probe_kmer_freqs) / n
    if mean_freq == 0:
        return 0.0

    variance = sum((f - mean_freq) ** 2 for f in probe_kmer_freqs) / n
    std = math.sqrt(variance)
    cv = std / mean_freq

    # Convert CV to [0, 1] score: CV=0 → 1.0, CV>=2 → ~0
    score = 1.0 / (1.0 + cv)
    return score


def generate_degenerate_consensus(variants: List[str]) -> str:
    """Generate IUPAC degenerate consensus from a list of sequence variants.

    All variants must be the same length. Positions where all variants agree
    get the unambiguous base; positions with variation get the IUPAC code.
    """
    if not variants:
        return ""
    if len(variants) == 1:
        return variants[0].upper()

    length = len(variants[0])
    consensus = []
    for i in range(length):
        bases = frozenset(v[i].upper() for v in variants if i < len(v))
        code = _IUPAC.get(bases, 'N')
        consensus.append(code)
    return ''.join(consensus)


def design_probes_for_monomer(
    monomer: str,
    monomer_id: str,
    kmer2tf,
    k: int = 23,
    probe_length: int = 40,
    min_gc: float = 0.35,
    max_gc: float = 0.65,
    min_tm: float = 70.0,
    max_tm: float = 95.0,
    top_n: int = 3,
) -> List[ProbeCandidate]:
    """Design FISH probes for a single monomer.

    Algorithm:
    1. Compute k-mer frequency profile for the circular monomer
    2. Slide a window of probe_length across the circular monomer
    3. For each window position: compute mean k-mer frequency and specificity
    4. Filter by GC% and Tm
    5. Rank by composite_score = normalized_freq * specificity
    6. Remove overlapping probes, return top-N
    """
    n = len(monomer)
    if n == 0:
        return []

    # If monomer is shorter than probe_length, tile it
    if n < probe_length:
        repeats_needed = (probe_length // n) + 2
        extended = monomer * repeats_needed
    else:
        extended = monomer + monomer[:probe_length - 1]  # circular wrap

    # Build frequency profile for the monomer
    freq_profile = compute_kmer_frequency_profile(monomer, kmer2tf, k)
    all_freqs = [f for _, _, f in freq_profile]
    if not all_freqs:
        return []

    sorted_freqs = sorted(all_freqs)
    monomer_median_freq = sorted_freqs[len(sorted_freqs) // 2]
    max_freq = max(all_freqs) if all_freqs else 1

    # Generate candidate probes
    candidates = []
    n_windows = min(n, len(extended) - probe_length + 1)

    for pos in range(n_windows):
        probe_seq = extended[pos:pos + probe_length]

        gc = compute_gc_content(probe_seq)
        if gc < min_gc or gc > max_gc:
            continue

        tm = estimate_melting_temp(probe_seq)
        if tm < min_tm or tm > max_tm:
            continue

        # Collect k-mer frequencies within this probe window
        probe_kmers_count = probe_length - k + 1
        if probe_kmers_count <= 0:
            continue

        # For probes spanning the circular boundary, compute freqs directly
        probe_freqs = []
        circular_monomer = monomer + monomer[:k - 1]
        for j in range(probe_kmers_count):
            kmer_pos = (pos + j) % n
            kmer = circular_monomer[kmer_pos:kmer_pos + k]
            if len(kmer) == k:
                probe_freqs.append(kmer2tf[kmer])

        if not probe_freqs:
            continue

        mean_freq = sum(probe_freqs) / len(probe_freqs)
        freq_score = mean_freq / max_freq if max_freq > 0 else 0.0
        specificity = compute_specificity_score(probe_freqs, monomer_median_freq)
        composite = freq_score * specificity

        candidates.append(ProbeCandidate(
            sequence=probe_seq,
            gc_content=gc,
            melting_temp=tm,
            frequency_score=freq_score,
            specificity_score=specificity,
            composite_score=composite,
            source_monomer_id=monomer_id,
            position_in_monomer=pos,
        ))

    # Sort by composite score descending
    candidates.sort(key=lambda p: p.composite_score, reverse=True)

    # Remove overlapping probes (greedy)
    selected = []
    used_positions = set()
    for probe in candidates:
        pos = probe.position_in_monomer
        probe_range = set(range(pos, pos + probe_length))
        # Also check circular overlap
        probe_range = set(p % n for p in probe_range)
        if probe_range & used_positions:
            continue
        selected.append(probe)
        used_positions |= probe_range
        if len(selected) >= top_n:
            break

    return selected


def design_probes(
    monomers: List[str],
    kmer2tf,
    k: int = 23,
    probe_length: int = 40,
    min_gc: float = 0.35,
    max_gc: float = 0.65,
    top_n: int = 3,
) -> List[ProbeCandidate]:
    """Design FISH probes for a list of monomers.

    Returns all probes across all monomers, sorted by composite_score.
    """
    all_probes = []
    for i, monomer in enumerate(monomers):
        monomer_id = f"monomer_{i}_{len(monomer)}bp"
        probes = design_probes_for_monomer(
            monomer, monomer_id, kmer2tf, k=k,
            probe_length=probe_length,
            min_gc=min_gc, max_gc=max_gc,
            top_n=top_n,
        )
        all_probes.extend(probes)

    all_probes.sort(key=lambda p: p.composite_score, reverse=True)
    return all_probes


def write_probe_fasta(probes: List[ProbeCandidate], path: str) -> None:
    """Write probes to FASTA file."""
    with open(path, 'w') as fh:
        for i, probe in enumerate(probes):
            header = (
                f">probe_{i}_{probe.source_monomer_id}"
                f"_pos{probe.position_in_monomer}"
                f"_gc{probe.gc_content:.2f}"
                f"_tm{probe.melting_temp:.1f}"
                f"_score{probe.composite_score:.3f}"
            )
            fh.write(f"{header}\n{probe.sequence}\n")


def write_probe_tsv(probes: List[ProbeCandidate], path: str) -> None:
    """Write probes to TSV file with all metrics."""
    header = "\t".join([
        "probe_id", "source_monomer", "position", "length",
        "sequence", "gc_content", "melting_temp",
        "frequency_score", "specificity_score", "composite_score",
    ])
    with open(path, 'w') as fh:
        fh.write(header + "\n")
        for i, probe in enumerate(probes):
            row = "\t".join([
                f"probe_{i}",
                probe.source_monomer_id,
                str(probe.position_in_monomer),
                str(len(probe.sequence)),
                probe.sequence,
                f"{probe.gc_content:.4f}",
                f"{probe.melting_temp:.1f}",
                f"{probe.frequency_score:.4f}",
                f"{probe.specificity_score:.4f}",
                f"{probe.composite_score:.4f}",
            ])
            fh.write(row + "\n")
