#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 20.02.2024
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import argparse
import logging
import shutil
import sys
import time
from .core_functions.index_tools import compute_and_get_index, compute_and_get_index_for_fasta, get_index
from .core_functions.tr_finder import tr_greedy_finder_bidirectional
from .core_functions.probe_design import (
    design_probes, write_probe_fasta, write_probe_tsv,
    generate_degenerate_consensus,
)
from tqdm import tqdm

from typing import Dict, List, Set
from collections import defaultdict

try:
    from extractr_rs import (
        tr_greedy_finder_bidirectional as _rs_bidirectional,
        find_monomer_variants as _rs_find_monomer_variants,
        design_probes_for_monomer as _rs_design_probes_for_monomer,
    )
    _HAS_RUST = True
except ImportError:
    _HAS_RUST = False

log = logging.getLogger("extracTR")


def _setup_logging(debug=False):
    level = logging.DEBUG if debug else logging.INFO
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    ))
    log.setLevel(level)
    log.addHandler(handler)


def _elapsed(t0):
    """Format elapsed time since t0."""
    dt = time.time() - t0
    if dt < 60:
        return f"{dt:.1f}s"
    return f"{dt / 60:.1f}min"


def check_dependencies(need_index_tools=True):
    """Check that required external tools are available."""
    missing = []
    if need_index_tools:
        if not shutil.which("compute_aindex.py"):
            missing.append("  compute_aindex.py — install: pip install aindex2")
        if not shutil.which("jellyfish"):
            missing.append("  jellyfish — install: conda install -c bioconda jellyfish")
    if missing:
        log.error("Missing dependencies:")
        for m in missing:
            log.error(m)
        sys.exit(1)

class KmerPathFinder:
    def __init__(self, kmer2tf, k=23, sdat=None):
        """
        Инициализация искателя путей с использованием индекса kmer2tf.

        Args:
            kmer2tf: Индекс, который возвращает частоту k-мера или 0, если его нет
            k: Длина k-мера (default: 23)
            sdat: Optional list of (kmer, tf) tuples for Rust acceleration
        """
        self.kmer2tf = kmer2tf
        self.k = k
        self.sdat = sdat
        
    def _is_valid_kmer(self, kmer: str) -> bool:
        """
        Проверяет, существует ли k-мер в индексе.
        """
        return self.kmer2tf[kmer] > 0
        
    def _get_sequence_from_path(self, kmers: List[str]) -> str:
        """
        Преобразует путь из k-меров в последовательность нуклеотидов.
        
        Args:
            kmers: Список k-меров в пути
            
        Returns:
            str: Полная последовательность нуклеотидов
        """
        if not kmers:
            return ""
        
        sequence = kmers[0] + ''.join(kmer[-1] for kmer in kmers[1:])
        return sequence
    
    def find_all_sequences(self, start: str, end: str, max_length: int) -> Set[str]:
        """
        Поиск всех возможных последовательностей между двумя k-мерами.
        
        Args:
            start (str): Начальный k-мер
            end (str): Конечный k-мер
            max_length (int): Максимальная длина последовательности
            
        Returns:
            Set[str]: Множество всех возможных последовательностей
        """
        if not (self._is_valid_kmer(start) and self._is_valid_kmer(end)):
            return set()
            
        # Проверяем, что длина k-меров совпадает
        if len(end) != self.k or len(start) != self.k:
            return set()
            
        sequences = set()
        paths_explored = 0
        
        # Создаем progress bar с примерной оценкой возможных путей
        # Используем 4^(max_length/k) как грубую оценку максимального числа путей
        estimated_paths = min(10000, int(4 ** (max_length/self.k)))
        pbar = tqdm(total=estimated_paths, desc="Finding sequences")
        
        def dfs(current_kmer: str, path: List[str], visited: Set[str]):
            nonlocal paths_explored
            current_seq = self._get_sequence_from_path(path)
            
            if paths_explored % 100 == 0:  # Обновляем progress bar каждые 100 путей
                pbar.update(100)
                pbar.refresh()
            
            # Проверяем длину текущей последовательности
            if len(current_seq) > max_length:
                return
                
            # Если достигли конечного k-мера, добавляем последовательность
            if current_kmer == end:
                paths_explored += 1
                sequences.add(current_seq)
                return
            
            # Генерируем все возможные следующие нуклеотиды
            suffix = current_kmer[1:]
            for base in 'ACGT':
                next_kmer = suffix + base
                if self._is_valid_kmer(next_kmer) and next_kmer not in visited:
                    visited.add(next_kmer)
                    path.append(next_kmer)
                    dfs(next_kmer, path, visited)
                    path.pop()
                    visited.remove(next_kmer)
        
        # Начинаем поиск
        visited = {start}
        try:
            dfs(start, [start], visited)
        finally:
            pbar.close()
            
        return sequences

    
    def get_path_frequencies(self, path: List[str]) -> List[int]:
        """
        Получить частоты k-меров для заданного пути.

        Args:
            path (List[str]): Путь из k-меров

        Returns:
            List[int]: Список частот для каждого k-мера в пути
        """
        return [self.kmer2tf[kmer] for kmer in path]

    def find_monomer_variants(self, monomer: str, max_variants: int = 10, length_tolerance: float = 0.15) -> List[str]:
        """
        Find sequence variants of a monomer via iterative DFS cycle search in the de Bruijn graph.

        Looks for cycles of length approximately len(monomer) that start from the
        same k-mer as the consensus monomer.

        Args:
            monomer: Consensus monomer sequence
            max_variants: Maximum number of variants to collect
            length_tolerance: Fractional tolerance for cycle length (default 15%)

        Returns:
            List of variant sequences (including the original if found as a cycle)
        """
        if _HAS_RUST and self.sdat is not None:
            return _rs_find_monomer_variants(monomer, self.sdat, k=self.k, max_variants=max_variants, length_tolerance=length_tolerance)
        k = self.k
        if len(monomer) < k:
            return [monomer]

        # Build a tandem of monomer to ensure we have valid k-mers at boundaries
        tandem = monomer + monomer[:k - 1]
        start_kmer = tandem[:k]
        if not self._is_valid_kmer(start_kmer):
            return [monomer]

        target_len = len(monomer)
        min_len = int(target_len * (1 - length_tolerance))
        max_len = int(target_len * (1 + length_tolerance))

        variants = set()
        variants.add(monomer)  # always include original

        # Iterative DFS with explicit stack
        # Stack items: (current_kmer, seq_len, base_index)
        # base_index tracks which of ACGT we try next when we revisit this frame
        stack = [(start_kmer, k, 0)]
        path = [start_kmer]
        visited = {start_kmer}

        while stack and len(variants) < max_variants:
            current_kmer, seq_len, base_idx = stack[-1]

            found_next = False
            for bi in range(base_idx, 4):
                if len(variants) >= max_variants:
                    break
                base = 'ACGT'[bi]
                next_kmer = current_kmer[1:] + base
                if not self._is_valid_kmer(next_kmer):
                    continue

                new_len = seq_len + 1

                # Check if we completed a cycle back to start
                if next_kmer == start_kmer and min_len <= new_len <= max_len:
                    seq = self._get_sequence_from_path(path + [next_kmer])
                    variant = seq[:new_len]
                    variants.add(variant)
                    # Update base_index so we resume from next base
                    stack[-1] = (current_kmer, seq_len, bi + 1)
                    found_next = True
                    break

                if next_kmer not in visited and new_len <= max_len:
                    # Save resume point: next time revisit this frame, start from bi+1
                    stack[-1] = (current_kmer, seq_len, bi + 1)
                    # Push new frame
                    visited.add(next_kmer)
                    path.append(next_kmer)
                    stack.append((next_kmer, new_len, 0))
                    found_next = True
                    break

            if not found_next:
                # Backtrack
                stack.pop()
                if path:
                    removed = path.pop()
                    visited.discard(removed)

        return list(variants)



def run_it():

    parser = argparse.ArgumentParser(description="Extract and analyze tandem repeats from raw DNA sequences.")
    parser.add_argument("-1", "--fastq1", help="Input file with DNA sequences in FASTQ format.", default=None, required=False)
    parser.add_argument("-2", "--fastq2", help="Input file with DNA sequences in FASTQ format (skip for SE).", default=None, required=False)
    parser.add_argument("-f", "--fasta", help="Input genome fasta file", required=False, default=None)
    parser.add_argument("--aindex", help="Prefix for precomputed index", required=False, default=None)
    parser.add_argument("-o", "--output", help="Output file with tandem repeats in CSV format.", required=True)
    parser.add_argument("-t", "--threads", help="Number of threads to use.", default=32, type=int, required=False)
    parser.add_argument("-c", "--coverage", help="Data coverage, set 1 for genome assembly", type=float, required=True)
    parser.add_argument("--lu", help="Minimal repeat kmers coverage [100 * coverage].", default=None, type=int, required=False)
    parser.add_argument("-k", "--k", help="K-mer size to use for aindex.", default=23, type=int, required=False)
    parser.add_argument("--probe-length", help="FISH probe length in bp.", default=40, type=int, required=False)
    parser.add_argument("--top-probes", help="Number of top probes per monomer.", default=3, type=int, required=False)
    parser.add_argument("--min-gc", help="Minimum GC content for probes.", default=0.35, type=float, required=False)
    parser.add_argument("--max-gc", help="Maximum GC content for probes.", default=0.65, type=float, required=False)
    parser.add_argument("--skip-probes", help="Skip FISH probe design step.", action="store_true", default=False)
    parser.add_argument("--skip-variants", help="Skip variant enrichment step.", action="store_true", default=False)
    parser.add_argument("--debug", help="Show verbose diagnostic output.", action="store_true", default=False)
    args = parser.parse_args()
    
    settings = {
        "fastq1": args.fastq1,
        "fastq2": args.fastq2,
        "fasta": args.fasta,
        "output": args.output,
        "aindex": args.aindex,
        "threads": args.threads,
        "coverage": args.coverage,
        "lu": args.lu,
        "k": args.k,
        "min_fraction_to_continue": 30,
        "probe_length": args.probe_length,
        "top_probes": args.top_probes,
        "min_gc": args.min_gc,
        "max_gc": args.max_gc,
        "skip_probes": args.skip_probes,
        "skip_variants": args.skip_variants,
        "debug": args.debug,
    }
    
    fastq1 = settings.get("fastq1", None)
    fastq2 = settings.get("fastq2", None)
    fasta = settings.get("fasta", None)
    threads = settings.get("threads", 32)
    coverage = settings.get("coverage", 1.0)
    debug = settings.get("debug", False)
    if settings["lu"] is None:
        settings["lu"] = int(100 * settings["coverage"])
    lu = int(settings.get("lu"))
    if lu <= 1:
        lu = 2
    prefix = settings.get("output", "test")
    min_fraction_to_continue = settings.get("min_fraction_to_continue", 30)
    k = settings.get("k", 23)

    _setup_logging(debug)
    pipeline_t0 = time.time()

    log.info("extracTR v0.3.0 | k=%d, coverage=%.1f, lu=%d, backend=%s",
             k, coverage, lu, "rust" if _HAS_RUST else "python")

    # Check dependencies (skip when precomputed index is provided)
    if not settings["aindex"]:
        check_dependencies(need_index_tools=True)

    ### step 1. Compute aindex for reads (precomputed index takes priority)
    t0 = time.time()
    log.info("[1/6] Building k-mer index...")
    if settings["aindex"]:
        log.info("  Loading precomputed index: %s", settings["aindex"])
        kmer2tf, sdat = get_index(settings["aindex"], lu)
    elif fastq1 and fastq2:
        log.info("  Input: PE reads %s, %s", fastq1, fastq2)
        kmer2tf, sdat = compute_and_get_index(fastq1, fastq2, prefix, threads, lu=lu, debug=debug)
    elif fastq1 and not fastq2:
        log.info("  Input: SE reads %s", fastq1)
        kmer2tf, sdat = compute_and_get_index(fastq1, None, prefix, threads, lu=lu, debug=debug)
    elif fasta:
        log.info("  Input: FASTA %s", fasta)
        kmer2tf, sdat = compute_and_get_index_for_fasta(fasta, prefix, threads, lu=lu, debug=debug)
    else:
        raise Exception("No input data")
    log.info("  Loaded %d k-mers with tf >= %d (%s)", len(sdat), lu, _elapsed(t0))

    ### step 2. Find tandem repeats using circular path in de bruijn graph
    t0 = time.time()
    log.info("[2/6] Detecting tandem repeats (bidirectional greedy search)...")

    if _HAS_RUST:
        repeats = _rs_bidirectional(sdat, max_depth=30_000, coverage=coverage, min_fraction_to_continue=min_fraction_to_continue, k=k, lu=lu)
    else:
        repeats = tr_greedy_finder_bidirectional(sdat, kmer2tf, max_depth=30_000, coverage=coverage, min_fraction_to_continue=min_fraction_to_continue, k=k, lu=lu)

    all_predicted_trs = []
    all_predicted_te = []
    status_counts = defaultdict(int)
    for i, (status, second_status, next_rid, next_i, seq) in enumerate(repeats):
        status_counts[status] += 1
        if status == "tr":
            seq = seq[:-k]
            log.debug("  TR #%d: len=%d seq=%s", len(all_predicted_trs), len(seq), seq[:60])
            all_predicted_trs.append(seq)
        elif status == "frag":
            pass
        elif status == "zero":
            all_predicted_te.append(seq)
        elif status == "long":
            pass
        elif status == "extended":
            all_predicted_te.append(seq)
        else:
            raise Exception(f"Unknown status: {status}")

    log.info("  Found %d TRs, %d dispersed elements (tr=%d frag=%d zero=%d long=%d extended=%d) (%s)",
             len(all_predicted_trs), len(all_predicted_te),
             status_counts.get("tr", 0), status_counts.get("frag", 0),
             status_counts.get("zero", 0), status_counts.get("long", 0),
             status_counts.get("extended", 0), _elapsed(t0))

    ### step 3. Save results
    t0 = time.time()
    log.info("[3/6] Saving results...")

    output_file = f"{prefix}.fa"
    with open(output_file, "w") as fh:
        for i, seq in enumerate(all_predicted_trs):
            fh.write(f">{i}_{len(seq)}bp\n{seq}\n")
    log.info("  TRs → %s (%d sequences)", output_file, len(all_predicted_trs))

    output_file = f"{prefix}_te.fa"
    with open(output_file, "w") as fh:
        for i, seq in enumerate(all_predicted_te):
            fh.write(f">{i}_{len(seq)}bp\n{seq}\n")
    log.info("  Dispersed → %s (%d sequences, %s)", output_file, len(all_predicted_te), _elapsed(t0))

    ### step 4. Analyze repeat borders
    log.info("[4/6] Repeat border analysis — skipped (not implemented)")

    ### step 5. Enrich repeat variants
    all_variants = {}  # monomer_id -> list of variant sequences
    if not settings["skip_variants"] and all_predicted_trs:
        t0 = time.time()
        log.info("[5/6] Enriching monomer variants (%d monomers)...", len(all_predicted_trs))
        finder = KmerPathFinder(kmer2tf, k=k, sdat=sdat)
        variants_file = f"{prefix}_variants.fa"
        total_variants = 0
        with open(variants_file, "w") as fh:
            for i, monomer in enumerate(tqdm(all_predicted_trs, desc="Variant enrichment")):
                variants = finder.find_monomer_variants(monomer, max_variants=10)
                monomer_id = f"monomer_{i}_{len(monomer)}bp"
                all_variants[monomer_id] = variants
                total_variants += len(variants)
                for vi, var in enumerate(variants):
                    fh.write(f">{monomer_id}_var{vi}_{len(var)}bp\n{var}\n")
        log.info("  Variants → %s (%d total variants, %s)", variants_file, total_variants, _elapsed(t0))

        # Write IUPAC consensus for each monomer's variants
        consensus_file = f"{prefix}_consensus.fa"
        n_consensus = 0
        with open(consensus_file, "w") as fh:
            for monomer_id, variants in all_variants.items():
                if len(variants) > 1:
                    # Align variants by length (only same-length for degenerate consensus)
                    by_len = defaultdict(list)
                    for v in variants:
                        by_len[len(v)].append(v)
                    for length, group in by_len.items():
                        if len(group) > 1:
                            consensus = generate_degenerate_consensus(group)
                            fh.write(f">{monomer_id}_consensus_{length}bp_n{len(group)}\n{consensus}\n")
                            n_consensus += 1
        log.info("  Consensus → %s (%d sequences)", consensus_file, n_consensus)
    else:
        log.info("[5/6] Variant enrichment — skipped")

    ### step 6. FISH probe design
    if not settings["skip_probes"] and all_predicted_trs:
        t0 = time.time()
        log.info("[6/6] Designing FISH probes (probe_length=%d)...", settings["probe_length"])
        if _HAS_RUST:
            from .core_functions.probe_design import ProbeCandidate
            all_probes = []
            for i, monomer in enumerate(all_predicted_trs):
                monomer_id = f"monomer_{i}_{len(monomer)}bp"
                raw = _rs_design_probes_for_monomer(
                    monomer, monomer_id, sdat, k=k,
                    probe_length=settings["probe_length"],
                    min_gc=settings["min_gc"], max_gc=settings["max_gc"],
                    top_n=settings["top_probes"],
                )
                for seq, gc, tm, fs, ss, cs, mid, pos in raw:
                    all_probes.append(ProbeCandidate(
                        sequence=seq, gc_content=gc, melting_temp=tm,
                        frequency_score=fs, specificity_score=ss,
                        composite_score=cs, source_monomer_id=mid, position_in_monomer=pos,
                    ))
            all_probes.sort(key=lambda p: p.composite_score, reverse=True)
            probes = all_probes
        else:
            probes = design_probes(
                all_predicted_trs,
                kmer2tf,
                k=k,
                probe_length=settings["probe_length"],
                min_gc=settings["min_gc"],
                max_gc=settings["max_gc"],
                top_n=settings["top_probes"],
            )
        probes_fasta = f"{prefix}_probes.fa"
        probes_tsv = f"{prefix}_probes.tsv"
        write_probe_fasta(probes, probes_fasta)
        write_probe_tsv(probes, probes_tsv)
        log.info("  Probes → %s, %s (%d probes, %s)", probes_fasta, probes_tsv, len(probes), _elapsed(t0))
    else:
        log.info("[6/6] FISH probe design — skipped")

    log.info("Done. Total time: %s", _elapsed(pipeline_t0))


if __name__ == "__main__":
    run_it()