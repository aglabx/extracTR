#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 20.02.2024
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import argparse
from .core_functions.index_tools import compute_and_get_index, compute_and_get_index_for_fasta, get_index
from .core_functions.tr_finder import tr_greedy_finder_bidirectional
from tqdm import tqdm

from typing import Dict, List, Set
from collections import defaultdict

class KmerPathFinder:
    def __init__(self, kmer2tf):
        """
        Инициализация искателя путей с использованием индекса kmer2tf.
        
        Args:
            kmer2tf: Индекс, который возвращает частоту k-мера или 0, если его нет
        """
        self.kmer2tf = kmer2tf
        self.k = 23
        
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
    }
    
    fastq1 = settings.get("fastq1", None)
    fastq2 = settings.get("fastq2", None)
    fasta = settings.get("fasta", None)
    threads = settings.get("threads", 32)
    coverage = settings.get("coverage", 1.0)
    if settings["lu"] is None:
        settings["lu"] = 100 * settings["coverage"]
    lu = settings.get("lu")
    if lu <= 1:
        lu = 2
    prefix = settings.get("output", "test")
    min_fraction_to_continue = settings.get("min_fraction_to_continue", 30)
    k = settings.get("k", 23)

    ### step 1. Compute aindex for reads
    if fastq1 and fastq2:
        kmer2tf, sdat = compute_and_get_index(fastq1, fastq2, prefix, threads, lu=lu)
    elif fastq1 and not fastq2:
        ### SE fastq case
        kmer2tf, sdat = compute_and_get_index(fastq1, None, prefix, threads, lu=lu)
    elif fasta:
        kmer2tf, sdat = compute_and_get_index_for_fasta(fasta, prefix, threads, lu=lu)
    elif settings["aindex"]:
        kmer2tf, sdat = get_index(settings["aindex"], lu)
    else:
        raise Exception("No input data")

    ### step 2. Find tandem repeats using circular path in de bruijn graph

    repeats = tr_greedy_finder_bidirectional(sdat, kmer2tf, max_depth=30_000, coverage=coverage, min_fraction_to_continue=min_fraction_to_continue, k=k, lu=lu)

    all_predicted_trs = []
    all_predicted_te = []
    for i, (status, second_status, next_rid, next_i, seq) in enumerate(repeats):
        if status == "tr":
            seq = seq[:-k]
            # print(status, second_status, next_rid, next_i, len(seq), seq)
            all_predicted_trs.append(seq)
        elif status == "frag":
            pass
        elif status == "zero":
            all_predicted_te.append(seq)
        elif status == "long":
            pass
        else:
            # print(status, second_status, next_rid, next_i, len(seq), seq)
            raise Exception("Unknown status")
        
    

    ### step 3. Save results to CSV

    print(f"Predicted {len(all_predicted_trs)} tandem repeats.")

    output_file = f"{prefix}.fa"
    with open(output_file, "w") as fh:
        for i, seq in enumerate(all_predicted_trs):
            fh.write(f">{i}_{len(seq)}bp\n{seq}\n")

    print(f"Predicted {len(all_predicted_te)} dispersed elements.")

    output_file = f"{prefix}_te.fa"
    with open(output_file, "w") as fh:
        for i, seq in enumerate(all_predicted_te):
            fh.write(f">{i}_{len(seq)}bp\n{seq}\n")


    ### step 4. Analyze repeat borders

    ### step 5. Enrich repeats variants
    # k = 23
    # for i, seq in enumerate(all_predicted_trs):
    #     if len(seq) > 10:
    #         continue
    #     monomer_length = max(2 * k + len(seq), 2 * len(seq))
    #     monomer_n = monomer_length // len(seq)
    #     consensus = seq * monomer_n
    #     print(monomer_length, monomer_n, consensus, len(consensus))

    #     finder = KmerPathFinder(kmer2tf)
    
    #     # Найти все пути между ACGT и TACG с максимальной длиной 4
    #     paths = finder.find_all_sequences(consensus[:k], consensus[-k:], 2 * len(consensus))
    
    #     print("Найденные пути:")
    #     for path in paths:
    #         frequencies = finder.get_path_frequencies(path)
    #         print(f"Patj: {' -> '.join(path)}")
    #         print(f"Freqs: {frequencies}")
    #         print()
    #     input("Press Enter to continue...")


if __name__ == "__main__":
    run_it()