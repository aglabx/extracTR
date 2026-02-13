#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Unit tests for KmerPathFinder class in extractr.py.

Uses mock kmer2tf — no aindex dependency.
"""

import sys
import unittest
from unittest.mock import MagicMock
from collections import defaultdict

# Mock the aindex C extension before importing extractr
sys.modules.setdefault("aindex", MagicMock())

from extractr.extractr import KmerPathFinder


def build_kmer_index(sequences, k, base_tf=500):
    """Build a mock kmer2tf defaultdict from sequences."""
    kmer2tf = defaultdict(int)
    for seq in sequences:
        for i in range(len(seq) - k + 1):
            kmer2tf[seq[i:i + k]] = base_tf
    return kmer2tf


class TestKmerPathFinderInit(unittest.TestCase):
    def test_default_k(self):
        kmer2tf = defaultdict(int)
        finder = KmerPathFinder(kmer2tf)
        self.assertEqual(finder.k, 23)

    def test_custom_k(self):
        kmer2tf = defaultdict(int)
        finder = KmerPathFinder(kmer2tf, k=31)
        self.assertEqual(finder.k, 31)

    def test_stores_kmer2tf(self):
        kmer2tf = defaultdict(int)
        kmer2tf["AAAAA"] = 100
        finder = KmerPathFinder(kmer2tf, k=5)
        self.assertEqual(finder.kmer2tf["AAAAA"], 100)


class TestIsValidKmer(unittest.TestCase):
    def test_valid(self):
        kmer2tf = defaultdict(int)
        kmer2tf["ATGCG"] = 500
        finder = KmerPathFinder(kmer2tf, k=5)
        self.assertTrue(finder._is_valid_kmer("ATGCG"))

    def test_invalid(self):
        kmer2tf = defaultdict(int)
        finder = KmerPathFinder(kmer2tf, k=5)
        self.assertFalse(finder._is_valid_kmer("ZZZZZ"))

    def test_zero_frequency(self):
        kmer2tf = defaultdict(int)
        kmer2tf["AAAAA"] = 0
        finder = KmerPathFinder(kmer2tf, k=5)
        self.assertFalse(finder._is_valid_kmer("AAAAA"))


class TestGetSequenceFromPath(unittest.TestCase):
    def test_empty_path(self):
        finder = KmerPathFinder(defaultdict(int), k=3)
        self.assertEqual(finder._get_sequence_from_path([]), "")

    def test_single_kmer(self):
        finder = KmerPathFinder(defaultdict(int), k=3)
        self.assertEqual(finder._get_sequence_from_path(["ATG"]), "ATG")

    def test_two_kmers(self):
        finder = KmerPathFinder(defaultdict(int), k=3)
        # ATG -> TGC → sequence = ATGC
        self.assertEqual(finder._get_sequence_from_path(["ATG", "TGC"]), "ATGC")

    def test_three_kmers(self):
        finder = KmerPathFinder(defaultdict(int), k=3)
        # ATG -> TGC -> GCA → sequence = ATGCA
        self.assertEqual(
            finder._get_sequence_from_path(["ATG", "TGC", "GCA"]), "ATGCA"
        )

    def test_longer_k(self):
        finder = KmerPathFinder(defaultdict(int), k=5)
        # ATGCA -> TGCAT → sequence = ATGCAT
        self.assertEqual(
            finder._get_sequence_from_path(["ATGCA", "TGCAT"]), "ATGCAT"
        )


class TestGetPathFrequencies(unittest.TestCase):
    def test_basic(self):
        kmer2tf = defaultdict(int)
        kmer2tf["ATG"] = 100
        kmer2tf["TGC"] = 200
        kmer2tf["GCA"] = 50
        finder = KmerPathFinder(kmer2tf, k=3)

        freqs = finder.get_path_frequencies(["ATG", "TGC", "GCA"])
        self.assertEqual(freqs, [100, 200, 50])

    def test_missing_kmer(self):
        kmer2tf = defaultdict(int)
        finder = KmerPathFinder(kmer2tf, k=3)
        freqs = finder.get_path_frequencies(["AAA", "BBB"])
        self.assertEqual(freqs, [0, 0])

    def test_empty_path(self):
        finder = KmerPathFinder(defaultdict(int), k=3)
        self.assertEqual(finder.get_path_frequencies([]), [])


class TestFindAllSequences(unittest.TestCase):
    def test_invalid_start(self):
        kmer2tf = defaultdict(int)
        finder = KmerPathFinder(kmer2tf, k=3)
        result = finder.find_all_sequences("AAA", "CCC", 10)
        self.assertEqual(result, set())

    def test_wrong_kmer_length(self):
        kmer2tf = defaultdict(int)
        kmer2tf["AA"] = 100
        kmer2tf["CCC"] = 100
        finder = KmerPathFinder(kmer2tf, k=3)
        result = finder.find_all_sequences("AA", "CCC", 10)
        self.assertEqual(result, set())

    def test_simple_path(self):
        """Find path in a simple linear graph: ATG -> TGC -> GCA."""
        k = 3
        kmer2tf = defaultdict(int)
        kmer2tf["ATG"] = 100
        kmer2tf["TGC"] = 100
        kmer2tf["GCA"] = 100
        finder = KmerPathFinder(kmer2tf, k=k)

        result = finder.find_all_sequences("ATG", "GCA", 10)
        self.assertGreater(len(result), 0)
        # The sequence should be ATGCA
        self.assertIn("ATGCA", result)

    def test_start_equals_end(self):
        """When start == end, the DFS should find it immediately."""
        k = 3
        kmer2tf = defaultdict(int)
        kmer2tf["ATG"] = 100
        finder = KmerPathFinder(kmer2tf, k=k)

        result = finder.find_all_sequences("ATG", "ATG", 10)
        # Should find the trivial path (just the kmer itself)
        self.assertIn("ATG", result)


class TestFindMonomerVariants(unittest.TestCase):
    def test_short_monomer(self):
        """Monomer shorter than k should return just the monomer."""
        kmer2tf = defaultdict(int)
        finder = KmerPathFinder(kmer2tf, k=5)
        result = finder.find_monomer_variants("ATG", max_variants=10)
        self.assertEqual(result, ["ATG"])

    def test_invalid_start_kmer(self):
        """If the start k-mer has zero frequency, return just the monomer."""
        kmer2tf = defaultdict(int)
        finder = KmerPathFinder(kmer2tf, k=3)
        monomer = "ATGCATGC"
        result = finder.find_monomer_variants(monomer)
        self.assertEqual(result, [monomer])

    def test_always_includes_original(self):
        """The original monomer should always be in the result."""
        k = 3
        monomer = "ATGCAT"
        kmer2tf = defaultdict(int)
        # Add all k-mers from the circular monomer
        tandem = monomer * 3
        for i in range(len(tandem) - k + 1):
            kmer2tf[tandem[i:i + k]] = 500
        finder = KmerPathFinder(kmer2tf, k=k)

        result = finder.find_monomer_variants(monomer, max_variants=10)
        self.assertIn(monomer, result)

    def test_max_variants_limit(self):
        """Should not return more than max_variants."""
        k = 3
        monomer = "ATGC"
        # Build a rich graph where many cycles exist
        kmer2tf = defaultdict(int)
        for a in "ACGT":
            for b in "ACGT":
                for c in "ACGT":
                    kmer2tf[a + b + c] = 500
        finder = KmerPathFinder(kmer2tf, k=k)

        result = finder.find_monomer_variants(monomer, max_variants=5)
        self.assertLessEqual(len(result), 5)

    def test_returns_list(self):
        kmer2tf = defaultdict(int)
        finder = KmerPathFinder(kmer2tf, k=5)
        result = finder.find_monomer_variants("ATGCATGCATGCATGCATGCATGCA")
        self.assertIsInstance(result, list)

    def test_circular_monomer_finds_variants(self):
        """For a monomer with a clear circular path and alternative paths,
        should find at least the original."""
        k = 3
        monomer = "ATGCG"  # 5bp
        tandem = monomer * 3
        kmer2tf = defaultdict(int)
        for i in range(len(tandem) - k + 1):
            kmer2tf[tandem[i:i + k]] = 500
        finder = KmerPathFinder(kmer2tf, k=k)

        result = finder.find_monomer_variants(monomer, max_variants=10)
        self.assertGreaterEqual(len(result), 1)
        self.assertIn(monomer, result)


if __name__ == "__main__":
    unittest.main()
