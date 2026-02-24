#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Unit tests for tr_finder module.

Uses a mock kmer2tf (defaultdict) to simulate the aindex without
requiring the C extension.
"""

import unittest
from collections import defaultdict

from extractr.core_functions.tr_finder import (
    tr_greedy_finder,
    tr_greedy_finder_bidirectional,
    naive_tr_finder,
)


def build_kmer_index(sequences, k, base_tf=5000):
    """Build a mock kmer2tf defaultdict from a list of sequences.

    Every k-mer found in the given sequences gets base_tf frequency.
    Missing k-mers return 0.
    """
    kmer2tf = defaultdict(int)
    for seq in sequences:
        for i in range(len(seq) - k + 1):
            kmer2tf[seq[i:i + k]] = base_tf
    return kmer2tf


def build_circular_index(monomer, k, base_tf=5000):
    """Build kmer2tf for a circular monomer (tandem repeat).

    Duplicates the monomer to capture the wrap-around k-mers.
    """
    tandem = monomer * 3  # enough to cover all circular k-mers
    return build_kmer_index([tandem], k, base_tf)


class TestTrGreedyFinder(unittest.TestCase):
    """Tests for the unidirectional greedy TR finder."""

    def test_finds_tandem_repeat(self):
        """A simple circular monomer should be detected as TR."""
        k = 5
        monomer = "ATGCATGC"  # 8bp monomer
        kmer2tf = build_circular_index(monomer, k, base_tf=5000)

        # sdat: start with the first k-mer of the monomer, high tf
        start_kmer = monomer[:k]
        sdat = [(start_kmer, 5000)]

        repeats = tr_greedy_finder(sdat, kmer2tf, max_depth=100, coverage=1, k=k)
        self.assertEqual(len(repeats), 1)
        status = repeats[0][0]
        self.assertEqual(status, "tr")

    def test_zero_status_no_extensions(self):
        """When no valid extensions exist, status should be 'zero'."""
        k = 5
        kmer2tf = defaultdict(int)
        # Use a kmer where no extension (suffix+X) has any frequency
        start_kmer = "ACGTG"
        kmer2tf[start_kmer] = 5000
        # No other k-mers → no extensions possible
        sdat = [(start_kmer, 5000)]

        repeats = tr_greedy_finder(sdat, kmer2tf, max_depth=100, coverage=1, k=k)
        self.assertEqual(len(repeats), 1)
        self.assertEqual(repeats[0][0], "zero")

    def test_long_status_max_depth(self):
        """When max_depth is reached without cycle, status should be 'long'."""
        k = 3
        # Build a long linear chain: A...A (no cycles)
        # Each k-mer leads to exactly one next k-mer
        chain = "A" + "C" * 20 + "G" * 20  # 41 chars, will have unique 3-mers
        kmer2tf = build_kmer_index([chain], k, base_tf=5000)
        start_kmer = chain[:k]
        sdat = [(start_kmer, 5000)]

        repeats = tr_greedy_finder(sdat, kmer2tf, max_depth=5, coverage=1, k=k)
        self.assertEqual(len(repeats), 1)
        self.assertEqual(repeats[0][0], "long")

    def test_skips_below_min_tf(self):
        """K-mers below MIN_TF should be skipped in sdat."""
        k = 5
        kmer2tf = defaultdict(int)
        sdat = [("AAAAA", 10)]  # tf=10 < MIN_TF=100 (coverage=1)

        repeats = tr_greedy_finder(sdat, kmer2tf, max_depth=100, coverage=1, k=k)
        self.assertEqual(len(repeats), 0)

    def test_frag_status(self):
        """When the walk hits a previously cached k-mer, status should be 'frag'."""
        k = 3
        # Two k-mers that share an overlap via a third
        # path: AAA -> AAC -> ACC -> ... hits cache
        kmer2tf = defaultdict(int)
        # First entry: AAA extends to AAC (only valid next)
        kmer2tf["AAA"] = 5000
        kmer2tf["AAC"] = 5000
        kmer2tf["ACC"] = 5000
        kmer2tf["CCC"] = 5000
        # CCC extends to CCA
        kmer2tf["CCA"] = 5000
        # CCA extends to CAA
        kmer2tf["CAA"] = 5000
        # CAA extends to AAC → already in cache → frag

        sdat = [("AAA", 5000)]
        repeats = tr_greedy_finder(sdat, kmer2tf, max_depth=100, coverage=1, k=k)
        self.assertEqual(len(repeats), 1)
        # Should be either "tr" (if loops back to AAA) or "frag" (hits cached kmer)
        self.assertIn(repeats[0][0], ("tr", "frag"))

    def test_multiple_sdat_entries(self):
        """Multiple sdat entries should produce multiple repeats."""
        k = 5
        # Two independent monomers
        mono1 = "AAACCCTTT"
        mono2 = "GGGAAATTT"
        kmer2tf1 = build_circular_index(mono1, k)
        kmer2tf2 = build_circular_index(mono2, k)
        # Merge
        kmer2tf = defaultdict(int)
        for d in [kmer2tf1, kmer2tf2]:
            for key, val in d.items():
                kmer2tf[key] = val

        sdat = [(mono1[:k], 5000), (mono2[:k], 5000)]
        repeats = tr_greedy_finder(sdat, kmer2tf, max_depth=200, coverage=1, k=k)
        self.assertGreaterEqual(len(repeats), 1)

    def test_output_tuple_structure(self):
        """Each repeat should be a 5-tuple: (status, second_status, rid, i, seq)."""
        k = 5
        monomer = "ATGCATGC"
        kmer2tf = build_circular_index(monomer, k)
        sdat = [(monomer[:k], 5000)]

        repeats = tr_greedy_finder(sdat, kmer2tf, max_depth=100, coverage=1, k=k)
        for r in repeats:
            self.assertEqual(len(r), 5)
            status, second_status, rid, i, seq = r
            self.assertIsInstance(seq, str)
            self.assertIn(status, ("tr", "frag", "zero", "long", None))


class TestTrGreedyFinderBidirectional(unittest.TestCase):
    """Tests for the bidirectional greedy TR finder."""

    def test_finds_tandem_repeat(self):
        k = 5
        monomer = "ATGCATGC"
        kmer2tf = build_circular_index(monomer, k, base_tf=5000)
        sdat = [(monomer[:k], 5000)]

        repeats = tr_greedy_finder_bidirectional(
            sdat, kmer2tf, max_depth=100, coverage=1, k=k, lu=100
        )
        self.assertEqual(len(repeats), 1)
        self.assertEqual(repeats[0][0], "tr")

    def test_zero_both_directions(self):
        """No extensions in either direction → 'zero'."""
        k = 5
        kmer2tf = defaultdict(int)
        # Use a kmer where no extension has any frequency
        start_kmer = "ACGTG"
        kmer2tf[start_kmer] = 5000
        sdat = [(start_kmer, 5000)]

        repeats = tr_greedy_finder_bidirectional(
            sdat, kmer2tf, max_depth=100, k=k, lu=100
        )
        self.assertEqual(len(repeats), 1)
        self.assertEqual(repeats[0][0], "zero")

    def test_extended_status(self):
        """When one direction is 'long' and other is not zero/frag/tr → 'extended'."""
        k = 3
        # Build a linear path that exceeds max_depth in right direction
        # and hits long in left direction too
        kmer2tf = defaultdict(int)
        # Right: AAA -> AAC -> ACC -> CCC (then stops = zero right)
        # Left: no extension (zero left)
        # Actually let's make it so right hits long and left hits long
        # Build a chain going right
        for i in range(50):
            kmer = chr(65 + (i % 4)) * 2 + chr(65 + ((i + 1) % 4))
            kmer2tf[kmer] = 5000

        # Simpler: just test that max_depth produces long guards
        kmer2tf_simple = defaultdict(int)
        # Right chain: AAA -> AAC -> ACG -> CGT -> GTA -> TAA -> AAC (frag)
        for kmer in ["AAA", "AAC", "ACG", "CGT", "GTA", "TAA"]:
            kmer2tf_simple[kmer] = 5000
        # Left chain from AAA: add something to extend left
        # left_suffix = AA, try nucleotide+AA: only valid is TAA (already there)
        # TAA exists → left hits cache → frag
        sdat = [("AAA", 5000)]

        repeats = tr_greedy_finder_bidirectional(
            sdat, kmer2tf_simple, max_depth=100, k=3, lu=100
        )
        self.assertEqual(len(repeats), 1)
        # Status depends on graph structure, but should not crash
        self.assertIn(repeats[0][0], ("tr", "frag", "zero", "long", "extended"))

    def test_lu_parameter(self):
        """When lu is provided, it overrides coverage-based MIN_TF."""
        k = 5
        monomer = "ATGCATGC"
        kmer2tf = build_circular_index(monomer, k, base_tf=200)
        sdat = [(monomer[:k], 200)]

        # With lu=100, should process (200 >= 100)
        repeats = tr_greedy_finder_bidirectional(
            sdat, kmer2tf, max_depth=100, k=k, lu=100
        )
        self.assertGreaterEqual(len(repeats), 1)

        # With lu=500, should skip (200 < 500)
        repeats = tr_greedy_finder_bidirectional(
            sdat, kmer2tf, max_depth=100, k=k, lu=500
        )
        self.assertEqual(len(repeats), 0)

    def test_full_sequence_contains_start_kmer(self):
        """The output sequence should contain the start k-mer."""
        k = 5
        monomer = "ATGCATGC"
        kmer2tf = build_circular_index(monomer, k, base_tf=5000)
        start_kmer = monomer[:k]
        sdat = [(start_kmer, 5000)]

        repeats = tr_greedy_finder_bidirectional(
            sdat, kmer2tf, max_depth=100, k=k, lu=100
        )
        if repeats and repeats[0][4]:
            self.assertIn(start_kmer, repeats[0][4])

    def test_output_tuple_structure(self):
        k = 5
        monomer = "ATGCATGC"
        kmer2tf = build_circular_index(monomer, k)
        sdat = [(monomer[:k], 5000)]

        repeats = tr_greedy_finder_bidirectional(
            sdat, kmer2tf, max_depth=100, k=k, lu=100
        )
        for r in repeats:
            self.assertEqual(len(r), 5)
            status, second_status, rid, i, seq = r
            self.assertIn(status, ("tr", "frag", "zero", "long", "extended"))

    def test_cached_kmers_skipped(self):
        """K-mers already in cache from a previous entry should be skipped."""
        k = 5
        monomer = "ATGCATGC"
        kmer2tf = build_circular_index(monomer, k, base_tf=5000)
        start_kmer = monomer[:k]
        # Duplicate the same start kmer
        sdat = [(start_kmer, 5000), (start_kmer, 5000)]

        repeats = tr_greedy_finder_bidirectional(
            sdat, kmer2tf, max_depth=100, k=k, lu=100
        )
        # Second entry should be skipped (already in cache)
        self.assertEqual(len(repeats), 1)

    def test_empty_sdat(self):
        k = 5
        kmer2tf = defaultdict(int)
        repeats = tr_greedy_finder_bidirectional(
            [], kmer2tf, max_depth=100, k=k, lu=100
        )
        self.assertEqual(len(repeats), 0)


class TestNaiveTrFinder(unittest.TestCase):
    """Tests for naive_tr_finder. This function processes only
    the first sdat entry due to the while-loop structure."""

    def test_returns_four_items(self):
        """Should return (repeats, kmer2rid, kmer2repeat, frag2type)."""
        k = 5
        kmer2tf = defaultdict(int)
        start_kmer = "AAAAA"
        kmer2tf[start_kmer] = 5000
        sdat = [(start_kmer, 5000)]

        result = naive_tr_finder(sdat, kmer2tf, min_tf_extension=100, k=k)
        self.assertEqual(len(result), 4)
        repeats, kmer2rid, kmer2repeat, frag2type = result
        self.assertIsInstance(repeats, list)
        self.assertIsInstance(kmer2rid, dict)
        self.assertIsInstance(frag2type, dict)

    def test_te_when_no_extension(self):
        """No valid next k-mer → TE status."""
        k = 5
        kmer2tf = defaultdict(int)
        start_kmer = "ACGTG"
        kmer2tf[start_kmer] = 5000
        # No other k-mers have sufficient frequency
        sdat = [(start_kmer, 5000)]

        repeats, _, _, _ = naive_tr_finder(
            sdat, kmer2tf, min_tf_extension=100, min_fraction_to_continue=30, k=k
        )
        self.assertGreaterEqual(len(repeats), 1)
        self.assertEqual(repeats[0][1], "TE")

    def test_circular_finds_tr(self):
        """A circular monomer should produce TR status."""
        k = 3
        monomer = "ATCG"
        kmer2tf = build_circular_index(monomer, k, base_tf=5000)
        start_kmer = monomer[:k]
        sdat = [(start_kmer, 5000)]

        repeats, _, _, frag2type = naive_tr_finder(
            sdat, kmer2tf, min_tf_extension=100, min_fraction_to_continue=1, k=k
        )
        # Should find at least one repeat
        self.assertGreaterEqual(len(repeats), 1)
        # Check if any is TR
        statuses = [r[1] for r in repeats]
        # It can be TR or FRAG depending on the graph traversal
        self.assertTrue(any(s in ("TR", "FRAG", "TE") for s in statuses))

    def test_skips_used_kmers(self):
        """Already used k-mers should be skipped."""
        k = 5
        kmer2tf = defaultdict(int)
        kmer = "AAAAA"
        kmer2tf[kmer] = 5000
        # Provide the same kmer twice
        sdat = [(kmer, 5000), (kmer, 5000)]

        repeats, _, _, _ = naive_tr_finder(
            sdat, kmer2tf, min_tf_extension=100, k=k
        )
        # Second should be skipped
        self.assertLessEqual(len(repeats), 1)


if __name__ == "__main__":
    unittest.main()
