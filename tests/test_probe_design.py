#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Unit tests for probe_design module.

Uses a mock kmer2tf dict — no aindex dependency required.
"""

import unittest
from collections import defaultdict

from extractr.core_functions.probe_design import (
    ProbeCandidate,
    compute_gc_content,
    estimate_melting_temp,
    compute_kmer_frequency_profile,
    compute_specificity_score,
    generate_degenerate_consensus,
    design_probes_for_monomer,
    design_probes,
    write_probe_fasta,
    write_probe_tsv,
)


class MockKmer2TF:
    """Dict-like mock that returns a default frequency for any k-mer."""

    def __init__(self, data=None, default=100):
        self._data = data or {}
        self._default = default

    def __getitem__(self, key):
        return self._data.get(key, self._default)


class TestGCContent(unittest.TestCase):
    def test_empty(self):
        self.assertEqual(compute_gc_content(""), 0.0)

    def test_all_gc(self):
        self.assertAlmostEqual(compute_gc_content("GCGCGC"), 1.0)

    def test_all_at(self):
        self.assertAlmostEqual(compute_gc_content("ATATAT"), 0.0)

    def test_mixed(self):
        self.assertAlmostEqual(compute_gc_content("ATGC"), 0.5)

    def test_case_insensitive(self):
        self.assertAlmostEqual(compute_gc_content("atgc"), 0.5)


class TestMeltingTemp(unittest.TestCase):
    def test_empty(self):
        self.assertEqual(estimate_melting_temp(""), 0.0)

    def test_short_wallace(self):
        # 10 bases: 5 AT + 5 GC → Tm = 2*5 + 4*5 = 30
        seq = "ATGCAATGCG"  # 4 GC, 6 AT for len=10
        tm = estimate_melting_temp(seq)
        gc = sum(1 for c in seq.upper() if c in 'GC')
        at = len(seq) - gc
        expected = 2 * at + 4 * gc
        self.assertAlmostEqual(tm, expected)

    def test_long_salt_adjusted(self):
        seq = "A" * 20  # all AT, 20bp
        tm = estimate_melting_temp(seq)
        # Tm = 64.9 + 41*(0 - 16.4)/20 = 64.9 - 33.62 = 31.28
        self.assertAlmostEqual(tm, 64.9 + 41.0 * (0 - 16.4) / 20, places=1)

    def test_50gc_long(self):
        seq = "ATGC" * 10  # 40bp, 50% GC
        tm = estimate_melting_temp(seq)
        gc = 20
        expected = 64.9 + 41.0 * (gc - 16.4) / 40
        self.assertAlmostEqual(tm, expected, places=1)


class TestKmerFrequencyProfile(unittest.TestCase):
    def test_profile_length(self):
        monomer = "ATGCATGCATGCATGCATGCATGCATG"  # 27bp
        k = 5
        mock = MockKmer2TF(default=50)
        profile = compute_kmer_frequency_profile(monomer, mock, k)
        # Should have len(monomer) entries (circular)
        self.assertEqual(len(profile), len(monomer))

    def test_profile_values(self):
        monomer = "AAAAA"
        k = 3
        data = {"AAA": 200, "AAA": 200}
        mock = MockKmer2TF(data, default=100)
        profile = compute_kmer_frequency_profile(monomer, mock, k)
        # All k-mers wrap around the circular monomer
        for kmer, pos, freq in profile:
            self.assertEqual(len(kmer), k)


class TestSpecificityScore(unittest.TestCase):
    def test_uniform_frequencies(self):
        # All same frequency → CV=0 → score=1.0
        freqs = [100, 100, 100, 100]
        score = compute_specificity_score(freqs, 100.0)
        self.assertAlmostEqual(score, 1.0)

    def test_variable_frequencies(self):
        # High variance → low score
        freqs = [10, 1000, 10, 1000]
        score = compute_specificity_score(freqs, 100.0)
        self.assertLess(score, 0.6)

    def test_empty(self):
        self.assertEqual(compute_specificity_score([], 100.0), 0.0)

    def test_zero_median(self):
        self.assertEqual(compute_specificity_score([1, 2, 3], 0.0), 0.0)


class TestDegenerateConsensus(unittest.TestCase):
    def test_empty(self):
        self.assertEqual(generate_degenerate_consensus([]), "")

    def test_single(self):
        self.assertEqual(generate_degenerate_consensus(["ATGC"]), "ATGC")

    def test_no_variation(self):
        self.assertEqual(generate_degenerate_consensus(["ATG", "ATG"]), "ATG")

    def test_one_position_varies(self):
        # Position 1: A vs G → R
        result = generate_degenerate_consensus(["AAT", "GAT"])
        self.assertEqual(result[0], "R")  # A/G = R
        self.assertEqual(result[1], "A")
        self.assertEqual(result[2], "T")

    def test_all_four_bases(self):
        result = generate_degenerate_consensus(["A", "C", "G", "T"])
        self.assertEqual(result, "N")

    def test_two_bases_ct(self):
        result = generate_degenerate_consensus(["C", "T"])
        self.assertEqual(result, "Y")  # C/T = Y


class TestDesignProbesForMonomer(unittest.TestCase):
    def test_returns_probes(self):
        # 60bp monomer with uniform k-mer frequencies
        monomer = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
        mock = MockKmer2TF(default=500)
        probes = design_probes_for_monomer(
            monomer, "test_monomer", mock,
            k=5, probe_length=20, top_n=3,
            min_gc=0.0, max_gc=1.0,
            min_tm=0.0, max_tm=200.0,
        )
        self.assertGreater(len(probes), 0)
        self.assertLessEqual(len(probes), 3)

    def test_probe_fields(self):
        monomer = "GCGCATATATGCGCATATAT" * 3  # 60bp
        mock = MockKmer2TF(default=200)
        probes = design_probes_for_monomer(
            monomer, "m0", mock,
            k=5, probe_length=15, top_n=1,
            min_gc=0.0, max_gc=1.0,
            min_tm=0.0, max_tm=200.0,
        )
        if probes:
            p = probes[0]
            self.assertIsInstance(p, ProbeCandidate)
            self.assertEqual(len(p.sequence), 15)
            self.assertEqual(p.source_monomer_id, "m0")
            self.assertGreaterEqual(p.gc_content, 0.0)
            self.assertLessEqual(p.gc_content, 1.0)
            self.assertGreaterEqual(p.composite_score, 0.0)

    def test_gc_filter(self):
        # All-AT monomer → GC=0 → should be filtered if min_gc=0.3
        monomer = "ATATATATATATATATATATAT" * 3
        mock = MockKmer2TF(default=500)
        probes = design_probes_for_monomer(
            monomer, "at_mono", mock,
            k=5, probe_length=15, top_n=5,
            min_gc=0.3, max_gc=0.7,
        )
        self.assertEqual(len(probes), 0)

    def test_empty_monomer(self):
        mock = MockKmer2TF(default=100)
        probes = design_probes_for_monomer("", "empty", mock, k=5, probe_length=10)
        self.assertEqual(len(probes), 0)

    def test_short_monomer(self):
        # Monomer shorter than probe_length — should still work via tiling
        monomer = "ATGCATGC"  # 8bp
        mock = MockKmer2TF(default=300)
        probes = design_probes_for_monomer(
            monomer, "short", mock,
            k=5, probe_length=20, top_n=2,
            min_gc=0.0, max_gc=1.0,
            min_tm=0.0, max_tm=200.0,
        )
        # May or may not find probes depending on overlap removal, but shouldn't crash
        self.assertIsInstance(probes, list)


class TestDesignProbes(unittest.TestCase):
    def test_multiple_monomers(self):
        monomers = [
            "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC",
            "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA",
        ]
        mock = MockKmer2TF(default=400)
        probes = design_probes(
            monomers, mock, k=5, probe_length=20, top_n=2,
            min_gc=0.0, max_gc=1.0,
        )
        self.assertGreater(len(probes), 0)
        # Should be sorted by composite_score descending
        for i in range(len(probes) - 1):
            self.assertGreaterEqual(probes[i].composite_score, probes[i + 1].composite_score)


class TestWriteOutput(unittest.TestCase):
    def _make_probe(self):
        return ProbeCandidate(
            sequence="ATGCATGCATGCATGCATGC",
            gc_content=0.5,
            melting_temp=60.0,
            frequency_score=0.8,
            specificity_score=0.9,
            composite_score=0.72,
            source_monomer_id="monomer_0_100bp",
            position_in_monomer=5,
        )

    def test_write_fasta(self, tmp_path=None):
        import tempfile, os
        probe = self._make_probe()
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            path = f.name
        try:
            write_probe_fasta([probe], path)
            with open(path) as fh:
                lines = fh.readlines()
            self.assertTrue(lines[0].startswith(">"))
            self.assertIn("probe_0", lines[0])
            self.assertEqual(lines[1].strip(), probe.sequence)
        finally:
            os.unlink(path)

    def test_write_tsv(self):
        import tempfile, os
        probe = self._make_probe()
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            path = f.name
        try:
            write_probe_tsv([probe], path)
            with open(path) as fh:
                lines = fh.readlines()
            self.assertEqual(len(lines), 2)  # header + 1 data row
            header = lines[0].strip().split("\t")
            self.assertIn("probe_id", header)
            self.assertIn("composite_score", header)
            data = lines[1].strip().split("\t")
            self.assertEqual(len(data), len(header))
        finally:
            os.unlink(path)


if __name__ == "__main__":
    unittest.main()
