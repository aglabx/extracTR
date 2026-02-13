#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Unit tests for sdat_tools module."""

import unittest
import tempfile
import os

from extractr.core_functions.sdat_tools import (
    load_sdat_as_dict,
    load_sdat_as_list,
    compute_abundace_anomaly,
)


class TestLoadSdatAsDict(unittest.TestCase):
    def _write_sdat(self, lines):
        fd, path = tempfile.mkstemp(suffix=".sdat")
        os.close(fd)
        with open(path, 'w') as fh:
            for kmer, tf in lines:
                fh.write(f"{kmer}\t{tf}\n")
        return path

    def test_basic_loading(self):
        data_lines = [
            ("AAAAA", 500),
            ("CCCCC", 300),
            ("GGGGG", 100),
        ]
        path = self._write_sdat(data_lines)
        try:
            data = load_sdat_as_dict(path, minimal_tf=50)
            self.assertEqual(data["AAAAA"], 500)
            self.assertEqual(data["CCCCC"], 300)
            self.assertEqual(data["GGGGG"], 100)
        finally:
            os.unlink(path)

    def test_cutoff(self):
        data_lines = [
            ("AAAAA", 500),
            ("CCCCC", 300),
            ("GGGGG", 50),  # below 100 cutoff → loaded but stops after
            ("TTTTT", 10),  # should NOT be loaded
        ]
        path = self._write_sdat(data_lines)
        try:
            data = load_sdat_as_dict(path, minimal_tf=100)
            self.assertIn("AAAAA", data)
            self.assertIn("CCCCC", data)
            # GGGGG is loaded (line is read before break) but TTTTT is not
            self.assertIn("GGGGG", data)
            self.assertNotIn("TTTTT", data)
        finally:
            os.unlink(path)

    def test_empty_file(self):
        path = self._write_sdat([])
        try:
            data = load_sdat_as_dict(path, minimal_tf=1)
            self.assertEqual(len(data), 0)
        finally:
            os.unlink(path)


class TestLoadSdatAsList(unittest.TestCase):
    def _write_sdat(self, lines):
        fd, path = tempfile.mkstemp(suffix=".sdat")
        os.close(fd)
        with open(path, 'w') as fh:
            for kmer, tf in lines:
                fh.write(f"{kmer}\t{tf}\n")
        return path

    def test_basic_loading(self):
        data_lines = [
            ("AAAAA", 500),
            ("CCCCC", 300),
        ]
        path = self._write_sdat(data_lines)
        try:
            data = load_sdat_as_list(path, minimal_tf=100)
            self.assertEqual(len(data), 2)
            self.assertEqual(data[0], ("AAAAA", 500))
            self.assertEqual(data[1], ("CCCCC", 300))
        finally:
            os.unlink(path)

    def test_returns_list_of_tuples(self):
        data_lines = [("AAAAA", 1000)]
        path = self._write_sdat(data_lines)
        try:
            data = load_sdat_as_list(path, minimal_tf=100)
            self.assertIsInstance(data, list)
            self.assertIsInstance(data[0], tuple)
            self.assertEqual(len(data[0]), 2)
        finally:
            os.unlink(path)

    def test_cutoff_stops_reading(self):
        data_lines = [
            ("AAAAA", 500),
            ("CCCCC", 50),   # below cutoff → included, then stops
            ("GGGGG", 10),   # not read
        ]
        path = self._write_sdat(data_lines)
        try:
            data = load_sdat_as_list(path, minimal_tf=100)
            self.assertEqual(len(data), 2)
        finally:
            os.unlink(path)


class TestComputeAbundanceAnomaly(unittest.TestCase):
    def test_basic(self):
        sdat = [("AAAAA", 3000), ("CCCCC", 1000)]
        ref_sdat = {"AAAAA": 100, "CCCCC": 50}
        coverage = 30

        all_rep, kmer2diff = compute_abundace_anomaly(sdat, ref_sdat, coverage)
        self.assertEqual(len(all_rep), 2)
        # Check sorted descending by anomaly value
        self.assertGreaterEqual(all_rep[0][0], all_rep[1][0])

    def test_absent_in_ref(self):
        sdat = [("AAAAA", 3000)]
        ref_sdat = {}
        coverage = 30

        all_rep, kmer2diff = compute_abundace_anomaly(sdat, ref_sdat, coverage)
        self.assertEqual(len(all_rep), 1)
        # b should be 0
        self.assertEqual(all_rep[0][3], 0)
        # v should be very large (a / 0.000001)
        self.assertGreater(all_rep[0][0], 1000)

    def test_revcomp_in_diff(self):
        sdat = [("AAAAA", 3000)]
        ref_sdat = {"AAAAA": 100}
        coverage = 30

        _, kmer2diff = compute_abundace_anomaly(sdat, ref_sdat, coverage)
        # Both kmer and revcomp should be in the dict
        self.assertIn("AAAAA", kmer2diff)
        self.assertIn("TTTTT", kmer2diff)

    def test_values_computation(self):
        sdat = [("ACGTG", 600)]
        ref_sdat = {"ACGTG": 200}
        coverage = 10

        all_rep, _ = compute_abundace_anomaly(sdat, ref_sdat, coverage)
        kmer_entry = all_rep[0]
        a = 600 // 10  # 60
        b = 200 // 1   # 200
        self.assertEqual(kmer_entry[2], a)
        self.assertEqual(kmer_entry[3], b)

    def test_empty_sdat(self):
        all_rep, kmer2diff = compute_abundace_anomaly([], {}, 30)
        self.assertEqual(len(all_rep), 0)
        self.assertEqual(len(kmer2diff), 0)


if __name__ == "__main__":
    unittest.main()
