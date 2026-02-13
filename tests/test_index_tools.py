#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Unit tests for index_tools module.

Tests the pure-python logic. Functions that call aindex.load_aindex or
os.system are tested with mocks.
"""

import sys
import unittest
from unittest.mock import patch, MagicMock
import tempfile
import os

# Mock the aindex C extension before importing
sys.modules.setdefault("aindex", MagicMock())

from extractr.core_functions.index_tools import (
    get_index,
    compute_and_get_index,
    compute_and_get_index_for_fasta,
    get_kmer_right_read_fragments,
    get_kmer_left_read_fragments,
)


class TestGetIndex(unittest.TestCase):
    """Test get_index function with mocked aindex."""

    @patch("extractr.core_functions.index_tools.aindex")
    @patch("extractr.core_functions.index_tools.load_sdat_as_list")
    def test_builds_correct_settings(self, mock_load_sdat, mock_aindex):
        mock_load_sdat.return_value = [("AAAAA", 500)]
        mock_aindex.load_aindex.return_value = MagicMock()

        kmer2tf, sdat = get_index("/path/to/prefix", lu=100)

        mock_load_sdat.assert_called_once_with("/path/to/prefix.23.sdat", minimal_tf=100)
        mock_aindex.load_aindex.assert_called_once()
        call_args = mock_aindex.load_aindex.call_args
        settings = call_args[0][0]
        self.assertEqual(settings["index_prefix"], "/path/to/prefix.23")
        self.assertEqual(settings["aindex_prefix"], "/path/to/prefix.23")
        self.assertEqual(settings["reads_file"], "/path/to/prefix.reads")

    @patch("extractr.core_functions.index_tools.aindex")
    @patch("extractr.core_functions.index_tools.load_sdat_as_list")
    def test_returns_kmer2tf_and_sdat(self, mock_load_sdat, mock_aindex):
        expected_sdat = [("AAAAA", 500), ("CCCCC", 300)]
        mock_load_sdat.return_value = expected_sdat
        mock_kmer2tf = MagicMock()
        mock_aindex.load_aindex.return_value = mock_kmer2tf

        kmer2tf, sdat = get_index("/some/prefix", lu=50)
        self.assertEqual(sdat, expected_sdat)
        self.assertIs(kmer2tf, mock_kmer2tf)


class TestComputeAndGetIndex(unittest.TestCase):
    """Test compute_and_get_index with mocked aindex and os.system."""

    @patch("extractr.core_functions.index_tools.aindex")
    @patch("extractr.core_functions.index_tools.load_sdat_as_list")
    @patch("extractr.core_functions.index_tools.os.system")
    @patch("extractr.core_functions.index_tools.os.path.isfile")
    def test_computes_when_files_missing(self, mock_isfile, mock_system, mock_load_sdat, mock_aindex):
        mock_isfile.return_value = False
        mock_load_sdat.return_value = [("AAAAA", 500)]
        mock_aindex.load_aindex.return_value = MagicMock()

        compute_and_get_index("reads1.fq", "reads2.fq", "out", 32, lu=100)
        # Should call os.system with compute_aindex.py
        mock_system.assert_called_once()
        cmd = mock_system.call_args[0][0]
        self.assertIn("compute_aindex.py", cmd)
        self.assertIn("reads1.fq,reads2.fq", cmd)

    @patch("extractr.core_functions.index_tools.aindex")
    @patch("extractr.core_functions.index_tools.load_sdat_as_list")
    @patch("extractr.core_functions.index_tools.os.path.isfile")
    def test_skips_computation_when_files_exist(self, mock_isfile, mock_load_sdat, mock_aindex):
        mock_isfile.return_value = True
        mock_load_sdat.return_value = [("AAAAA", 500)]
        mock_aindex.load_aindex.return_value = MagicMock()

        compute_and_get_index("reads1.fq", "reads2.fq", "out", 32, lu=100)
        # os.system should NOT be called
        # (no mock_system, would fail if called)

    @patch("extractr.core_functions.index_tools.aindex")
    @patch("extractr.core_functions.index_tools.load_sdat_as_list")
    @patch("extractr.core_functions.index_tools.os.system")
    @patch("extractr.core_functions.index_tools.os.path.isfile")
    def test_se_mode(self, mock_isfile, mock_system, mock_load_sdat, mock_aindex):
        mock_isfile.return_value = False
        mock_load_sdat.return_value = []
        mock_aindex.load_aindex.return_value = MagicMock()

        compute_and_get_index("reads1.fq", None, "out", 32, lu=100)
        cmd = mock_system.call_args[0][0]
        self.assertIn("-t se", cmd)
        self.assertIn("reads1.fq", cmd)


class TestComputeAndGetIndexForFasta(unittest.TestCase):

    @patch("extractr.core_functions.index_tools.aindex")
    @patch("extractr.core_functions.index_tools.load_sdat_as_list")
    @patch("extractr.core_functions.index_tools.os.system")
    @patch("extractr.core_functions.index_tools.os.path.isfile")
    def test_fasta_mode(self, mock_isfile, mock_system, mock_load_sdat, mock_aindex):
        mock_isfile.return_value = False
        mock_load_sdat.return_value = []
        mock_aindex.load_aindex.return_value = MagicMock()

        compute_and_get_index_for_fasta("genome.fa", "out", 32, lu=100)
        cmd = mock_system.call_args[0][0]
        self.assertIn("-t fasta", cmd)
        self.assertIn("genome.fa", cmd)


class TestGetKmerReadFragments(unittest.TestCase):
    """Test get_kmer_right_read_fragments and get_kmer_left_read_fragments."""

    def _make_mock_kmer2tf(self, positions, reads_data):
        mock = MagicMock()
        mock.pos.return_value = positions
        mock.reads = reads_data
        return mock

    def test_right_fragments_split_springs(self):
        reads_data = b"AAAATTTTT~CCCC\nother_read"
        mock = self._make_mock_kmer2tf([0], reads_data)

        result = get_kmer_right_read_fragments("AAA", mock, read_length=20, topk=10)
        self.assertEqual(len(result), 1)
        # Should split by \n first, then by ~
        self.assertEqual(result[0], b"AAAATTTTT")

    def test_right_fragments_no_split(self):
        reads_data = b"AAAATTTTT~CCCC\nother_read"
        mock = self._make_mock_kmer2tf([0], reads_data)

        result = get_kmer_right_read_fragments(
            "AAA", mock, read_length=20, topk=10, split_springs=False
        )
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0], b"AAAATTTTT~CCCC")

    def test_right_fragments_topk(self):
        reads_data = b"AAAA\nBBBB\nCCCC"
        mock = self._make_mock_kmer2tf([0, 5, 10], reads_data)

        result = get_kmer_right_read_fragments("A", mock, read_length=4, topk=2)
        self.assertEqual(len(result), 2)

    def test_left_fragments(self):
        reads_data = b"PREFIX~SUFFIX\nDATA"
        mock = self._make_mock_kmer2tf([10], reads_data)

        result = get_kmer_left_read_fragments(
            "AAA", mock, read_length=5, topk=10, k=3
        )
        self.assertEqual(len(result), 1)

    def test_empty_positions(self):
        mock = self._make_mock_kmer2tf([], b"")
        result = get_kmer_right_read_fragments("AAA", mock, topk=10)
        self.assertEqual(len(result), 0)


if __name__ == "__main__":
    unittest.main()
