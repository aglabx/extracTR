#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Unit tests for helpers module."""

import unittest
import tempfile
import gzip
import os

from extractr.core_functions.helpers import get_revcomp, sc_iter_fasta_brute


class TestGetRevcomp(unittest.TestCase):
    def test_simple(self):
        self.assertEqual(get_revcomp("ATCG"), "CGAT")

    def test_single_base(self):
        self.assertEqual(get_revcomp("A"), "T")
        self.assertEqual(get_revcomp("T"), "A")
        self.assertEqual(get_revcomp("C"), "G")
        self.assertEqual(get_revcomp("G"), "C")

    def test_palindrome(self):
        # ATAT → revcomp = ATAT
        self.assertEqual(get_revcomp("ATAT"), "ATAT")

    def test_lowercase(self):
        self.assertEqual(get_revcomp("atcg"), "cgat")

    def test_mixed_case(self):
        self.assertEqual(get_revcomp("AtCg"), "cGaT")

    def test_n_base(self):
        self.assertEqual(get_revcomp("ANCG"), "CGNT")

    def test_empty(self):
        self.assertEqual(get_revcomp(""), "")

    def test_tilde(self):
        # Tilde maps to tilde
        self.assertEqual(get_revcomp("~"), "~")

    def test_brackets(self):
        # [ maps to ], ] maps to [
        self.assertEqual(get_revcomp("["), "]")
        self.assertEqual(get_revcomp("]"), "[")

    def test_longer_sequence(self):
        seq = "ATCGATCGATCG"
        rc = get_revcomp(seq)
        # Double revcomp should give original
        self.assertEqual(get_revcomp(rc), seq)

    def test_double_revcomp_identity(self):
        for seq in ["ACGT", "AAAA", "GCGCGC", "ATCGNNN"]:
            self.assertEqual(get_revcomp(get_revcomp(seq)), seq)


class TestScIterFastaBrute(unittest.TestCase):
    def _write_fasta(self, content, suffix=".fa", gz=False):
        """Helper to write temp FASTA file."""
        if gz:
            suffix = ".fa.gz"
        fd, path = tempfile.mkstemp(suffix=suffix)
        os.close(fd)
        if gz:
            with gzip.open(path, 'wt') as fh:
                fh.write(content)
        else:
            with open(path, 'w') as fh:
                fh.write(content)
        return path

    def test_single_record(self):
        path = self._write_fasta(">seq1\nATGCATGC\n")
        try:
            records = list(sc_iter_fasta_brute(path))
            self.assertEqual(len(records), 1)
            self.assertEqual(records[0][0], ">seq1")
            self.assertEqual(records[0][1], "ATGCATGC")
        finally:
            os.unlink(path)

    def test_multiple_records(self):
        content = ">seq1\nATGC\n>seq2\nGCTA\n>seq3\nAAAA\n"
        path = self._write_fasta(content)
        try:
            records = list(sc_iter_fasta_brute(path))
            self.assertEqual(len(records), 3)
            self.assertEqual(records[0], (">seq1", "ATGC"))
            self.assertEqual(records[1], (">seq2", "GCTA"))
            self.assertEqual(records[2], (">seq3", "AAAA"))
        finally:
            os.unlink(path)

    def test_multiline_sequence(self):
        content = ">seq1\nATGC\nGCTA\nAAAA\n"
        path = self._write_fasta(content)
        try:
            records = list(sc_iter_fasta_brute(path))
            self.assertEqual(len(records), 1)
            self.assertEqual(records[0][1], "ATGCGCTAAAAA")
        finally:
            os.unlink(path)

    def test_lower_flag(self):
        content = ">seq1\nATGC\n"
        path = self._write_fasta(content)
        try:
            records = list(sc_iter_fasta_brute(path, lower=True))
            self.assertEqual(records[0][1], "atgc")
        finally:
            os.unlink(path)

    def test_inmem_flag(self):
        content = ">seq1\nATGC\n>seq2\nGCTA\n"
        path = self._write_fasta(content)
        try:
            records = list(sc_iter_fasta_brute(path, inmem=True))
            self.assertEqual(len(records), 2)
            self.assertEqual(records[0][1], "ATGC")
        finally:
            os.unlink(path)

    def test_gzipped_file(self):
        content = ">seq1\nATGCATGC\n>seq2\nTTTT\n"
        path = self._write_fasta(content, gz=True)
        try:
            records = list(sc_iter_fasta_brute(path))
            self.assertEqual(len(records), 2)
            self.assertEqual(records[0][1], "ATGCATGC")
            self.assertEqual(records[1][1], "TTTT")
        finally:
            os.unlink(path)

    def test_empty_file(self):
        path = self._write_fasta("")
        try:
            records = list(sc_iter_fasta_brute(path))
            self.assertEqual(len(records), 0)
        finally:
            os.unlink(path)

    def test_no_trailing_newline(self):
        content = ">seq1\nATGC"
        path = self._write_fasta(content)
        try:
            records = list(sc_iter_fasta_brute(path))
            self.assertEqual(len(records), 1)
            self.assertEqual(records[0][1], "ATGC")
        finally:
            os.unlink(path)

    def test_header_with_description(self):
        content = ">seq1 some description here\nATGC\n"
        path = self._write_fasta(content)
        try:
            records = list(sc_iter_fasta_brute(path))
            self.assertEqual(records[0][0], ">seq1 some description here")
        finally:
            os.unlink(path)


if __name__ == "__main__":
    unittest.main()
