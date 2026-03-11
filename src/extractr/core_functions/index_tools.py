#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 20.02.2024
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import logging
import os
import subprocess
import time
from .sdat_tools import load_sdat_as_list
import aindex

log = logging.getLogger("extracTR")


def _run_command(command, debug=False):
    """Run a shell command, suppressing output unless debug is enabled."""
    log.debug("  Running: %s", command)
    t0 = time.time()
    if debug:
        result = subprocess.run(command, shell=True)
    else:
        result = subprocess.run(command, shell=True,
                                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed (exit {result.returncode}): {command}")
    dt = time.time() - t0
    log.debug("  Command finished in %.1fs", dt)


def _needs_gz_workaround(fastq1, fastq2):
    """Check if any input files are gzipped (jellyfish can't read .gz directly)."""
    for f in (fastq1, fastq2):
        if f and f.endswith(".gz"):
            return True
    return False


def _decompress_for_index(fastq1, fastq2, prefix, threads, lu, reads_type, debug=False):
    """Use jellyfish generator to handle .gz files without decompressing to disk."""
    gen_file = f"{prefix}.generator.sh"
    files = [f for f in (fastq1, fastq2) if f]
    with open(gen_file, "w") as fh:
        for f in files:
            if f.endswith(".gz"):
                fh.write(f"zcat {f}\n")
            else:
                fh.write(f"cat {f}\n")
    n_generators = min(len(files), 2)
    jf2_file = f"{prefix}.23.jf2"
    command = (
        f"jellyfish count -m 23 -t {threads} -s 5G -C -L {lu} "
        f"-o {jf2_file} -g {gen_file} -G {n_generators}"
    )
    _run_command(command, debug=debug)
    os.remove(gen_file)
    # dump + sort via compute_aindex.py with existing jf2
    # Pass a real file as -i to satisfy existence check (it won't be used since -j is set)
    real_file = files[0]
    command = f"compute_aindex.py -i {real_file} -j {jf2_file} -t {reads_type} -o {prefix} --lu {lu} --sort 1 -P {threads} --onlyindex 1"
    _run_command(command, debug=debug)


def compute_and_get_index(fastq1, fastq2, prefix, threads, lu=2, debug=False):

    lu = int(lu)
    sdat_file = f"{prefix}.23.sdat"
    index_prefix_file = f"{prefix}.23"

    if not os.path.isfile(sdat_file) or not os.path.isfile(index_prefix_file):
      log.info("  Computing aindex (threads=%d)...", threads)
      if _needs_gz_workaround(fastq1, fastq2):
          log.info("  Using jellyfish generator for .gz files")
          reads_type = "fastq" if fastq1 and fastq2 else "se"
          _decompress_for_index(fastq1, fastq2, prefix, threads, lu, reads_type, debug=debug)
      elif fastq1 and fastq2:
          command = f"compute_aindex.py -i {fastq1},{fastq2} -t fastq -o {prefix} --lu {lu} --sort 1 -P {threads} --onlyindex 1"
          _run_command(command, debug=debug)
      elif fastq1 and fastq2 is None:
          command = f"compute_aindex.py -i {fastq1} -t se -o {prefix} --lu {lu} --sort 1 -P {threads} --onlyindex 1"
          _run_command(command, debug=debug)
      elif fastq2 and fastq1 is None:
          command = f"compute_aindex.py -i {fastq2} -t se -o {prefix} --lu {lu} --sort 1 -P {threads} --onlyindex 1"
          _run_command(command, debug=debug)
    else:
      log.info("  Index files found, skipping computation")

    log.info("  Loading sdat: %s", sdat_file)
    sdat = load_sdat_as_list(sdat_file, minimal_tf=lu)

    log.info("  Loading aindex: %s", f"{prefix}.23")
    settings = {
        "index_prefix": f"{prefix}.23",
        "aindex_prefix": f"{prefix}.23",
        "reads_file": f"{prefix}.reads",
        "max_tf": 10000000,
    }

    kmer2tf = aindex.load_aindex(settings, skip_reads=True, skip_aindex=True)
    return kmer2tf, sdat

def compute_and_get_index_for_fasta(fasta_file, prefix, threads, lu=2, debug=False):

    lu = int(lu)
    sdat_file = f"{prefix}.23.sdat"
    index_prefix_file = f"{prefix}.23"

    if not os.path.isfile(sdat_file) or not os.path.isfile(index_prefix_file):
        log.info("  Computing aindex for FASTA (threads=%d)...", threads)
        command = f"compute_aindex.py -i {fasta_file} -t fasta -o {prefix} --lu {lu} --sort 1 -P {threads} --onlyindex 1"
        _run_command(command, debug=debug)
    else:
        log.info("  Index files found, skipping computation")

    log.info("  Loading sdat: %s", sdat_file)
    sdat = load_sdat_as_list(sdat_file, minimal_tf=lu)

    log.info("  Loading aindex: %s", f"{prefix}.23")
    settings = {
        "index_prefix": f"{prefix}.23",
        "aindex_prefix": f"{prefix}.23",
        "reads_file": f"{prefix}.reads",
        "max_tf": 10000000,
    }

    kmer2tf = aindex.load_aindex(settings, skip_reads=True, skip_aindex=True)
    return kmer2tf, sdat

def get_index(prefix, lu=2):
    """Load a precomputed aindex without recomputing."""
    sdat_file = f"{prefix}.23.sdat"
    log.info("  Loading sdat: %s", sdat_file)
    sdat = load_sdat_as_list(sdat_file, minimal_tf=lu)

    log.info("  Loading aindex: %s", f"{prefix}.23")
    settings = {
        "index_prefix": f"{prefix}.23",
        "aindex_prefix": f"{prefix}.23",
        "reads_file": f"{prefix}.reads",
        "max_tf": 10000000,
    }

    kmer2tf = aindex.load_aindex(settings, skip_reads=True, skip_aindex=True)
    return kmer2tf, sdat


def print_kmer_right_read_fragments(kmer, kmer2tf, read_length=310, topk=100, split_springs=True):
  pos = kmer2tf.pos(kmer)
  for p in pos[:topk]:
    if split_springs:
      print(kmer2tf.reads[p:p+read_length].split(b"\n")[0].split("~")[0])
    else:
      print(kmer2tf.reads[p:p+read_length].split(b"\n")[0])

def print_kmer_left_read_fragments(kmer, kmer2tf, topk=100, read_length=310, k=23, split_springs=True):
  pos = kmer2tf.pos(kmer)
  for p in pos[:topk]:
    if split_springs:
      print(kmer2tf.reads[p-read_length:p+k].split(b"\n")[-1].split(b"~")[-1])
    else:
      print(kmer2tf.reads[p-read_length:p+k].split(b"\n")[-1])


def get_kmer_right_read_fragments(kmer, kmer2tf, read_length=310, topk=100, split_springs=True):
  pos = kmer2tf.pos(kmer)
  reads = []
  for p in pos[:topk]:
    if split_springs:
      reads.append(kmer2tf.reads[p:p+read_length].split(b"\n")[0].split(b"~")[0])
    else:
      reads.append(kmer2tf.reads[p:p+read_length].split(b"\n")[0])
  return reads

def get_kmer_left_read_fragments(kmer, kmer2tf, topk=100, read_length=310, k=23, split_springs=True):
  pos = kmer2tf.pos(kmer)
  reads = []
  for p in pos[:topk]:
    if split_springs:
      reads.append(kmer2tf.reads[p-read_length:p+k].split(b"\n")[-1].split(b"~")[-1])
    else:
      reads.append(kmer2tf.reads[p-read_length:p+k].split(b"\n")[-1])
  return reads
