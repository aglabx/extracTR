use pyo3::prelude::*;

pub mod kmer;
pub mod kmer_index;
pub mod tr_finder;
pub mod monomer_variants;
pub mod probe_design;

use kmer_index::KmerIndex;

/// Bidirectional TR finder with DFS backtracking and two-tier threshold.
/// - lu: seed selection threshold (high). Only k-mers with tf >= lu are used as DFS roots.
/// - ext_lu: extension threshold (low). During DFS, extensions with tf >= ext_lu are followed.
/// Returns Vec<(status, second_status, next_rid, next_i, sequence)>.
#[pyfunction]
#[pyo3(signature = (sdat, max_depth=30_000, coverage=30.0, min_fraction_to_continue=30, k=23, lu=None, max_backtracks=1000, ext_lu=None))]
fn tr_greedy_finder_bidirectional(
    sdat: Vec<(String, u32)>,
    max_depth: usize,
    coverage: f64,
    min_fraction_to_continue: u32,
    k: usize,
    lu: Option<u32>,
    max_backtracks: usize,
    ext_lu: Option<u32>,
) -> PyResult<Vec<(String, Option<String>, Option<usize>, Option<usize>, Option<String>)>> {
    let results = tr_finder::tr_greedy_finder_bidirectional(
        &sdat,
        max_depth,
        coverage,
        min_fraction_to_continue,
        k,
        lu,
        max_backtracks,
        ext_lu,
    );
    Ok(results
        .into_iter()
        .map(|r| {
            (
                r.status,
                r.second_status,
                r.next_rid,
                r.next_i,
                r.sequence,
            )
        })
        .collect())
}

/// Unidirectional greedy TR finder.
#[pyfunction]
#[pyo3(signature = (sdat, max_depth=30_000, coverage=30.0, min_fraction_to_continue=30, k=23, lu=None))]
fn tr_greedy_finder(
    sdat: Vec<(String, u32)>,
    max_depth: usize,
    coverage: f64,
    min_fraction_to_continue: u32,
    k: usize,
    lu: Option<u32>,
) -> PyResult<Vec<(String, Option<String>, Option<usize>, Option<usize>, Option<String>)>> {
    let results = tr_finder::tr_greedy_finder(
        &sdat,
        max_depth,
        coverage,
        min_fraction_to_continue,
        k,
        lu,
    );
    Ok(results
        .into_iter()
        .map(|r| {
            (
                r.status,
                r.second_status,
                r.next_rid,
                r.next_i,
                r.sequence,
            )
        })
        .collect())
}

/// Naive TR finder.
/// Returns Vec<(rid, status, length, start_kmer, sequence)>.
#[pyfunction]
#[pyo3(signature = (sdat, min_tf_extension=3000, min_fraction_to_continue=30, k=23))]
fn naive_tr_finder(
    sdat: Vec<(String, u32)>,
    min_tf_extension: u32,
    min_fraction_to_continue: u32,
    k: usize,
) -> PyResult<Vec<(usize, String, usize, String, String)>> {
    Ok(tr_finder::naive_tr_finder(
        &sdat,
        min_tf_extension,
        min_fraction_to_continue,
        k,
    ))
}

/// Find monomer variants via DFS cycle search.
/// sdat is used to build the k-mer index.
#[pyfunction]
#[pyo3(signature = (monomer, sdat, k=23, max_variants=10, length_tolerance=0.15))]
fn find_monomer_variants(
    monomer: &str,
    sdat: Vec<(String, u32)>,
    k: usize,
    max_variants: usize,
    length_tolerance: f64,
) -> PyResult<Vec<String>> {
    let index = KmerIndex::from_sdat(&sdat, k);
    Ok(monomer_variants::find_monomer_variants(
        monomer,
        &index,
        k,
        max_variants,
        length_tolerance,
    ))
}

/// Design FISH probes for a single monomer.
/// Returns Vec<(sequence, gc, tm, freq_score, specificity, composite, monomer_id, position)>.
#[pyfunction]
#[pyo3(signature = (monomer, monomer_id, sdat, k=23, probe_length=40, min_gc=0.35, max_gc=0.65, min_tm=70.0, max_tm=95.0, top_n=3))]
fn design_probes_for_monomer(
    monomer: &str,
    monomer_id: &str,
    sdat: Vec<(String, u32)>,
    k: usize,
    probe_length: usize,
    min_gc: f64,
    max_gc: f64,
    min_tm: f64,
    max_tm: f64,
    top_n: usize,
) -> PyResult<Vec<(String, f64, f64, f64, f64, f64, String, usize)>> {
    let index = KmerIndex::from_sdat(&sdat, k);
    let probes = probe_design::design_probes_for_monomer(
        monomer,
        monomer_id,
        &index,
        k,
        probe_length,
        min_gc,
        max_gc,
        min_tm,
        max_tm,
        top_n,
    );
    Ok(probes
        .into_iter()
        .map(|p| {
            (
                p.sequence,
                p.gc_content,
                p.melting_temp,
                p.frequency_score,
                p.specificity_score,
                p.composite_score,
                p.source_monomer_id,
                p.position_in_monomer,
            )
        })
        .collect())
}

/// Reverse complement of a DNA string.
#[pyfunction]
fn get_revcomp(sequence: &str) -> String {
    kmer::revcomp_string(sequence)
}

/// A PyO3 Python module implemented in Rust.
#[pymodule]
fn extractr_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(tr_greedy_finder_bidirectional, m)?)?;
    m.add_function(wrap_pyfunction!(tr_greedy_finder, m)?)?;
    m.add_function(wrap_pyfunction!(naive_tr_finder, m)?)?;
    m.add_function(wrap_pyfunction!(find_monomer_variants, m)?)?;
    m.add_function(wrap_pyfunction!(design_probes_for_monomer, m)?)?;
    m.add_function(wrap_pyfunction!(get_revcomp, m)?)?;
    Ok(())
}
