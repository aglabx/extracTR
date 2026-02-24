/// FISH probe design for satellite DNA monomers.

use crate::kmer::encode_kmer;
use crate::kmer_index::KmerIndex;

/// Probe candidate with quality metrics.
pub struct ProbeCandidate {
    pub sequence: String,
    pub gc_content: f64,
    pub melting_temp: f64,
    pub frequency_score: f64,
    pub specificity_score: f64,
    pub composite_score: f64,
    pub source_monomer_id: String,
    pub position_in_monomer: usize,
}

/// Compute GC fraction of a DNA sequence.
#[inline]
pub fn compute_gc_content(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let gc = seq
        .iter()
        .filter(|&&b| b == b'G' || b == b'C' || b == b'g' || b == b'c')
        .count();
    gc as f64 / seq.len() as f64
}

/// Estimate melting temperature.
/// For len <= 14: Tm = 2*(A+T) + 4*(G+C)  (Wallace rule)
/// For len > 14:  Tm = 81.5 + 41*(G+C)/N - 675/N  (salt-adjusted)
#[inline]
pub fn estimate_melting_temp(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let n = seq.len() as f64;
    let gc = seq
        .iter()
        .filter(|&&b| b == b'G' || b == b'C' || b == b'g' || b == b'c')
        .count() as f64;
    let at = n - gc;

    if seq.len() <= 14 {
        2.0 * at + 4.0 * gc
    } else {
        81.5 + 41.0 * gc / n - 675.0 / n
    }
}

/// Compute specificity score based on coefficient of variation.
/// Low CV = specific. Returns [0, 1] where 1 is maximally specific.
#[inline]
pub fn compute_specificity_score(
    probe_kmer_freqs: &[u32],
    monomer_median_freq: f64,
) -> f64 {
    if probe_kmer_freqs.is_empty() || monomer_median_freq <= 0.0 {
        return 0.0;
    }

    let n = probe_kmer_freqs.len() as f64;
    let mean_freq: f64 = probe_kmer_freqs.iter().map(|&f| f as f64).sum::<f64>() / n;
    if mean_freq == 0.0 {
        return 0.0;
    }

    let variance: f64 = probe_kmer_freqs
        .iter()
        .map(|&f| {
            let diff = f as f64 - mean_freq;
            diff * diff
        })
        .sum::<f64>()
        / n;
    let std = variance.sqrt();
    let cv = std / mean_freq;

    1.0 / (1.0 + cv)
}

/// Design FISH probes for a single monomer.
/// Returns Vec of ProbeCandidate sorted by composite_score descending.
pub fn design_probes_for_monomer(
    monomer: &str,
    monomer_id: &str,
    index: &KmerIndex,
    k: usize,
    probe_length: usize,
    min_gc: f64,
    max_gc: f64,
    min_tm: f64,
    max_tm: f64,
    top_n: usize,
) -> Vec<ProbeCandidate> {
    let n = monomer.len();
    if n == 0 {
        return Vec::new();
    }

    let monomer_bytes = monomer.as_bytes();

    // Build extended monomer for circular wrapping
    let extended: Vec<u8> = if n < probe_length {
        let repeats_needed = (probe_length / n) + 2;
        let mut ext = Vec::with_capacity(n * repeats_needed);
        for _ in 0..repeats_needed {
            ext.extend_from_slice(monomer_bytes);
        }
        ext
    } else {
        let mut ext = Vec::with_capacity(n + probe_length);
        ext.extend_from_slice(monomer_bytes);
        ext.extend_from_slice(&monomer_bytes[..probe_length - 1]);
        ext
    };

    // Build circular monomer for k-mer lookups
    let circular: Vec<u8> = {
        let mut c = Vec::with_capacity(n + k);
        c.extend_from_slice(monomer_bytes);
        c.extend_from_slice(&monomer_bytes[..k.min(n).saturating_sub(1)]);
        c
    };

    // Build frequency profile
    let mut all_freqs: Vec<u32> = Vec::with_capacity(n);
    for i in 0..n {
        if i + k <= circular.len() {
            let kmer = &circular[i..i + k];
            let freq = match encode_kmer(kmer) {
                Some(e) => index.get(e),
                None => 0,
            };
            all_freqs.push(freq);
        }
    }

    if all_freqs.is_empty() {
        return Vec::new();
    }

    let mut sorted_freqs = all_freqs.clone();
    sorted_freqs.sort();
    let monomer_median_freq = sorted_freqs[sorted_freqs.len() / 2] as f64;
    let max_freq = *all_freqs.iter().max().unwrap_or(&1) as f64;
    if max_freq == 0.0 {
        return Vec::new();
    }

    // Generate candidate probes
    let n_windows = n.min(extended.len().saturating_sub(probe_length) + 1);
    let mut candidates: Vec<ProbeCandidate> = Vec::new();

    for pos in 0..n_windows {
        if pos + probe_length > extended.len() {
            break;
        }
        let probe_seq = &extended[pos..pos + probe_length];

        let gc = compute_gc_content(probe_seq);
        if gc < min_gc || gc > max_gc {
            continue;
        }

        let tm = estimate_melting_temp(probe_seq);
        if tm < min_tm || tm > max_tm {
            continue;
        }

        let probe_kmers_count = probe_length.saturating_sub(k) + 1;
        if probe_kmers_count == 0 {
            continue;
        }

        // Collect k-mer frequencies within this probe window
        let mut probe_freqs: Vec<u32> = Vec::with_capacity(probe_kmers_count);
        for j in 0..probe_kmers_count {
            let kmer_pos = (pos + j) % n;
            if kmer_pos + k <= circular.len() {
                let kmer = &circular[kmer_pos..kmer_pos + k];
                if kmer.len() == k {
                    let freq = match encode_kmer(kmer) {
                        Some(e) => index.get(e),
                        None => 0,
                    };
                    probe_freqs.push(freq);
                }
            }
        }

        if probe_freqs.is_empty() {
            continue;
        }

        let mean_freq: f64 =
            probe_freqs.iter().map(|&f| f as f64).sum::<f64>() / probe_freqs.len() as f64;
        let freq_score = mean_freq / max_freq;
        let specificity = compute_specificity_score(&probe_freqs, monomer_median_freq);
        let composite = freq_score * specificity;

        candidates.push(ProbeCandidate {
            sequence: String::from_utf8_lossy(probe_seq).to_string(),
            gc_content: gc,
            melting_temp: tm,
            frequency_score: freq_score,
            specificity_score: specificity,
            composite_score: composite,
            source_monomer_id: monomer_id.to_string(),
            position_in_monomer: pos,
        });
    }

    // Sort by composite score descending
    candidates.sort_by(|a, b| {
        b.composite_score
            .partial_cmp(&a.composite_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // Remove overlapping probes (greedy)
    let mut selected: Vec<ProbeCandidate> = Vec::new();
    let mut used_positions: std::collections::HashSet<usize> = std::collections::HashSet::new();

    for probe in candidates {
        let pos = probe.position_in_monomer;
        let probe_positions: std::collections::HashSet<usize> =
            (pos..pos + probe_length).map(|p| p % n).collect();
        if probe_positions.intersection(&used_positions).next().is_some() {
            continue;
        }
        used_positions.extend(&probe_positions);
        selected.push(probe);
        if selected.len() >= top_n {
            break;
        }
    }

    selected
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gc_content() {
        assert!((compute_gc_content(b"ACGT") - 0.5).abs() < 1e-10);
        assert!((compute_gc_content(b"AAAA") - 0.0).abs() < 1e-10);
        assert!((compute_gc_content(b"GCGC") - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_melting_temp_short() {
        // ACGT: A+T=2, G+C=2 => Tm = 2*2 + 4*2 = 12
        assert!((estimate_melting_temp(b"ACGT") - 12.0).abs() < 1e-10);
    }

    #[test]
    fn test_melting_temp_long() {
        // 20-mer all G: gc=20, at=0
        let seq = b"GGGGGGGGGGGGGGGGGGGG";
        let expected = 81.5 + 41.0 * 20.0 / 20.0 - 675.0 / 20.0;
        assert!((estimate_melting_temp(seq) - expected).abs() < 1e-10);
    }

    #[test]
    fn test_specificity_uniform() {
        // All same frequency -> CV=0 -> score=1.0
        let freqs = vec![100, 100, 100, 100];
        assert!((compute_specificity_score(&freqs, 100.0) - 1.0).abs() < 1e-10);
    }
}
