/// Tandem repeat finders ported from tr_finder.py.

use std::collections::HashMap;
use crate::kmer::{encode_kmer, decode_kmer};
use crate::kmer_index::KmerIndex;

/// Bases in order matching Python's alphabet = ["A", "C", "T", "G"]
const ALPHABET: [u8; 4] = [0, 1, 3, 2]; // A=0, C=1, T=3, G=2

/// Result of a greedy TR search for one seed.
pub struct FinderResult {
    pub status: String,
    pub second_status: Option<String>,
    pub next_rid: Option<usize>,
    pub next_i: Option<usize>,
    pub sequence: Option<String>,
}

/// Bidirectional greedy TR finder (mirrors Python lines 173-351).
pub fn tr_greedy_finder_bidirectional(
    sdat: &[(String, u32)],
    max_depth: usize,
    coverage: f64,
    _min_fraction_to_continue: u32,
    k: usize,
    lu: Option<u32>,
) -> Vec<FinderResult> {
    let min_tf: u32 = match lu {
        Some(val) => {
            let v = val;
            if v <= 1 { 2 } else { v }
        }
        None => {
            let v = (coverage * 100.0) as u32;
            if v <= 1 { 2 } else { v }
        }
    };

    // Build encoded index from sdat
    let index = KmerIndex::from_sdat(sdat, k);

    // Pre-encode sdat k-mers for iteration
    let encoded_sdat: Vec<(u64, u32)> = sdat
        .iter()
        .filter_map(|(s, tf)| encode_kmer(s.as_bytes()).map(|e| (e, *tf)))
        .collect();

    let mut repeats: Vec<FinderResult> = Vec::new();
    let mut rid: usize = 0;

    // cache: encoded_kmer -> (rid, strand, length)
    let mut cache: HashMap<u64, (usize, u8, usize)> = HashMap::new();

    for &(start_encoded, tf) in &encoded_sdat {
        if tf < min_tf {
            break;
        }
        if cache.contains_key(&start_encoded) {
            continue;
        }

        let mut final_status: Option<&str> = None;
        let mut second_status: Option<String> = None;
        let mut next_rid: Option<usize> = None;
        let mut next_i: Option<usize> = None;
        let mut total_length = k;

        cache.insert(start_encoded, (rid, 0, total_length));

        // === Right extension ===
        // Compute prefix: drop leftmost base of start_kmer (shift left, mask to k-1 bases)
        let k_minus_1_mask = if k - 1 == 32 {
            u64::MAX
        } else {
            (1u64 << ((k - 1) * 2)) - 1
        };
        let mut right_prefix = start_encoded & k_minus_1_mask;
        let mut right_seq: Vec<u8> = Vec::new(); // stores base codes (0-3)
        let mut right_length: usize = 0;
        let mut right_status: Option<&str> = None;

        while right_length < max_depth {
            // Try all 4 bases, collect valid extensions
            let mut best_tf: u32 = 0;
            let mut best_base: u8 = 0;
            let mut has_solution = false;

            for &base in &ALPHABET {
                let new_kmer_full = ((right_prefix << 2) | (base as u64))
                    & if k == 32 { u64::MAX } else { (1u64 << (k * 2)) - 1 };
                let ctf = index.get(new_kmer_full);
                if ctf < min_tf {
                    continue;
                }
                if !has_solution || ctf > best_tf {
                    best_tf = ctf;
                    best_base = base;
                    has_solution = true;
                }
            }

            if !has_solution {
                right_status = Some("zero");
                break;
            }

            // Python sorts ascending and takes [-1] (max), but solutions.sort(reverse=True) and [0]
            // is equivalent — we already picked max above.
            let new_kmer_full = ((right_prefix << 2) | (best_base as u64))
                & if k == 32 { u64::MAX } else { (1u64 << (k * 2)) - 1 };

            right_seq.push(best_base);
            total_length += 1;

            if new_kmer_full == start_encoded {
                right_status = Some("tr");
                break;
            }

            if let Some(&(existing_rid, _strand, i)) = cache.get(&new_kmer_full) {
                // Bug-compatible: `existing_rid not in repeats` where repeats is Vec<tuple>
                // int-in-list-of-tuples is always False, so second_status is always "self"
                second_status = Some("self".to_string());
                next_rid = Some(existing_rid);
                next_i = Some(i);
                right_status = Some("frag");
                break;
            }

            cache.insert(new_kmer_full, (rid, 0, total_length));
            right_prefix = new_kmer_full & k_minus_1_mask;
            right_length += 1;
        }

        if right_status.is_none() {
            right_status = Some("long");
        }

        let full_seq: Option<String>;

        if right_status == Some("tr") {
            final_status = Some("tr");
            // Reconstruct: start_kmer + right_seq bases
            let start_str = decode_kmer(start_encoded, k);
            let right_str: String = right_seq
                .iter()
                .map(|&b| match b {
                    0 => 'A',
                    1 => 'C',
                    2 => 'G',
                    3 => 'T',
                    _ => unreachable!(),
                })
                .collect();
            full_seq = Some(format!("{}{}", start_str, right_str));
        } else {
            // === Left extension ===
            // suffix = start_kmer[:-1] = start_encoded >> 2 (upper k-1 bases)
            let mut left_suffix = start_encoded >> 2;
            let mut left_seq: Vec<u8> = Vec::new();
            let mut left_length: usize = 0;
            let mut left_status: Option<&str> = None;

            while left_length < max_depth {
                let mut best_tf: u32 = 0;
                let mut best_base: u8 = 0;
                let mut has_solution = false;

                for &base in &ALPHABET {
                    // new_kmer = base + left_suffix (k-1 bases) = base at position k-1 (high bits)
                    let new_kmer_full =
                        left_suffix | ((base as u64) << ((k - 1) * 2));
                    let ctf = index.get(new_kmer_full);
                    if ctf < min_tf {
                        continue;
                    }
                    if !has_solution || ctf > best_tf {
                        best_tf = ctf;
                        best_base = base;
                        has_solution = true;
                    }
                }

                if !has_solution {
                    left_status = Some("zero");
                    break;
                }

                let new_kmer_full =
                    left_suffix | ((best_base as u64) << ((k - 1) * 2));

                left_seq.push(best_base);
                total_length += 1;

                if new_kmer_full == start_encoded {
                    left_status = Some("tr");
                    break;
                }

                if let Some(&(existing_rid, _strand, i)) = cache.get(&new_kmer_full) {
                    second_status = Some("self".to_string());
                    next_rid = Some(existing_rid);
                    next_i = Some(i);
                    left_status = Some("frag");
                    break;
                }

                cache.insert(new_kmer_full, (rid, 0, total_length));
                // new suffix = new_kmer[:-1] = new_kmer >> 2
                left_suffix = new_kmer_full >> 2;
                left_length += 1;
            }

            if left_status.is_none() {
                left_status = Some("long");
            }

            // Determine overall status
            let rs = right_status.unwrap();
            let ls = left_status.unwrap();
            if ls == "tr" || rs == "tr" {
                final_status = Some("tr");
            } else if ls == "frag" || rs == "frag" {
                final_status = Some("frag");
            } else if ls == "zero" && rs == "zero" {
                final_status = Some("zero");
            } else {
                final_status = Some("extended");
            }

            // Reconstruct sequence: reversed(left_seq) + start_kmer + right_seq
            let start_str = decode_kmer(start_encoded, k);
            let left_str: String = left_seq
                .iter()
                .rev()
                .map(|&b| match b {
                    0 => 'A',
                    1 => 'C',
                    2 => 'G',
                    3 => 'T',
                    _ => unreachable!(),
                })
                .collect();
            let right_str: String = right_seq
                .iter()
                .map(|&b| match b {
                    0 => 'A',
                    1 => 'C',
                    2 => 'G',
                    3 => 'T',
                    _ => unreachable!(),
                })
                .collect();
            full_seq = Some(format!("{}{}{}", left_str, start_str, right_str));
        }

        let seq_out = match &full_seq {
            Some(s) if s.is_empty() => None,
            other => other.clone(),
        };

        repeats.push(FinderResult {
            status: final_status.unwrap_or("zero").to_string(),
            second_status,
            next_rid,
            next_i,
            sequence: seq_out,
        });
        rid += 1;
    }

    repeats
}

/// Unidirectional (right-only) greedy TR finder (mirrors Python lines 99-171).
pub fn tr_greedy_finder(
    sdat: &[(String, u32)],
    max_depth: usize,
    coverage: f64,
    _min_fraction_to_continue: u32,
    k: usize,
    lu: Option<u32>,
) -> Vec<FinderResult> {
    let min_tf: u32 = match lu {
        Some(val) => {
            if val <= 1 { 2 } else { val }
        }
        None => {
            let v = (coverage * 100.0) as u32;
            if v <= 1 { 2 } else { v }
        }
    };

    let index = KmerIndex::from_sdat(sdat, k);

    let encoded_sdat: Vec<(u64, u32)> = sdat
        .iter()
        .filter_map(|(s, tf)| encode_kmer(s.as_bytes()).map(|e| (e, *tf)))
        .collect();

    let mut repeats: Vec<FinderResult> = Vec::new();
    let mut rid: usize = 0;
    let mut cache: HashMap<u64, (usize, u8, usize)> = HashMap::new();

    let k_minus_1_mask = if k - 1 == 32 {
        u64::MAX
    } else {
        (1u64 << ((k - 1) * 2)) - 1
    };

    for &(start_encoded, tf) in &encoded_sdat {
        if tf < min_tf {
            break;
        }
        if cache.contains_key(&start_encoded) {
            continue;
        }

        let mut status: Option<&str> = None;
        let mut second_status: Option<String> = None;
        let mut next_rid: Option<usize> = None;
        let mut next_i: Option<usize> = None;
        let mut length = k;

        let mut seq: Vec<u8> = Vec::new(); // base codes for extensions
        let mut prefix = start_encoded & k_minus_1_mask;

        cache.insert(start_encoded, (rid, 0, length));

        loop {
            let mut best_tf: u32 = 0;
            let mut best_base: u8 = 0;
            let mut has_solution = false;

            for &base in &ALPHABET {
                let new_kmer = ((prefix << 2) | (base as u64))
                    & if k == 32 { u64::MAX } else { (1u64 << (k * 2)) - 1 };
                let ctf = index.get(new_kmer);
                if ctf < min_tf {
                    continue;
                }
                if !has_solution || ctf > best_tf {
                    best_tf = ctf;
                    best_base = base;
                    has_solution = true;
                }
            }

            if !has_solution {
                status = Some("zero");
                break;
            }

            let kmer = ((prefix << 2) | (best_base as u64))
                & if k == 32 { u64::MAX } else { (1u64 << (k * 2)) - 1 };
            seq.push(best_base);

            if kmer == start_encoded {
                status = Some("tr");
                break;
            }

            prefix = kmer & k_minus_1_mask;
            length += 1;

            if let Some(&(existing_rid, _strand, i)) = cache.get(&kmer) {
                // Bug-compatible: always "self"
                second_status = Some("self".to_string());
                next_rid = Some(existing_rid);
                next_i = Some(i);
                status = Some("frag");
                break;
            }

            cache.insert(kmer, (rid, 0, length));

            if length == max_depth {
                status = Some("long");
                break;
            }
        }

        // Build sequence: start_kmer + extension bases
        let start_str = decode_kmer(start_encoded, k);
        let ext_str: String = seq
            .iter()
            .map(|&b| match b {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                3 => 'T',
                _ => unreachable!(),
            })
            .collect();
        let full_seq = format!("{}{}", start_str, ext_str);

        repeats.push(FinderResult {
            status: status.unwrap_or("zero").to_string(),
            second_status,
            next_rid,
            next_i,
            sequence: Some(full_seq),
        });
        rid += 1;
    }

    repeats
}

/// Naive TR finder — returns only repeats list (other return values unused in benchmark).
/// Mirrors Python lines 13-97.
pub fn naive_tr_finder(
    sdat: &[(String, u32)],
    min_tf_extension: u32,
    min_fraction_to_continue: u32,
    k: usize,
) -> Vec<(usize, String, usize, String, String)> {
    // (rid, status, length, start_kmer, sequence)
    use std::collections::HashSet;
    use crate::kmer::revcomp_string;

    let mut repeats: Vec<(usize, String, usize, String, String)> = Vec::new();
    let mut used_kmers: HashSet<String> = HashSet::new();
    let mut rid: usize = 0;

    // Build a simple string->u32 map from sdat for lookups
    let mut kmer2tf: HashMap<String, u32> = HashMap::new();
    for (kmer, tf) in sdat {
        kmer2tf.insert(kmer.clone(), *tf);
    }

    let alphabet = ['A', 'C', 'T', 'G'];

    for (start_kmer, _tf) in sdat {
        let mut step: usize = 0;
        if used_kmers.contains(start_kmer) {
            continue;
        }
        used_kmers.insert(start_kmer.clone());
        used_kmers.insert(revcomp_string(start_kmer));

        let mut seq = start_kmer.clone();

        loop {
            let kmer = &seq[seq.len() - k..];
            let ori_tf = *kmer2tf.get(kmer).unwrap_or(&0);

            let mut solutions: Vec<(u32, char)> = Vec::new();
            let prefix = &kmer[1..];
            for &nucleotide in &alphabet {
                let candidate = format!("{}{}", prefix, nucleotide);
                let ctf = *kmer2tf.get(&candidate).unwrap_or(&0);
                if ctf < min_tf_extension {
                    continue;
                }
                let fr = (100 * ctf as u64 / (ori_tf as u64 + 1)) as u32;
                solutions.push((fr, nucleotide));
            }
            solutions.sort();

            if !solutions.is_empty()
                && solutions.last().unwrap().0 >= min_fraction_to_continue
            {
                let best = solutions.last().unwrap().1;
                seq.push(best);
                step += 1;

                let last_kmer = &seq[seq.len() - k..];

                if used_kmers.contains(last_kmer) && last_kmer != start_kmer.as_str() {
                    // FRAG
                    for iii in 0..=(seq.len() - k) {
                        let km = &seq[iii..iii + k];
                        used_kmers.insert(km.to_string());
                        used_kmers.insert(revcomp_string(km));
                    }
                    repeats.push((
                        rid,
                        "FRAG".to_string(),
                        seq.len(),
                        start_kmer.clone(),
                        seq.clone(),
                    ));
                    rid += 1;
                    break;
                }

                used_kmers.insert(last_kmer.to_string());
                used_kmers.insert(revcomp_string(last_kmer));

                if last_kmer == start_kmer.as_str() {
                    // TR found
                    for iii in 0..=(seq.len() - k) {
                        let km = &seq[iii..iii + k];
                        used_kmers.insert(km.to_string());
                        used_kmers.insert(revcomp_string(km));
                    }
                    repeats.push((
                        rid,
                        "TR".to_string(),
                        step,
                        start_kmer.clone(),
                        seq[..step].to_string(),
                    ));
                    rid += 1;
                    break;
                }

                continue;
            }

            // Terminal — no valid extension
            for iii in 0..=(seq.len() - k) {
                let km = &seq[iii..iii + k];
                used_kmers.insert(km.to_string());
                used_kmers.insert(revcomp_string(km));
            }
            repeats.push((
                rid,
                "TE".to_string(),
                seq.len(),
                start_kmer.clone(),
                seq.clone(),
            ));
            rid += 1;
            break;
        }
    }

    repeats
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_circular_sdat(monomer: &str, tf: u32, k: usize) -> Vec<(String, u32)> {
        // Generate all k-mers from a circular monomer, all with the same tf
        let extended = format!("{}{}", monomer, &monomer[..k - 1]);
        let mut sdat = Vec::new();
        for i in 0..monomer.len() {
            sdat.push((extended[i..i + k].to_string(), tf));
        }
        sdat
    }

    #[test]
    fn test_bidirectional_finds_tr() {
        let sdat = make_circular_sdat("AATGG", 1000, 3);
        let results = tr_greedy_finder_bidirectional(&sdat, 30000, 1.0, 30, 3, Some(100));
        // The first seed should find a TR
        assert!(!results.is_empty());
        assert_eq!(results[0].status, "tr");
    }

    #[test]
    fn test_unidirectional_finds_tr() {
        let sdat = make_circular_sdat("AATGG", 1000, 3);
        let results = tr_greedy_finder(&sdat, 30000, 1.0, 30, 3, Some(100));
        assert!(!results.is_empty());
        assert_eq!(results[0].status, "tr");
    }
}
