/// Tandem repeat finders ported from tr_finder.py.

use std::collections::{HashMap, HashSet};
use crate::kmer::{encode_kmer, decode_kmer};
use crate::kmer_index::KmerIndex;

/// Bases in order matching Python's alphabet = ["A", "C", "T", "G"]
const ALPHABET: [u8; 4] = [0, 1, 3, 2]; // A=0, C=1, T=3, G=2

/// Result of a TR search for one seed.
pub struct FinderResult {
    pub status: String,
    pub second_status: Option<String>,
    pub next_rid: Option<usize>,
    pub next_i: Option<usize>,
    pub sequence: Option<String>,
}

/// Result of a DFS extension in one direction.
struct DfsResult {
    found_tr: bool,
    path: Vec<u8>,       // base codes in extension order
    status: &'static str,
    backtracks: usize,
}

/// Stack frame for iterative DFS.
struct DfsFrame {
    kmer: u64,
    extensions: Vec<(u32, u8)>, // sorted by tf desc
    ext_index: usize,
}

/// Get valid right-extensions sorted by tf descending.
fn get_right_extensions(prefix: u64, index: &KmerIndex, min_tf: u32, k_mask: u64) -> Vec<(u32, u8)> {
    let mut exts = Vec::with_capacity(4);
    for &base in &ALPHABET {
        let kmer = ((prefix << 2) | (base as u64)) & k_mask;
        let tf = index.get(kmer);
        if tf >= min_tf {
            exts.push((tf, base));
        }
    }
    exts.sort_by(|a, b| b.0.cmp(&a.0));
    exts
}

/// Get valid left-extensions sorted by tf descending.
fn get_left_extensions(suffix: u64, index: &KmerIndex, min_tf: u32, k: usize) -> Vec<(u32, u8)> {
    let mut exts = Vec::with_capacity(4);
    for &base in &ALPHABET {
        let kmer = suffix | ((base as u64) << ((k - 1) * 2));
        let tf = index.get(kmer);
        if tf >= min_tf {
            exts.push((tf, base));
        }
    }
    exts.sort_by(|a, b| b.0.cmp(&a.0));
    exts
}

/// Decode a slice of base codes (0-3) to a DNA string.
fn decode_bases(bases: &[u8]) -> String {
    bases.iter().map(|&b| match b {
        0 => 'A', 1 => 'C', 2 => 'G', 3 => 'T', _ => unreachable!(),
    }).collect()
}

/// DFS right extension looking for a cycle back to start_encoded.
///
/// Uses iterative DFS with explicit stack. Extensions are tried in tf-descending
/// order (best first). Backtracking is bounded by max_backtracks.
fn dfs_extend_right(
    start_encoded: u64,
    index: &KmerIndex,
    _k: usize,
    min_tf: u32,
    max_depth: usize,
    max_backtracks: usize,
    k_mask: u64,
    k_minus_1_mask: u64,
) -> DfsResult {
    let start_prefix = start_encoded & k_minus_1_mask;
    let root_exts = get_right_extensions(start_prefix, index, min_tf, k_mask);
    if root_exts.is_empty() {
        return DfsResult { found_tr: false, path: Vec::new(), status: "zero", backtracks: 0 };
    }

    let mut visited = HashSet::new();
    visited.insert(start_encoded);

    let mut path_bases: Vec<u8> = Vec::new();
    let mut best_path: Vec<u8> = Vec::new();
    let mut best_status: &str = "zero";
    let mut backtracks: usize = 0;

    let mut stack: Vec<DfsFrame> = vec![DfsFrame {
        kmer: start_encoded,
        extensions: root_exts,
        ext_index: 0,
    }];

    while let Some(frame) = stack.last_mut() {
        if frame.ext_index >= frame.extensions.len() {
            // All extensions exhausted at this node — backtrack
            let popped = stack.pop().unwrap();
            if stack.is_empty() {
                // Root exhausted, done
                break;
            }
            visited.remove(&popped.kmer);
            path_bases.pop();
            backtracks += 1;
            if backtracks >= max_backtracks {
                break;
            }
            continue;
        }

        let (_ctf, base) = frame.extensions[frame.ext_index];
        frame.ext_index += 1;

        let parent_prefix = frame.kmer & k_minus_1_mask;
        let new_kmer = ((parent_prefix << 2) | (base as u64)) & k_mask;

        // TR check: cycle back to start
        if new_kmer == start_encoded {
            path_bases.push(base);
            return DfsResult { found_tr: true, path: path_bases, status: "tr", backtracks };
        }

        // Skip if already on current path (avoid non-TR self-loops)
        if visited.contains(&new_kmer) {
            continue;
        }

        // Max depth check
        if path_bases.len() + 1 >= max_depth {
            if path_bases.len() + 1 > best_path.len() {
                best_path = path_bases.clone();
                best_path.push(base);
                best_status = "long";
            }
            continue;
        }

        // Extend into new_kmer
        path_bases.push(base);
        visited.insert(new_kmer);

        let next_prefix = new_kmer & k_minus_1_mask;
        let next_exts = get_right_extensions(next_prefix, index, min_tf, k_mask);

        if next_exts.is_empty() {
            // Dead end — update best path, immediate backtrack
            if path_bases.len() > best_path.len() {
                best_path = path_bases.clone();
                best_status = "zero";
            }
            path_bases.pop();
            visited.remove(&new_kmer);
            backtracks += 1;
            if backtracks >= max_backtracks {
                break;
            }
        } else {
            stack.push(DfsFrame {
                kmer: new_kmer,
                extensions: next_exts,
                ext_index: 0,
            });
        }
    }

    DfsResult { found_tr: false, path: best_path, status: best_status, backtracks }
}

/// DFS left extension looking for a cycle back to start_encoded.
fn dfs_extend_left(
    start_encoded: u64,
    index: &KmerIndex,
    k: usize,
    min_tf: u32,
    max_depth: usize,
    max_backtracks: usize,
) -> DfsResult {
    let start_suffix = start_encoded >> 2;
    let root_exts = get_left_extensions(start_suffix, index, min_tf, k);
    if root_exts.is_empty() {
        return DfsResult { found_tr: false, path: Vec::new(), status: "zero", backtracks: 0 };
    }

    let mut visited = HashSet::new();
    visited.insert(start_encoded);

    let mut path_bases: Vec<u8> = Vec::new();
    let mut best_path: Vec<u8> = Vec::new();
    let mut best_status: &str = "zero";
    let mut backtracks: usize = 0;

    let mut stack: Vec<DfsFrame> = vec![DfsFrame {
        kmer: start_encoded,
        extensions: root_exts,
        ext_index: 0,
    }];

    while let Some(frame) = stack.last_mut() {
        if frame.ext_index >= frame.extensions.len() {
            let popped = stack.pop().unwrap();
            if stack.is_empty() {
                break;
            }
            visited.remove(&popped.kmer);
            path_bases.pop();
            backtracks += 1;
            if backtracks >= max_backtracks {
                break;
            }
            continue;
        }

        let (_ctf, base) = frame.extensions[frame.ext_index];
        frame.ext_index += 1;

        let parent_suffix = frame.kmer >> 2;
        let new_kmer = parent_suffix | ((base as u64) << ((k - 1) * 2));

        if new_kmer == start_encoded {
            path_bases.push(base);
            return DfsResult { found_tr: true, path: path_bases, status: "tr", backtracks };
        }

        if visited.contains(&new_kmer) {
            continue;
        }

        if path_bases.len() + 1 >= max_depth {
            if path_bases.len() + 1 > best_path.len() {
                best_path = path_bases.clone();
                best_path.push(base);
                best_status = "long";
            }
            continue;
        }

        path_bases.push(base);
        visited.insert(new_kmer);

        let next_suffix = new_kmer >> 2;
        let next_exts = get_left_extensions(next_suffix, index, min_tf, k);

        if next_exts.is_empty() {
            if path_bases.len() > best_path.len() {
                best_path = path_bases.clone();
                best_status = "zero";
            }
            path_bases.pop();
            visited.remove(&new_kmer);
            backtracks += 1;
            if backtracks >= max_backtracks {
                break;
            }
        } else {
            stack.push(DfsFrame {
                kmer: new_kmer,
                extensions: next_exts,
                ext_index: 0,
            });
        }
    }

    DfsResult { found_tr: false, path: best_path, status: best_status, backtracks }
}

/// Cache all k-mers along a right extension path.
fn cache_right_path(
    start: u64, path: &[u8], k: usize, rid: usize,
    cache: &mut HashMap<u64, (usize, u8, usize)>,
    k_mask: u64, k_minus_1_mask: u64,
) {
    let mut kmer = start;
    for (i, &base) in path.iter().enumerate() {
        let prefix = kmer & k_minus_1_mask;
        kmer = ((prefix << 2) | (base as u64)) & k_mask;
        cache.entry(kmer).or_insert((rid, 0, k + i + 1));
    }
}

/// Cache all k-mers along a left extension path.
fn cache_left_path(
    start: u64, path: &[u8], k: usize, rid: usize,
    cache: &mut HashMap<u64, (usize, u8, usize)>,
) {
    let mut kmer = start;
    for (i, &base) in path.iter().enumerate() {
        let suffix = kmer >> 2;
        kmer = suffix | ((base as u64) << ((k - 1) * 2));
        cache.entry(kmer).or_insert((rid, 0, k + i + 1));
    }
}

/// Bidirectional TR finder with DFS backtracking.
///
/// For each seed k-mer (from sdat, sorted by tf desc):
/// 1. Try DFS right extension to find a cycle back to seed
/// 2. If TR found right → done
/// 3. Otherwise try DFS left extension with remaining backtrack budget
/// 4. If TR found left → done
/// 5. Otherwise report best non-TR path with status zero/long/extended
pub fn tr_greedy_finder_bidirectional(
    sdat: &[(String, u32)],
    max_depth: usize,
    coverage: f64,
    _min_fraction_to_continue: u32,
    k: usize,
    lu: Option<u32>,
    max_backtracks: usize,
) -> Vec<FinderResult> {
    let min_tf: u32 = match lu {
        Some(val) => if val <= 1 { 2 } else { val },
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

    let k_mask = if k == 32 { u64::MAX } else { (1u64 << (k * 2)) - 1 };
    let k_minus_1_mask = if k - 1 == 32 { u64::MAX } else { (1u64 << ((k - 1) * 2)) - 1 };

    let mut repeats: Vec<FinderResult> = Vec::new();
    let mut rid: usize = 0;
    let mut cache: HashMap<u64, (usize, u8, usize)> = HashMap::new();

    for &(start_encoded, tf) in &encoded_sdat {
        if tf < min_tf {
            break;
        }
        if cache.contains_key(&start_encoded) {
            continue;
        }

        cache.insert(start_encoded, (rid, 0, k));

        // === Right DFS ===
        let right_result = dfs_extend_right(
            start_encoded, &index, k, min_tf, max_depth, max_backtracks,
            k_mask, k_minus_1_mask,
        );

        if right_result.found_tr {
            // Cache k-mers on the cycle path
            cache_right_path(start_encoded, &right_result.path, k, rid,
                             &mut cache, k_mask, k_minus_1_mask);
            let start_str = decode_kmer(start_encoded, k);
            let path_str = decode_bases(&right_result.path);
            repeats.push(FinderResult {
                status: "tr".to_string(),
                second_status: None,
                next_rid: None,
                next_i: None,
                sequence: Some(format!("{}{}", start_str, path_str)),
            });
        } else {
            // === Left DFS with remaining budget ===
            let remaining_bt = max_backtracks.saturating_sub(right_result.backtracks);
            let left_result = dfs_extend_left(
                start_encoded, &index, k, min_tf, max_depth, remaining_bt,
            );

            if left_result.found_tr {
                cache_left_path(start_encoded, &left_result.path, k, rid, &mut cache);
                let mut left_bases = left_result.path;
                left_bases.reverse();
                let left_str = decode_bases(&left_bases);
                let start_str = decode_kmer(start_encoded, k);
                repeats.push(FinderResult {
                    status: "tr".to_string(),
                    second_status: None,
                    next_rid: None,
                    next_i: None,
                    sequence: Some(format!("{}{}", left_str, start_str)),
                });
            } else {
                // No TR found — report best non-TR path
                let r_status = right_result.status;
                let l_status = left_result.status;
                let final_status = if r_status == "zero" && l_status == "zero" {
                    "zero"
                } else if r_status == "long" || l_status == "long" {
                    "long"
                } else {
                    "extended"
                };

                // Cache both paths
                cache_right_path(start_encoded, &right_result.path, k, rid,
                                 &mut cache, k_mask, k_minus_1_mask);
                cache_left_path(start_encoded, &left_result.path, k, rid, &mut cache);

                // Build sequence: reversed(left) + start + right
                let mut left_bases = left_result.path;
                left_bases.reverse();
                let left_str = decode_bases(&left_bases);
                let start_str = decode_kmer(start_encoded, k);
                let right_str = decode_bases(&right_result.path);
                let full_seq = format!("{}{}{}", left_str, start_str, right_str);
                let seq_out = if full_seq.is_empty() { None } else { Some(full_seq) };

                repeats.push(FinderResult {
                    status: final_status.to_string(),
                    second_status: None,
                    next_rid: None,
                    next_i: None,
                    sequence: seq_out,
                });
            }
        }

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

        let status: Option<&str>;
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
        let ext_str = decode_bases(&seq);
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
        let results = tr_greedy_finder_bidirectional(&sdat, 30000, 1.0, 30, 3, Some(100), 1000);
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

    #[test]
    fn test_backtracking_finds_tr_at_fork() {
        // Monomer "ATCGA" (5bp) with k=3.
        // Cycle: ATC→TCG→CGA→GAA→AAT→ATC (all tf=100)
        // Dead-end fork: TCT has tf=200 (higher than TCG=100)
        // Greedy from ATC would pick TCT and fail.
        // DFS backtracks from the TCT dead end and finds the cycle via TCG.
        let sdat = vec![
            ("TCT".to_string(), 200),  // dead-end fork (highest tf)
            ("ATC".to_string(), 100),
            ("TCG".to_string(), 100),
            ("CGA".to_string(), 100),
            ("GAA".to_string(), 100),
            ("AAT".to_string(), 100),
        ];
        let results = tr_greedy_finder_bidirectional(&sdat, 30000, 1.0, 30, 3, Some(50), 1000);

        // TCT is processed first as seed (highest tf). No cycle → not TR.
        // ATC is processed next. DFS right: tries TCT(200) first → dead end →
        // backtracks → tries TCG(100) → CGA → GAA → AAT → ATC = TR!
        let tr_results: Vec<_> = results.iter().filter(|r| r.status == "tr").collect();
        assert!(!tr_results.is_empty(), "DFS should find TR via backtracking");
    }

    #[test]
    fn test_backtracking_with_multiple_forks() {
        // Monomer "ATCGTA" (6bp) with k=3 and two forks
        // Cycle: ATC→TCG→CGT→GTA→TAA→AAT→ATC (all tf=100)
        // Fork 1: TCT tf=200 (dead end)
        // Fork 2: GTG tf=150 (dead end)
        let sdat = vec![
            ("TCT".to_string(), 200),
            ("GTG".to_string(), 150),
            ("ATC".to_string(), 100),
            ("TCG".to_string(), 100),
            ("CGT".to_string(), 100),
            ("GTA".to_string(), 100),
            ("TAA".to_string(), 100),
            ("AAT".to_string(), 100),
        ];
        let results = tr_greedy_finder_bidirectional(&sdat, 30000, 1.0, 30, 3, Some(50), 1000);
        let tr_results: Vec<_> = results.iter().filter(|r| r.status == "tr").collect();
        assert!(!tr_results.is_empty(), "DFS should find TR despite two dead-end forks");
    }

    #[test]
    fn test_no_backtrack_needed_for_simple_repeat() {
        // Simple repeat: no forks, greedy path = correct path
        let sdat = make_circular_sdat("AATGG", 1000, 3);
        let results = tr_greedy_finder_bidirectional(&sdat, 30000, 1.0, 30, 3, Some(100), 0);
        // Even with max_backtracks=0 (pure greedy), should find TR on first try
        assert!(!results.is_empty());
        assert_eq!(results[0].status, "tr");
    }
}
