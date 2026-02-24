/// Tandem repeat finders ported from tr_finder.py.
///
/// Key design: greedy extension (highest tf first) with DFS backtracking at
/// forks. The global cache is checked and written during DFS to match the
/// Python code's behavior: each k-mer is explored at most once across all seeds.

use std::collections::HashMap;
use crate::kmer::{encode_kmer, decode_kmer};
use crate::kmer_index::KmerIndex;

/// Bases in order matching Python's alphabet = ["A", "C", "T", "G"]
const ALPHABET: [u8; 4] = [0, 1, 3, 2]; // A=0, C=1, T=3, G=2

/// Hard limit on total DFS extension attempts per seed direction.
const MAX_DFS_STEPS: usize = 100_000;

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
}

/// Stack frame for iterative DFS.
struct DfsFrame {
    kmer: u64,
    extensions: Vec<(u32, u8)>,  // sorted by descending tf
    ext_index: usize,
}

/// Get valid right-extensions sorted by descending tf (highest first = greedy).
fn get_right_extensions(prefix: u64, index: &KmerIndex, min_tf: u32, k_mask: u64) -> Vec<(u32, u8)> {
    let mut exts = Vec::with_capacity(4);
    for &base in &ALPHABET {
        let kmer = ((prefix << 2) | (base as u64)) & k_mask;
        let tf = index.get(kmer);
        if tf >= min_tf {
            exts.push((tf, base));
        }
    }
    // Highest tf first — greedy, follow the strongest signal
    exts.sort_by(|a, b| b.0.cmp(&a.0));
    exts
}

/// Get valid left-extensions sorted by descending tf (highest first = greedy).
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

/// DFS right extension with greedy ordering and global cache interaction.
///
/// At each node, extensions are tried in descending tf order (greedy = follow
/// strongest signal). If the greedy choice hits a dead end or cached k-mer,
/// the DFS backtracks and tries the next-best extension.
///
/// The global cache is checked and updated during DFS:
/// - Before extending into a k-mer, check if it's in cache → skip (like Python's "frag" stop)
/// - After extending into a k-mer, add it to cache → prevents re-exploration by future seeds
/// - On backtrack, k-mers are NOT removed from cache (once explored, always cached)
fn dfs_extend_right(
    start_encoded: u64,
    index: &KmerIndex,
    k: usize,
    min_tf: u32,
    max_depth: usize,
    max_backtracks: usize,
    k_mask: u64,
    k_minus_1_mask: u64,
    cache: &mut HashMap<u64, (usize, u8, usize)>,
    rid: usize,
) -> DfsResult {
    let start_prefix = start_encoded & k_minus_1_mask;
    let root_exts = get_right_extensions(start_prefix, index, min_tf, k_mask);
    if root_exts.is_empty() {
        return DfsResult { found_tr: false, path: Vec::new(), status: "zero" };
    }

    let mut path_bases: Vec<u8> = Vec::new();
    let mut best_path: Vec<u8> = Vec::new();
    let mut best_status: &str = "zero";
    let mut backtracks: usize = 0;
    let mut steps: usize = 0;

    let mut stack: Vec<DfsFrame> = vec![DfsFrame {
        kmer: start_encoded,
        extensions: root_exts,
        ext_index: 0,
    }];

    while let Some(frame) = stack.last_mut() {
        if frame.ext_index >= frame.extensions.len() {
            // All extensions exhausted — unwind
            stack.pop();
            if stack.is_empty() {
                break;
            }
            // Don't remove from cache — once explored, always cached
            path_bases.pop();
            continue;
        }

        let ext_idx = frame.ext_index;
        let (_ext_tf, base) = frame.extensions[ext_idx];
        frame.ext_index += 1;

        // Count fork decisions: trying 2nd+ extension at a node
        if ext_idx > 0 {
            backtracks += 1;
            if backtracks >= max_backtracks {
                break;
            }
        }

        steps += 1;
        if steps >= MAX_DFS_STEPS {
            break;
        }

        let parent_prefix = frame.kmer & k_minus_1_mask;
        let new_kmer = ((parent_prefix << 2) | (base as u64)) & k_mask;

        // TR check: cycle back to start
        if new_kmer == start_encoded {
            path_bases.push(base);
            return DfsResult { found_tr: true, path: path_bases, status: "tr" };
        }

        // Skip if already in global cache (explored by this or previous seed)
        if cache.contains_key(&new_kmer) {
            if path_bases.len() + 1 > best_path.len() {
                best_path = path_bases.clone();
                best_path.push(base);
                best_status = "frag";
            }
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

        // Extend into new_kmer — add to global cache immediately
        path_bases.push(base);
        cache.insert(new_kmer, (rid, 0, k + path_bases.len()));

        let next_prefix = new_kmer & k_minus_1_mask;
        let next_exts = get_right_extensions(next_prefix, index, min_tf, k_mask);

        if next_exts.is_empty() {
            // Dead end
            if path_bases.len() > best_path.len() {
                best_path = path_bases.clone();
                best_status = "zero";
            }
            path_bases.pop();
            continue;
        }

        stack.push(DfsFrame {
            kmer: new_kmer,
            extensions: next_exts,
            ext_index: 0,
        });
    }

    DfsResult { found_tr: false, path: best_path, status: best_status }
}

/// DFS left extension with greedy ordering and global cache interaction.
fn dfs_extend_left(
    start_encoded: u64,
    index: &KmerIndex,
    k: usize,
    min_tf: u32,
    max_depth: usize,
    max_backtracks: usize,
    cache: &mut HashMap<u64, (usize, u8, usize)>,
    rid: usize,
) -> DfsResult {
    let start_suffix = start_encoded >> 2;
    let root_exts = get_left_extensions(start_suffix, index, min_tf, k);
    if root_exts.is_empty() {
        return DfsResult { found_tr: false, path: Vec::new(), status: "zero" };
    }

    let mut path_bases: Vec<u8> = Vec::new();
    let mut best_path: Vec<u8> = Vec::new();
    let mut best_status: &str = "zero";
    let mut backtracks: usize = 0;
    let mut steps: usize = 0;

    let mut stack: Vec<DfsFrame> = vec![DfsFrame {
        kmer: start_encoded,
        extensions: root_exts,
        ext_index: 0,
    }];

    while let Some(frame) = stack.last_mut() {
        if frame.ext_index >= frame.extensions.len() {
            stack.pop();
            if stack.is_empty() {
                break;
            }
            path_bases.pop();
            continue;
        }

        let ext_idx = frame.ext_index;
        let (_ext_tf, base) = frame.extensions[ext_idx];
        frame.ext_index += 1;

        if ext_idx > 0 {
            backtracks += 1;
            if backtracks >= max_backtracks {
                break;
            }
        }

        steps += 1;
        if steps >= MAX_DFS_STEPS {
            break;
        }

        let parent_suffix = frame.kmer >> 2;
        let new_kmer = parent_suffix | ((base as u64) << ((k - 1) * 2));

        if new_kmer == start_encoded {
            path_bases.push(base);
            return DfsResult { found_tr: true, path: path_bases, status: "tr" };
        }

        if cache.contains_key(&new_kmer) {
            if path_bases.len() + 1 > best_path.len() {
                best_path = path_bases.clone();
                best_path.push(base);
                best_status = "frag";
            }
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
        cache.insert(new_kmer, (rid, 0, k + path_bases.len()));

        let next_suffix = new_kmer >> 2;
        let next_exts = get_left_extensions(next_suffix, index, min_tf, k);

        if next_exts.is_empty() {
            if path_bases.len() > best_path.len() {
                best_path = path_bases.clone();
                best_status = "zero";
            }
            path_bases.pop();
            continue;
        }

        stack.push(DfsFrame {
            kmer: new_kmer,
            extensions: next_exts,
            ext_index: 0,
        });
    }

    DfsResult { found_tr: false, path: best_path, status: best_status }
}

/// Bidirectional TR finder: greedy with DFS backtracking at forks.
///
/// Matches the Python bidirectional finder's behavior:
/// 1. Extensions sorted by descending tf (greedy = follow strongest signal)
/// 2. Global cache checked during extension (stop at previously-explored territory)
/// 3. New k-mers added to global cache during extension (once explored, always cached)
///
/// Two-tier threshold:
/// - `lu` (seed_min_tf): high threshold for selecting starting k-mers from sdat
/// - `ext_lu` (ext_min_tf): low threshold for DFS extension filtering
pub fn tr_greedy_finder_bidirectional(
    sdat: &[(String, u32)],
    max_depth: usize,
    coverage: f64,
    _min_fraction_to_continue: u32,
    k: usize,
    lu: Option<u32>,
    max_backtracks: usize,
    ext_lu: Option<u32>,
) -> Vec<FinderResult> {
    let seed_min_tf: u32 = match lu {
        Some(val) => if val <= 1 { 2 } else { val },
        None => {
            let v = (coverage * 100.0) as u32;
            if v <= 1 { 2 } else { v }
        }
    };

    let ext_min_tf: u32 = match ext_lu {
        Some(val) => if val <= 1 { 2 } else { val },
        None => seed_min_tf,
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
        if tf < seed_min_tf {
            break;
        }
        if cache.contains_key(&start_encoded) {
            continue;
        }

        cache.insert(start_encoded, (rid, 0, k));

        // === Right extension (greedy + DFS backtracking) ===
        let right_result = dfs_extend_right(
            start_encoded, &index, k, ext_min_tf, max_depth, max_backtracks,
            k_mask, k_minus_1_mask, &mut cache, rid,
        );

        if right_result.found_tr {
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
            // === Left extension with remaining budget ===
            let left_result = dfs_extend_left(
                start_encoded, &index, k, ext_min_tf, max_depth, max_backtracks,
                &mut cache, rid,
            );

            if left_result.found_tr {
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
                // No TR found — determine status
                let r_status = right_result.status;
                let l_status = left_result.status;
                let final_status = if r_status == "frag" || l_status == "frag" {
                    "frag"
                } else if r_status == "zero" && l_status == "zero" {
                    "zero"
                } else if r_status == "long" || l_status == "long" {
                    "long"
                } else {
                    "extended"
                };

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

        let mut seq: Vec<u8> = Vec::new();
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
    use std::collections::HashSet;
    use crate::kmer::revcomp_string;

    let mut repeats: Vec<(usize, String, usize, String, String)> = Vec::new();
    let mut used_kmers: HashSet<String> = HashSet::new();
    let mut rid: usize = 0;

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
        let results = tr_greedy_finder_bidirectional(&sdat, 30000, 1.0, 30, 3, Some(100), 1000, None);
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
        // Cycle k-mers at tf=1000 (high, processed as seeds first).
        // Dead-end fork: TCT has tf=200 (lower than cycle, but still above min_tf).
        // When DFS from ATC seed tries right extensions:
        //   TCT(200) and TCG(1000) are both valid.
        //   Greedy picks TCG(1000) first → finds cycle immediately.
        // DFS with max_backtracks=0 also works since greedy first choice is correct.
        let sdat = vec![
            ("ATC".to_string(), 1000),
            ("TCG".to_string(), 1000),
            ("CGA".to_string(), 1000),
            ("GAA".to_string(), 1000),
            ("AAT".to_string(), 1000),
            ("TCT".to_string(), 200),  // dead-end fork (lower tf than cycle)
        ];
        let results = tr_greedy_finder_bidirectional(&sdat, 30000, 1.0, 30, 3, Some(50), 1000, None);
        let tr_results: Vec<_> = results.iter().filter(|r| r.status == "tr").collect();
        assert!(!tr_results.is_empty(), "Should find TR — greedy picks correct extension first");
    }

    #[test]
    fn test_backtracking_from_dead_end() {
        // Cycle ATCGA at tf=1000. Dead-end fork TCT at tf=2000 (higher tf).
        // From ATC seed, greedy tries TCT(2000) first → dead end.
        // DFS backtracks and tries TCG(1000) → finds cycle.
        // Key: ATC is the first seed (tf=1000 in sdat), TCT is NOT a separate
        // seed because it only appears as an extension candidate (not in sdat
        // at higher position than ATC).
        let sdat = vec![
            ("ATC".to_string(), 1000),
            ("TCG".to_string(), 1000),
            ("CGA".to_string(), 1000),
            ("GAA".to_string(), 1000),
            ("AAT".to_string(), 1000),
            ("TCT".to_string(), 2000),  // dead-end: higher tf but lower position in sdat
        ];
        // Note: sdat must be sorted by tf descending for realistic behavior
        let mut sorted = sdat;
        sorted.sort_by(|a, b| b.1.cmp(&a.1));
        // Now TCT(2000) is first seed. But TCT has no right extensions and
        // its left extension is ATC which will be cached.
        // To properly test backtracking, we need TCT to NOT be a seed.
        // Set seed_lu=1500 so only TCT(2000) is above seed threshold but
        // ext_lu=50 so all k-mers are available for extension.
        // Wait - that means only TCT is a seed, which doesn't help.
        //
        // Better approach: put dead-end inside the cycle's tf range but
        // make it NOT a separate k-mer in sdat.
        // Actually, the simplest test: all cycle k-mers at same tf,
        // dead-end also at same tf. First seed is cycle k-mer (deterministic
        // because sdat order is preserved).
        let sdat2 = vec![
            ("ATC".to_string(), 1000),  // first in sdat → first seed
            ("TCG".to_string(), 1000),
            ("CGA".to_string(), 1000),
            ("GAA".to_string(), 1000),
            ("AAT".to_string(), 1000),
            ("TCT".to_string(), 1000),  // dead-end, same tf, later in sdat
        ];
        let results = tr_greedy_finder_bidirectional(&sdat2, 30000, 1.0, 30, 3, Some(50), 1000, None);
        let tr_results: Vec<_> = results.iter().filter(|r| r.status == "tr").collect();
        assert!(!tr_results.is_empty(), "DFS should find TR even with dead-end fork at same tf");
    }

    #[test]
    fn test_backtracking_with_multiple_forks() {
        // Cycle ATCGTA (6bp) with k=3. Multiple dead-end forks available
        // but all at lower tf than cycle k-mers.
        let sdat = vec![
            ("ATC".to_string(), 1000),
            ("TCG".to_string(), 1000),
            ("CGT".to_string(), 1000),
            ("GTA".to_string(), 1000),
            ("TAA".to_string(), 1000),
            ("AAT".to_string(), 1000),
            ("TCT".to_string(), 200),   // dead-end fork
            ("GTG".to_string(), 150),   // dead-end fork
        ];
        let results = tr_greedy_finder_bidirectional(&sdat, 30000, 1.0, 30, 3, Some(50), 1000, None);
        let tr_results: Vec<_> = results.iter().filter(|r| r.status == "tr").collect();
        assert!(!tr_results.is_empty(), "DFS should find TR despite multiple dead-end forks");
    }

    #[test]
    fn test_no_backtrack_needed_for_simple_repeat() {
        let sdat = make_circular_sdat("AATGG", 1000, 3);
        let results = tr_greedy_finder_bidirectional(&sdat, 30000, 1.0, 30, 3, Some(100), 0, None);
        assert!(!results.is_empty());
        assert_eq!(results[0].status, "tr");
    }

    #[test]
    fn test_two_tier_threshold_finds_tr_with_low_ext_lu() {
        let sdat = vec![
            ("ATC".to_string(), 200),  // seed (tf >= seed_lu)
            ("TCG".to_string(), 30),   // cycle (tf < seed_lu, tf >= ext_lu)
            ("CGA".to_string(), 30),
            ("GAA".to_string(), 30),
            ("AAT".to_string(), 30),
        ];
        // Without ext_lu → cycle breaks (tf=30 < seed_lu=100)
        let results_no_ext = tr_greedy_finder_bidirectional(&sdat, 30000, 1.0, 30, 3, Some(100), 1000, None);
        let tr_no_ext: Vec<_> = results_no_ext.iter().filter(|r| r.status == "tr").collect();
        assert!(tr_no_ext.is_empty(), "Without ext_lu, cycle should break (tf=30 < seed_lu=100)");

        // With ext_lu=10 → DFS traverses tf=30 k-mers
        let results_ext = tr_greedy_finder_bidirectional(&sdat, 30000, 1.0, 30, 3, Some(100), 1000, Some(10));
        let tr_ext: Vec<_> = results_ext.iter().filter(|r| r.status == "tr").collect();
        assert!(!tr_ext.is_empty(), "With ext_lu=10, DFS should find TR through tf=30 k-mers");
    }

    #[test]
    fn test_cache_prevents_reexploration() {
        // Two overlapping cycles sharing k-mer TCG:
        // Cycle1: ATC→TCG→CGA→GAA→AAT→ATC (tf=1000)
        // Cycle2: GTC→TCG→CGA→GAG→AGG→GGT→GTC (tf=500)
        //
        // Seed1 = ATC (tf=1000) → finds cycle1, caches all its k-mers
        // When seed2 = GTC (tf=500) tries to extend, TCG is in cache → frag
        // This is correct: TCG is already "claimed" by cycle1
        let sdat = vec![
            ("ATC".to_string(), 1000),
            ("TCG".to_string(), 1000),
            ("CGA".to_string(), 1000),
            ("GAA".to_string(), 1000),
            ("AAT".to_string(), 1000),
            ("GTC".to_string(), 500),
            ("GAG".to_string(), 500),
            ("AGG".to_string(), 500),
            ("GGT".to_string(), 500),
        ];
        let results = tr_greedy_finder_bidirectional(&sdat, 30000, 1.0, 30, 3, Some(100), 1000, None);
        let tr_results: Vec<_> = results.iter().filter(|r| r.status == "tr").collect();
        assert_eq!(tr_results.len(), 1, "Should find exactly 1 TR (cycle1); cycle2 blocked by cache");
    }
}
