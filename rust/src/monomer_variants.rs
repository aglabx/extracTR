/// Monomer variant discovery via iterative DFS cycle search in de Bruijn graph.

use std::collections::HashSet;
use crate::kmer::{encode_kmer, decode_kmer};
use crate::kmer_index::KmerIndex;

/// Find sequence variants of a monomer via iterative DFS cycle search.
///
/// Looks for cycles of length approximately len(monomer) that start from the
/// same k-mer as the consensus monomer.
///
/// Returns variant sequences including the original.
pub fn find_monomer_variants(
    monomer: &str,
    index: &KmerIndex,
    k: usize,
    max_variants: usize,
    length_tolerance: f64,
) -> Vec<String> {
    if monomer.len() < k {
        return vec![monomer.to_string()];
    }

    // Build tandem to ensure valid k-mers at boundaries
    let tandem = format!("{}{}", monomer, &monomer[..k - 1]);
    let start_kmer = &tandem[..k];
    let start_encoded = match encode_kmer(start_kmer.as_bytes()) {
        Some(e) => e,
        None => return vec![monomer.to_string()],
    };

    // Check if start k-mer exists in index
    if index.get(start_encoded) == 0 {
        return vec![monomer.to_string()];
    }

    let target_len = monomer.len();
    let min_len = (target_len as f64 * (1.0 - length_tolerance)) as usize;
    let max_len = (target_len as f64 * (1.0 + length_tolerance)) as usize;

    let mut variants: HashSet<String> = HashSet::new();
    variants.insert(monomer.to_string());

    let k_mask = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (k * 2)) - 1
    };
    let k_minus_1_mask = if k - 1 == 32 {
        u64::MAX
    } else {
        (1u64 << ((k - 1) * 2)) - 1
    };

    // Iterative DFS with explicit stack
    // Stack items: (current_kmer_encoded, seq_len, base_index)
    let mut stack: Vec<(u64, usize, u8)> = vec![(start_encoded, k, 0)];
    // Path stores encoded k-mers for sequence reconstruction
    let mut path: Vec<u64> = vec![start_encoded];
    let mut visited: HashSet<u64> = HashSet::new();
    visited.insert(start_encoded);

    while !stack.is_empty() && variants.len() < max_variants {
        let frame = stack.last_mut().unwrap();
        let current_kmer = frame.0;
        let seq_len = frame.1;
        let base_idx = frame.2;

        let mut found_next = false;

        for bi in base_idx..4 {
            if variants.len() >= max_variants {
                break;
            }
            let base = bi; // 0=A, 1=C, 2=G, 3=T
            let prefix = current_kmer & k_minus_1_mask;
            let next_kmer = ((prefix << 2) | (base as u64)) & k_mask;

            if index.get(next_kmer) == 0 {
                continue;
            }

            let new_len = seq_len + 1;

            // Check if we completed a cycle back to start
            if next_kmer == start_encoded && min_len <= new_len && new_len <= max_len {
                // Reconstruct sequence from path
                let seq = reconstruct_sequence(&path, next_kmer, k);
                let variant = &seq[..new_len];
                variants.insert(variant.to_string());
                // Update base_index so we resume from next base
                frame.2 = bi + 1;
                found_next = true;
                break;
            }

            if !visited.contains(&next_kmer) && new_len <= max_len {
                // Save resume point
                frame.2 = bi + 1;
                // Push new frame
                visited.insert(next_kmer);
                path.push(next_kmer);
                stack.push((next_kmer, new_len, 0));
                found_next = true;
                break;
            }
        }

        if !found_next {
            // Backtrack
            stack.pop();
            if let Some(removed) = path.pop() {
                visited.remove(&removed);
            }
        }
    }

    variants.into_iter().collect()
}

/// Reconstruct sequence from a path of encoded k-mers.
fn reconstruct_sequence(path: &[u64], last_kmer: u64, k: usize) -> String {
    if path.is_empty() {
        return decode_kmer(last_kmer, k);
    }

    let mut seq = decode_kmer(path[0], k);
    for &kmer in &path[1..] {
        // Each subsequent k-mer adds one base (the last base)
        let last_base = (kmer & 3) as u8;
        seq.push(match last_base {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            _ => unreachable!(),
        });
    }
    // Add last base from the closing k-mer
    let last_base = (last_kmer & 3) as u8;
    seq.push(match last_base {
        0 => 'A',
        1 => 'C',
        2 => 'G',
        3 => 'T',
        _ => unreachable!(),
    });
    seq
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_variants_returns_original() {
        // With an empty index, should return just the monomer
        let index = KmerIndex::from_encoded(vec![], 3);
        let variants = find_monomer_variants("AATGG", &index, 3, 10, 0.15);
        assert_eq!(variants.len(), 1);
        assert!(variants.contains(&"AATGG".to_string()));
    }

    #[test]
    fn test_find_variants_with_cycle() {
        // Build index with all k-mers from circular AATGG
        let monomer = "AATGG";
        let k = 3;
        let extended = format!("{}{}", monomer, &monomer[..k - 1]);
        let mut pairs = Vec::new();
        for i in 0..monomer.len() {
            let kmer = &extended[i..i + k];
            if let Some(e) = encode_kmer(kmer.as_bytes()) {
                pairs.push((e, 100u32));
            }
        }
        let index = KmerIndex::from_encoded(pairs, k);
        let variants = find_monomer_variants(monomer, &index, k, 10, 0.15);
        assert!(variants.contains(&monomer.to_string()));
    }
}
