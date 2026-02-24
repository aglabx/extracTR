/// HashMap<u64, u32> wrapper for k-mer frequency lookup.

use std::collections::HashMap;
use crate::kmer::encode_kmer;

pub struct KmerIndex {
    map: HashMap<u64, u32>,
    k: usize,
}

impl KmerIndex {
    /// Build index from sdat: Vec<(String, u32)> sorted by tf descending.
    pub fn from_sdat(sdat: &[(String, u32)], k: usize) -> Self {
        let mut map = HashMap::with_capacity(sdat.len());
        for (kmer_str, tf) in sdat {
            if let Some(encoded) = encode_kmer(kmer_str.as_bytes()) {
                map.insert(encoded, *tf);
            }
        }
        KmerIndex { map, k }
    }

    /// Build index from pre-encoded pairs.
    pub fn from_encoded(pairs: Vec<(u64, u32)>, k: usize) -> Self {
        let mut map = HashMap::with_capacity(pairs.len());
        for (encoded, tf) in pairs {
            map.insert(encoded, tf);
        }
        KmerIndex { map, k }
    }

    /// Get frequency of a k-mer, returns 0 if not found.
    #[inline]
    pub fn get(&self, kmer: u64) -> u32 {
        *self.map.get(&kmer).unwrap_or(&0)
    }

    /// Get k-mer size.
    #[inline]
    pub fn k(&self) -> usize {
        self.k
    }

    /// Number of entries.
    pub fn len(&self) -> usize {
        self.map.len()
    }

    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_sdat() {
        let sdat = vec![
            ("ACGTACGTACGTACGTACGTACG".to_string(), 100),
            ("TTTTTTTTTTTTTTTTTTTTTTT".to_string(), 50),
        ];
        let index = KmerIndex::from_sdat(&sdat, 23);
        assert_eq!(index.len(), 2);

        let key = encode_kmer(b"ACGTACGTACGTACGTACGTACG").unwrap();
        assert_eq!(index.get(key), 100);

        // Unknown k-mer returns 0
        let unknown = encode_kmer(b"AAAAAAAAAAAAAAAAAAAAAAA").unwrap();
        assert_eq!(index.get(unknown), 0);
    }
}
