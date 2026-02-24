/// K-mer encoding/decoding using u64 (2 bits per base, max k=32).
/// A=0b00, C=0b01, G=0b10, T=0b11

const BASE_TO_BITS: [u8; 256] = {
    let mut table = [255u8; 256];
    table[b'A' as usize] = 0;
    table[b'a' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'c' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'g' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b't' as usize] = 3;
    table
};

const BITS_TO_BASE: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// Encode a k-mer string to u64. Returns None if any base is invalid or k > 32.
#[inline]
pub fn encode_kmer(kmer: &[u8]) -> Option<u64> {
    let k = kmer.len();
    if k == 0 || k > 32 {
        return None;
    }
    let mut encoded: u64 = 0;
    for &base in kmer {
        let bits = BASE_TO_BITS[base as usize];
        if bits == 255 {
            return None;
        }
        encoded = (encoded << 2) | (bits as u64);
    }
    Some(encoded)
}

/// Decode a u64-encoded k-mer back to a String.
#[inline]
pub fn decode_kmer(encoded: u64, k: usize) -> String {
    let mut result = vec![0u8; k];
    let mut val = encoded;
    for i in (0..k).rev() {
        result[i] = BITS_TO_BASE[(val & 3) as usize];
        val >>= 2;
    }
    // SAFETY: all bytes are valid ASCII
    unsafe { String::from_utf8_unchecked(result) }
}

/// Compute reverse complement of an encoded k-mer.
/// Complement each 2-bit pair (XOR 0b11), then reverse the order of pairs.
#[inline]
pub fn revcomp_encoded(encoded: u64, k: usize) -> u64 {
    // Complement all bits (XOR with mask of k*2 bits all set)
    let mask = if k == 32 { u64::MAX } else { (1u64 << (k * 2)) - 1 };
    let complemented = encoded ^ mask;
    // Reverse 2-bit pairs
    reverse_2bit_pairs(complemented, k)
}

/// Reverse 2-bit pairs in the lower k*2 bits of val.
#[inline]
fn reverse_2bit_pairs(val: u64, k: usize) -> u64 {
    let mut result: u64 = 0;
    let mut v = val;
    for _ in 0..k {
        result = (result << 2) | (v & 3);
        v >>= 2;
    }
    result
}

/// Extend a k-mer to the right: drop leftmost base, add new base on right.
/// prefix is the current k-mer, base is 0-3 (A/C/G/T), k is k-mer length.
#[inline]
pub fn extend_right(prefix: u64, base: u8, k: usize) -> u64 {
    let mask = if k == 32 { u64::MAX } else { (1u64 << (k * 2)) - 1 };
    ((prefix << 2) | (base as u64)) & mask
}

/// Extend a k-mer to the left: drop rightmost base, add new base on left.
#[inline]
pub fn extend_left(suffix: u64, base: u8, k: usize) -> u64 {
    (suffix >> 2) | ((base as u64) << ((k - 1) * 2))
}

/// Get reverse complement of a DNA string (matching Python get_revcomp).
pub fn revcomp_string(seq: &str) -> String {
    seq.bytes()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'a' => b't',
            b't' => b'a',
            b'c' => b'g',
            b'g' => b'c',
            b'N' => b'N',
            b'n' => b'n',
            b'~' => b'~',
            b'[' => b']',
            b']' => b'[',
            other => other,
        })
        .map(|b| b as char)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_decode_roundtrip() {
        let kmer = "ACGTACGTACGTACGTACGTACG"; // 23-mer
        let encoded = encode_kmer(kmer.as_bytes()).unwrap();
        let decoded = decode_kmer(encoded, 23);
        assert_eq!(kmer, decoded);
    }

    #[test]
    fn test_encode_single_bases() {
        assert_eq!(encode_kmer(b"A"), Some(0));
        assert_eq!(encode_kmer(b"C"), Some(1));
        assert_eq!(encode_kmer(b"G"), Some(2));
        assert_eq!(encode_kmer(b"T"), Some(3));
    }

    #[test]
    fn test_encode_invalid() {
        assert_eq!(encode_kmer(b"N"), None);
        assert_eq!(encode_kmer(b""), None);
    }

    #[test]
    fn test_revcomp() {
        // ACGT -> revcomp = ACGT
        let kmer = "ACGT";
        let encoded = encode_kmer(kmer.as_bytes()).unwrap();
        let rc = revcomp_encoded(encoded, 4);
        assert_eq!(decode_kmer(rc, 4), "ACGT");

        // AAA -> TTT
        let encoded = encode_kmer(b"AAA").unwrap();
        let rc = revcomp_encoded(encoded, 3);
        assert_eq!(decode_kmer(rc, 3), "TTT");
    }

    #[test]
    fn test_extend_right() {
        // ACG + T -> CGT
        let kmer = encode_kmer(b"ACG").unwrap();
        let extended = extend_right(kmer, 3, 3); // T=3
        assert_eq!(decode_kmer(extended, 3), "CGT");
    }

    #[test]
    fn test_extend_left() {
        // ACG + T on left -> TAC
        let kmer = encode_kmer(b"ACG").unwrap();
        let extended = extend_left(kmer, 3, 3); // T=3
        assert_eq!(decode_kmer(extended, 3), "TAC");
    }

    #[test]
    fn test_revcomp_string() {
        assert_eq!(revcomp_string("ACGT"), "ACGT");
        assert_eq!(revcomp_string("AAA"), "TTT");
        assert_eq!(revcomp_string("AACG"), "CGTT");
    }
}
