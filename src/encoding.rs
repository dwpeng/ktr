/// Lookup table: ASCII base to 2-bit value (0=A, 1=C, 2=G, 3=T), 0xFF for invalid.
const BASE_TO_BITS: [u8; 256] = {
    let mut table = [0xFFu8; 256];
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

/// Convert a single ASCII base to 2-bit value.
/// Returns None for non-ACTG characters.
#[inline]
pub fn base_to_bits(b: u8) -> Option<u8> {
    let v = BASE_TO_BITS[b as usize];
    if v == 0xFF { None } else { Some(v) }
}

/// Check if a base is valid (A/C/G/T, case-insensitive).
#[inline]
pub fn is_valid_base(b: u8) -> bool {
    BASE_TO_BITS[b as usize] != 0xFF
}

/// Encode a kmer starting at position `start` in `seq` into a u64 integer.
/// Returns None if any base is invalid.
///
/// Layout: first base in the highest 2 bits, last base in the lowest 2 bits.
#[inline]
pub fn encode_kmer(seq: &[u8], start: usize, k: usize) -> Option<u64> {
    let mut code = 0u64;
    for i in 0..k {
        let b = seq[start + i];
        let bits = BASE_TO_BITS[b as usize];
        if bits == 0xFF {
            return None;
        }
        code = (code << 2) | bits as u64;
    }
    Some(code)
}

/// Mask for kmer of length k: keeps the lowest 2*k bits.
#[inline]
pub fn kmer_mask(k: usize) -> u64 {
    if k >= 32 {
        !0u64
    } else {
        (1u64 << (2 * k)) - 1
    }
}

/// Incrementally update a kmer code when sliding the window by one position.
/// `next_base` is the newly entered base at position i+k.
/// Returns None if the new base is invalid.
#[inline]
pub fn update_kmer(code: u64, next_base: u8, k: usize) -> Option<u64> {
    let bits = BASE_TO_BITS[next_base as usize];
    if bits == 0xFF {
        return None;
    }
    let mask = kmer_mask(k);
    Some(((code << 2) | bits as u64) & mask)
}

/// Decode a u64 kmer code back to ASCII bases for debugging/stats.
pub fn decode_kmer(code: u64, k: usize) -> Vec<u8> {
    const BITS_TO_BASE: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let mut result = Vec::with_capacity(k);
    for i in (0..k).rev() {
        let bits = ((code >> (2 * i)) & 0b11) as usize;
        result.push(BITS_TO_BASE[bits]);
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_base_to_bits() {
        assert_eq!(base_to_bits(b'A'), Some(0));
        assert_eq!(base_to_bits(b'c'), Some(1));
        assert_eq!(base_to_bits(b'N'), None);
    }

    #[test]
    fn test_encode_kmer() {
        let seq = b"ACGTA";
        let code = encode_kmer(seq, 0, 5).unwrap();
        // A=00, C=01, G=10, T=11, A=00
        assert_eq!(code, 0b00_01_10_11_00u64);
    }

    #[test]
    fn test_encode_kmer_with_n() {
        let seq = b"ACNGT";
        assert!(encode_kmer(seq, 0, 5).is_none());
    }

    #[test]
    fn test_update_kmer() {
        let seq = b"ACGTAA";
        let code = encode_kmer(seq, 0, 5).unwrap();
        // Slide: ACGTA -> CGTAA
        let new_code = update_kmer(code, b'A', 5).unwrap();
        let expected = encode_kmer(seq, 1, 5).unwrap();
        assert_eq!(new_code, expected);
    }

    #[test]
    fn test_decode_roundtrip() {
        let seq = b"ACGTA";
        let code = encode_kmer(seq, 0, 5).unwrap();
        assert_eq!(decode_kmer(code, 5), b"ACGTA");
    }
}
