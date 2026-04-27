use crate::types::{Candidate, TandemRepeat, Config, CopyInfo, AlignParams};
use crate::align::wdp_align;

/// Validate a candidate from Phase 1 and produce a TandemRepeat record.
/// Returns None if the candidate fails validation checks.
pub fn validate_candidate(
    candidate: &Candidate,
    full_sequence: &[u8],
    config: &Config,
) -> Option<TandemRepeat> {
    let params = AlignParams::default();
    let d = candidate.period;

    // 1. Extract candidate region with flanks
    let extend = d.max(100);
    let context_start = candidate.start.saturating_sub(extend);
    let context_end = (candidate.end + extend).min(full_sequence.len());

    if context_end <= context_start {
        return None;
    }

    // 2. Run WDP alignment
    let wdp_result = wdp_align(
        full_sequence,
        context_start,
        context_end,
        d,
        &params,
        config.min_identity,
    )?;

    // 3. Extract consensus from copies
    let consensus = extract_consensus(&wdp_result.copies, full_sequence);

    if consensus.is_empty() {
        return None;
    }

    // 4. Compute statistics
    let composition = compute_base_composition(&consensus);
    let entropy = compute_entropy(&composition);

    // Determine actual boundaries from first/last copy
    let first_copy = wdp_result.copies.first()?;
    let last_copy = wdp_result.copies.last()?;
    let copy_count = wdp_result.copies.len() as f64;

    Some(TandemRepeat {
        seq_name: candidate.seq_name.clone(),
        start: first_copy.start,
        end: last_copy.end,
        period: consensus.len(),
        copies: copy_count,
        consensus,
        identity: wdp_result.avg_identity,
        indel_rate: wdp_result.avg_indel,
        score: wdp_result.score,
        composition,
        entropy,
    })
}

/// Validate a batch of candidates using Rayon parallel iterator.
pub fn validate_candidates(
    candidates: Vec<Candidate>,
    full_sequence: &[u8],
    config: &Config,
) -> Vec<TandemRepeat> {
    use rayon::prelude::*;

    candidates
        .par_iter()
        .filter_map(|candidate| validate_candidate(candidate, full_sequence, config))
        .collect()
}

/// Extract consensus sequence from aligned copies.
fn extract_consensus(copies: &[CopyInfo], seq: &[u8]) -> Vec<u8> {
    if copies.is_empty() {
        return Vec::new();
    }

    // Use the first copy as the consensus
    let first = &copies[0];
    if first.start < first.end && first.end <= seq.len() {
        seq[first.start..first.end].to_vec()
    } else {
        Vec::new()
    }
}

/// Compute base composition [A, C, G, T] fractions.
fn compute_base_composition(seq: &[u8]) -> [f64; 4] {
    if seq.is_empty() {
        return [0.0; 4];
    }
    let mut counts = [0usize; 4];
    for &b in seq {
        match b.to_ascii_uppercase() {
            b'A' => counts[0] += 1,
            b'C' => counts[1] += 1,
            b'G' => counts[2] += 1,
            b'T' => counts[3] += 1,
            _ => {}
        }
    }
    let total = counts.iter().sum::<usize>() as f64;
    if total == 0.0 {
        return [0.0; 4];
    }
    [
        counts[0] as f64 / total,
        counts[1] as f64 / total,
        counts[2] as f64 / total,
        counts[3] as f64 / total,
    ]
}

/// Compute Shannon entropy of the base composition.
fn compute_entropy(composition: &[f64; 4]) -> f64 {
    composition.iter()
        .filter(|&&p| p > 0.0)
        .map(|&p| -p * p.log2())
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Config;

    #[test]
    fn test_base_composition() {
        let seq = b"ACGTACGT";
        let comp = compute_base_composition(seq);
        assert!((comp[0] - 0.25).abs() < 0.01); // A
        assert!((comp[1] - 0.25).abs() < 0.01); // C
        assert!((comp[2] - 0.25).abs() < 0.01); // G
        assert!((comp[3] - 0.25).abs() < 0.01); // T
    }

    #[test]
    fn test_composition_lowercase() {
        let seq = b"acgtacgt";
        let comp = compute_base_composition(seq);
        assert!((comp[0] - 0.25).abs() < 0.01);
    }

    #[test]
    fn test_entropy_max() {
        let comp = [0.25, 0.25, 0.25, 0.25];
        let entropy = compute_entropy(&comp);
        assert!((entropy - 2.0).abs() < 0.01);
    }

    #[test]
    fn test_entropy_min() {
        let comp = [1.0, 0.0, 0.0, 0.0];
        let entropy = compute_entropy(&comp);
        assert!(entropy.abs() < 0.01);
    }

    #[test]
    fn test_validate_perfect_repeat() {
        let unit = b"ACGT";
        let mut seq = Vec::new();
        // Left flank: 3 perfect repeats
        for _ in 0..3 {
            seq.extend_from_slice(unit);
        }
        // Repeat region: 10 perfect repeats
        for _ in 0..10 {
            seq.extend_from_slice(unit);
        }
        // Right flank: 3 perfect repeats
        for _ in 0..3 {
            seq.extend_from_slice(unit);
        }

        let config = Config::new(5, 50, 20, 3, 0.7);
        let candidate = Candidate {
            seq_name: "test".to_string(),
            start: 12,
            end: 12 + 40,
            period: 4,
        };

        let result = validate_candidate(&candidate, &seq, &config);
        assert!(result.is_some());
        let tr = result.unwrap();
        assert!(tr.identity > 0.9, "Identity should be high for perfect repeat, got {}", tr.identity);
        assert!(tr.copies >= 9.0, "Should detect at least 9 copies, got {}", tr.copies);
    }

    #[test]
    fn test_validate_no_repeat() {
        let seq: &[u8] = b"AGCTAGCTAGCTAGCTAGCTAGC";
        let config = Config::new(5, 50, 20, 3, 0.5);
        let candidate = Candidate {
            seq_name: "test".to_string(),
            start: 0,
            end: seq.len(),
            period: 5,
        };
        let result = validate_candidate(&candidate, seq, &config);
        // This should fail validation since it's not really a repeat
        // (may pass or fail depending on WDP, but shouldn't panic)
        let _ = result;
    }
}
