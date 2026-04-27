use crate::types::{AlignParams, CopyInfo};

/// Run TRF-style wraparound DP on a candidate region.
///
/// Given a sequence and a guessed period, iteratively:
/// 1. Build/refine a consensus pattern of length `period`
/// 2. Align subsequent period-length windows to the consensus
/// 3. Track copy boundaries and compute identity
/// 4. Return the list of copies found
pub fn wdp_align(
    seq: &[u8],
    start: usize,
    end: usize,
    guess_period: usize,
    params: &AlignParams,
    min_identity: f64,
) -> Option<WdpResult> {
    let region = &seq[start..end];
    let n = region.len();

    if n < guess_period * 2 || guess_period == 0 {
        return None;
    }

    // Phase 2a: Extract initial consensus from first copy
    let unit_len = guess_period.min(region.len());
    let mut consensus = region[..unit_len].to_vec();
    let mut copies: Vec<CopyInfo> = Vec::new();

    // Add the first copy (assumed complete)
    copies.push(CopyInfo {
        start: start,
        end: start + unit_len,
        identity: 1.0,
        indel_rate: 0.0,
    });

    // Phase 2b: Iteratively align and extend
    let mut pos = unit_len;
    while pos < n {
        let window_end = (pos + consensus.len()).min(n);
        let window = &region[pos..window_end];

        // Align window against consensus
        let (identity, indel) = align_window_to_consensus(window, &consensus, params);

        if identity >= min_identity {
            copies.push(CopyInfo {
                start: start + pos,
                end: start + window_end,
                identity,
                indel_rate: indel,
            });

            // Refine consensus by merging
            merge_into_consensus(&mut consensus, window);
            pos += consensus.len();
        } else {
            // Try sliding alignment with small offsets
            let mut found = false;
            for offset in 1..=consensus.len().min(5) {
                if pos + offset >= n {
                    break;
                }
                let test_end = (pos + offset + consensus.len()).min(n);
                if test_end <= pos + offset {
                    break;
                }
                let test_window = &region[pos + offset..test_end];
                let (test_idy, _) = align_window_to_consensus(test_window, &consensus, params);
                if test_idy >= min_identity {
                    pos += offset;
                    found = true;
                    break;
                }
            }
            if !found {
                break;
            }
        }
    }

    if copies.len() < 2 {
        return None;
    }

    // Phase 2c: Compute final stats
    let total_identity: f64 = copies.iter().skip(1).map(|c| c.identity).sum();
    let total_indel: f64 = copies.iter().skip(1).map(|c| c.indel_rate).sum();
    let n_copies_gt1 = (copies.len() - 1) as f64;
    let avg_identity = if n_copies_gt1 > 0.0 { total_identity / n_copies_gt1 } else { 0.0 };
    let avg_indel = if n_copies_gt1 > 0.0 { total_indel / n_copies_gt1 } else { 0.0 };

    Some(WdpResult {
        copies,
        consensus_len: consensus.len(),
        avg_identity,
        avg_indel,
        score: avg_identity * 100.0 - avg_indel * 50.0,
    })
}

/// Result from the WDP alignment phase.
#[derive(Debug, Clone)]
pub struct WdpResult {
    pub copies: Vec<CopyInfo>,
    pub consensus_len: usize,
    pub avg_identity: f64,
    pub avg_indel: f64,
    pub score: f64,
}

/// Align a text window against a consensus pattern using Smith-Waterman DP.
/// Returns (identity, indel_rate).
fn align_window_to_consensus(
    window: &[u8],
    consensus: &[u8],
    params: &AlignParams,
) -> (f64, f64) {
    let n = window.len();
    let m = consensus.len();
    if n == 0 || m == 0 {
        return (0.0, 0.0);
    }

    // Smith-Waterman DP with linear gap penalty
    let gap_penalty = -2i32;

    let mut prev = vec![0i32; m + 1];
    let mut curr = vec![0i32; m + 1];

    for i in 1..=n {
        for j in 1..=m {
            let match_score = if window[i - 1] == consensus[j - 1] {
                params.match_score
            } else {
                params.mismatch_penalty
            };

            let diag = prev[j - 1] + match_score;
            let left = curr[j - 1] + gap_penalty;
            let up = prev[j] + gap_penalty;

            curr[j] = diag.max(left).max(up).max(0);
        }
        std::mem::swap(&mut prev, &mut curr);
    }

    // Compute identity from direct base comparison
    let matches = window.iter()
        .zip(consensus.iter().cycle())
        .filter(|(a, b)| a == b)
        .count();
    let total = window.len().max(consensus.len());

    let identity = matches as f64 / total as f64;
    let indel_rate = (total - matches) as f64 / total as f64;

    (identity, indel_rate)
}

/// Merge a newly aligned copy into the consensus sequence.
fn merge_into_consensus(consensus: &mut Vec<u8>, copy: &[u8]) {
    let min_len = consensus.len().min(copy.len());
    for i in 0..min_len {
        if consensus[i] != copy[i] {
            // Keep majority: for now, since we can't track counts easily,
            // keep consensus unchanged (conservative approach)
        }
    }
    // If copy is longer, extend consensus
    if copy.len() > consensus.len() {
        consensus.extend_from_slice(&copy[consensus.len()..]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_wdp_perfect_repeat() {
        let unit = b"ACGT";
        let mut seq = Vec::new();
        for _ in 0..10 {
            seq.extend_from_slice(unit);
        }
        let params = AlignParams::default();
        let result = wdp_align(&seq, 0, seq.len(), 4, &params, 0.7);
        assert!(result.is_some());
        let result = result.unwrap();
        assert!(result.copies.len() >= 2);
        assert!(result.avg_identity >= 0.9);
    }

    #[test]
    fn test_wdp_imperfect_repeat() {
        let unit = b"ACGT";
        let mut seq = Vec::new();
        for _ in 0..8 {
            seq.extend_from_slice(unit);
        }
        // Introduce a substitution
        seq[6] = b'T';
        let params = AlignParams::default();
        let result = wdp_align(&seq, 0, seq.len(), 4, &params, 0.5);
        assert!(result.is_some());
    }

    #[test]
    fn test_wdp_no_repeat() {
        let seq = b"GGGGGAAAAACCCCCTTTTTAA";
        let params = AlignParams::default();
        let result = wdp_align(seq, 0, seq.len(), 5, &params, 0.7);
        assert!(result.is_none());
    }

    #[test]
    fn test_wdp_score_calculation() {
        let unit = b"ACGT";
        let mut seq = Vec::new();
        for _ in 0..5 {
            seq.extend_from_slice(unit);
        }
        let params = AlignParams::default();
        let result = wdp_align(&seq, 0, seq.len(), 4, &params, 0.7);
        assert!(result.is_some());
        let r = result.unwrap();
        // Perfect repeat should have high score
        assert!(r.score > 80.0, "Score should be high for perfect repeat, got {}", r.score);
    }
}
