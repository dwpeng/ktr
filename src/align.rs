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

            // Refine consensus by merging (cap to 2x the guessed period)
            merge_into_consensus(&mut consensus, window, Some(guess_period * 2));
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
///
/// Uses a full DP + traceback matrix to compute identity from the actual
/// alignment path, not from naive zipping.
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

    // Use gap_open_penalty as the linear gap penalty (approximates TRF behaviour).
    let gap_penalty = params.gap_open_penalty;

    // DP score matrix and traceback directions.
    // trace[i][j]: 0 = stop (score hit 0), 1 = diagonal, 2 = up, 3 = left
    let mut dp = vec![vec![0i32; m + 1]; n + 1];
    let mut trace = vec![vec![0u8; m + 1]; n + 1];

    let mut max_score = 0i32;
    let mut max_i = 0;
    let mut max_j = 0;

    for i in 1..=n {
        for j in 1..=m {
            let match_score = if window[i - 1] == consensus[j - 1] {
                params.match_score
            } else {
                params.mismatch_penalty
            };

            let diag = dp[i - 1][j - 1] + match_score;
            let up = dp[i - 1][j] + gap_penalty;
            let left = dp[i][j - 1] + gap_penalty;

            // Smith-Waterman recurrence: pick the best of the three or 0.
            if diag >= up && diag >= left && diag > 0 {
                dp[i][j] = diag;
                trace[i][j] = 1;
            } else if up >= left && up > 0 {
                dp[i][j] = up;
                trace[i][j] = 2;
            } else if left > 0 {
                dp[i][j] = left;
                trace[i][j] = 3;
            } else {
                dp[i][j] = 0;
                trace[i][j] = 0;
            }

            if dp[i][j] > max_score {
                max_score = dp[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }

    // Trace back from the cell with the maximum score to compute alignment stats.
    let mut matches = 0usize;
    let mut aligned_len = 0usize;
    let mut i = max_i;
    let mut j = max_j;

    while i > 0 && j > 0 && trace[i][j] != 0 {
        match trace[i][j] {
            1 => {
                // Diagonal: aligned bases
                if window[i - 1] == consensus[j - 1] {
                    matches += 1;
                }
                aligned_len += 1;
                i -= 1;
                j -= 1;
            }
            2 => {
                // Up: gap in consensus (insertion in window)
                aligned_len += 1;
                i -= 1;
            }
            3 => {
                // Left: gap in window (deletion)
                aligned_len += 1;
                j -= 1;
            }
            _ => break,
        }
    }

    let identity = if aligned_len > 0 {
        matches as f64 / aligned_len as f64
    } else {
        0.0
    };
    let indel_rate = if aligned_len > 0 {
        (aligned_len - matches) as f64 / aligned_len as f64
    } else {
        0.0
    };

    (identity, indel_rate)
}

/// Merge a newly aligned copy into the consensus sequence.
///
/// Conservative approach: mismatches keep the original consensus base.
/// Never extends on a longer copy; optionally caps total length to prevent
/// unbounded growth.
fn merge_into_consensus(consensus: &mut Vec<u8>, copy: &[u8], max_len: Option<usize>) {
    let min_len = consensus.len().min(copy.len());
    for i in 0..min_len {
        if consensus[i] != copy[i] {
            // Keep majority: for now, since we can't track counts easily,
            // keep consensus unchanged (conservative approach)
        }
    }
    // Only ever truncate, never extend -- be conservative.
    if copy.len() < consensus.len() {
        consensus.truncate(copy.len());
    }
    // Cap length to prevent unbounded growth.
    if let Some(max) = max_len {
        if consensus.len() > max {
            consensus.truncate(max);
        }
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
