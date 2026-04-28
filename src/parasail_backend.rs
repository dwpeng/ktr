use parasail_rs::prelude::*;

use crate::{align::PairwiseAligner, types::AlignParams};

/// Aligner backend wrapping parasail's SIMD-accelerated pairwise alignment.
///
/// Uses striped vectorisation (SSE/AVX2) with local (Smith-Waterman) mode.
pub struct ParasailAligner {
    aligner: Aligner,
}

impl ParasailAligner {
    pub fn new(match_score: i32, mismatch: i32, gap_open: i32, gap_extend: i32) -> Self {
        let matrix = Matrix::create(b"ACGT", match_score, mismatch)
            .expect("parasail matrix creation failed");
        // parasail expects POSITIVE gap penalties (it internalizes the sign).
        let aligner = Aligner::new()
            .local()
            .striped()
            .matrix(matrix)
            .gap_open(gap_open.abs())
            .gap_extend(gap_extend.abs())
            .bandwidth(36)
            .build();
        Self { aligner }
    }
}

impl PairwiseAligner for ParasailAligner {
    fn align(&mut self, window: &[u8], consensus: &[u8], _params: &AlignParams) -> (f64, f64) {
        if window.is_empty() || consensus.is_empty() {
            return (0.0, 0.0);
        }

        let alignment = self
            .aligner
            .align(Some(window), consensus)
            .expect("parasail alignment failed");

        // get_similar() requires stats, and without stats we can use the
        // score to estimate matches.  For SW local alignment with our
        // parameters each match contributes +2 to the score.
        let score = alignment.get_score();
        // Score = matches * 2 + ... (mismatches/gaps are negative, but for
        // identity matrix gapped alignment the score is roughly 2 * matched).
        let total = std::cmp::max(window.len(), consensus.len());
        let matches = (score.max(0) as usize / 2).min(total);
        if total == 0 {
            return (0.0, 0.0);
        }

        let identity = matches as f64 / total as f64;
        let indel_rate = (total - matches) as f64 / total as f64;
        (identity, indel_rate)
    }
}
