/// CLI and algorithm configuration parameters.
#[derive(Debug, Clone)]
pub struct Config {
    pub k: usize,
    pub max_period: usize,
    pub min_run_length: usize,
    pub min_matches: usize,
    pub min_identity: f64,
    pub window_size: usize,
    pub chunk_size: usize,
    pub max_gap: usize,
    pub debug: bool,
}

impl Config {
    pub fn new(k: usize, max_period: usize, min_run_length: usize,
               min_matches: usize, min_identity: f64) -> Self {
        assert!(k > 0 && k <= 32, "kmer length must be between 1 and 32, got {}", k);
        Self {
            k,
            max_period,
            min_run_length,
            min_matches,
            min_identity,
            window_size: max_period * 2,
            chunk_size: 500_000,
            max_gap: 20,
            debug: false,
        }
    }
}

/// A candidate detected during Phase 1 streaming scan.
#[derive(Debug, Clone)]
pub struct Candidate {
    pub seq_name: String,
    pub start: usize,
    pub end: usize,
    pub period: usize,
}

/// Period voting info tracked during Phase 1.
#[derive(Debug, Clone)]
pub struct VoteInfo {
    pub total: usize,
    pub first_pos: usize,
    pub last_pos: usize,
    pub cluster_end: usize,
    pub window_end: usize,
}

impl VoteInfo {
    pub fn new(pos: usize) -> Self {
        Self { total: 1, first_pos: pos, last_pos: pos, cluster_end: pos, window_end: pos }
    }
}

/// Final tandem repeat record output.
#[derive(Debug, Clone)]
pub struct TandemRepeat {
    pub seq_name: String,
    pub start: usize,
    pub end: usize,
    pub period: usize,
    pub copies: f64,
    pub consensus: Vec<u8>,
    pub identity: f64,
    pub indel_rate: f64,
    pub score: f64,
    pub composition: [f64; 4],  // A, C, G, T fraction
    pub entropy: f64,
}

/// Information about one aligned copy.
#[derive(Debug, Clone)]
pub struct CopyInfo {
    pub start: usize,
    pub end: usize,
    pub identity: f64,
    pub indel_rate: f64,
}

/// WDP alignment scoring parameters.
#[derive(Debug, Clone, Copy)]
pub struct AlignParams {
    pub match_score: i32,
    pub mismatch_penalty: i32,
    pub gap_open_penalty: i32,
    pub gap_extend_penalty: i32,
}

impl AlignParams {
    /// Compute alignment parameters adapted to repeat period.
    /// Longer repeats inherently tolerate more divergence (per-base mutation
    /// rate accumulates over more bases), so penalties are softened.
    /// The min_identity threshold in Config still gates overall quality.
    pub fn for_period(period: usize) -> Self {
        // Scale: log₂-based so impact tapers off — p=10→1.0, p=100→0.7, p=500→0.55
        let scale = (2.0 / (period.max(2) as f64).log2()).clamp(0.5, 1.2);
        Self {
            match_score: 2,
            mismatch_penalty: -(5.0 * scale).round().max(2.0) as i32,
            gap_open_penalty: -(7.0 * scale).round().max(2.0) as i32,
            gap_extend_penalty: -(2.0 * scale).round().max(1.0) as i32,
        }
    }
}

impl Default for AlignParams {
    fn default() -> Self {
        Self::for_period(100)
    }
}
