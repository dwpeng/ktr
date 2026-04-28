use std::collections::{HashMap, VecDeque};
use crate::types::{Config, Candidate, VoteInfo};
use crate::encoding::{encode_kmer, update_kmer};

/// Internal state for the streaming Phase 1 scan.
struct ScanState<'a> {
    config: &'a Config,
    seq: &'a [u8],
    seq_name: &'a str,
    seq_len: usize,

    // Data structures
    last_seen: HashMap<u64, usize>,
    window: VecDeque<(usize, u64)>,

    // Run tracking
    run_start: Option<usize>,
    run_end: usize,
    run_periods: Vec<usize>,


    // Period voting
    period_votes: HashMap<usize, VoteInfo>,

    // Output
    candidates: Vec<Candidate>,
}

impl<'a> ScanState<'a> {
    fn new(config: &'a Config, seq: &'a [u8], seq_name: &'a str) -> Self {
        Self {
            config,
            seq,
            seq_name,
            seq_len: seq.len(),
            last_seen: HashMap::new(),
            window: VecDeque::with_capacity(config.window_size),
            run_start: None,
            run_end: 0,
            run_periods: Vec::new(),
            period_votes: HashMap::new(),
            candidates: Vec::new(),
        }
    }

    /// Process the entire sequence and return candidates.
    fn run(&mut self) -> Vec<Candidate> {
        let k = self.config.k;
        let last_pos = self.seq_len.saturating_sub(k);
        let mut code = 0u64;
        let mut valid = false;

        for i in 0..=last_pos {
            // Rebuild from scratch when the sliding window has no valid kmer
            // (first iteration or after an invalid base slides through).
            if !valid {
                if i + k > self.seq_len {
                    self.flush_run(i + k);
                    self.reset_run();
                    continue;
                }
                match encode_kmer(self.seq, i, k) {
                    Some(c) => { code = c; valid = true; }
                    None => { self.flush_run(i + k); self.reset_run(); continue; }
                }
            } else {
                let enter = i + k - 1;
                match update_kmer(code, self.seq[enter], k) {
                    Some(c) => code = c,
                    None => {
                        self.flush_run(i + k);
                        self.reset_run();
                        valid = false;
                        continue;
                    }
                }
            };

            match self.last_seen.get(&code).copied() {
                Some(prev_pos) => {
                    let d = i - prev_pos;
                    if d <= self.config.max_period {
                        // REPEAT within max_period
                        self.handle_repeat(i, code, d);
                        self.update_structures(code, i);
                        self.periodic_cleanup(i);
                        continue;
                    }
                    // d > max_period: treat as non-repeat
                    self.flush_run(i + k);
                    self.reset_run();
                }
                None => {
                    // NOT a repeat
                    self.flush_run(i + k);
                    self.reset_run();
                }
            }

            self.update_structures(code, i);
            self.periodic_cleanup(i);
        }

        // Flush final run at end of sequence
        self.flush_run(self.seq_len);
        std::mem::take(&mut self.candidates)
    }

    fn handle_repeat(&mut self, i: usize, _code: u64, d: usize) {
        // Update period votes
        self.period_votes.entry(d)
            .and_modify(|v| {
                v.total += 1;
                v.last_pos = i;
                // Track the dense cluster end: consecutive votes close together
                // indicate a true TR; isolated votes (gap > d*3) are likely noise.
                if i - v.cluster_end <= d * 3 {
                    v.cluster_end = i;
                }
                v.window_end = i;
            })
            .or_insert_with(|| VoteInfo::new(i));

        // Update run state
        if self.run_start.is_none() {
            self.run_start = Some(i);
        }
        self.run_end = i;
        if !self.run_periods.contains(&d) {
            self.run_periods.push(d);
        }
    }

    fn update_structures(&mut self, code: u64, pos: usize) {
        self.last_seen.insert(code, pos);
        self.window.push_back((pos, code));
        if self.window.len() > self.config.window_size {
            let (out_pos, out_code) = self.window.pop_front().unwrap();
            if self.last_seen.get(&out_code) == Some(&out_pos) {
                self.last_seen.remove(&out_code);
            }
        }
    }

    fn flush_run(&mut self, end_pos: usize) {
        let rs = match self.run_start {
            Some(s) => s,
            None => return,
        };
        let run_len = end_pos - rs;
        if run_len < self.config.min_run_length {
            return;
        }

        // Evaluate each period using its vote positions as the TR boundaries,
        // preventing non-TR flanks from diluting the concentration.
        let best = self.run_periods.iter()
            .filter_map(|&d| self.period_votes.get(&d).map(|v| (d, v)))
            .filter_map(|(d, v)| {
                // The first period vote occurs at the 2nd copy start.
                // Subtract one period to find the true TR start.
                let tr_start = v.first_pos.saturating_sub(d);
                let tr_end = (v.cluster_end + self.config.k).min(end_pos);
                // Use the overlap of this run with the period span so that
                // mutations (which fragment the run) don't pull in sequence
                // from before the current fragment.
                let eff_start = rs.max(tr_start);
                let eff_end = tr_end.max(eff_start);
                let eff_len = eff_end - eff_start;
                if eff_len >= self.config.min_run_length
                    && v.total >= self.config.min_matches
                    && (v.total as f64 / eff_len as f64) >= self.config.min_concentration
                {
                    Some((d, v.total, tr_start, tr_end))
                } else {
                    None
                }
            })
            .max_by_key(|&(_, total, _, _)| total);

        if let Some((d, _, tr_start, tr_end)) = best {
            self.candidates.push(Candidate {
                seq_name: self.seq_name.to_string(),
                start: tr_start,
                end: tr_end,
                period: d,
            });
        }
    }

    fn reset_run(&mut self) {
        self.run_start = None;
        self.run_end = 0;
        self.run_periods.clear();
        self.period_votes.clear();
    }

    fn periodic_cleanup(&mut self, i: usize) {
        if i % 1000 != 0 {
            return;
        }
        self.period_votes.retain(|&d, v| {
            i - v.window_end <= d * 10
        });
    }
}

/// Run Phase 1 streaming scan on a single sequence.
pub fn scan_sequence(
    seq_name: &str,
    seq: &[u8],
    config: &Config,
) -> Vec<Candidate> {
    let mut state = ScanState::new(config, seq, seq_name);
    state.run()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_tandem(unit: &[u8], copies: usize) -> Vec<u8> {
        let mut seq = Vec::with_capacity(unit.len() * copies + 50);
        // Add flanking unique sequence
        seq.extend_from_slice(b"GCTAGCTAGCTAGCTAGCTAGC");
        for _ in 0..copies {
            seq.extend_from_slice(unit);
        }
        seq.extend_from_slice(b"GCTAGCTAGCTAGCTAGCTAGC");
        seq
    }

    #[test]
    fn test_scan_tandem_repeat() {
        let seq = make_tandem(b"ACGT", 20);
        let config = Config::new(5, 50, 20, 3, 0.7);
        let candidates = scan_sequence("test", &seq, &config);
        assert!(!candidates.is_empty(), "Should detect at least one candidate");
        // Period should be 4 or a multiple/vote
        assert!(candidates.iter().any(|c| c.period == 4 || c.period == 8),
            "Expected period 4 or related, got: {:?}", candidates);
    }

    #[test]
    fn test_scan_longer_period() {
        let seq = make_tandem(b"ACGTACGT", 10);  // period 8
        let config = Config::new(5, 50, 30, 3, 0.7);
        let candidates = scan_sequence("test", &seq, &config);
        assert!(!candidates.is_empty(), "Should detect period-8 repeat");
    }

    #[test]
    fn test_scan_no_repeat() {
        let seq = b"GCTAGCTAGCTAGCTAGCTAGC";
        let config = Config::new(5, 50, 40, 3, 0.7);
        let candidates = scan_sequence("test", seq, &config);
        assert!(candidates.is_empty(), "Random-like sequence should not produce candidates");
    }

    #[test]
    fn test_scan_empty_seq() {
        let seq = b"";
        let config = Config::new(5, 50, 20, 3, 0.7);
        let candidates = scan_sequence("empty", seq, &config);
        assert!(candidates.is_empty());
    }

    #[test]
    fn test_scan_with_n() {
        let mut seq = Vec::new();
        seq.extend_from_slice(b"GCTAGCTAGC");
        seq.extend_from_slice(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        seq.push(b'N');
        seq.extend_from_slice(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        seq.extend_from_slice(b"GCTAGCTAGC");

        let config = Config::new(5, 50, 20, 3, 0.7);
        let candidates = scan_sequence("test", &seq, &config);
        assert_eq!(candidates.len(), 2, "N should split into two separate runs");
    }
}
