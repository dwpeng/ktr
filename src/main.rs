use std::{
    fs::File,
    io::{self, BufWriter, Write},
};

use clap::Parser;
use helicase::{HelicaseParser, ParserOptions, input::FromFile};
use ktr::{
    cli::Cli,
    output, scanner,
    types::{Candidate, Config},
    validate,
};
use rayon::prelude::*;

// Helicase parser config (different from types::Config)
const HELICASE_CONFIG: helicase::Config = ParserOptions::default().config();

fn main() -> io::Result<()> {
    let cli = Cli::parse();
    let config = Config::from(&cli);

    if cli.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build_global()
            .unwrap();
    }

    let path = &cli.fasta;
    let mut parser = helicase::FastxParser::<HELICASE_CONFIG>::from_file(path)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

    // Set up output writer
    let stdout = io::stdout();
    let mut out: Box<dyn Write> = match &cli.output {
        Some(path) => {
            let file = File::create(path)?;
            Box::new(BufWriter::new(file))
        }
        None => Box::new(BufWriter::new(stdout.lock())),
    };
    output::write_header(&mut out)?;

    // Process each FASTA record
    while let Some(_event) = parser.next() {
        let header = String::from_utf8_lossy(parser.get_header()).into_owned();
        let seq = parser.get_dna_string();

        if seq.is_empty() {
            continue;
        }

        // Parse sequence name (first word of header)
        let seq_name = header.split_whitespace().next().unwrap_or(&header);

        // Phase 1: stream scan for candidates (chunked, parallel, multi-k)
        let overlap = config.max_period * 2;
        let chunk_size = if config.chunk_size == 0 {
            seq.len()
        } else {
            config.chunk_size
        };

        // Build chunk ranges first (safe to share &seq across threads)
        let mut chunk_ranges: Vec<(usize, usize)> = Vec::new();
        let mut chunk_start: usize = 0;
        while chunk_start < seq.len() {
            let chunk_end = (chunk_start + chunk_size + overlap).min(seq.len());
            chunk_ranges.push((chunk_start, chunk_end));
            chunk_start += chunk_size;
        }

        // Run Phase 1 for each k-mer length and collect all candidates
        let all_candidates: Vec<Candidate> = config
            .k_values
            .iter()
            .flat_map(|&k| {
                let k_config = config.with_k(k);
                chunk_ranges
                    .par_iter()
                    .flat_map(|&(cs, ce)| {
                        let chunk = &seq[cs..ce];
                        let mut cands = scanner::scan_sequence(seq_name, chunk, &k_config);
                        for c in &mut cands {
                            c.start += cs;
                            c.end += cs;
                        }
                        cands
                    })
                    .collect::<Vec<_>>()
            })
            .collect();

        if all_candidates.is_empty() {
            continue;
        }

        // Dedup candidates across k values to avoid redundant Phase 2 validation
        let candidates = dedup_candidates(all_candidates);

        // Phase 2: parallel validation
        let repeats = validate::validate_candidates(candidates, seq, &config);

        // Filter by minimum score, then write
        for tr in repeats.iter().filter(|r| r.score >= config.min_score) {
            output::write_record(&mut out, tr, &config, seq)?;
        }
    }

    Ok(())
}

/// Dedup overlapping candidates from different k-mer scans.
/// When two candidates overlap significantly, keeps the better one.
/// Quality is measured as approximate copies = span / period.
fn dedup_candidates(mut candidates: Vec<Candidate>) -> Vec<Candidate> {
    if candidates.len() < 2 {
        return candidates;
    }
    candidates.sort_by(|a, b| a.seq_name.cmp(&b.seq_name).then(a.start.cmp(&b.start)));

    let mut deduped: Vec<Candidate> = Vec::with_capacity(candidates.len());
    for c in candidates {
        if let Some(last) = deduped.last_mut() {
            if last.seq_name == c.seq_name {
                let overlap = last.end.min(c.end).saturating_sub(last.start.max(c.start));
                let min_span = (last.end - last.start).min(c.end - c.start);
                if overlap > 0 && (overlap as f64 / min_span as f64) > 0.5 {
                    // Overlap > 50%: same region, keep candidate with higher
                    // quality (span / period = approximate copy count).
                    let last_q = (last.end - last.start) as f64 / last.period.max(1) as f64;
                    let c_q = (c.end - c.start) as f64 / c.period.max(1) as f64;
                    if c_q > last_q && c.period > 0 {
                        *last = c;
                    }
                    continue;
                }
            }
        }
        deduped.push(c);
    }
    deduped
}
