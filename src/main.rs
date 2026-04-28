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

        // Phase 1: stream scan for candidates (chunked, parallel)
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

        // Scan chunks in parallel via Rayon
        let candidates: Vec<Candidate> = chunk_ranges
            .par_iter()
            .flat_map(|&(cs, ce)| {
                let chunk = &seq[cs..ce];
                let mut cands = scanner::scan_sequence(seq_name, chunk, &config);
                for c in &mut cands {
                    c.start += cs;
                    c.end += cs;
                }
                cands
            })
            .collect();

        if candidates.is_empty() {
            continue;
        }

        // Phase 2: parallel validation
        let repeats = validate::validate_candidates(candidates, seq, &config);

        // Filter by minimum score, then write
        for tr in repeats.iter().filter(|r| r.score >= config.min_score) {
            output::write_record(&mut out, tr, &config, seq)?;
        }
    }

    Ok(())
}
