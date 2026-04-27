pub mod cli;
pub mod types;
pub mod encoding;
pub mod scanner;
pub mod align;
pub mod validate;
pub mod output;

use std::fs::File;
use std::io::{self, BufWriter, Write};
use clap::Parser;
use helicase::input::FromFile;
use helicase::HelicaseParser;
use helicase::ParserOptions;
use cli::Cli;

// Helicase parser config (different from types::Config)
const HELICASE_CONFIG: helicase::Config = ParserOptions::default()
    .config();

fn main() -> io::Result<()> {
    let cli = Cli::parse();
    let config = types::Config::from(&cli);

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
        None => Box::new(stdout.lock()),
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

        // Phase 1: stream scan for candidates
        let candidates = scanner::scan_sequence(seq_name, seq, &config);
        if candidates.is_empty() {
            continue;
        }

        eprintln!("{}: {} candidates from Phase 1", seq_name, candidates.len());

        // Phase 2: parallel validation
        let repeats = validate::validate_candidates(candidates, seq, &config);

        // Write results
        for tr in &repeats {
            output::write_record(&mut out, tr)?;
        }

        eprintln!("{}: {} tandem repeats validated", seq_name, repeats.len());
    }

    Ok(())
}
