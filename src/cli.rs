use clap::Parser;

use crate::types::Config;

#[derive(Parser, Debug)]
#[command(
    name = "ktr",
    version,
    about = "Streaming tandem repeat finder for FASTA genomes"
)]
pub struct Cli {
    /// Input FASTA file
    pub fasta: String,

    /// k-mer length(s), comma-separated for multi-k (e.g. 3,5,7)
    #[arg(short, long, default_value = "5")]
    pub k: String,

    /// Maximum tandem repeat period (also sets sliding window = 2 * max_period)
    #[arg(short, long, default_value = "5000")]
    pub period: usize,

    /// Output file (default: stdout)
    #[arg(short, long)]
    pub output: Option<String>,

    /// Number of threads for Phase 2 parallel validation
    #[arg(short, long, default_value = "0")]
    pub threads: usize,

    /// Debug mode: output repeat sequence and consensus unit in TSV
    #[arg(short, long)]
    pub debug: bool,

    /// Chunk size for splitting long sequences (0 = no chunking, default: 500000)
    #[arg(short, long, default_value = "500000")]
    pub chunk_size: usize,

    /// Minimum concentration of period votes in a run (0.0 to 1.0, default: 0.20)
    #[arg(long, default_value = "0.20")]
    pub min_concentration: f64,

    /// Minimum run length of consecutive repeat kmers (in bases)
    #[arg(long, default_value = "40")]
    pub min_run_length: usize,

    /// Minimum kmer matches for a candidate period
    #[arg(long, default_value = "3")]
    pub min_matches: usize,

    /// Minimum sequence identity between copies (0.0 to 1.0)
    #[arg(long, default_value = "0.70")]
    pub min_identity: f64,

    /// Minimum score for output filtering (default: 0, no filter)
    #[arg(long, default_value = "0")]
    pub min_score: f64,
}

impl From<&Cli> for Config {
    fn from(cli: &Cli) -> Self {
        // Parse comma-separated k values
        let k_values: Vec<usize> = cli
            .k
            .split(',')
            .map(|s| {
                s.trim()
                    .parse()
                    .expect("Invalid k-mer length, expected positive integers")
            })
            .collect();
        assert!(
            !k_values.is_empty(),
            "At least one k-mer length must be specified"
        );
        let primary_k = k_values[0];

        let mut config = Config::new(
            primary_k,
            cli.period,
            cli.min_run_length,
            cli.min_matches,
            cli.min_identity,
        );
        config.k_values = k_values;
        config.debug = cli.debug;
        config.chunk_size = cli.chunk_size;
        config.min_concentration = cli.min_concentration;
        config.min_score = cli.min_score;
        config
    }
}
