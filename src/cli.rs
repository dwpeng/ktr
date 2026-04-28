use clap::Parser;
use crate::types::Config;

#[derive(Parser, Debug)]
#[command(name = "ktr", version, about = "Streaming tandem repeat finder for FASTA genomes")]
pub struct Cli {
    /// Input FASTA file
    pub fasta: String,

    /// kmer length (default: 5)
    #[arg(long, default_value = "5")]
    pub k: usize,

    /// Maximum tandem repeat period (also sets sliding window = 2 * max_period)
    #[arg(long, default_value = "5000")]
    pub max_period: usize,

    /// Minimum run length of consecutive repeat kmers (in bases)
    #[arg(long, default_value = "40")]
    pub min_run_length: usize,

    /// Minimum kmer matches for a candidate period
    #[arg(long, default_value = "3")]
    pub min_matches: usize,

    /// Minimum sequence identity between copies (0.0 to 1.0)
    #[arg(long, default_value = "0.70")]
    pub min_identity: f64,

    /// Output file (default: stdout)
    #[arg(long)]
    pub output: Option<String>,

    /// Number of threads for Phase 2 parallel validation
    #[arg(long, default_value = "0")]
    pub threads: usize,

    /// Debug mode: output repeat sequence and consensus unit in TSV
    #[arg(long)]
    pub debug: bool,

    /// Minimum concentration of period votes in a run (0.0 to 1.0, default: 0.20)
    #[arg(long, default_value = "0.20")]
    pub min_concentration: f64,

    /// Chunk size for splitting long sequences (0 = no chunking, default: 500000)
    #[arg(long, default_value = "500000")]
    pub chunk_size: usize,

    /// Minimum score for output filtering (default: 0, no filter)
    #[arg(long, default_value = "0")]
    pub min_score: f64,
}

impl From<&Cli> for Config {
    fn from(cli: &Cli) -> Self {
        let mut config = Config::new(
            cli.k,
            cli.max_period,
            cli.min_run_length,
            cli.min_matches,
            cli.min_identity,
        );
        config.debug = cli.debug;
        config.chunk_size = cli.chunk_size;
        config.min_concentration = cli.min_concentration;
        config.min_score = cli.min_score;
        config
    }
}
