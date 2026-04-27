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
}

impl From<&Cli> for Config {
    fn from(cli: &Cli) -> Self {
        Config::new(
            cli.k,
            cli.max_period,
            cli.min_run_length,
            cli.min_matches,
            cli.min_identity,
        )
    }
}
