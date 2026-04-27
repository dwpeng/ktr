pub mod align;
pub mod cli;
pub mod encoding;
pub mod scanner;
pub mod types;

use clap::Parser;
use cli::Cli;

fn main() {
    let cli = Cli::parse();
    let config = types::Config::from(&cli);

    // Set thread count if specified
    if cli.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build_global()
            .unwrap();
    }

    println!("ktr v{} - k={}, max_period={}, min_run_length={}, min_identity={}",
        env!("CARGO_PKG_VERSION"), config.k, config.max_period,
        config.min_run_length, config.min_identity);
}
