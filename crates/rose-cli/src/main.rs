//! Command-line entry point for `rose-rs`.

mod args;
mod pipeline;
mod plot;

use std::process::ExitCode;

use clap::Parser;

use crate::args::Cli;

fn main() -> ExitCode {
    let raw_args: Vec<String> = std::env::args().collect();
    let cli = Cli::parse();

    if cli.threads > 0 {
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build_global();
    }

    match pipeline::dispatch(cli.command, raw_args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("rose error: {err:#}");
            ExitCode::from(1)
        }
    }
}
