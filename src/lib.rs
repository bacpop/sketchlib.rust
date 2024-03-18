//! DOCS
//!

#![warn(missing_docs)]
use std::time::Instant;

extern crate num_cpus;
use ska::io_utils::*;

pub mod cli;
use crate::cli::*;

pub mod sketch;
use crate::sketch::sketch_files;

pub mod hashing;

#[doc(hidden)]
pub fn main() {
    let args = cli_args();
    if args.verbose {
        simple_logger::init_with_level(log::Level::Info).unwrap();
    } else {
        simple_logger::init_with_level(log::Level::Warn).unwrap();
    }

    eprintln!("sketchlib.rust: fast biological distances at multiple k-mer lengths");
    let start = Instant::now();
    match &args.command {
        Commands::Sketch {
            seq_files,
            file_list,
            output,
            k,
            single_strand,
            min_count,
            min_qual,
            threads,
        } => {
            check_threads(*threads);

            // Read input
            let input_files = get_input_list(file_list, seq_files);
            // Build, merge
            let rc = !*single_strand;

            let sketches = sketch_files();
        }
        Commands::Dist {
            db1,
            db2,
            output,
            threads,
        } => {
            check_threads(*threads);
            todo!();
        }
    }
    let end = Instant::now();

    eprintln!("sketchlib done in {}s", end.duration_since(start).as_secs());
    eprintln!("✍️🧬");
    log::info!("Complete");
}
