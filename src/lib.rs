//! DOCS
//!

#![warn(missing_docs)]
use std::time::Instant;

extern crate num_cpus;

pub mod cli;
use crate::cli::*;

pub mod sketch;
use crate::sketch::sketch_files;

pub mod multisketch;
use crate::multisketch::MultiSketch;
pub mod sketch_datafile;

pub mod io;
use crate::io::get_input_list;

pub mod bloom_filter;
pub mod hashing;

/// Default k-mer size for sketching
pub const DEFAULT_KMER: usize = 17;

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
            k_vals,
            mut sketch_size,
            single_strand,
            min_count,
            min_qual,
            threads,
        } => {
            check_threads(*threads);

            // Read input
            log::info!("Getting input files");
            let input_files = get_input_list(file_list, seq_files);
            // TODO this is very clunky, better replace fastx type
            let names: Vec<String> = input_files.iter().map(|x| x.0.to_string()).collect();
            // Build, merge
            let rc = !*single_strand;
            // Set expected sketchsize
            sketch_size /= u64::BITS as u64;

            log::info!(
                "Running sketching: k:{:?}; sketch_size:{}; threads:{}",
                k_vals,
                sketch_size * u64::BITS as u64,
                threads
            );
            let mut sketches = sketch_files(&output, &input_files, k_vals, sketch_size, rc, *min_count, *min_qual);
            log::info!("Saving sketch metadata");
            let sketch_vec = MultiSketch::new(&names, &mut sketches, sketch_size, &k_vals, &format!("{output}.skm"));
            sketch_vec.save_metadata().expect("Error saving metadata");
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
    eprintln!("🧬🖋️");
    log::info!("Complete");
}
