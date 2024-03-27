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
            // Build, merge
            let rc = !*single_strand;
            // Set expected sketchsize
            sketch_size /= u64::BITS as u64;

            log::info!(
                "Running sketching: k:{:?}; sketch_size:{}; threads{}",
                k_vals,
                sketch_size * u64::BITS as u64,
                threads
            );
            let mut sketches = MultiSketch::new(input_files.len(), sketch_size, &k_vals, &output);
            sketch_files(&mut sketches, &input_files, k_vals, sketch_size, rc, *min_count, *min_qual);
            sketches.save();

            // TODO FIRST
            // Maybe the thing to do is not have usigs in the sketch class at all
            // Write straight to a memmapped buffer.
            // Use split_at_mut to safely create chunks of memory for each thread
            // The multisketch class will then read parts of this buffer
            // This does commit to having all in memory
            // Alternative is to create a mutex on file and write out with index in the order
            // sketched
            // Then when reading, create a byte interval to read from and read
            // the file in large chunks (so don't use memmap at all)
            // memmap:
            // basically would be an advantage with a random access pattern
            // if we can write serially then no point

            // TODO saving
            // make a multisketch class which has a vec of sketchsize * n_samples * k-mer
            // Can stride as you like, probably having a single sample's bins layed out next to each other
            // then k-mers next to each other. So when reading can set a range sketchsize * n_kmers to read in a block
            // Will want to make a nice class around access which deals with reading the right part of the [u8] (see u64::from_le_bytes())

            // TODO metadata
            // Probably serde with a no serialise option on the k-mer field would be ok, but maybe check performance for
            // e.g. 2M of these (as otherwise will want to be able to get them via random access)

            // TODO save on the fly
            // Once this is done, just write the sketches out into that vec
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
