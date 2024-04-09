//! DOCS
//!

#![warn(missing_docs)]
use std::io::Write;
use std::time::Instant;

#[macro_use]
extern crate arrayref;
extern crate num_cpus;
use indicatif::ProgressBar;

pub mod cli;
use crate::cli::*;

pub mod sketch;
use crate::sketch::sketch_files;

pub mod multisketch;
use crate::multisketch::MultiSketch;

pub mod sketch_datafile;

pub mod distances;
use crate::distances::DistanceMatrix;

pub mod io;
use crate::io::{get_input_list, parse_kmers, set_ostream};

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
            k_seq,
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
            let kmers = parse_kmers(k_vals, k_seq);
            // TODO this is very clunky, better replace fastx type
            let mut names: Vec<String> = input_files.iter().map(|x| x.0.to_string()).collect();
            // Build, merge
            let rc = !*single_strand;
            // Set expected sketchsize
            sketch_size /= u64::BITS as u64;

            log::info!(
                "Running sketching: k:{:?}; sketch_size:{}; threads:{}",
                kmers,
                sketch_size * u64::BITS as u64,
                threads
            );
            let mut sketches = sketch_files(
                &output,
                &input_files,
                &kmers,
                sketch_size,
                rc,
                *min_count,
                *min_qual,
            );
            log::info!("Saving sketch metadata");
            let sketch_vec = MultiSketch::new(&mut names, &mut sketches, sketch_size, &kmers);
            sketch_vec
                .save_metadata(output)
                .expect("Error saving metadata");
        }
        Commands::Dist {
            ref_db,
            query_db,
            output,
            subset,
            kmer,
            threads,
        } => {
            check_threads(*threads);
            let mut output_file = set_ostream(output);

            let ref_db_name = if ref_db.ends_with(".skm") {
                &ref_db[0..ref_db.len() - 4]
            } else {
                ref_db.as_str()
            };
            log::info!("Loading sketch metadata from {}.skm", ref_db_name);
            let mut references = MultiSketch::load(ref_db_name)
                .expect(&format!("Could not read sketch metadata from {ref_db}.skm"));
            log::info!("Read sketches:\n{references:?}");

            // TODO deal with subsetting
            log::info!("Loading sketch data from {}.skd", ref_db_name);
            references.read_sketch_data(ref_db_name);

            match query_db {
                None => {
                    // TODO parallelise
                    // Self mode
                    log::info!("Calculating all ref vs ref distances");
                    let k_idx = if let Some(k) = kmer {
                        references.get_k_idx(*k)
                    } else {
                        None
                    };
                    let mut distances = DistanceMatrix::new(&references, None, k_idx.is_some());
                    let bar = ProgressBar::new(distances.n_distances as u64);
                    for i in 0..references.number_samples_loaded() {
                        for j in (i + 1)..references.number_samples_loaded() {
                            if let Some(k) = k_idx {
                                let dist = references.jaccard_dist(i, j, k);
                                distances.add_jaccard_dist(dist);
                            } else {
                                let dist = references.core_acc_dist(i, j);
                                distances.add_core_acc_dist(dist.0, dist.1);
                            }
                            bar.inc(1);
                        }
                    }

                    log::info!("Writing out in long matrix form");
                    write!(output_file, "{distances}").expect("Error writing output distances");
                }
                Some(db2) => {
                    // TODO Ref v query mode
                    log::info!("Calculating all ref vs query distances");
                    todo!()
                }
            }
        }
        Commands::Info {
            skm_file,
            sample_info,
        } => {
            let ref_db_name = if skm_file.ends_with(".skm") {
                &skm_file[0..skm_file.len() - 4]
            } else {
                skm_file.as_str()
            };
            let sketches = MultiSketch::load(ref_db_name).expect(&format!(
                "Could not read sketch metadata from {ref_db_name}.skm"
            ));
            println!("{sketches:?}");
            if *sample_info {
                log::info!("Printing sample info");
                println!("{sketches}");
            }
        }
    }
    let end = Instant::now();

    eprintln!(
        "🧬🖋️ sketchlib done in {}s",
        end.duration_since(start).as_secs()
    );
    log::info!("Complete");
}
