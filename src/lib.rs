//! Fast distance calculations between biological sequences (DNA, AA or structures
//! via the 3di alphabet). Distances are based on bindash approximations of the Jaccard
//! distance, with the PopPUNK method to calculate core and accessory distances. nthash/aahash
//! are used for hash functions to create the sketches
//!
//! This package is a work in progress, but is mature enough for research use. See README.md
//! for current CLI usage.

// #![warn(missing_docs)]

use std::io::Write;
use std::sync::mpsc;
use std::time::Instant;

#[macro_use]
extern crate arrayref;
extern crate num_cpus;
use anyhow::Error;
use indicatif::ParallelProgressIterator;
use rayon::prelude::*;
use sketch::num_bins;
use utils::strip_sketch_extension;

pub mod cli;
use crate::cli::*;

pub mod sketch;
use crate::hashing::HashType;
use crate::sketch::sketch_files;

pub mod multisketch;
use crate::multisketch::MultiSketch;

pub mod inverted;
use crate::inverted::Inverted;

pub mod sketch_datafile;

pub mod distances;
use crate::distances::*;

pub mod io;
use crate::io::{get_input_list, parse_kmers, read_subset_names, reorder_input_files, set_ostream};
pub mod structures;

pub mod bloom_filter;
pub mod hashing;

pub mod utils;
use crate::utils::get_progress_bar;

use std::fs::{File, OpenOptions};
use std::io::copy;

use std::io::BufRead;
use std::path::Path;

/// Default k-mer size for (genome) sketching
pub const DEFAULT_KMER: usize = 21;

#[doc(hidden)]
pub fn main() -> Result<(), Error> {
    let args = cli_args();
    if args.quiet {
        simple_logger::init_with_level(log::Level::Error).unwrap();
    } else if args.verbose {
        simple_logger::init_with_level(log::Level::Info).unwrap();
        // simple_logger::init_with_level(log::Level::Trace).unwrap();
    } else {
        simple_logger::init_with_level(log::Level::Warn).unwrap();
    }

    let mut print_success = true;
    let start = Instant::now();
    let result = match &args.command {
        Commands::Sketch {
            seq_files,
            file_list,
            concat_fasta,
            #[cfg(feature = "3di")]
            convert_pdb,
            output,
            kmers,
            sketch_size,
            seq_type,
            level,
            single_strand,
            min_count,
            min_qual,
            threads,
        } => {
            if *concat_fasta && matches!(*seq_type, HashType::DNA | HashType::PDB) {
                panic!("--concat-fasta currently only supported with --seq-type aa");
            }

            // An extra thread is needed for the writer. This doesn't 'overuse' CPU
            check_and_set_threads(*threads + 1);

            // Read input
            log::info!("Getting input files");
            let input_files = get_input_list(file_list, seq_files);
            log::info!("Parsed {} samples in input list", input_files.len());
            let kmers = parse_kmers(kmers);
            // Build, merge
            let rc = !*single_strand;
            // Set aa level
            let seq_type = if let HashType::AA(_) = seq_type {
                HashType::AA(level.clone())
            } else {
                seq_type.clone()
            };

            let (_, sketch_bins, _) = num_bins(*sketch_size);
            log::info!(
                "Running sketching: k:{:?}; sketch_size:{}; seq:{:?}; threads:{}",
                kmers,
                sketch_bins,
                seq_type,
                threads
            );
            let mut sketches = sketch_files(
                output,
                &input_files,
                *concat_fasta,
                #[cfg(feature = "3di")]
                *convert_pdb,
                &kmers,
                sketch_bins,
                &seq_type,
                rc,
                *min_count,
                *min_qual,
                args.quiet,
            );
            let sketch_vec = MultiSketch::new(&mut sketches, sketch_bins, &kmers, seq_type);
            sketch_vec
                .save_metadata(output)
                .expect("Error saving metadata");
            Ok(())
        }
        Commands::Dist {
            ref_db,
            query_db,
            output,
            mut knn,
            subset,
            kmer,
            ani,
            threads,
        } => {
            check_and_set_threads(*threads);

            let mut output_file = set_ostream(output);

            let ref_db_name = utils::strip_sketch_extension(ref_db);

            let mut references = MultiSketch::load(ref_db_name)
                .unwrap_or_else(|_| panic!("Could not read sketch metadata from {ref_db}.skm"));

            log::info!("Loading sketch data from {}.skd", ref_db_name);
            if let Some(subset_file) = subset {
                let subset_names = read_subset_names(subset_file);
                references.read_sketch_data_block(ref_db_name, &subset_names);
            } else {
                references.read_sketch_data(ref_db_name);
            }
            log::info!("Read reference sketches:\n{references:?}");
            let n = references.number_samples_loaded();
            if let Some(nn) = knn {
                if nn >= n {
                    log::warn!("knn={nn} is higher than number of samples={n}");
                    knn = Some(n - 1);
                }
            }
            let (dist_type, k_idx, k_f32) = set_k(&references, *kmer, *ani);

            // Read queries if supplied. Note no subsetting here
            let queries = if let Some(query_db_name) = query_db {
                let mut queries = MultiSketch::load(query_db_name).unwrap_or_else(|_| {
                    panic!("Could not read sketch metadata from {query_db_name}.skm")
                });
                log::info!("Loading query sketch data from {}.skd", query_db_name);
                queries.read_sketch_data(query_db_name);
                log::info!("Read query sketches:\n{queries:?}");
                Some(queries)
            } else {
                None
            };

            match queries {
                None => {
                    // Ref v ref functions
                    match knn {
                        None => {
                            // Self mode (dense)
                            log::info!("Calculating all ref vs ref distances");
                            let distances = self_dists_all(
                                &references,
                                n,
                                k_idx,
                                k_f32,
                                dist_type,
                                *ani,
                                args.quiet,
                            );
                            log::info!("Writing out in long matrix form");
                            write!(output_file, "{distances}")
                                .expect("Error writing output distances");
                        }
                        Some(nn) => {
                            // Self mode (sparse)
                            log::info!("Calculating sparse ref vs ref distances with {nn} nearest neighbours");
                            let distances = self_dists_knn(&references, n, nn, k_idx, k_f32, dist_type, *ani, args.quiet);

                            log::info!("Writing out in sparse matrix form");
                            write!(output_file, "{distances}")
                                .expect("Error writing output distances");
                        }
                    }
                }
                Some(query_db) => {
                    // Ref v query mode
                    log::info!("Calculating all ref vs query distances");

                    let nq = query_db.number_samples_loaded();
                    let distances = self_query_dists_all(&references, &query_db, n, nq, k_idx, k_f32, dist_type, *ani, args.quiet);

                    log::info!("Writing out in long matrix form");
                    write!(output_file, "{distances}").expect("Error writing output distances");
                }
            }
            Ok(())
        }
        Commands::Merge { db1, db2, output } => {
            let ref_db_name1 = utils::strip_sketch_extension(db1);
            let ref_db_name2 = utils::strip_sketch_extension(db2);

            log::info!("Reading input metadata");
            let mut sketches1: MultiSketch = MultiSketch::load(ref_db_name1).unwrap_or_else(|_| {
                panic!("Could not read sketch metadata from {}.skm", ref_db_name1)
            });

            let sketches2: MultiSketch = MultiSketch::load(ref_db_name2).unwrap_or_else(|_| {
                panic!("Could not read sketch metadata from {}.skm", ref_db_name2)
            });
            // check compatibility
            if !sketches1.is_compatible_with(&sketches2) {
                panic!("Databases are not compatible for merging.")
            }

            log::info!("Merging metadata to {}.skm", output);
            let merged_sketch = sketches1.merge_sketches(&sketches2);
            // merge metadata
            merged_sketch
                .save_metadata(output)
                .unwrap_or_else(|_| panic!("Couldn't save metadata to {}", output));

            // merge actual sketch data
            log::info!("Merging and saving sketch data to {}.skd", output);
            utils::save_sketch_data(ref_db_name1, ref_db_name2, output)
        }
        Commands::Inverted { command } => match command {
            InvertedCommands::Build {
                seq_files,
                file_list,
                output,
                species_names,
                single_strand,
                min_count,
                min_qual,
                threads,
                sketch_size,
                kmer_length,
            } => {
                // An extra thread is needed for the writer
                check_and_set_threads(*threads + 1);

                // Get input files
                log::info!("Getting input files");
                let input_files: Vec<(String, String, Option<String>)> =
                    get_input_list(file_list, seq_files);
                log::info!("Parsed {} samples in input list", input_files.len());

                // Reordering by species, or default
                let file_order = if let Some(species_name_file) = species_names {
                    reorder_input_files(&input_files, species_name_file)
                } else {
                    (0..input_files.len()).collect()
                };

                let rc = !*single_strand;
                let seq_type = &HashType::DNA;
                let inverted = Inverted::new(
                    &input_files,
                    &file_order,
                    *kmer_length,
                    *sketch_size, // unconstrained, equals the number of bins here, doesn't need to be a multiple of 64
                    seq_type,
                    rc,
                    *min_count,
                    *min_qual,
                    args.quiet,
                );
                inverted.save(output)?;
                log::info!("Index info:\n{inverted:?}");
                Ok(())
            }
            InvertedCommands::Query {
                ski,
                seq_files,
                file_list,
                output,
                query_type,
                min_count,
                min_qual,
                threads,
            } => {
                let mut output_file = set_ostream(output);
                let inverted_index = Inverted::load(strip_sketch_extension(ski))?;

                // Get input files
                log::info!("Getting input queries");
                let input_files: Vec<(String, String, Option<String>)> =
                    get_input_list(file_list, seq_files);
                log::info!("Parsed {} samples in input query list", input_files.len());

                log::info!("Sketching input queries");
                check_and_set_threads(*threads + 1);
                let (queries, query_names) =
                    inverted_index.sketch_queries(&input_files, *min_count, *min_qual, args.quiet);

                log::info!("Running queries in mode: {query_type}");
                if *query_type == InvertedQueryType::MatchCount {
                    // Header
                    write!(output_file, "Query")?;
                    for name in inverted_index.sample_names() {
                        write!(output_file, ",{name}")?;
                    }
                    writeln!(output_file)?;
                }
                let (tx, rx) = mpsc::channel();
                let percent = false;
                let progress_bar = get_progress_bar(queries.len(), percent, args.quiet);
                rayon::scope(|s| {
                    s.spawn(|_| {
                        queries
                            .par_iter()
                            .progress_with(progress_bar)
                            .zip(query_names)
                            .map(|(q, q_name)| match query_type {
                                InvertedQueryType::MatchCount => {
                                    (q_name, inverted_index.query_against_inverted_index(q))
                                }
                                InvertedQueryType::AllBins => {
                                    (q_name, inverted_index.all_shared_bins(q))
                                }
                                InvertedQueryType::AnyBin => {
                                    (q_name, inverted_index.any_shared_bins(q))
                                }
                            })
                            .for_each_with(tx, |tx, dists| {
                                let _ = tx.send(dists);
                            });
                    });
                });
                for (q_name, dist) in rx {
                    write!(output_file, "{q_name}")?;
                    if *query_type == InvertedQueryType::MatchCount {
                        for distance in dist {
                            write!(output_file, ",{distance}")?;
                        }
                    } else {
                        for r_name in dist
                            .iter()
                            .map(|idx| inverted_index.sample_at(*idx as usize))
                        {
                            write!(output_file, ",{r_name}")?;
                        }
                    }
                    writeln!(output_file)?;
                }
                Ok(())
            }
            InvertedCommands::Precluster {
                ski,
                ref_db,
                output,
                knn,
                ani,
                threads,
            } => {
                // TODO
                // 1. Check ski and ref_db samples and k-mer length match. Will also need to
                // check order/make a mapping between the two
                // 2. Get all the comparison pairs from the inverted index
                // 3. Run the jaccard/ani and knn comparison from dists as above, will always be sparse
                // 3a. Ideally should be made into a function

                let mut output_file = set_ostream(output);
                let inverted_index = Inverted::load(strip_sketch_extension(ski))?;

                let prefilter_pairs = inverted_index.any_shared_bin_list(args.quiet);
                log::info!(
                    "Identified {} prefilter pairs from a max of {}",
                    prefilter_pairs.len(),
                    inverted_index.sample_names().len() * (inverted_index.sample_names().len() - 1)
                        / 2
                );

                Ok(())
            }
        },

        Commands::Append {
            db,
            seq_files,
            file_list,
            output,
            single_strand,
            min_count,
            min_qual,
            concat_fasta,
            threads,
            level,
        } => {
            // An extra thread is needed for the writer. This doesn't 'overuse' CPU
            check_and_set_threads(*threads + 1);
            //get input files
            log::info!("Getting input files");
            let input_files: Vec<(String, String, Option<String>)> =
                get_input_list(file_list, seq_files);
            log::info!("Parsed {} samples in input list", input_files.len());

            //check if any of the new files are already existant in the db
            let db_metadata: MultiSketch = MultiSketch::load(db)?;

            if !db_metadata.append_compatibility(&input_files) {
                panic!("Databases are not compatible for merging.")
            }
            log::info!("Passed concat check");

            // read out sketching information needed to sketch the new files
            let kmers = db_metadata.kmer_lengths();
            let rc = !*single_strand;
            let sketch_size = db_metadata.sketch_size;
            let seq_type = db_metadata.get_hash_type();
            if *concat_fasta && matches!(*seq_type, HashType::DNA | HashType::PDB) {
                panic!("--concat-fasta currently only supported with --seq-type aa");
            }

            log::info!(
                "Running sketching: k:{:?}; sketch_size:{}; seq:{:?}; threads:{}",
                kmers,
                sketch_size * u64::BITS as u64,
                seq_type,
                threads,
            );

            let seq_type = if let HashType::AA(_) = seq_type {
                HashType::AA(level.clone())
            } else {
                seq_type.clone()
            };
            // sketch genomes and save them to concat output file
            let mut db2_sketches = sketch_files(
                output,
                &input_files,
                *concat_fasta,
                #[cfg(feature = "3di")]
                false,
                kmers,
                sketch_size,
                &seq_type,
                rc,
                *min_count,
                *min_qual,
                args.quiet,
            );
            let mut db2_metadata =
                MultiSketch::new(&mut db2_sketches, sketch_size, kmers, seq_type);

            // save skd data from db1 and from freshly sketched input files
            log::info!("Merging and saving sketch data to {}.skd", output);

            let mut output_file = OpenOptions::new()
                .create(true)
                .append(true)
                .open(format!("{}.skd", output))?;
            // stream sketch data directly to concat output file
            let mut db_sketch = File::open(format!("{}.skd", db))?;
            copy(&mut db_sketch, &mut output_file)?;

            // merge and update skm from db1 and the new just sketched sketch
            let concat_metadata = db2_metadata.merge_sketches(&db_metadata);
            concat_metadata
                .save_metadata(output)
                .unwrap_or_else(|_| panic!("Could not save metadata to {}", output));
            Ok(())
        }

        Commands::Delete {
            db,
            samples,
            output_file,
        } => {
            let ref_db = utils::strip_sketch_extension(db);

            log::info!("Reading input genomes");
            let path = Path::new(samples);
            let file = File::open(path)?;
            let reader = std::io::BufReader::new(file);

            // Read in genomes to
            let ids: Vec<String> = reader.lines().map_while(Result::ok).collect();

            log::info!("Reading input metadata");
            let mut sketches: MultiSketch = MultiSketch::load(ref_db)
                .unwrap_or_else(|_| panic!("Could not read sketch metadata from {}.skm", ref_db));

            // write new .skm
            sketches.remove_metadata(output_file, &ids)?;

            // remove samples from .skd file
            log::info!("Remove genomes and writing output");
            sketches.remove_genomes(ref_db, output_file, &ids)?;

            log::info!("Finished writing filtered sketch data to {}", output_file);

            Ok(())
        }

        Commands::Info {
            skm_file,
            sample_info,
        } => {
            if skm_file.ends_with(".ski") {
                let ski_file = &skm_file[0..skm_file.len() - 4];
                let index = Inverted::load(ski_file).unwrap_or_else(|err| {
                    println!("Read error: {err}");
                    panic!("Could not read inverted index from {ski_file}.ski")
                });
                if *sample_info {
                    log::info!("Printing sample info");
                    println!("{index}");
                } else {
                    log::info!("Printing inverted index info");
                    println!("{index:?}");
                }
            } else {
                let ref_db_name = if skm_file.ends_with(".skm") || skm_file.ends_with(".skd") {
                    &skm_file[0..skm_file.len() - 4]
                } else {
                    skm_file.as_str()
                };
                let sketches = MultiSketch::load(ref_db_name).unwrap_or_else(|_| {
                    panic!("Could not read sketch metadata from {ref_db_name}.skm")
                });
                if *sample_info {
                    log::info!("Printing sample info");
                    println!("{sketches}");
                } else {
                    log::info!("Printing database info");
                    println!("{sketches:?}");
                }
            }

            print_success = false; // Turn the final message off
            Ok(())
        }
    };
    let end = Instant::now();

    log::info!("Complete");
    if print_success && !args.quiet {
        eprintln!(
            "🧬🖋️ sketchlib done in {}s",
            end.duration_since(start).as_secs()
        );
    }
    result
}
