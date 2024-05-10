//! DOCS
//!

// #![warn(missing_docs)]

use std::collections::BinaryHeap;
use std::io::Write;
use std::time::Instant;

#[macro_use]
extern crate arrayref;
extern crate num_cpus;
use anyhow::Error;
use indicatif::{ParallelProgressIterator, ProgressStyle};
use rayon::prelude::*;

pub mod cli;
use crate::cli::*;

pub mod sketch;
use crate::hashing::HashType;
use crate::sketch::sketch_files;

pub mod multisketch;
use crate::multisketch::MultiSketch;

pub mod jaccard;
use crate::jaccard::{ani_pois, core_acc_dist, jaccard_dist};

pub mod sketch_datafile;

pub mod distances;
use crate::distances::*;

pub mod io;
use crate::io::{get_input_list, parse_kmers, read_subset_names, set_ostream};
pub mod structures;

pub mod bloom_filter;
pub mod hashing;

pub mod utils;
use std::fs::{File, OpenOptions};
use std::io::copy;

/// Default k-mer size for (genome) sketching
pub const DEFAULT_KMER: usize = 17;
/// Chunk size in parallel distance calculations
pub const CHUNK_SIZE: usize = 1000;

#[doc(hidden)]
pub fn main() -> Result<(), Error> {
    let args = cli_args();
    if args.verbose {
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
            k_vals,
            k_seq,
            mut sketch_size,
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
            check_threads(*threads + 1);

            // Read input
            log::info!("Getting input files");
            let input_files = get_input_list(file_list, seq_files);
            log::info!("Parsed {} samples in input list", input_files.len());
            let kmers = parse_kmers(k_vals, k_seq);
            // Build, merge
            let rc = !*single_strand;
            // Set expected sketchsize
            sketch_size = sketch_size.div_ceil(u64::BITS as u64);
            // Set aa level
            let seq_type = if let HashType::AA(_) = seq_type {
                HashType::AA(level.clone())
            } else {
                seq_type.clone()
            };

            log::info!(
                "Running sketching: k:{:?}; sketch_size:{}; seq:{:?}; threads:{}",
                kmers,
                sketch_size * u64::BITS as u64,
                seq_type,
                threads
            );
            let mut sketches = sketch_files(
                output,
                &input_files,
                *concat_fasta,
                #[cfg(feature = "3di")]
                convert_pdb,
                &kmers,
                sketch_size,
                &seq_type,
                rc,
                *min_count,
                *min_qual,
            );
            let sketch_vec = MultiSketch::new(&mut sketches, sketch_size, &kmers, seq_type);
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
            check_threads(*threads);

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

            // Set type of distances to use
            let k_idx;
            let mut k_f32 = 0.0;
            let dist_type = if let Some(k) = kmer {
                k_idx = references.get_k_idx(*k);
                k_f32 = *k as f32;
                DistType::Jaccard(*k, *ani)
            } else {
                k_idx = None;
                DistType::CoreAcc
            };
            log::info!("{dist_type}");

            let bar_style =
                ProgressStyle::with_template("{percent}% {bar:80.cyan/blue} eta:{eta}").unwrap();
            // TODO: possible improvement would be to load sketch slices when i, j change
            // This would require a change to core_acc where multiple k-mer lengths are loaded at once
            // Overall this would be nicer I think (not sure about speed)
            match queries {
                None => {
                    // Ref v ref functions
                    match knn {
                        None => {
                            // Self mode (dense)
                            log::info!("Calculating all ref vs ref distances");
                            let mut distances = DistanceMatrix::new(&references, None, dist_type);
                            let par_chunk = CHUNK_SIZE * distances.n_dist_cols();
                            distances
                                .dists_mut()
                                .par_chunks_mut(par_chunk)
                                .progress_with_style(bar_style)
                                .enumerate()
                                .for_each(|(chunk_idx, dist_slice)| {
                                    // Get first i, j index for the chunk
                                    let start_dist_idx = chunk_idx * CHUNK_SIZE;
                                    let mut i = calc_row_idx(start_dist_idx, n);
                                    let mut j = calc_col_idx(start_dist_idx, i, n);
                                    for dist_idx in 0..CHUNK_SIZE {
                                        // TODO might be good to try and move this if out of the loop... but may not matter
                                        if let Some(k) = k_idx {
                                            let mut dist = jaccard_dist(
                                                references.get_sketch_slice(i, k),
                                                references.get_sketch_slice(j, k),
                                                references.sketch_size,
                                            );
                                            dist = if *ani { ani_pois(dist, k_f32) } else { dist };
                                            dist_slice[dist_idx] = dist;
                                        } else {
                                            let dist =
                                                core_acc_dist(&references, &references, i, j);
                                            dist_slice[dist_idx * 2] = dist.0;
                                            dist_slice[dist_idx * 2 + 1] = dist.1;
                                        }

                                        // Move to next index in upper triangle
                                        j += 1;
                                        if j >= n {
                                            i += 1;
                                            j = i + 1;
                                            // End of all dists reached (final chunk)
                                            if i >= (n - 1) {
                                                break;
                                            }
                                        }
                                    }
                                });

                            log::info!("Writing out in long matrix form");
                            write!(output_file, "{distances}")
                                .expect("Error writing output distances");
                        }
                        Some(nn) => {
                            // Self mode (sparse)
                            log::info!("Calculating sparse ref vs ref distances with {nn} nearest neighbours");
                            let mut sp_distances =
                                SparseDistanceMatrix::new(&references, nn, dist_type);
                            // TODO is it possible to add a template to the trait so this code is only written once? Maybe not
                            match sp_distances.dists_mut() {
                                DistVec::Jaccard(distances) => {
                                    let k = k_idx.unwrap();
                                    distances
                                        .par_chunks_mut(nn)
                                        .progress_with_style(bar_style)
                                        .enumerate()
                                        .for_each(|(i, row_dist_slice)| {
                                            let mut heap = BinaryHeap::with_capacity(nn);
                                            let i_sketch = references.get_sketch_slice(i, k);
                                            for j in 0..n {
                                                if i == j {
                                                    continue;
                                                }
                                                let mut dist = jaccard_dist(
                                                    i_sketch,
                                                    references.get_sketch_slice(j, k),
                                                    references.sketch_size,
                                                );
                                                dist =
                                                    if *ani { ani_pois(dist, k_f32) } else { dist };
                                                let dist_item = SparseJaccard(j, dist);
                                                if heap.len() < nn
                                                    || dist_item < *heap.peek().unwrap()
                                                {
                                                    heap.push(dist_item);
                                                    if heap.len() > nn {
                                                        heap.pop();
                                                    }
                                                }
                                            }
                                            debug_assert_eq!(row_dist_slice.len(), heap.len());
                                            row_dist_slice
                                                .clone_from_slice(&heap.into_sorted_vec());
                                        });
                                }
                                DistVec::CoreAcc(distances) => {
                                    distances
                                        .par_chunks_mut(nn)
                                        .progress_with_style(bar_style)
                                        .enumerate()
                                        .for_each(|(i, row_dist_slice)| {
                                            let mut heap = BinaryHeap::with_capacity(nn);
                                            for j in 0..n {
                                                if i == j {
                                                    continue;
                                                }
                                                let dists =
                                                    core_acc_dist(&references, &references, i, j);
                                                let dist_item = SparseCoreAcc(j, dists.0, dists.1);
                                                if heap.len() < nn
                                                    || dist_item < *heap.peek().unwrap()
                                                {
                                                    heap.push(dist_item);
                                                    if heap.len() > nn {
                                                        heap.pop();
                                                    }
                                                }
                                            }
                                            debug_assert_eq!(row_dist_slice.len(), heap.len());
                                            row_dist_slice
                                                .clone_from_slice(&heap.into_sorted_vec());
                                        });
                                }
                            }

                            log::info!("Writing out in sparse matrix form");
                            write!(output_file, "{sp_distances}")
                                .expect("Error writing output distances");
                        }
                    }
                }
                Some(query_db) => {
                    // Ref v query mode
                    log::info!("Calculating all ref vs query distances");
                    let mut distances =
                        DistanceMatrix::new(&references, Some(&query_db), dist_type);
                    let par_chunk = CHUNK_SIZE * distances.n_dist_cols();
                    distances
                        .dists_mut()
                        .par_chunks_mut(par_chunk)
                        .progress_with_style(bar_style)
                        .enumerate()
                        .for_each(|(chunk_idx, dist_slice)| {
                            // Get first i, j index for the chunk
                            let start_dist_idx = chunk_idx * CHUNK_SIZE;
                            let (mut i, mut j) = calc_query_indices(start_dist_idx, n);
                            for dist_idx in 0..CHUNK_SIZE {
                                if let Some(k) = k_idx {
                                    let mut dist = jaccard_dist(
                                        references.get_sketch_slice(i, k),
                                        query_db.get_sketch_slice(j, k),
                                        references.sketch_size,
                                    );
                                    dist = if *ani { ani_pois(dist, k_f32) } else { dist };
                                    dist_slice[dist_idx] = dist;
                                } else {
                                    let dist = core_acc_dist(&references, &query_db, i, j);
                                    dist_slice[dist_idx * 2] = dist.0;
                                    dist_slice[dist_idx * 2 + 1] = dist.1;
                                }

                                // Move to next index
                                j += 1;
                                if j >= n {
                                    i += 1;
                                    j = 0;
                                    // End of all dists reached (final chunk)
                                    if i >= n {
                                        break;
                                    }
                                }
                            }
                        });

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
            check_threads(*threads + 1);
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
                kmers,
                sketch_size,
                &seq_type,
                rc,
                *min_count,
                *min_qual,
            );
            let mut db2_metadata =
                MultiSketch::new(&mut db2_sketches, sketch_size, kmers, seq_type);

            // save skd data from db1 and from freshly sketched input files
            log::info!("Merging and saving sketch data to {}.skd", output);
            // let mut output_file = File::create(format!("{}.skd", output))?;
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

        Commands::Info {
            skm_file,
            sample_info,
        } => {
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
            print_success = false; // Turn the final message off
            Ok(())
        }
    };
    let end = Instant::now();

    log::info!("Complete");
    if print_success {
        eprintln!(
            "🧬🖋️ sketchlib done in {}s",
            end.duration_since(start).as_secs()
        );
    }
    result
}
