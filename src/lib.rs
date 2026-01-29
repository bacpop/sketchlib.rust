//! Fast distance calculations between biological sequences (DNA, AA or structures
//! via the 3di alphabet). Distances are based on bindash approximations of the Jaccard
//! distance, with the [PopPUNK method](https://poppunk.bacpop.org/index.html) to calculate core and accessory distances. nthash/aahash
//! are used for hash functions to create the sketches
//!
//! ## Files/databases
//!
//! Sketch databases have two files: `.skm` which is the metadata (samples names, base counts etc)
//! and `.skd` which is the actual sketch data. These must have the same prefix.
//!
//! Inverted indexes are `.ski` files, and should be specified using their full name (not just the prefix).
//! These optionally include and `.skq` file which is needed if used for preclustering.
//!
//! ## Usage
//! With all options we typically recommend using `-v` to see all progress during the run.
//!
//! ### Sketching
//!
//! Using input fasta/fastq files, create a sketch database. Run `sketchlib sketch -h` to see the help.
//!
//! - List .fasta files on the command line, or use `-f` to provide a file(s). Inputs can be gzipped or not, this is automatically detected.
//!   From file, these are one line per sample listing:
//!     - One column (fasta input): file name, which is also used as the sample name
//!     - Two columns (fasta input): sample name and file name
//!     - Three columns (fastq input): sample name and two read files
//! - To set the k-mer size in the sketch database you can either give a list of sizes with `--k-vals`
//!   or a sequence `--k-seq` with start,stop,step. e.g. `--k-seq 17,29,4` would sketch at k=17, 21, 25 and 29.
//! - Set the sketch size with `-s`. Typically 1000 is enough for species level resolution, 10000 for within-species/strain
//!   resolution and 100000-1000000 for SNP level resolution.
//! - To sketch amino acid sequences use `--seq-type aa --concat-fasta` if you have the typical case
//!   of each fasta file being a multifasta with many aa sequences. Each one will then be its own sample.
//! - You can also sketch structures with .pdb input, see 'Enabling PDB->3Di' below. This is experimental.
//!
//! ### Distances
//!
//! To compute internal all-vs-all core and accessory distances use:
//! ```bash
//! sketchlib dist db_name
//! ```
//! Note the database names can be the prefix, or the full path to the .skm file. The output
//! is in pairwise 'long' format, which lists the upper triangle of the distance matrix row-by-row.
//!
//! To calculate distances between two different sample sets, each in their own sketch database, use:
//! ```bash
//! sketchlib dist db1 db2
//! ```
//! For example, if you want to query distances of a new sample against an existing database,
//! first sketch the new sample with e.g. `sketchlib sketch -o db2 new_sample.fasta`, then
//! run the above command.
//!
//! Modifiers:
//! - Use `-k` to calculate Jaccard distance at the given k. Otherwise the default is to
//!   calculate across multiple k and output core and accessory distances.
//! - Use `--ani` with `-k` to transform the Jaccard distance into average nucleotide identity.
//! - Use `--subset` to provide a list of sample names to include in the distance calculations,
//!   only these sample will be loaded from the `.skd` file.
//! - Use `-o` to write the distances to a file. The default it to write to stdout, so you can also
//!   use `>` to redirect to a file (progress messages are written to stderr).
//! - Use `--knn` to only keep this many nearest neighbour distances. For very large databases
//!   it may be useful to keep only ~50 distances. This makes the memory use manageable. This sparse output
//!   can be used with e.g. [mandrake](https://github.com/bacpop/mandrake).
//!
//! ### Inverted indexes
//!
//! Inverted indexes can be used for:
//!
//! - Compressed storage of large numbers of sketches.
//! - Fast querying of new samples against large numbers of sketches.
//! - Preclustering to speed up distance operations.
//!
//! #### Building
//!
//! Similar to a normal sketch:
//! ```bash
//! sketchlib inverted build -o inverted -v -k 21 -s 10 -f rfile.txt
//! ```
//! Provide sample labels (for example species, or clusters) with `--species-names`,
//! tab separated sample and label. These do not need to totally overlap with the samples in
//! the database. Samples will be reordered so that clustered samples are next to each
//! other in the index, reducing size and increasing efficiency.
//!
//! #### Querying
//!
//! Query samples can be provided as a list or with `-f`:
//! ```bash
//! sketchlib inverted query -v -f qfile.txt --query-type match-count inverted.ski
//! ```
//! Queries will be sketched anew each time, we do not yet support saving these sketches.
//!
//! Three query types are supported:
//! - `match-count` (default). Gives the count of bins matching between samples and queries.
//! - `all-bins`. Give samples which have identical sketches to the query.
//! - `any-bins`. Gives samples which have at least one bin matching with the query.
//!
//! To convert from counts to a Jaccard index, you can use the count (intersection, c) from
//! the first mode using the sketch size (s) by J = c / (2s - c).
//!
//! All bins will (rapidly) use AND operations to find very close neighbours, any bins
//! will use OR operations to rule out very distant neighbours.
//!
//! #### Preclustering
//! This is an accelerated nearest neighbour query reducing the total number of comparisons, that requires:
//! - An inverted index file, and corresponding `.skq`, generated with the `--write-skq`
//!   flag to `inverted build`. The inverted index should use a small sketch size (e.g. ~10).
//! - A standard sketch database with `.skd` and `.skm`.
//!
//! So with `inverted.ski`, `inverted.skq`, `standard.skd` and `standard.skm` one can run:
//! ```bash
//! sketchlib inverted precluster -v --knn 10 inverted.ski --skd standard --ani
//! ```
//!
//! ### Other operations
//!
//! - `merge` joins two existing sketch databases.
//! - `append` sketches new input samples, and adds them to an existing database.
//! - `delete` removes samples from a sketch database.
//!
//! ## Enabling PDB->3Di
//! conda doesn't work, so make sure it is deactivated
//! ```bash
//! export PYO3_PYTHON=python3
//! python3 -m venv 3di_venv
//! source 3di_venv/bin/activate
//! python3 -m pip install numpy biopython mini3di
//! cargo run -F 3di
//! export PYTHONPATH=${PYTHONPATH}:$(realpath ./)/3di_venv/lib/python3.12/site-packages
//! ```

#![warn(missing_docs)]
#![allow(clippy::too_many_arguments)]

#[cfg(not(target_arch = "wasm32"))]
use std::io::Write;
#[cfg(not(target_arch = "wasm32"))]
use std::sync::mpsc;
#[cfg(not(target_arch = "wasm32"))]
use std::time::Instant;

#[macro_use]
extern crate arrayref;
extern crate num_cpus;
use anyhow::Error;
#[cfg(not(target_arch = "wasm32"))]
use indicatif::ParallelProgressIterator;
#[cfg(not(target_arch = "wasm32"))]
use rayon::prelude::*;

pub mod cli;
#[cfg(not(target_arch = "wasm32"))]
use crate::cli::*;
#[cfg(target_arch = "wasm32")]
use crate::cli::{InvertedQueryType, DEFAULT_MINCOUNT, DEFAULT_MINQUAL};

#[cfg(not(target_arch = "wasm32"))]
use crate::hashing::HashType;

#[cfg(not(target_arch = "wasm32"))]
use hashbrown::{HashMap, HashSet};

pub mod sketch;
#[cfg(not(target_arch = "wasm32"))]
use crate::sketch::multisketch::MultiSketch;
#[cfg(not(target_arch = "wasm32"))]
use crate::sketch::sketch_datafile::SketchArrayReader;
#[cfg(not(target_arch = "wasm32"))]
use crate::sketch::{num_bins, sketch_files};

pub mod inverted;
use crate::inverted::Inverted;

pub mod distances;
#[cfg(not(target_arch = "wasm32"))]
use crate::distances::*;

pub mod io;

#[cfg(not(target_arch = "wasm32"))]
use crate::io::{
    get_input_list, parse_kmers, parse_metadata_info, read_completeness_file, read_subset_names,
    reorder_input_files, set_ostream,
};

pub mod structures;

pub mod hashing;

pub mod utils;
use crate::utils::get_progress_bar;
#[cfg(not(target_arch = "wasm32"))]
use crate::utils::strip_sketch_extension;

#[cfg(target_arch = "wasm32")]
pub mod fastx_wasm;

#[cfg(not(target_arch = "wasm32"))]
use std::fs::{File, OpenOptions};
#[cfg(not(target_arch = "wasm32"))]
use std::io::copy;

#[cfg(not(target_arch = "wasm32"))]
use std::io::BufRead;
#[cfg(not(target_arch = "wasm32"))]
use std::path::Path;

/// Default k-mer size for (genome) sketching
pub const DEFAULT_KMER: usize = 21;

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;
#[cfg(target_arch = "wasm32")]
extern crate console_error_panic_hook;

#[doc(hidden)]
#[cfg(not(target_arch = "wasm32"))]
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
                "Running sketching: k:{kmers:?}; sketch_size:{sketch_bins}; seq:{seq_type:?}; threads:{threads}"
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
            completeness_file,
            completeness_cutoff,
        } => {
            check_and_set_threads(*threads);

            let mut output_file = set_ostream(output);

            let ref_db_name = utils::strip_sketch_extension(ref_db);

            let mut references = MultiSketch::load_metadata(ref_db_name)
                .unwrap_or_else(|_| panic!("Could not read sketch metadata from {ref_db}.skm"));

            log::info!("Loading sketch data from {ref_db_name}.skd");
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
            // Read in completeness (parallel implementation)
            let completeness_vec: Option<Vec<f64>> = if let Some(file_path) = completeness_file {
                Some(read_completeness_file(file_path, &references)?)
            } else {
                None
            };

            let dist_type = set_k(&references, *kmer, *ani).unwrap_or_else(|e| {
                panic!("Error setting k size: {e}");
            });

            // Read queries if supplied. Note no subsetting here
            let queries = if let Some(query_db_name) = query_db {
                let mut queries = MultiSketch::load_metadata(query_db_name).unwrap_or_else(|_| {
                    panic!("Could not read sketch metadata from {query_db_name}.skm")
                });
                log::info!("Loading query sketch data from {query_db_name}.skd");
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
                                dist_type,
                                args.quiet,
                                completeness_vec.as_ref(),
                                *completeness_cutoff,
                            );
                            log::info!("Writing out in long matrix form");
                            write!(output_file, "{distances}")
                                .expect("Error writing output distances");
                        }
                        Some(nn) => {
                            // Self mode (sparse)
                            log::info!("Calculating sparse ref vs ref distances with {nn} nearest neighbours");
                            let distances = self_dists_knn(
                                &references,
                                n,
                                nn,
                                dist_type,
                                args.quiet,
                                completeness_vec.as_ref(),
                                *completeness_cutoff,
                            );

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
                    let distances = self_query_dists_all(
                        &references,
                        &query_db,
                        n,
                        nq,
                        dist_type,
                        args.quiet,
                        completeness_vec.as_ref(),
                        *completeness_cutoff,
                    );

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
            let mut sketches1: MultiSketch = MultiSketch::load_metadata(ref_db_name1)
                .unwrap_or_else(|_| {
                    panic!("Could not read sketch metadata from {ref_db_name1}.skm")
                });

            let sketches2: MultiSketch =
                MultiSketch::load_metadata(ref_db_name2).unwrap_or_else(|_| {
                    panic!("Could not read sketch metadata from {ref_db_name2}.skm")
                });
            // check compatibility
            if !sketches1.is_compatible_with(&sketches2) {
                panic!("Databases are not compatible for merging.")
            }

            log::info!("Merging metadata to {output}.skm");
            let merged_sketch = sketches1.merge_sketches(&sketches2);
            // merge metadata
            merged_sketch
                .save_metadata(output)
                .unwrap_or_else(|_| panic!("Couldn't save metadata to {output}"));

            // merge actual sketch data
            log::info!("Merging and saving sketch data to {output}.skd");
            utils::save_sketch_data(ref_db_name1, ref_db_name2, output)
        }
        Commands::Inverted { command } => match command {
            InvertedCommands::Build {
                seq_files,
                file_list,
                output,
                write_skq,
                species_names,
                metadata,
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
                let input_files: Vec<(String, Vec<String>)> = get_input_list(file_list, seq_files);
                log::info!("Parsed {} samples in input list", input_files.len());

                let mut differentsamples: HashSet<String> = HashSet::new();

                for i in input_files.iter() {
                    differentsamples.insert(i.0.clone());
                }

                // Reordering by species, or default
                let (file_order, map_names_labels) = if let Some(species_name_file) = species_names
                {
                    reorder_input_files(&input_files, species_name_file)
                } else {
                    // Check first if there are repeated samples

                    let tmpnamesset = input_files
                        .iter()
                        .map(|x| x.0.clone())
                        .collect::<HashSet<String>>();
                    if tmpnamesset.len() == input_files.len() {
                        ((0..input_files.len()).collect(), None)
                    } else {
                        let mut tmpoutvec: Vec<usize> = vec![0; input_files.len()];
                        let mut tmpmap: HashMap<String, usize> = HashMap::new();

                        for (i, name) in tmpnamesset.iter().enumerate() {
                            tmpmap.insert(name.clone(), i);
                        }
                        for i in 0..tmpoutvec.len() {
                            tmpoutvec[i] = tmpmap[&input_files[i].0];
                        }

                        (tmpoutvec, None)
                    }
                };

                // If species labels were provided, create the list of them
                let species_labels_vec = if let Some(themaplabels) = map_names_labels {
                    let mut tmpvec: Vec<String> = vec!["".to_string(); differentsamples.len()];
                    file_order
                        .iter()
                        .zip(&input_files)
                        .for_each(|(idx, (name, _))| {
                            // log::info!("{:?} {:?}", name, idx);
                            tmpvec[*idx] = themaplabels.get(name).unwrap_or(&"".to_owned()).clone();
                        });
                    Some(tmpvec)
                } else {
                    None
                };

                // Parse metadata, if any
                let metadata_vec;
                if let Some(metadata_file) = metadata {
                    let tmpdict = parse_metadata_info(metadata_file);
                    let mut tmpvec: Vec<String> = vec!["".to_string(); differentsamples.len()];
                    file_order
                        .iter()
                        .zip(&input_files)
                        .for_each(|(idx, (name, _))| tmpvec[*idx] = tmpdict[name].clone());
                    metadata_vec = Some(tmpvec);
                } else {
                    metadata_vec = None;
                };

                let skq_file = if *write_skq {
                    Some(format!("{output}.skq"))
                } else {
                    None
                };

                let rc = !*single_strand;
                let seq_type = &HashType::DNA;
                let inverted = Inverted::new(
                    &input_files,
                    skq_file,
                    &file_order,
                    *kmer_length,
                    *sketch_size, // unconstrained, equals the number of bins here, doesn't need to be a multiple of 64
                    seq_type,
                    rc,
                    *min_count,
                    *min_qual,
                    args.quiet,
                    &metadata_vec,
                    &species_labels_vec,
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
                log::info!("Read inverted index:\n{inverted_index:?}");

                // Get input files
                log::info!("Getting input queries");
                let input_files: Vec<(String, Vec<String>)> = get_input_list(file_list, seq_files);
                log::info!("Parsed {} samples in input query list", input_files.len());

                log::info!("Sketching input queries");
                check_and_set_threads(*threads + 1); // Writer thread
                let (queries, query_names) =
                    inverted_index.sketch_queries(&input_files, *min_count, *min_qual, args.quiet);

                log::info!("Running queries in mode: {query_type}");
                // Header
                write!(output_file, "Query")?;
                if *query_type == InvertedQueryType::MatchCount {
                    for name in inverted_index.sample_names() {
                        write!(output_file, "\t{name}")?;
                    }
                    writeln!(output_file)?;
                } else {
                    writeln!(output_file, "\tMatches")?;
                }

                // Query loop (parallelised)
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
                                InvertedQueryType::AnyBins => {
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
                            write!(output_file, "\t{distance}")?;
                        }
                    } else if !dist.is_empty() {
                        write!(
                            output_file,
                            "\t{}",
                            inverted_index.sample_at(dist[0] as usize)
                        )?;
                        for r_name in dist
                            .iter()
                            .skip(1)
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
                skd,
                count,
                output,
                mut knn,
                ani,
                threads,
                completeness_file,
                completeness_cutoff,
            } => {
                check_and_set_threads(*threads);

                // Load the inverted index
                let input_prefix = strip_sketch_extension(ski);
                let inverted_index = Inverted::load(input_prefix)?;

                // Two mutually exclusive modes
                if *count {
                    // For count, count the total number of pairs prefilter yields
                    // Note this can be high memory and relatively long running (~90m and 50Gb for 661k samples, 32 threads)
                    let prefilter_pairs = inverted_index.any_shared_bin_list(args.quiet);
                    println!(
                        "Identified {} prefilter pairs from a max of {}",
                        prefilter_pairs.len(),
                        inverted_index.sample_names().len()
                            * (inverted_index.sample_names().len() - 1)
                            / 2
                    );
                } else if let Some(ref_db_input) = skd {
                    let mut output_file = set_ostream(output);

                    // Open the .skq
                    let skq_filename = &format!("{input_prefix}.skq");
                    log::info!("Loading queries from {skq_filename}");
                    let (mmap, bin_stride, kmer_stride, sample_stride) =
                        (false, 1, 1, inverted_index.sketch_size());
                    let mut skq_reader = SketchArrayReader::open(
                        skq_filename,
                        mmap,
                        bin_stride,
                        kmer_stride,
                        sample_stride,
                    );
                    let skq_bins =
                        skq_reader.read_all_from_skq(sample_stride * inverted_index.sketch_size());

                    // Load the .skd/.skm
                    let ref_db_name = utils::strip_sketch_extension(ref_db_input);
                    let mut references =
                        MultiSketch::load_metadata(ref_db_name).unwrap_or_else(|_| {
                            panic!("Could not read sketch metadata from {ref_db_name}.skm")
                        });
                    log::info!("Loading sketch data from {ref_db_name}.skd");
                    references.read_sketch_data(ref_db_name);
                    log::info!("Read reference sketches:\n{references:?}");
                    let n = references.number_samples_loaded();
                    if knn >= n {
                        log::warn!("knn={knn} is higher than number of samples={n}");
                        knn = n - 1;
                    }

                    // Check that k-mer exists in the .skd, and find its index
                    let kmer = inverted_index.kmer();
                    // This panics if k not found. Maybe more graceful error if this happens
                    let dist_type = set_k(&references, Some(kmer), *ani).unwrap_or_else(|e| {
                        panic!("K-mer size {kmer} used for .ski not found in .skd: {e}");
                    });

                    // Read in completeness (parallel implementation)
                    let completeness_vec: Option<Vec<f64>> =
                        if let Some(file_path) = completeness_file {
                            Some(read_completeness_file(file_path, &references)?)
                        } else {
                            None
                        };
                    // Run the distances with both indexes
                    log::info!(
                        "Calculating sparse ref vs ref distances with {knn} nearest neighbours"
                    );
                    log::info!(
                        "Preclustering with k={} and s={}",
                        kmer,
                        inverted_index.sketch_size()
                    );
                    let distances = self_dists_knn_precluster(
                        &references,
                        &inverted_index,
                        &skq_bins,
                        skq_reader.sample_stride,
                        n,
                        knn,
                        dist_type,
                        args.quiet,
                        completeness_vec.as_ref(),
                        *completeness_cutoff,
                    );

                    // Write the results
                    log::info!("Writing out in sparse matrix form");
                    write!(output_file, "{distances}").expect("Error writing output distances");
                }

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
            let input_files: Vec<(String, Vec<String>)> = get_input_list(file_list, seq_files);
            log::info!("Parsed {} samples in input list", input_files.len());

            //check if any of the new files are already existant in the db
            let db_metadata: MultiSketch = MultiSketch::load_metadata(db)?;

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
            log::info!("Merging and saving sketch data to {output}.skd");

            let mut output_file = OpenOptions::new()
                .create(true)
                .append(true)
                .open(format!("{output}.skd"))?;
            // stream sketch data directly to concat output file
            let mut db_sketch = File::open(format!("{db}.skd"))?;
            copy(&mut db_sketch, &mut output_file)?;

            // merge and update skm from db1 and the new just sketched sketch
            let concat_metadata = db2_metadata.merge_sketches(&db_metadata);
            concat_metadata
                .save_metadata(output)
                .unwrap_or_else(|_| panic!("Could not save metadata to {output}"));
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
            let mut sketches: MultiSketch = MultiSketch::load_metadata(ref_db)
                .unwrap_or_else(|_| panic!("Could not read sketch metadata from {ref_db}.skm"));

            // write new .skm
            sketches.remove_metadata(output_file, &ids)?;

            // remove samples from .skd file
            log::info!("Remove genomes and writing output");
            sketches.remove_genomes(ref_db, output_file, &ids)?;

            log::info!("Finished writing filtered sketch data to {output_file}");

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
                let sketches = MultiSketch::load_metadata(ref_db_name).unwrap_or_else(|_| {
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

// WASM implementation
#[cfg(target_arch = "wasm32")]
#[doc(hidden)]
pub fn main() {
    panic!("You've compiled sketchlib.rust for WebAssembly support, you cannot use it as a normal binary anymore!");
}

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(js_namespace = console)]
    fn log(s: &str);
}

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
/// Function that allows to propagate panic error messages when compiling to wasm, see https://github.com/rustwasm/console_error_panic_hook
pub fn init_panic_hook() {
    console_error_panic_hook::set_once();
}

#[cfg(target_arch = "wasm32")]
/// Logging wrapper function for the WebAssembly version
pub fn logw(text: &str, typ: Option<&str>) {
    if let Some(thetyp) = typ {
        log((String::from("sketchlib.rust::") + thetyp + "::" + text).as_str());
    } else {
        log(text);
    }
}

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
/// Struct to interact with JS when working with WebAssembly
pub struct SketchlibData {
    out_probs: Vec<(f64, usize)>,
    index: Inverted,
}

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
impl SketchlibData {
    /// Constructor of the SketchlibData struct
    pub fn new(skifile: web_sys::File) -> Self {
        let inverted_index = Inverted::load(&skifile).expect("Failed loading Sketchlib index");

        logw(
            format!("Read inverted index:\n{inverted_index:?}").as_str(),
            Some("info"),
        );

        Self {
            out_probs: Vec::new(),
            index: inverted_index,
        }
    }

    /// Query some files against an inverted index
    pub fn query(&mut self, file1: web_sys::File, file2: Option<web_sys::File>) {
        // TEMPORAL BEGIN
        let min_count = &DEFAULT_MINCOUNT;
        let min_qual = &DEFAULT_MINQUAL;
        let query_type = &InvertedQueryType::MatchCount;
        // TEMPORAL END

        // Get input files
        let (queries, _query_names) =
            self.index
                .sketch_queries((&file1, file2.as_ref()), *min_count, *min_qual, false);

        logw(
            format!("Running query in mode: {query_type}").as_str(),
            Some("info"),
        );

        // Query loop (parallelised)
        let dist = match query_type {
            InvertedQueryType::MatchCount => self
                .index
                .query_against_inverted_index(queries[0].as_slice()),
            InvertedQueryType::AllBins => self.index.all_shared_bins(queries[0].as_slice()),
            InvertedQueryType::AnyBins => self.index.any_shared_bins(queries[0].as_slice()),
        };

        let mut outvec: Vec<(f64, usize)> = Vec::with_capacity(dist.len());

        for (i, d) in dist.iter().enumerate() {
            outvec.push((
                (*d as f64) / ((2 * self.index.sketch_size()) as f64 - *d as f64),
                i,
            ));
        }

        outvec.sort_by(|a, b| a.0.partial_cmp(&b.0).expect("NaN obtained!"));
        outvec.reverse();

        self.out_probs = outvec;
    }

    /// Mapping function.
    pub fn get_probs(&self, nouts: usize) -> String {
        if self.out_probs.is_empty() {
            panic!("No probabilities calculated!");
        }

        let mut results = json::JsonValue::new_array();

        logw(
            format!("Probabilities: {:?}", self.out_probs).as_str(),
            Some("info"),
        );

        results["probs"] = json::JsonValue::Array(
            self.out_probs
                .iter()
                .take(nouts)
                .map(|x| json::JsonValue::Number(x.0.into()))
                .collect(),
        );
        results["names"] = json::JsonValue::Array(
            self.out_probs
                .iter()
                .take(nouts)
                .map(|x| {
                    if let Some(labelsvec) = self.index.get_sample_labels() {
                        json::JsonValue::String(labelsvec[x.1].clone())
                    } else {
                        json::JsonValue::String("".to_string())
                    }
                })
                .collect(),
        );
        results["metadata"] = json::JsonValue::Array(
            self.out_probs
                .iter()
                .take(nouts)
                .map(|x| {
                    if let Some(metadatavec) = self.index.get_metadata() {
                        json::JsonValue::String(metadatavec[x.1].clone())
                    } else {
                        json::JsonValue::String("".to_string())
                    }
                })
                .collect(),
        );

        logw(results.dump().as_str(), Some("debug"));

        results.dump()
    }
}
