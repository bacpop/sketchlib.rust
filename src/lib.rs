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
//! - `any-bin`. Gives samples which have at least one bin matching with the query.
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

use std::io::Write;
use std::sync::mpsc;
use std::time::Instant;

#[macro_use]
extern crate arrayref;
extern crate num_cpus;
use anyhow::Error;
use hashbrown::HashSet;
use indicatif::ParallelProgressIterator;
use rayon::prelude::*;

pub mod cli;
use crate::cli::*;

use crate::hashing::HashType;

pub mod sketch;
use crate::sketch::Sketch;
use crate::sketch::multisketch::MultiSketch;
use crate::sketch::sketch_datafile::SketchArrayReader;
use crate::sketch::{num_bins, sketch_files, SIGN_MOD};

pub mod inverted;
use crate::inverted::Inverted;

pub mod distances;
use crate::distances::*;

pub mod io;
use crate::io::{get_input_list, parse_kmers, read_subset_names, reorder_input_files, set_ostream};
pub mod structures;

pub mod hashing;
// TODO: for containment. Maybe containment.rs
use crate::hashing::{nthash_iterator::NtHashIterator, bloom_filter::KmerFilter, RollHash};

pub mod utils;
use crate::utils::{get_progress_bar, strip_sketch_extension};

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
        simple_logger::init_with_level(log::Level::Debug).unwrap();
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
            let transposed = true;
            let sketch_vec = MultiSketch::new(&mut sketches, sketch_bins, &kmers, seq_type, transposed);
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
                            let distances = self_dists_all(&references, n, dist_type, args.quiet);
                            log::info!("Writing out in long matrix form");
                            write!(output_file, "{distances}")
                                .expect("Error writing output distances");
                        }
                        Some(nn) => {
                            // Self mode (sparse)
                            log::info!("Calculating sparse ref vs ref distances with {nn} nearest neighbours");
                            let distances =
                                self_dists_knn(&references, n, nn, dist_type, args.quiet);

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
                    let distances =
                        self_query_dists_all(&references, &query_db, n, nq, dist_type, args.quiet);

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
        Commands::Containment { command } => match command {
            ContainmentCommands::Index {
                seq_files,
                file_list,
                output,
                kmers,
                n_downsample,
                min_count,
                min_qual,
                threads,
            } => {
                todo!();
                // TODO: in read filter, bloom filter size should be based on downsampling rate
                // TODO: use/edit the inverted interface to get signs. Pack 4xu16 in a u64
                // in fact, just use inverted interface rather than messing with multisketch -- can turn off inverted part and just write the .skq
                // TODO: support streaming of reads (probably just on by default). Also consider how to downsample
                // TODO: could add abundance, for angular similarity support
            }
            ContainmentCommands::Query {
                skq_file,
                ref_inverted,
                output,
                single_strand,
                min_count,
                min_qual,
                threads,
            } => {
                check_and_set_threads(*threads);

                let mut output_file = set_ostream(output);
                let rc = !*single_strand;

                // Open the reference inverted index
                let inverted_index = Inverted::load(strip_sketch_extension(ref_inverted))?;
                log::info!("Read inverted index:\n{inverted_index:?}");
                let inverted_binsize: u64 = SIGN_MOD.div_ceil(inverted_index.sketch_size() as u64);

                // Open the metagenome .skq
                log::info!("Loading queries from {skq_file}");
                let sketch_size = 28423956;
                let kmer = 21;
                log::warn!("skq sketch size hard-coded as {sketch_size}; kmer as {kmer}");
                let (mmap, bin_stride, kmer_stride, sample_stride) =
                    (false, 1, 1, sketch_size);
                let mut skq_reader = SketchArrayReader::open(
                    skq_file,
                    mmap,
                    bin_stride,
                    kmer_stride,
                    sample_stride,
                );
                let skq_bins =
                    skq_reader.read_all_from_skq(sample_stride);
                let skq_binsize: u64 = SIGN_MOD.div_ceil(sample_stride as u64);
                log::debug!("inverted_binsize:{inverted_binsize} skq_binsize:{skq_binsize}");

                log::info!("Running queries");
                // Header
                write!(output_file, "Query")?;
                for name in inverted_index.sample_names() {
                    write!(output_file, "\t{name}")?;
                }
                writeln!(output_file)?;

                let mut match_counts = vec![0; inverted_index.n_samples()];
                let index = inverted_index.inverted_index_ref();

                for (query_bin_idx, query_bin_hash) in skq_bins.iter().enumerate() {
                    let inverted_bin_idx = ((query_bin_idx as u64 * skq_binsize) / inverted_binsize) as usize;
                    if let Some(matching_samples) = index[inverted_bin_idx].get(query_bin_hash) {
                        for sample_idx in matching_samples {
                            match_counts[sample_idx as usize] += 1;
                        }
                    }
                }

                write!(output_file, "SRR6269135")?;
                for distance in match_counts {
                    write!(output_file, "\t{}", distance as f64 / inverted_index.sketch_size() as f64)?;
                }
                writeln!(output_file)?;

                /*
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

                // TODO Need to loop over k-mers in loop below
                let k_idx = references.get_k_idx(*kmer).expect(&format!("K-mer size {kmer} not found in file"));
                let rc = !*single_strand;

                log::info!("Getting input files");
                let input_files = get_input_list(query_file_list, query_seq_files);
                log::info!("Parsed {} samples in input list", input_files.len());

                let percent = false;
                let progress_bar = get_progress_bar(input_files.len(), percent, args.quiet);
                let containment_vec: Vec<Vec<f64>> = input_files
                    .par_iter()
                    .progress_with(progress_bar)
                    .map(|(name, fastx1, fastx2)| {
                        // TODO: eventually add a query pre-sketch command (same as inverted write-skq)
                        let mut hash_it = NtHashIterator::new((fastx1, fastx2.as_ref()), rc, *min_qual);
                        let hash_it = hash_it.first_mut().expect("Empty hash iterator");

                        if hash_it.seq_len() == 0 {
                            panic!("{name} has no valid sequence");
                        }

                        let mut read_filter = if hash_it.reads() {
                            let mut filter = KmerFilter::new(*min_count);
                            filter.init();
                            Some(filter)
                        } else {
                            None
                        };

                        let signs =
                            Sketch::get_all_signs(hash_it, *kmer, &mut read_filter);

                        // TODO: fix this, sketches are transposed and 14 bits, new sketches are full u64
                        // TODO: use the distance matrix class to store and print results
                        let mut containments = Vec::with_capacity(references.number_samples_loaded());
                        for reference_idx in 0..references.number_samples_loaded() {
                            let mut intersection = 0;
                            let ref_bins = references.get_sketch_slice(reference_idx, k_idx);
                            for bin in ref_bins {
                                if signs.contains(bin) {
                                    intersection += 1;
                                }
                            }
                            containments.push(intersection as f64 / signs.len() as f64); // |A ∩ B| / |A|
                        }
                        containments

                    }).collect();

                containment_vec.iter().zip(input_files).for_each(|(c_vals, (name, _, _))| {
                    c_vals.iter().enumerate().for_each(|(query_idx, containment)| {
                        println!("{name}\t{}\t{containment}", references.sketch_name(query_idx));
                    });
                });
                */
                Ok(())
            }
        }
        Commands::Inverted { command } => match command {
            InvertedCommands::Build {
                seq_files,
                file_list,
                output,
                write_skq,
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
                let input_files: Vec<(String, String, Option<String>)> =
                    get_input_list(file_list, seq_files);
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
                            write!(output_file, "\t{distance}")?;
                        }
                    } else {
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
            let input_files: Vec<(String, String, Option<String>)> =
                get_input_list(file_list, seq_files);
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
                MultiSketch::new(&mut db2_sketches, sketch_size, kmers, seq_type, db_metadata.transposed);

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
