//! Fast distance calculations between biological sequences (DNA or AA).
//! Distances are based on bindash approximations of the Jaccard
//! distance, with the [PopPUNK method](https://poppunk.bacpop.org/index.html) to calculate core and accessory distances. nthash/aahash
//! are used for hash functions to create the sketches.
//!
//! ## Important biological considerations
//!
//! - Core/accessory distances are only tested within-species (>95% ANI). Using input
//!   above these distances is unsupported and may lead to poor estimation of distances
//!   without clear warning.
//! - Short k-mer lengths are likely to match at random, see [PopPUNK's docs](https://poppunk-docs.bacpop.org/sketching.html#choosing-the-right-k-mer-lengths)
//!   for information on how to select good lengths. Note that this library does not support
//!   random match correction.
//! - ANI distance resolution is highly affected by sketch size at higher mismatch
//!   levels, so note that if you see lots of samples at around 80% they may be much lower than
//!   this. We recommend checking the Jaccard values in this case, if they are close to 0
//!   you should increase the sketch size.
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
#![warn(missing_docs)]
#![allow(clippy::too_many_arguments)]

#[cfg(not(target_arch = "wasm32"))]
use std::time::Instant;

#[macro_use]
extern crate arrayref;
extern crate num_cpus;
use anyhow::Error;
#[cfg(not(target_arch = "wasm32"))]
pub mod api;
pub mod cli;
#[cfg(all(not(target_arch = "wasm32"), feature = "python"))]
pub mod python;
#[cfg(not(target_arch = "wasm32"))]
use crate::cli::*;
#[cfg(target_arch = "wasm32")]
use crate::cli::{InvertedQueryType, DEFAULT_MINCOUNT, DEFAULT_MINQUAL};

#[cfg(not(target_arch = "wasm32"))]
use crate::hashing::HashType;

pub mod sketch;
#[cfg(not(target_arch = "wasm32"))]
use crate::sketch::multisketch::MultiSketch;
#[cfg(not(target_arch = "wasm32"))]
use crate::sketch::sketch_files;

pub mod inverted;
use crate::inverted::Inverted;

pub mod distances;

pub mod io;

#[cfg(not(target_arch = "wasm32"))]
use crate::io::{get_input_list, parse_kmers};

pub mod hashing;

pub mod utils;
pub use crate::utils::get_progress_bar;
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
            api::sketch(api::SketchOptions {
                seq_files: seq_files.clone(),
                file_list: file_list.clone(),
                output: output.clone(),
                kmers: parse_kmers(kmers),
                sketch_size: *sketch_size,
                seq_type: seq_type.clone(),
                level: level.clone(),
                concat_fasta: *concat_fasta,
                single_strand: *single_strand,
                min_count: *min_count,
                min_qual: *min_qual,
                threads: *threads,
                quiet: args.quiet,
            })?;
            Ok(())
        }
        Commands::Dist {
            ref_db,
            query_db,
            output,
            knn,
            subset,
            kmer,
            ani,
            threads,
            ref_completeness_file,
            query_completeness_file,
            completeness_cutoff,
        } => {
            api::dist(api::DistOptions {
                ref_db: ref_db.clone(),
                query_db: query_db.clone(),
                output: output.clone(),
                knn: *knn,
                subset: subset.clone(),
                kmer: *kmer,
                ani: *ani,
                threads: *threads,
                ref_completeness_file: ref_completeness_file.clone(),
                query_completeness_file: query_completeness_file.clone(),
                completeness_cutoff: *completeness_cutoff,
                quiet: args.quiet,
            })?;
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
                api::inverted_build(api::InvertedBuildOptions {
                    seq_files: seq_files.clone(),
                    file_list: file_list.clone(),
                    output: output.clone(),
                    write_skq: *write_skq,
                    species_names: species_names.clone(),
                    metadata: metadata.clone(),
                    sketch_size: *sketch_size,
                    kmer_length: *kmer_length,
                    single_strand: *single_strand,
                    min_count: *min_count,
                    min_qual: *min_qual,
                    threads: *threads,
                    quiet: args.quiet,
                })?;
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
                api::inverted_query(
                    ski, seq_files, file_list, output, query_type, *min_count, *min_qual, *threads,
                    args.quiet,
                )?;
                Ok(())
            }
            InvertedCommands::Precluster {
                ski,
                skd,
                count,
                output,
                knn,
                ani,
                threads,
                ref_completeness_file,
                completeness_cutoff,
                retain_unmatched,
            } => {
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
                    api::inverted_precluster(api::InvertedPreclusterOptions {
                        ski: ski.clone(),
                        skd: ref_db_input.clone(),
                        output: output.clone(),
                        knn: *knn,
                        ani: *ani,
                        threads: *threads,
                        ref_completeness_file: ref_completeness_file.clone(),
                        completeness_cutoff: *completeness_cutoff,
                        retain_unmatched: retain_unmatched.clone(),
                        quiet: args.quiet,
                    })?;
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
            if *concat_fasta && matches!(*seq_type, HashType::DNA) {
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
            println!("{}", api::db_info(skm_file, *sample_info)?);
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
