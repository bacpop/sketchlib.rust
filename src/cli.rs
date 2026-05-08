//! Command line interface, built using [`crate::clap` with `Derive`](https://docs.rs/clap/latest/clap/_derive/_tutorial/index.html)
use core::fmt;

use clap::{ArgGroup, Args, Parser, Subcommand, ValueEnum};

use crate::DEFAULT_KMER;

use super::hashing::{AaLevel, HashType, DEFAULT_LEVEL};

/// Default single strand (which is equivalent to !rc)
pub const DEFAULT_STRAND: bool = false;
/// Default minimum k-mer count for FASTQ files
pub const DEFAULT_MINCOUNT: u16 = 5;
/// Default minimum base quality (PHRED score) for FASTQ files
pub const DEFAULT_MINQUAL: u8 = 20;
/// Default sketch size
pub const DEFAULT_SKETCHSIZE: u64 = 1000;
/// Default nearest neighbours
pub const DEFAULT_KNN: usize = 50;

/// Query types supported by bitvec operations
#[derive(Clone, Debug, PartialEq, PartialOrd, ValueEnum, Default)]
pub enum InvertedQueryType {
    #[default]
    /// Count the number of matching bins for each sample in the index
    MatchCount,
    /// Return samples which match at every bin
    AllBins,
    /// Return samples which match at least one bin
    AnyBins,
}

impl fmt::Display for InvertedQueryType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            InvertedQueryType::MatchCount => write!(f, "Count of matching bins")?,
            InvertedQueryType::AllBins => write!(f, "All bins matching")?,
            InvertedQueryType::AnyBins => write!(f, "At least one bin matching")?,
        }
        Ok(())
    }
}

#[doc(hidden)]
fn valid_cpus(s: &str) -> Result<usize, String> {
    let threads: usize = s
        .parse()
        .map_err(|_| format!("`{s}` isn't a valid number of cores"))?;
    if threads < 1 {
        Err("Threads must be one or higher".to_string())
    } else {
        Ok(threads)
    }
}

/// Prints a warning if more threads than available have been requested
pub fn check_and_set_threads(threads: usize) {
    let max_threads = num_cpus::get();
    if threads > max_threads {
        log::warn!("{threads} threads is greater than available cores {max_threads}");
    } else {
        log::info!("Using {threads} threads");
    }
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();
}

/// Options that apply to all subcommands
#[derive(Parser)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
pub struct MainArgs {
    #[doc(hidden)]
    #[command(subcommand)]
    pub command: Commands,

    /// Show progress messages
    #[arg(short, long, global = true)]
    pub verbose: bool,

    /// Don't show any messages
    #[arg(long, global = true)]
    pub quiet: bool,
}

/// K-mer sequence, or single k-mer value, definitions from the CLI
#[derive(Args)]
#[group(required = true, multiple = false)]
pub struct Kmers {
    /// K-mer list (comma separated k-mer values to sketch at)
    #[arg(short, long, required = true, value_delimiter = ',')]
    pub k_vals: Option<Vec<usize>>,

    /// K-mer linear sequence (start,end,step)
    #[arg(long, required = true, value_delimiter = ',')]
    pub k_seq: Option<Vec<usize>>,
}

/// Subcommands and their specific options
#[derive(Subcommand)]
pub enum Commands {
    #[command(group(
        ArgGroup::new("input")
            .required(true)
            .args(["seq_files", "file_list"]),
    ))]

    /// Create sketches from input data
    Sketch {
        /// List of input FASTA files
        #[arg(group = "input")]
        seq_files: Option<Vec<String>>,

        /// File listing input files (tab separated name, sequences, see README)
        #[arg(short, group = "input")]
        file_list: Option<String>,

        /// Treat every sequence in an input file as a new sample (aa only)
        // TODO: for now, could be extended to dna, but probably no need
        #[arg(long, default_value_t = false)]
        concat_fasta: bool,

        /// Input files are .pdb, convert them to 3Di first
        #[cfg(feature = "3di")]
        #[arg(long, default_value_t = false)]
        convert_pdb: bool,

        /// Output prefix
        #[arg(short)]
        output: String,

        /// K-mers to sketch
        #[command(flatten)]
        kmers: Kmers,

        /// Sketch size
        #[arg(short, long, default_value_t = DEFAULT_SKETCHSIZE)]
        sketch_size: u64,

        /// Type of sequence to hash
        #[arg(long, value_enum, default_value_t = HashType::DNA)]
        seq_type: HashType,

        /// aaHash 'level'
        #[arg(long, value_enum, default_value_t = DEFAULT_LEVEL)]
        level: AaLevel,

        /// Ignore reverse complement (all contigs are oriented along same strand)
        #[arg(long, default_value_t = DEFAULT_STRAND)]
        single_strand: bool,

        /// Minimum k-mer count (with reads)
        #[arg(long, default_value_t = DEFAULT_MINCOUNT)]
        min_count: u16,

        /// Minimum k-mer quality (with reads)
        #[arg(long, default_value_t = DEFAULT_MINQUAL)]
        min_qual: u8,

        /// Number of CPU threads
        #[arg(long, value_parser = valid_cpus, default_value_t = 1)]
        threads: usize,
    },

    /// Calculate pairwise distances using sketches
    Dist {
        /// The .skm file used as the reference
        #[arg(required = true)]
        ref_db: String,

        /// The .skm file used as the query (omit for ref v ref)
        #[arg(group = "query")]
        query_db: Option<String>,

        /// Output filename (omit to output to stdout)
        #[arg(short)]
        output: Option<String>,

        /// Calculate sparse distances with k nearest-neighbours
        #[arg(long, group = "query")]
        knn: Option<usize>,

        /// Sample names to analyse
        #[arg(long)]
        subset: Option<String>,

        /// K-mer length (if provided only calculate Jaccard distance)
        #[arg(short)]
        kmer: Option<usize>,

        /// Calculate ANI rather than Jaccard dists, using Poisson model
        #[arg(long, requires("kmer"), default_value_t = false)]
        ani: bool,

        /// Number of CPU threads
        #[arg(long, value_parser = valid_cpus, default_value_t = 1)]
        threads: usize,

        /// File listing sample and completeness estimate 0.0-1.0 (tab separated)
        #[arg(long)]
        completeness_file: Option<String>,

        /// minimum completeness for a sample to be corrected but the completeness correction
        #[arg(long, default_value_t = 0.64)]
        completeness_cutoff: f64,
    },

    /// Building and querying with inverted indices (.ski)
    Inverted {
        /// Interactions with inverted indices
        #[command(subcommand)]
        command: InvertedCommands,
    },

    /// Merge two sketch files (.skm and .skd pair)
    Merge {
        /// The first .skd (sketch data) file
        #[arg(required = true)]
        db1: String,

        /// The second .skd (sketch data) file
        #[arg(required = true)]
        db2: String,

        /// Output filename for the merged sketch
        #[arg(required = true, short)]
        output: String,
    },

    /// Append new genomes to be sketched to an existing sketch database
    Append {
        /// Sketching database basename (so without .skm or .skd)
        #[arg(required = true)]
        db: String,

        /// List of input FASTA files
        #[arg(group = "input")]
        seq_files: Option<Vec<String>>,

        /// File listing input files (tab separated name, sequences)
        #[arg(short, group = "input")]
        file_list: Option<String>,

        /// Output filename for the merged sketch
        #[arg(required = true, short)]
        output: String,

        /// Ignore reverse complement (all contigs are oriented along same strand)
        #[arg(long, default_value_t = DEFAULT_STRAND)]
        single_strand: bool,

        /// Minimum k-mer count (with reads)
        #[arg(long, default_value_t = DEFAULT_MINCOUNT)]
        min_count: u16,

        /// Minimum k-mer quality (with reads)
        #[arg(long, default_value_t = DEFAULT_MINQUAL)]
        min_qual: u8,

        /// Treat every sequence in an input file as a new sample (aa only)
        // TODO: for now, could be extended to dna, but probably no need
        #[arg(long, default_value_t = false)]
        concat_fasta: bool,

        /// Number of CPU threads
        #[arg(long, value_parser = valid_cpus, default_value_t = 1)]
        threads: usize,

        /// aaHash 'level'
        #[arg(long, value_enum, default_value_t = DEFAULT_LEVEL)]
        level: AaLevel,
    },

    /// Delete genome(s) from a database (input: one id per line)
    Delete {
        /// Sketching database
        #[arg(required = true)]
        db: String,

        /// Input file with IDs to delete (one ID per line)
        #[arg(required = true)]
        samples: String,

        /// output file name
        #[arg(required = true)]
        output_file: String,
    },

    /// Print information about a .skm file
    Info {
        /// Sketch metadata file (.skm) to describe
        skm_file: String,

        /// Write out the information for every sample contained
        #[arg(long, default_value_t = false)]
        sample_info: bool,
    },
}

/// Commands to support building and querying with an inverted index
#[derive(Subcommand)]
pub enum InvertedCommands {
    /// Create sketches from input data and store in an inverted index structure
    Build {
        /// List of input FASTA files
        #[arg(group = "input")]
        seq_files: Option<Vec<String>>,

        /// File listing input files (tab separated name, sequences, see README)
        #[arg(short, group = "input")]
        file_list: Option<String>,

        /// Output filename for the merged sketch
        #[arg(required = true, short)]
        output: String,

        /// Also write an .skq file, which is needed by 'precluster' mode
        #[arg(long)]
        write_skq: bool,

        /// File listing species names, or clusters, for phylogenetic ordering
        #[arg(long)]
        species_names: Option<String>,

        /// File listing species names, or clusters, for phylogenetic ordering
        #[arg(long)]
        metadata: Option<String>,

        /// Sketch size
        #[arg(short, long, default_value_t = DEFAULT_SKETCHSIZE)]
        sketch_size: u64,

        /// K-mer size
        #[arg(short, long, default_value_t = DEFAULT_KMER)]
        kmer_length: usize,

        /// Ignore reverse complement (all contigs are oriented along same strand)
        #[arg(long, default_value_t = DEFAULT_STRAND)]
        single_strand: bool,

        /// Minimum k-mer count (with reads)
        #[arg(long, default_value_t = DEFAULT_MINCOUNT)]
        min_count: u16,

        /// Minimum k-mer quality (with reads)
        #[arg(long, default_value_t = DEFAULT_MINQUAL)]
        min_qual: u8,

        /// Number of CPU threads
        #[arg(long, value_parser = valid_cpus, default_value_t = 1)]
        threads: usize,
    },

    /// Find distances against an inverted index
    Query {
        /// The inverted index (.ski) file used as the reference
        #[arg(required = true)]
        ski: String,

        /// List of input FASTA files
        #[arg(group = "input")]
        seq_files: Option<Vec<String>>,

        /// File listing input files (tab separated name, sequences, see README)
        #[arg(short, group = "input")]
        file_list: Option<String>,

        /// Output filename (omit to output to stdout)
        #[arg(short)]
        output: Option<String>,

        /// Type of query to perform
        #[arg(long, value_enum, default_value_t = InvertedQueryType::MatchCount)]
        query_type: InvertedQueryType,

        /// Minimum k-mer count (with reads)
        #[arg(long, default_value_t = DEFAULT_MINCOUNT)]
        min_count: u16,

        /// Minimum k-mer quality (with reads)
        #[arg(long, default_value_t = DEFAULT_MINQUAL)]
        min_qual: u8,

        /// Number of CPU threads
        #[arg(long, value_parser = valid_cpus, default_value_t = 1)]
        threads: usize,
    },

    #[command(group(
        ArgGroup::new("mode")
            .required(true)
            .args(["skd", "count"]),
    ))]
    /// Use an inverted index to reduce query comparisons
    Precluster {
        /// The inverted index (.ski) file used as the reference
        #[arg(required = true)]
        ski: String,

        /// The .skd/.skm file prefix, which must have the same samples and k-mer length as the .ski file
        #[arg(long, group = "mode")]
        skd: Option<String>,

        /// Output filename (omit to output to stdout)
        #[arg(short)]
        output: Option<String>,

        /// Do not run analysis, only return the number of comparisons
        #[arg(long, group = "mode")]
        count: bool,

        /// Reduce to a maximum of k nearest-neighbours
        #[arg(long, default_value_t = DEFAULT_KNN)]
        knn: usize,

        /// Calculate ANI rather than Jaccard dists, using Poisson model
        #[arg(long, default_value_t = false)]
        ani: bool,

        /// Number of CPU threads
        #[arg(long, value_parser = valid_cpus, default_value_t = 1)]
        threads: usize,

        /// Completeness file
        #[arg(long)]
        completeness_file: Option<String>,

        /// minimum completeness for a sample to be corrected but the completeness correction
        #[arg(long, default_value_t = 0.64)]
        completeness_cutoff: f64,
    },
}

/// Function to parse command line args into [`MainArgs`] struct
pub fn cli_args() -> MainArgs {
    MainArgs::parse()
}
