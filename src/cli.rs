//! Command line interface, built using [`crate::clap` with `Derive`](https://docs.rs/clap/latest/clap/_derive/_tutorial/index.html)
use clap::{ArgGroup, Args, Parser, Subcommand};

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
pub fn check_threads(threads: usize) {
    let max_threads = num_cpus::get();
    if threads > max_threads {
        log::warn!("{threads} threads is greater than available cores {max_threads}");
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
    },
    /// Create sketches from input data and store in an inverted index structure
    Inverted {
        /// List of input FASTA files
        #[arg(group = "input")]
        seq_files: Option<Vec<String>>,

        /// File listing input files (tab separated name, sequences, see README)
        #[arg(short, group = "input")]
        file_list: Option<String>,

        /// Output filename for the merged sketch
        #[arg(required = true, short)]
        output: String,

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

        /// aaHash 'level'
        #[arg(long, value_enum, default_value_t = DEFAULT_LEVEL)]
        level: AaLevel,
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

/// Function to parse command line args into [`MainArgs`] struct
pub fn cli_args() -> MainArgs {
    MainArgs::parse()
}
