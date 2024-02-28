//! Command line interface, built using [`crate::clap` with `Derive`](https://docs.rs/clap/latest/clap/_derive/_tutorial/index.html)
use clap::{ArgGroup, Parser, Subcommand};

/* C++ interface

Usage:
  sketchlib sketch <files>... -o <output> [-k <kseq>|--kmer <k>] [-s <size>] [--single-strand] [--codon-phased] [--min-count <count>] [--exact-counter] [--cpus <cpus>] [--gpu <gpu>]
  sketchlib sketch -l <file-list> -o <output> [-k <kseq>|--kmer <k>] [-s <size>] [--single-strand] [--codon-phased] [--min-count <count>] [--exact-counter] [--cpus <cpus>] [--gpu <gpu>]
  sketchlib query dist <db1> [<db2>] [-o <output>] [--adj-random] [--subset <file>] [--cpus <cpus>] [--gpu <gpu>]
  sketchlib query jaccard <db1> [<db2>] [-o <output>] [--kmer <k>] [--adj-random] [--subset <file>] [--cpus <cpus>]
  sketchlib query sparse <db1> (--kNN <k>|--threshold <max>) [-o <output>] [--accessory] [--adj-random] [--subset <file>] [--cpus <cpus>] [--gpu <gpu>]
  sketchlib query sparse jaccard <db1> --kNN <k> --kmer <k> [-o <output>] [--adj-random] [--subset <file>] [--cpus <cpus>]
  sketchlib join <db1> <db2> -o <output>
  sketchlib (add|remove) random <db1> [--single-strand] [--cpus <cpus>]
  sketchlib (-h | --help)
  sketchlib (--version)

Options:
  -h --help     Show this help.
  --version     Show version.

  -o <output>    Output prefix.
  -l <file-list> File with a list of input files.
  --cpus <cpus>  Number of CPU threads to use [default: 1].
  --gpu <gpu>    Use GPU with specified device ID [default: -1].

  -k <kseq>     Sequence of k-mers to sketch (min,max,step) [default: 15,31,4].
  --kmer <k>    Sketch (or distance) at a single k-mer length k.
  -s <size>     Sketch size [default: 10000].
  --single-strand  Ignore the reverse complement (e.g. in RNA viruses).
  --codon-phased  Use codon phased seeds X--X--X
  --min-count <count>  Minimum coverage count for k-mers from reads to be sketched [default: 20].
  --exact-counter  Use an exact k-mer count filter for reads (for genomes >10Mb)

  --adj-random  Adjust query matches for their chance of occurring at random
  --subset <file>  Only query samples matching names in file

  --kNN <k>  Use k nearest neighbours to sparsify
  --threshold <max>  Remove distances over max to sparsify
  --accessory  Use accessory distances rather than core to sparsify
*/

/// Default split k-mer size
pub const DEFAULT_KMER: usize = 17;
/// Default single strand (which is equivalent to !rc)
pub const DEFAULT_STRAND: bool = false;
/// Default minimum k-mer count for FASTQ files
pub const DEFAULT_MINCOUNT: u16 = 5;
/// Default minimum base quality (PHRED score) for FASTQ files
pub const DEFAULT_MINQUAL: u8 = 20;

#[doc(hidden)]
fn valid_kmer(s: &str) -> Result<usize, String> {
    let k: usize = s
        .parse()
        .map_err(|_| format!("`{s}` isn't a valid k-mer"))?;
    if k < 5 {
        Err("K-mer must a number greater than five".to_string())
    } else {
        Ok(k)
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
pub fn check_threads(threads: usize) {
    let max_threads = num_cpus::get();
    if threads > max_threads {
        log::warn!("{threads} threads is greater than available cores {max_threads}");
    }
}

/// Options that apply to all subcommands
#[derive(Parser)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
pub struct Args {
    #[doc(hidden)]
    #[command(subcommand)]
    pub command: Commands,

    /// Show progress messages
    #[arg(short, long, global = true)]
    pub verbose: bool,
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

        /// File listing input files (tab separated name, sequences)
        #[arg(short, group = "input")]
        file_list: Option<String>,

        /// Output prefix
        #[arg(short)]
        output: String,

        /// K-mer size
        #[arg(short, value_parser = valid_kmer, default_value_t = DEFAULT_KMER)]
        k: usize,

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
    /// Write an unordered alignment
    Dist {
        /// A .h5 file as the reference
        #[arg(required = true)]
        db1: String,

        /// A .h5 file as the query (omit for ref v ref)
        #[arg(required = true)]
        db2: Option<String>,

        /// Output filename (omit to output to stdout)
        #[arg(short)]
        output: Option<String>,

        /// Number of CPU threads
        #[arg(long, value_parser = valid_cpus, default_value_t = 1)]
        threads: usize,
    },
}

/// Function to parse command line args into [`Args`] struct
pub fn cli_args() -> Args {
    Args::parse()
}
