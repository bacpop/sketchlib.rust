//! Reusable API functions shared by the CLI and Python bindings.

use std::fmt;
use std::io::Write;
use std::sync::mpsc;

use anyhow::{bail, Result};
use hashbrown::{HashMap, HashSet};
use indicatif::ParallelProgressIterator;
use rayon::prelude::*;

use crate::cli::{
    check_and_set_threads, InvertedQueryType, RetainUnmatched, DEFAULT_MINCOUNT, DEFAULT_MINQUAL,
};
use crate::distances::{
    cross_dists_all, cross_dists_knn, self_dists_all, self_dists_knn, self_dists_knn_precluster,
    set_k,
};
use crate::hashing::{AaLevel, HashType, DEFAULT_LEVEL};
use crate::inverted::Inverted;
use crate::io::{
    get_input_list, parse_metadata_info, read_completeness_file, read_subset_names,
    reorder_input_files, set_ostream,
};
use crate::sketch::multisketch::MultiSketch;
use crate::sketch::sketch_datafile::SketchArrayReader;
use crate::sketch::{num_bins, sketch_files};
use crate::utils::{get_progress_bar, strip_sketch_extension};

/// File paths produced by sketching.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SketchPaths {
    /// Metadata file path.
    pub skm: String,
    /// Sketch data file path.
    pub skd: String,
}

/// File paths produced by inverted index construction.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct InvertedPaths {
    /// Inverted index file path.
    pub ski: String,
    /// Optional query sketch file path.
    pub skq: Option<String>,
}

/// Basic database metadata.
#[derive(Clone, Debug, PartialEq)]
pub struct DbInfo {
    /// Human-readable metadata in the same format as the CLI `info` command.
    pub text: String,
    /// Whether this is an inverted index.
    pub inverted: bool,
    /// Sketch size.
    pub sketch_size: usize,
    /// Number of samples.
    pub n_samples: usize,
    /// K-mer sizes.
    pub kmers: Vec<usize>,
}

/// Options for [`sketch`].
#[derive(Clone, Debug)]
pub struct SketchOptions {
    /// Input sequence files, mutually exclusive with `file_list`.
    pub seq_files: Option<Vec<String>>,
    /// File listing input samples.
    pub file_list: Option<String>,
    /// Output file prefix.
    pub output: String,
    /// K-mer sizes.
    pub kmers: Vec<usize>,
    /// Requested sketch size.
    pub sketch_size: u64,
    /// Sequence type.
    pub seq_type: HashType,
    /// aaHash level.
    pub level: AaLevel,
    /// Treat each FASTA record as a sample.
    pub concat_fasta: bool,
    /// Ignore reverse complement.
    pub single_strand: bool,
    /// Minimum k-mer count.
    pub min_count: u16,
    /// Minimum base quality.
    pub min_qual: u8,
    /// Number of worker threads.
    pub threads: usize,
    /// Suppress progress output.
    pub quiet: bool,
}

impl SketchOptions {
    /// Construct DNA sketching options for a file-list input.
    pub fn dna_file_list(
        file_list: String,
        output: String,
        sketch_size: u64,
        kmer_length: usize,
        threads: usize,
    ) -> Self {
        Self {
            seq_files: None,
            file_list: Some(file_list),
            output,
            kmers: vec![kmer_length],
            sketch_size,
            seq_type: HashType::DNA,
            level: DEFAULT_LEVEL,
            concat_fasta: false,
            single_strand: false,
            min_count: DEFAULT_MINCOUNT,
            min_qual: DEFAULT_MINQUAL,
            threads,
            quiet: true,
        }
    }
}

/// Options for [`dist`].
#[derive(Clone, Debug)]
pub struct DistOptions {
    /// Reference sketch database prefix or file.
    pub ref_db: String,
    /// Optional query sketch database prefix or file.
    pub query_db: Option<String>,
    /// Output distance filename. If omitted, writes to stdout.
    pub output: Option<String>,
    /// Optional sparse nearest-neighbor count.
    pub knn: Option<usize>,
    /// Optional subset file for reference samples.
    pub subset: Option<String>,
    /// Optional single k-mer length.
    pub kmer: Option<usize>,
    /// Return ANI rather than Jaccard distances.
    pub ani: bool,
    /// Number of worker threads.
    pub threads: usize,
    /// Optional reference completeness file.
    pub ref_completeness_file: Option<String>,
    /// Optional query completeness file.
    pub query_completeness_file: Option<String>,
    /// Completeness correction cutoff.
    pub completeness_cutoff: f64,
    /// Suppress progress output.
    pub quiet: bool,
}

/// Options for [`inverted_build`].
#[derive(Clone, Debug)]
pub struct InvertedBuildOptions {
    /// Input sequence files, mutually exclusive with `file_list`.
    pub seq_files: Option<Vec<String>>,
    /// File listing input samples.
    pub file_list: Option<String>,
    /// Output file prefix.
    pub output: String,
    /// Whether to write an `.skq` file.
    pub write_skq: bool,
    /// Optional species/label file.
    pub species_names: Option<String>,
    /// Optional metadata file.
    pub metadata: Option<String>,
    /// Inverted sketch size.
    pub sketch_size: u64,
    /// K-mer length.
    pub kmer_length: usize,
    /// Ignore reverse complement.
    pub single_strand: bool,
    /// Minimum k-mer count.
    pub min_count: u16,
    /// Minimum base quality.
    pub min_qual: u8,
    /// Number of worker threads.
    pub threads: usize,
    /// Suppress progress output.
    pub quiet: bool,
}

/// Options for [`inverted_precluster`].
#[derive(Clone, Debug)]
pub struct InvertedPreclusterOptions {
    /// Inverted index file.
    pub ski: String,
    /// Standard sketch database prefix or file.
    pub skd: String,
    /// Output distance filename. If omitted, writes to stdout.
    pub output: Option<String>,
    /// Sparse nearest-neighbor count.
    pub knn: usize,
    /// Return ANI rather than Jaccard distances.
    pub ani: bool,
    /// Number of worker threads.
    pub threads: usize,
    /// Optional reference completeness file.
    pub ref_completeness_file: Option<String>,
    /// Completeness correction cutoff.
    pub completeness_cutoff: f64,
    /// How to retain unmatched samples.
    pub retain_unmatched: Option<RetainUnmatched>,
    /// Suppress progress output.
    pub quiet: bool,
}

/// Create a sketch database.
pub fn sketch(options: SketchOptions) -> Result<SketchPaths> {
    if options.concat_fasta && matches!(options.seq_type, HashType::DNA) {
        bail!("--concat-fasta currently only supported with --seq-type aa");
    }
    check_and_set_threads(options.threads + 1);

    let input_files = get_input_list(&options.file_list, &options.seq_files);
    let rc = !options.single_strand;
    let seq_type = if matches!(options.seq_type, HashType::AA(_)) {
        HashType::AA(options.level)
    } else {
        options.seq_type
    };
    let (_, sketch_bins, _) = num_bins(options.sketch_size);
    let mut sketches = sketch_files(
        &options.output,
        &input_files,
        options.concat_fasta,
        &options.kmers,
        sketch_bins,
        &seq_type,
        rc,
        options.min_count,
        options.min_qual,
        options.quiet,
    );
    let sketch_vec = MultiSketch::new(&mut sketches, sketch_bins, &options.kmers, seq_type);
    sketch_vec.save_metadata(&options.output)?;

    Ok(SketchPaths {
        skm: format!("{}.skm", options.output),
        skd: format!("{}.skd", options.output),
    })
}

/// Calculate distances from one or two sketch databases.
pub fn dist(options: DistOptions) -> Result<Option<String>> {
    check_and_set_threads(options.threads);
    let mut output_file = set_ostream(&options.output);
    write_distances(&options, &mut output_file)?;
    Ok(options.output)
}

fn write_distances(options: &DistOptions, output_file: &mut impl Write) -> Result<()> {
    let ref_db_name = strip_sketch_extension(&options.ref_db);
    let mut references = MultiSketch::load_metadata(ref_db_name)?;

    if let Some(subset_file) = &options.subset {
        let subset_names = read_subset_names(subset_file);
        references.read_sketch_data_block(ref_db_name, &subset_names);
    } else {
        references.read_sketch_data(ref_db_name);
    }
    let n = references.number_samples_loaded();
    let ref_completeness_vec = if let Some(file_path) = &options.ref_completeness_file {
        Some(read_completeness_file(file_path, &references)?)
    } else {
        None
    };

    let dist_type = set_k(&references, options.kmer, options.ani)?;
    let queries = if let Some(query_db_name) = &options.query_db {
        let query_db_name = strip_sketch_extension(query_db_name);
        let mut queries = MultiSketch::load_metadata(query_db_name)?;
        queries.read_sketch_data(query_db_name);
        Some(queries)
    } else {
        None
    };

    match queries {
        None => match options.knn {
            None => {
                let distances = self_dists_all(
                    &references,
                    n,
                    dist_type,
                    options.quiet,
                    ref_completeness_vec.as_ref(),
                    options.completeness_cutoff,
                );
                write!(output_file, "{distances}")?;
            }
            Some(mut nn) => {
                if nn >= n {
                    log::warn!("knn={nn} is higher than number of samples={n}");
                    nn = n - 1;
                }
                let distances = self_dists_knn(
                    &references,
                    n,
                    nn,
                    dist_type,
                    options.quiet,
                    ref_completeness_vec.as_ref(),
                    options.completeness_cutoff,
                );
                write!(output_file, "{distances}")?;
            }
        },
        Some(query_db) => {
            let query_completeness_vec = if let Some(file_path) = &options.query_completeness_file {
                Some(read_completeness_file(file_path, &query_db)?)
            } else {
                None
            };
            let n_query = query_db.number_samples_loaded();
            match options.knn {
                Some(mut nn) => {
                    if nn > n {
                        log::warn!("knn={nn} is higher than number of reference samples={n}");
                        nn = n;
                    }
                    let distances = cross_dists_knn(
                        &references,
                        &query_db,
                        n,
                        n_query,
                        nn,
                        dist_type,
                        options.quiet,
                        ref_completeness_vec.as_ref(),
                        query_completeness_vec.as_ref(),
                        options.completeness_cutoff,
                    );
                    write!(output_file, "{distances}")?;
                }
                None => {
                    let distances = cross_dists_all(
                        &references,
                        &query_db,
                        n,
                        n_query,
                        dist_type,
                        options.quiet,
                        ref_completeness_vec.as_ref(),
                        query_completeness_vec.as_ref(),
                        options.completeness_cutoff,
                    );
                    write!(output_file, "{distances}")?;
                }
            }
        }
    }
    Ok(())
}

/// Build an inverted index.
pub fn inverted_build(options: InvertedBuildOptions) -> Result<InvertedPaths> {
    check_and_set_threads(options.threads + 1);
    let input_files = get_input_list(&options.file_list, &options.seq_files);
    let mut different_samples = HashSet::new();
    for input_file in input_files.iter() {
        different_samples.insert(input_file.0.clone());
    }

    let (file_order, map_names_labels) = if let Some(species_name_file) = &options.species_names {
        reorder_input_files(&input_files, species_name_file)
    } else {
        default_file_order(&input_files)
    };

    let species_labels_vec = map_names_labels.map(|labels| {
        let mut tmpvec = vec!["".to_string(); different_samples.len()];
        file_order
            .iter()
            .zip(&input_files)
            .for_each(|(idx, (name, _))| {
                tmpvec[*idx] = labels.get(name).unwrap_or(&"".to_owned()).clone();
            });
        tmpvec
    });

    let metadata_vec = if let Some(metadata_file) = &options.metadata {
        let tmpdict = parse_metadata_info(metadata_file);
        let mut tmpvec = vec!["".to_string(); different_samples.len()];
        file_order
            .iter()
            .zip(&input_files)
            .for_each(|(idx, (name, _))| tmpvec[*idx] = tmpdict[name].clone());
        Some(tmpvec)
    } else {
        None
    };

    let skq_file = options.write_skq.then(|| format!("{}.skq", options.output));
    let rc = !options.single_strand;
    let inverted = Inverted::new(
        &input_files,
        skq_file.clone(),
        &file_order,
        options.kmer_length,
        options.sketch_size,
        &HashType::DNA,
        rc,
        options.min_count,
        options.min_qual,
        options.quiet,
        &metadata_vec,
        &species_labels_vec,
    );
    inverted.save(&options.output)?;

    Ok(InvertedPaths {
        ski: format!("{}.ski", options.output),
        skq: skq_file,
    })
}

fn default_file_order(
    input_files: &[(String, Vec<String>)],
) -> (Vec<usize>, Option<HashMap<String, String>>) {
    let tmp_name_set = input_files
        .iter()
        .map(|x| x.0.clone())
        .collect::<HashSet<String>>();
    if tmp_name_set.len() == input_files.len() {
        ((0..input_files.len()).collect(), None)
    } else {
        let mut tmpoutvec = vec![0; input_files.len()];
        let mut tmpmap = HashMap::new();
        for (i, name) in tmp_name_set.iter().enumerate() {
            tmpmap.insert(name.clone(), i);
        }
        for i in 0..tmpoutvec.len() {
            tmpoutvec[i] = tmpmap[&input_files[i].0];
        }
        (tmpoutvec, None)
    }
}

/// Run inverted-index preclustering distances.
pub fn inverted_precluster(options: InvertedPreclusterOptions) -> Result<Option<String>> {
    check_and_set_threads(options.threads);
    let mut output_file = set_ostream(&options.output);
    write_inverted_precluster(&options, &mut output_file)?;
    Ok(options.output)
}

fn write_inverted_precluster(
    options: &InvertedPreclusterOptions,
    output_file: &mut impl Write,
) -> Result<()> {
    let input_prefix = strip_sketch_extension(&options.ski);
    let inverted_index = Inverted::load(input_prefix)?;
    let skq_filename = format!("{input_prefix}.skq");
    let (mmap, bin_stride, kmer_stride, sample_stride) =
        (false, 1, 1, inverted_index.sketch_size());
    let mut skq_reader =
        SketchArrayReader::open(&skq_filename, mmap, bin_stride, kmer_stride, sample_stride);
    let skq_bins = skq_reader.read_all_from_skq(sample_stride * inverted_index.sketch_size());

    let ref_db_name = strip_sketch_extension(&options.skd);
    let mut references = MultiSketch::load_metadata(ref_db_name)?;
    references.read_sketch_data(ref_db_name);
    let n = references.number_samples_loaded();
    let mut knn = options.knn;
    if knn >= n {
        log::warn!("knn={knn} is higher than number of samples={n}");
        knn = n - 1;
    }

    let kmer = inverted_index.kmer();
    let dist_type = set_k(&references, Some(kmer), options.ani)?;
    let ref_completeness_vec = if let Some(file_path) = &options.ref_completeness_file {
        Some(read_completeness_file(file_path, &references)?)
    } else {
        None
    };
    let distances = self_dists_knn_precluster(
        &references,
        &inverted_index,
        &skq_bins,
        skq_reader.sample_stride,
        n,
        knn,
        dist_type,
        options.quiet,
        ref_completeness_vec.as_ref(),
        options.completeness_cutoff,
        &options.retain_unmatched,
    );
    write!(output_file, "{distances}")?;
    Ok(())
}

/// Compute query-vs-reference sparse distances.
pub fn query_dist(
    reference_skm: String,
    query_skm: String,
    output: String,
    kmer_length: usize,
    knn: usize,
    threads: usize,
    ref_completeness_file: Option<String>,
    query_completeness_file: Option<String>,
    completeness_cutoff: f64,
    quiet: bool,
) -> Result<String> {
    dist(DistOptions {
        ref_db: reference_skm,
        query_db: Some(query_skm),
        output: Some(output.clone()),
        knn: Some(knn),
        subset: None,
        kmer: Some(kmer_length),
        ani: true,
        threads,
        ref_completeness_file,
        query_completeness_file,
        completeness_cutoff,
        quiet,
    })?;
    Ok(output)
}

/// Query an inverted index and write the result table.
pub fn inverted_query(
    ski: &str,
    seq_files: &Option<Vec<String>>,
    file_list: &Option<String>,
    output: &Option<String>,
    query_type: &InvertedQueryType,
    min_count: u16,
    min_qual: u8,
    threads: usize,
    quiet: bool,
) -> Result<()> {
    let mut output_file = set_ostream(output);
    let inverted_index = Inverted::load(strip_sketch_extension(ski))?;
    let input_files = get_input_list(file_list, seq_files);
    check_and_set_threads(threads + 1);
    let (queries, query_names) =
        inverted_index.sketch_queries(&input_files, min_count, min_qual, quiet);

    write!(output_file, "Query")?;
    if *query_type == InvertedQueryType::MatchCount {
        for name in inverted_index.sample_names() {
            write!(output_file, "\t{name}")?;
        }
        writeln!(output_file)?;
    } else {
        writeln!(output_file, "\tMatches")?;
    }

    let (tx, rx) = mpsc::channel();
    let progress_bar = get_progress_bar(queries.len(), false, quiet);
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
                    InvertedQueryType::AllBins => (q_name, inverted_index.all_shared_bins(q)),
                    InvertedQueryType::AnyBins => (q_name, inverted_index.any_shared_bins(q)),
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

/// Return database information.
pub fn db_info(db_file: &str, sample_info: bool) -> Result<DbInfo> {
    if db_file.ends_with(".ski") {
        let ski_file = strip_sketch_extension(db_file);
        let index = Inverted::load(ski_file)?;
        let text = if sample_info {
            format!("{index}")
        } else {
            format!("{index:?}")
        };
        Ok(DbInfo {
            text,
            inverted: true,
            sketch_size: index.sketch_size(),
            n_samples: index.n_samples(),
            kmers: vec![index.kmer()],
        })
    } else {
        let ref_db_name = strip_sketch_extension(db_file);
        let sketches = MultiSketch::load_metadata(ref_db_name)?;
        let text = if sample_info {
            format!("{sketches}")
        } else {
            format!("{sketches:?}")
        };
        Ok(DbInfo {
            text,
            inverted: false,
            sketch_size: sketches.sketch_size as usize,
            n_samples: sketches.number_samples_loaded(),
            kmers: sketches.kmer_lengths().to_vec(),
        })
    }
}

impl fmt::Display for DbInfo {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.text)
    }
}
