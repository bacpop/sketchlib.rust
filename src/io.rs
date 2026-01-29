//! Functions to read input fasta/fastq files
use crate::cli::Kmers;
use crate::sketch::multisketch::MultiSketch;

use std::fs::File;
use std::io::{stdout, BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::sync::Mutex;

use anyhow::Context;
use hashbrown::{HashMap, HashSet};
use rayon::prelude::*;
use regex::Regex;

/// Wrapper type for the three fields in an rfile
pub type InputFastx = (String, Vec<String>);

/// Given a list of input files, parses them into triples of name, filename and
/// [`None`] to be used as sketch input.
pub fn read_input_fastas(seq_files: &[String]) -> Vec<InputFastx> {
    let mut input_files = Vec::new();
    // matches the file name (no extension) in a full path
    let re_path =
        Regex::new(r"^.+/(.+\.(fa|fasta|fa\.gz|fasta\.gz|fastq|fastq\.gz|fq|fq\.gz))$").unwrap();
    // matches the file name (no extension) with no path
    let re_name =
        Regex::new(r"^(.+\.(fa|fasta|fa\.gz|fasta\.gz|fastq|fastq\.gz|fq|fq\.gz))$").unwrap();
    for file in seq_files {
        let caps = re_path.captures(file).or(re_name.captures(file));
        let name = match caps {
            Some(capture) => capture[1].to_string(),
            None => file.to_string(),
        };
        input_files.push((name, vec![file.to_string()]));
    }
    input_files
}

/// Give a reordering for input files given some labels, putting the same labels next to each other
pub fn reorder_input_files(
    input_files: &Vec<(String, Vec<String>)>,
    species_name_file: &str,
) -> (Vec<usize>, Option<HashMap<String, String>>) {
    // Set of names, so only these are read from the species order
    log::info!("Reordering samples using labels in {species_name_file}");
    let input_names: HashSet<String> = input_files.iter().map(|fastx| fastx.0.clone()).collect();

    let f = File::open(species_name_file).unwrap_or_else(|_| {
        panic!(
            "Unable to
        open species name file {species_name_file}"
        )
    });
    let f = BufReader::new(f);
    let mut species_labels: HashMap<String, usize> = HashMap::new(); // Stores [species, species order index]
    let mut map_names_labels: HashMap<String, String> = HashMap::new(); // Stores [sample name, species]
    let mut label_order: Vec<(String, usize)> = Vec::with_capacity(input_files.len()); // Stores [species, index-to-be]
    let mut order_idx = 0;
    // Read through labels, assign each name to a cluster in increasing order
    for line in f.lines() {
        let line = line.expect("Unable to read line in species_list");
        let fields: Vec<&str> = line.split_terminator('\t').collect();
        if input_names.contains(fields[0]) {
            if let Some(idx) = species_labels.get(fields[1]) {
                label_order.push((fields[0].to_string(), *idx));
            } else {
                species_labels.insert(fields[1].to_string(), order_idx);
                label_order.push((fields[0].to_string(), order_idx));
                order_idx += 1;
            }
        }
        map_names_labels.insert(fields[0].to_string(), fields[1].to_string());
    }
    log::info!(
        "{} samples with {} unique labels",
        label_order.len(),
        species_labels.len()
    );
    // Order the found labels by cluster they are associated with
    label_order.sort_unstable_by_key(|k| k.1);

    // Create a lookup table for name -> new index
    let mut reordered_dict = HashMap::with_capacity(label_order.len());
    for (new_idx, (reordered_name, _)) in label_order.iter().enumerate() {
        reordered_dict.insert(reordered_name, new_idx);
    }

    if reordered_dict.is_empty() {
        log::warn!("Could not find any sample names in {species_name_file}");
        ((0..input_files.len()).collect(), None)
    } else {
        // Use lookup table to create a list of new index for each input sample
        // This deals with missing labels
        let mut sample_order = Vec::with_capacity(input_files.len());
        let mut new_idx = reordered_dict.len() - 1;
        for sample_name in input_files {
            let sample_idx = if let Some(order) = reordered_dict.get(&sample_name.0) {
                *order
            } else {
                new_idx += 1;
                new_idx
            };
            sample_order.push(sample_idx);
        }
        log::info!(
            "Found {} of {} input samples with given labels",
            input_files
                .iter()
                .filter(|name| reordered_dict.contains_key(&name.0))
                .count(),
            input_files.len()
        );
        (sample_order, Some(map_names_labels))
    }
}

/// Parse the metadata for creating an inverted index, so that we can include them in it
pub fn parse_metadata_info(metadata_file: &str) -> HashMap<String, String> {
    let f = File::open(metadata_file)
        .unwrap_or_else(|_| panic!("Unable to open species name file {metadata_file}"));
    let f = BufReader::new(f);
    let mut out_dict: HashMap<String, String> = HashMap::new(); // Stores [label, metadata]
                                                                // Read through labels, saves all metadata in the dictionary
    for line in f.lines() {
        let line = line.expect("Unable to read line in metadata");
        let fields: Vec<&str> = line.split_terminator("\t").collect();
        if out_dict.contains_key(fields[0]) {
            panic!("Some entry in metadata is duplicated");
        } else {
            out_dict.insert(fields[0].to_string(), fields[1].to_string());
        }
    }

    log::info!("Got metadata for {} labels", out_dict.len(),);

    out_dict
}

/// Validate and sort k-mer lists provided via the CLI
pub fn parse_kmers(k: &Kmers) -> Vec<usize> {
    if k.k_vals.is_some() && k.k_seq.is_some() {
        panic!("Only one of --k-vals or --k-seq should be specified");
    }

    let mut kmers = if let Some(k) = &k.k_vals {
        k.clone().to_vec()
    } else if let Some(k) = &k.k_seq {
        (k[0]..=k[1]).step_by(k[2]).collect()
    } else {
        panic!("Must specify --k-vals or --k-seq");
    };

    kmers.sort_unstable();
    if !kmers.iter().all(|&k| k >= 3) {
        panic!("K-mers must be >=3");
    }

    kmers
}

/// Set a buffered stream to write to.
///
/// Either a file (if [`Some`]) or stdout otherwise (if [`None`]).
pub fn set_ostream(oprefix: &Option<String>) -> BufWriter<Box<dyn Write>> {
    let out_writer = match oprefix {
        Some(prefix) => {
            let path = Path::new(prefix);
            Box::new(File::create(path).unwrap()) as Box<dyn Write>
        }
        None => Box::new(stdout()) as Box<dyn Write>,
    };
    BufWriter::new(out_writer)
}

/// Obtain a list of input files and names from command line input.
///
/// If `file_list` is provided, read each line as `name\tseq1\tseq2`, where
/// `seq2` is optional, and if present the reverse fastqs. Otherwise, treat
/// as fasta.
///
/// If `seq_files` are provided use [`read_input_fastas`].
pub fn get_input_list(
    file_list: &Option<String>,
    seq_files: &Option<Vec<String>>,
) -> Vec<InputFastx> {
    if file_list.is_none() && seq_files.is_none() {
        panic!("No input files provided");
    }
    // Read input
    match file_list {
        Some(files) => {
            let mut input_files: Vec<InputFastx> = Vec::new();
            let f = File::open(files).unwrap_or_else(|_| {
                panic!(
                    "
            Unable to open file_list {files}"
                )
            });
            let f = BufReader::new(f);
            for line in f.lines() {
                let line = line.expect("Unable to read line in file_list");
                let fields: Vec<&str> = line.split_whitespace().collect();
                // 1 entry:           one fastX with name = file
                // 2 entries:         name, fastX
                // 3 or more entries: name, fastX1, fastX2, ...
                let parsed_input = match fields.len() {
                    0 => panic!("Unable to parse line in file_list"),
                    1 => ((fields[0].to_string()), vec![fields[0].to_string()]),
                    2 => ((fields[0].to_string()), vec![fields[1].to_string()]),
                    _ => {
                        let mut tmpvec = Vec::with_capacity(fields.len() - 1);
                        for file in fields.iter().skip(1) {
                            tmpvec.push(file.to_string());
                        }
                        ((fields[0].to_string()), tmpvec)
                    }
                };
                input_files.push(parsed_input);
            }
            input_files
        }
        None => read_input_fastas(seq_files.as_ref().unwrap()),
    }
}

/// Read the sample names provided in a `--subset`` file
pub fn read_subset_names(subset_file: &str) -> Vec<String> {
    let f = File::open(subset_file).unwrap_or_else(|_| panic!("Unable to open {subset_file}"));
    let f = BufReader::new(f);
    let mut subset_names = Vec::new();
    for line in f.lines() {
        subset_names.push(line.unwrap());
    }
    log::info!("Read {} names to subset", subset_names.len());
    subset_names
}

/// Read completeness values from a file and create a vector for each genome in the provided sketches.
/// If a genome is not found in the file, a default value of 1.0 is used.
pub fn read_completeness_file(
    completeness_file: &str,
    sketches: &MultiSketch,
) -> Result<Vec<f64>, crate::Error> {
    // Pre-allocate vector with default values (1.0 for missing genomes)
    let mut completeness_vec = vec![1.0_f64; sketches.number_samples_loaded()];
    let missing_genomes = Mutex::new(Vec::new());

    // Open file with buffered reader for streaming
    let f = File::open(completeness_file)
        .with_context(|| format!("Failed to open completeness file: {completeness_file}"))?;
    let f = BufReader::new(f);

    // Read lines and collect completeness values
    let lines: Vec<String> = f.lines().collect::<Result<Vec<_>, _>>().with_context(|| {
        format!("Failed to read lines from completeness file: {completeness_file}")
    })?;

    // Parse all lines in parallel and collect completeness values
    let updates: Vec<(usize, f64)> = lines
        .par_iter()
        .filter_map(|line| {
            if let Some((genome_id, completeness_str)) = line.split_once('\t') {
                if let Ok(completeness) = completeness_str.parse::<f64>() {
                    // Use MultiSketch's method to get the logical index
                    if let Some(index) = sketches.get_sample_index(genome_id) {
                        Some((index, completeness))
                    } else {
                        // Track missing genomes
                        missing_genomes.lock().unwrap().push(genome_id.to_string());
                        None
                    }
                } else {
                    None
                }
            } else {
                None
            }
        })
        .collect();

    // Apply updates to completeness_vec (thread safe)
    for (index, completeness) in updates {
        completeness_vec[index] = completeness;
    }

    // Report missing genomes
    let missing = missing_genomes.into_inner().unwrap();
    if !missing.is_empty() {
        log::warn!(
            "Found {} genomes not in completeness file, using default 1.0: {}",
            missing.len(),
            missing.join(", ")
        );
    }

    Ok(completeness_vec)
}

#[cfg(test)]
use pretty_assertions::assert_eq;
#[cfg(test)]
use tempfile::NamedTempFile;

#[test]
fn test_reorder_input_files() {
    // Create a temporary file for species labels
    let mut temp_file = NamedTempFile::new().expect("Failed to create temp file");
    writeln!(
        temp_file,
        "sample1\tspecies A\nsample2\tspeciesB\nsample3\tspecies A\nsample4\tspeciesC\nsample6\tspeciesD"
    )
    .expect("Failed to write to temp file");

    let input_files = vec![
        ("sample1".to_string(), vec!["assembly1.fa".to_string()]),
        ("sample2".to_string(), vec!["assembly2.fa".to_string()]),
        ("sample3".to_string(), vec!["assembly3.fa".to_string()]),
        ("sample4".to_string(), vec!["assembly4.fa".to_string()]),
        ("sample5".to_string(), vec!["assembly5.fa".to_string()]),
    ];

    let species_name_file = temp_file.path().to_str().unwrap();
    let reordered_indices = reorder_input_files(&input_files, species_name_file);

    assert_eq!(reordered_indices.0.len(), input_files.len());
    assert_eq!(reordered_indices.0, vec![0, 2, 1, 3, 4]) // 1(A), 3(A), 2(B), 4(C), 5(NA) (sample6 not included)
}
