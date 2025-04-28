//! Functions to read input fasta/fastq files
use crate::cli::Kmers;

use std::fs::File;
use std::io::{stdout, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use hashbrown::{HashMap, HashSet};
use regex::Regex;

pub type InputFastx = (String, String, Option<String>);

pub fn read_input_fastas(seq_files: &[String]) -> Vec<InputFastx> {
    let mut input_files = Vec::new();
    // matches the file name (no extension) in a full path
    let re_path = Regex::new(r"^.+/(.+\.?i:fa|fasta|fastq|fastq\.gz)$").unwrap();
    // matches the file name (no extension) with no path
    let re_name = Regex::new(r"^(.+\.?i:fa|fasta|fastq|fastq\.gz)$").unwrap();
    for file in seq_files {
        let caps = re_path.captures(file).or(re_name.captures(file));
        let name = match caps {
            Some(capture) => capture[1].to_string(),
            None => file.to_string(),
        };
        input_files.push((name, file.to_string(), None));
    }
    input_files
}

/// Give a reordering for input files given some labels, putting the same labels next to each other
pub fn reorder_input_files(
    input_files: &Vec<(String, String, Option<String>)>,
    species_name_file: &str,
) -> Vec<usize> {
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
    let mut label_order: Vec<(String, usize)> = Vec::with_capacity(input_files.len());
    let mut order_idx = 0;
    // Read through labels, assign each name to a cluster in increasing order
    for line in f.lines() {
        let line = line.expect("Unable to read line in species_list");
        let fields: Vec<&str> = line.split_terminator("\t").collect();
        if input_names.contains(fields[0]) {
            if let Some(idx) = species_labels.get(fields[1]) {
                label_order.push((fields[0].to_string(), *idx));
            } else {
                species_labels.insert(fields[1].to_string(), order_idx);
                label_order.push((fields[0].to_string(), order_idx));
                order_idx += 1;
            }
        }
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

    let mut sample_order = Vec::new();
    if reordered_dict.is_empty() {
        log::warn!("Could not find any sample names in {species_name_file}");
        sample_order = (0..input_files.len()).collect();
    } else {
        // Use lookup table to create a list of new index for each input sample
        // This deals with missing labels
        sample_order.reserve_exact(input_files.len());
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
                .len()
                .saturating_sub(new_idx - (reordered_dict.len() - 1)),
            input_files.len()
        );
    }
    sample_order
}

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
                // 1 entry: fasta with name = file
                // 2 entries: fasta name, file
                // 3 entries: fastq name, file1, file2
                let parsed_input = match fields.len() {
                    1 => ((fields[0].to_string()), fields[0].to_string(), None),
                    2 => ((fields[0].to_string()), fields[1].to_string(), None),
                    3 => (
                        (fields[0].to_string()),
                        fields[1].to_string(),
                        Some(fields[2].to_string()),
                    ),
                    _ => {
                        panic!("Unable to parse line in file_list")
                    }
                };
                input_files.push(parsed_input);
            }
            input_files
        }
        None => read_input_fastas(seq_files.as_ref().unwrap()),
    }
}

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
        ("sample1".to_string(), "assembly1.fa".to_string(), None),
        ("sample2".to_string(), "assembly2.fa".to_string(), None),
        ("sample3".to_string(), "assembly3.fa".to_string(), None),
        ("sample4".to_string(), "assembly4.fa".to_string(), None),
        ("sample5".to_string(), "assembly5.fa".to_string(), None),
    ];

    let species_name_file = temp_file.path().to_str().unwrap();
    let reordered_indices = reorder_input_files(&input_files, species_name_file);

    assert_eq!(reordered_indices.len(), input_files.len());
    assert_eq!(reordered_indices, vec![0, 2, 1, 3, 4]) // 1(A), 3(A), 2(B), 4(C), 5(NA) (sample6 not included)
}
