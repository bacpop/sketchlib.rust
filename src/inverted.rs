use std::sync::mpsc;

extern crate needletail;
use hashbrown::HashMap;
use indicatif::ParallelProgressIterator;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use super::hashing::{nthash_iterator::NtHashIterator, HashType, RollHash};
use crate::hashing::aahash_iterator::AaHashIterator;
use crate::io::InputFastx;
use crate::sketch::*;
use crate::utils::get_progress_bar;
use crate::bloom_filter::KmerFilter;
use anyhow::Error;
use std::fs::File;
use std::io::{BufReader, BufWriter};

#[derive(Serialize, Deserialize, Debug, Default, Clone, PartialEq)]
pub struct Inverted {
    index: Vec<HashMap<u64, Vec<u32>>>,
    sample_names: Vec<String>,
}

impl Inverted {
    pub fn new(
        input_files: &[InputFastx],
        k: usize,
        sketch_size: u64,
        seq_type: &HashType,
        rc: bool,
        min_count: u16,
        min_qual: u8,
        quiet: bool,
    ) -> Self {
        log::info!("Creating sketches");
        let (sketches, sample_names) = Self::sketch_files_inverted(
            input_files,
            k,
            sketch_size,
            seq_type,
            rc,
            min_count,
            min_qual,
            quiet,
        );
        log::info!("Inverting sketch order");
        Self {
            index: Self::build_inverted_index(&sketches, sketch_size),
            sample_names,
        }
    }

    pub fn save(&self, file_prefix: &str) -> Result<(), Error> {
        let filename = format!("{}.ski", file_prefix);
        log::info!("Saving inverted index to {filename}");
        let serial_file = BufWriter::new(File::create(filename)?);
        let mut compress_writer = snap::write::FrameEncoder::new(serial_file);
        ciborium::ser::into_writer(self, &mut compress_writer)?;
        Ok(())
    }

    pub fn load(file_prefix: &str) -> Result<Self, Error> {
        let filename = format!("{}.ski", file_prefix);
        log::info!("Loading inverted index from {filename}");
        let ski_file = BufReader::new(File::open(filename)?);
        let decompress_reader = snap::read::FrameDecoder::new(ski_file);
        let ski_obj: Self = ciborium::de::from_reader(decompress_reader)?;
        Ok(ski_obj)
    }

    pub fn query_against_inverted_index(&self, query_sigs: &[u64], n_samples: usize) -> Vec<u32> {
        let mut match_counts = vec![0; n_samples];

        for (bin_idx, query_bin_hash) in query_sigs.iter().enumerate() {
            if let Some(matching_samples) = self.index[bin_idx].get(query_bin_hash) {
                for &sample_idx in matching_samples {
                    match_counts[sample_idx as usize] += 1;
                }
            }
        }
        match_counts
    }
    // // example for how query against II might work
    // // iterate over bins
    // for query_bin, inverted_index_bins in query_sigs.zip(inverted_index) { // iterating over u64, HashMap<u64, Vec<String>>
    //     // look up bin value in the hash map
    //     if query in inverted_index_bins {
    //         let same_bin_samples = inverted_index_bins[query_bin]; // sample_bin_samples: Vec<String>
    //         for sample in same_bin_samples { // for each sample String
    //             dist_vec[sample] += 1;
    //         }
    //         // add one to all samples distances in the vec
    //         for
    //     }
    // }

    fn sketch_files_inverted(
        input_files: &[InputFastx],
        k: usize,
        sketch_size: u64,
        seq_type: &HashType,
        rc: bool,
        min_count: u16,
        min_qual: u8,
        quiet: bool,
    ) -> (Vec<Vec<u64>>, Vec<String>) {
        let (tx, rx) = mpsc::channel();

        let percent = false;
        let progress_bar = get_progress_bar(input_files.len(), percent, quiet);
        rayon::scope(|s| {
            s.spawn(move |_| {
                input_files
                    .par_iter()
                    .enumerate()
                    .progress_with(progress_bar)
                    .map(|(genome_idx, (name, fastx1, fastx2))| {
                        let mut hash_its: Vec<Box<dyn RollHash>> = match seq_type {
                            HashType::DNA => {
                                NtHashIterator::new((fastx1, fastx2.as_ref()), rc, min_qual)
                                    .into_iter()
                                    .map(|it| Box::new(it) as Box<dyn RollHash>)
                                    .collect()
                            }
                            HashType::AA(level) => {
                                AaHashIterator::new(fastx1, level.clone(), false)
                                    .into_iter()
                                    .map(|it| Box::new(it) as Box<dyn RollHash>)
                                    .collect()
                            }
                            _ => unimplemented!("Inverted index for PDB files not supported"),
                        };

                        if let Some(hash_it) = hash_its.first_mut() {
                            if hash_it.seq_len() == 0 {
                                panic!("Genome {} has no valid sequence", genome_idx);
                            }

                            let mut read_filter = if hash_it.reads() {
                                let mut filter = KmerFilter::new(min_count);
                                filter.init();
                                Some(filter)
                            } else {
                                None
                            };

                            let (signs, densified) =  Sketch::get_signs(&mut **hash_it, k, &mut read_filter, sketch_size);
                            if densified {
                                log::trace!("{name} was densified");
                            }
                            (genome_idx, signs)
                        } else {
                            panic!("Empty hash iterator for {name}");
                        }
                    })
                    .for_each_with(tx, |tx, result| {
                        let _ = tx.send(result);
                    });
            });
        });

        let mut sketch_results = Vec::new();
        while let Ok((genome_idx, sketch)) = rx.recv() {
            while sketch_results.len() <= genome_idx {
                sketch_results.push(Vec::new());
            }
            sketch_results[genome_idx] = sketch;
        }

        // TODO: this feels unnecessary -- can't we just extract from input files when needed?
        let sample_names: Vec<String> = input_files
            .iter()
            .map(|(name, _, _)| name.clone())
            .collect();

        (sketch_results, sample_names)
    }

    fn build_inverted_index(
        genome_sketches: &Vec<Vec<u64>>,
        sketch_size: u64,
    ) -> Vec<HashMap<u64, Vec<u32>>> {
        // initialize inverted index structure
        let mut inverted_index: Vec<HashMap<u64, Vec<u32>>> =
            vec![HashMap::new(); sketch_size as usize];

        // process each sketch
        for (genome_idx, genome_signs) in genome_sketches.into_iter().enumerate() {
            for (i, hash) in genome_signs.iter().enumerate() {
                // add current genome to the inverted index at the current position
                inverted_index[i]
                    .entry(*hash)
                    .and_modify(|genome_list| genome_list.push(genome_idx as u32))
                    .or_insert(vec![genome_idx as u32]);
            }
        }
        inverted_index
    }
}
