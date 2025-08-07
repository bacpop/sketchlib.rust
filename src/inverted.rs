//! The class to support .ski creation, reading and writing, containing an inverted
//! index of multiple sketches
use std::fmt;
use std::sync::mpsc;

extern crate needletail;
use hashbrown::HashMap;
use indicatif::{ParallelProgressIterator, ProgressIterator};
use rayon::prelude::*;
use roaring::{RoaringBitmap, RoaringTreemap};
use serde::{Deserialize, Serialize};

use super::hashing::{
    bloom_filter::KmerFilter, nthash_iterator::NtHashIterator, HashType, RollHash,
};
use crate::distances::distance_matrix::square_to_condensed;
use crate::io::InputFastx;
use crate::sketch::{sketch_datafile::SketchArrayWriter, Sketch};
use crate::utils::get_progress_bar;
use anyhow::Error;
use std::fs::File;
use std::io::{BufReader, BufWriter};

type InvSketches = (Vec<Vec<u16>>, Vec<String>);

/// An inverted index and associated metadata
#[derive(Serialize, Deserialize, Default, Clone, PartialEq)]
pub struct Inverted {
    index: Vec<HashMap<u16, RoaringBitmap>>,
    n_samples: usize,
    sample_names: Vec<String>,
    kmer_size: usize,
    sketch_version: String,
    rc: bool,
    hash_type: HashType,
}

impl Inverted {
    /// Sketch files without transposing bins, and invert the index. Files
    /// may be reordered via [`crate::io::reorder_input_files`]
    ///
    /// Optionally write untransposed samples to an skq file with `write_skq`.
    pub fn new(
        input_files: &[InputFastx],
        write_skq: Option<String>,
        file_order: &[usize],
        k: usize,
        sketch_size: u64,
        seq_type: &HashType,
        rc: bool,
        min_count: u16,
        min_qual: u8,
        quiet: bool,
    ) -> Self {
        log::info!("Creating sketches");
        let (sketches, names) = Self::sketch_files_inverted(
            input_files,
            file_order,
            k,
            sketch_size,
            seq_type,
            rc,
            min_count,
            min_qual,
            quiet,
        );
        if let Some(skq_file) = write_skq {
            log::info!("Writing bins for use with precluster as {skq_file}");
            let (bin_stride, kmer_stride, sample_stride) = (1, 1, sketch_size as usize);
            let mut skq_writer =
                SketchArrayWriter::new(&skq_file, bin_stride, kmer_stride, sample_stride);
            for sketch in &sketches {
                skq_writer.write_sketch(sketch);
            }
        }
        log::info!("Inverting sketch order");
        Self {
            index: Self::build_inverted_index(&sketches, sketch_size),
            n_samples: names.len(),
            sample_names: names,
            kmer_size: k,
            sketch_version: env!("CARGO_PKG_VERSION").to_string(),
            rc,
            hash_type: seq_type.clone(),
        }
    }

    /// Sketch query files in a compatible manner with the index,
    /// used when querying the index on the fly
    pub fn sketch_queries(
        &self,
        input_files: &[InputFastx],
        min_count: u16,
        min_qual: u8,
        quiet: bool,
    ) -> InvSketches {
        let file_order: Vec<usize> = (0..input_files.len()).collect();
        Self::sketch_files_inverted(
            input_files,
            &file_order,
            self.kmer_size,
            self.index.len() as u64,
            &self.hash_type,
            self.rc,
            min_count,
            min_qual,
            quiet,
        )
    }

    /// Sample names in the index
    pub fn sample_names(&self) -> &Vec<String> {
        &self.sample_names
    }

    /// Number of samples in the index
    pub fn n_samples(&self) -> usize {
        self.sample_names.len()
    }

    /// Sample name at the given index
    pub fn sample_at(&self, idx: usize) -> &str {
        &self.sample_names[idx]
    }

    /// K-mer length used for the index
    pub fn kmer(&self) -> usize {
        self.kmer_size
    }

    /// Sketch size used for the index
    pub fn sketch_size(&self) -> usize {
        self.index.len()
    }

    /// Saves to `file_prefix.ski`, using MessagePack as the serialisation format
    pub fn save(&self, file_prefix: &str) -> Result<(), Error> {
        let filename = format!("{file_prefix}.ski");
        log::info!("Saving inverted index to {filename}");
        let serial_file = BufWriter::new(File::create(filename)?);
        let mut compress_writer = snap::write::FrameEncoder::new(serial_file);
        rmp_serde::encode::write(&mut compress_writer, self)?;
        Ok(())
    }

    /// Loads from `file_prefix.ski`, using MessagePack as the serialisation format
    // NB MessagePack rather the CBOR uses here because of
    // https://github.com/enarx/ciborium/issues/96
    pub fn load(file_prefix: &str) -> Result<Self, Error> {
        let filename = format!("{file_prefix}.ski");
        log::info!("Loading inverted index from {filename}");
        let ski_file = BufReader::new(File::open(filename)?);
        let decompress_reader = snap::read::FrameDecoder::new(ski_file);
        let ski_obj: Self = rmp_serde::decode::from_read(decompress_reader)?;
        Ok(ski_obj)
    }

    /// Query a single sample against the index, returning the number of matches
    /// of bins against all samples in the index
    pub fn query_against_inverted_index(&self, query_sigs: &[u16]) -> Vec<u32> {
        let mut match_counts = vec![0; self.sample_names.len()];

        for (bin_idx, query_bin_hash) in query_sigs.iter().enumerate() {
            if let Some(matching_samples) = self.index[bin_idx].get(query_bin_hash) {
                for sample_idx in matching_samples {
                    match_counts[sample_idx as usize] += 1;
                }
            }
        }
        match_counts
    }

    pub fn inverted_index_ref(&self) -> &Vec<HashMap<u16, RoaringBitmap>> {
        &self.index
    }

    /// Return indexes of samples where all sketch bins are the same as a query
    pub fn all_shared_bins(&self, query_sigs: &[u16]) -> Vec<u32> {
        let mut matching_bits = RoaringBitmap::new();
        matching_bits.insert_range(0..self.sample_names.len() as u32);
        for (bin_idx, query_bin_hash) in query_sigs.iter().enumerate() {
            if let Some(matching_samples) = self.index[bin_idx].get(query_bin_hash) {
                matching_bits &= matching_samples;
            } else {
                matching_bits.clear();
                break;
            }
        }

        matching_bits.iter().collect()
    }

    /// Return indexes of samples where at least one sketch bins is the same as a query
    pub fn any_shared_bins(&self, query_sigs: &[u16]) -> Vec<u32> {
        let mut matching_bits = RoaringBitmap::new();
        for (bin_idx, query_bin_hash) in query_sigs.iter().enumerate() {
            if let Some(matching_samples) = self.index[bin_idx].get(query_bin_hash) {
                matching_bits |= matching_samples;
            }
        }

        matching_bits.iter().collect()
    }

    /// Get the pairwise comparison list used by prefiler
    pub fn any_shared_bin_list(&self, quiet: bool) -> RoaringTreemap {
        let percent = false;
        let progress_bar = get_progress_bar(self.index.len(), percent, quiet);
        self.index
            .iter()
            .progress_with(progress_bar)
            .map(|bin| {
                bin.par_values()
                    .map(|hash_pres| {
                        let mut pair_map_hash = RoaringTreemap::new();
                        let samples_together: Vec<u32> = hash_pres.iter().collect();
                        for (i, sample1_idx) in samples_together.iter().enumerate() {
                            for sample2_idx in samples_together.iter().skip(i + 1) {
                                pair_map_hash.insert(square_to_condensed(
                                    *sample1_idx as usize,
                                    *sample2_idx as usize,
                                    self.n_samples,
                                ) as u64);
                            }
                        }
                        pair_map_hash
                    })
                    .reduce(RoaringTreemap::new, |pair_map_hash_a, pair_map_hash_b| {
                        pair_map_hash_a | pair_map_hash_b
                    }) // Reduction from pair_map_hash produces pair_map_bin
            })
            .reduce(|pair_map_all, pair_map_bin| pair_map_all | pair_map_bin)
            .unwrap()
        // Reduction from pair_map_bin produces pair_map_all
    }

    fn sketch_files_inverted(
        input_files: &[InputFastx],
        file_order: &[usize],
        k: usize,
        sketch_size: u64,
        seq_type: &HashType,
        rc: bool,
        min_count: u16,
        min_qual: u8,
        quiet: bool,
    ) -> InvSketches {
        if *seq_type != HashType::DNA {
            unimplemented!("Inverted index only supported for DNA");
        }

        let (tx, rx) = mpsc::channel();
        let percent = false;
        let progress_bar = get_progress_bar(input_files.len(), percent, quiet);
        rayon::scope(|s| {
            s.spawn(move |_| {
                input_files
                    .par_iter()
                    .zip(file_order)
                    .progress_with(progress_bar)
                    .map(|((name, fastx1, fastx2), genome_idx)| {
                        let mut hash_it = NtHashIterator::new((fastx1, fastx2.as_ref()), rc, min_qual);
                        let hash_it = hash_it.first_mut().expect("Empty hash iterator");
                        if hash_it.seq_len() == 0 {
                            panic!("{name} has no valid sequence");
                        }

                        let mut read_filter = if hash_it.reads() {
                            let mut filter = KmerFilter::new(min_count);
                            filter.init();
                            Some(filter)
                        } else {
                            None
                        };

                        let (signs, densified) =
                            Sketch::get_signs(hash_it, k, &mut read_filter, sketch_size);
                        if densified {
                            log::trace!("{name} was densified");
                        }
                        (*genome_idx, signs)
                    })
                    .for_each_with(tx, |tx, result| {
                        let _ = tx.send(result);
                    });
            });
        });

        let mut sketch_results = vec![Vec::new(); input_files.len()];
        while let Ok((genome_idx, sketch)) = rx.recv() {
            let sketch_u16 = sketch.iter().map(|h| *h as u16).collect();
            sketch_results[genome_idx] = sketch_u16;
        }

        // Sample names in the correct order
        // (clones names, but reference would be annoying here)
        let mut sample_names: Vec<String> = vec!["".to_string(); input_files.len()];
        file_order
            .iter()
            .zip(input_files)
            .for_each(|(idx, (name, _, _))| sample_names[*idx] = name.to_string());

        (sketch_results, sample_names)
    }

    fn build_inverted_index(
        genome_sketches: &[Vec<u16>],
        sketch_size: u64,
    ) -> Vec<HashMap<u16, RoaringBitmap>> {
        // initialize inverted index structure
        let mut inverted_index: Vec<HashMap<u16, RoaringBitmap>> =
            vec![HashMap::new(); sketch_size as usize];

        // process each sketch
        // this could be parallelised over sketch bin, but probably not worth it
        // (NB doing par_iter on the below doesn't quite work as multiple mut borrows
        // of inverted_index needed)
        for (genome_idx, genome_signs) in genome_sketches.iter().enumerate() {
            for (i, hash) in genome_signs.iter().enumerate() {
                // add current genome to the inverted index at the current position
                inverted_index[i]
                    .entry(*hash)
                    .and_modify(|genome_list| {
                        genome_list.insert(genome_idx as u32);
                    })
                    .or_insert_with(|| {
                        let mut rb = RoaringBitmap::new();
                        rb.insert(genome_idx as u32);
                        rb
                    });
            }
        }
        inverted_index
    }
}

impl fmt::Debug for Inverted {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "sketch_version={}\nsequence_type={:?}\nsketch_size={}\nn_samples={}\nkmer={}\nrc={}\ninverted=true\n",
            self.sketch_version,
            self.hash_type,
            self.index.len(),
            self.sample_names.len(),
            self.kmer_size,
            self.rc,
        )?;
        let mut sizes = Vec::new();
        for bin in self.index.iter() {
            sizes.push(bin.len());
        }
        write!(
            f,
            "max_hashes_per_bin={}\nmin_hashes_per_bin={}\navg_hashes_per_bin={}",
            sizes.iter().max().unwrap(),
            sizes.iter().min().unwrap(),
            sizes.iter().sum::<usize>() as f64 / sizes.len() as f64
        )
    }
}

impl fmt::Display for Inverted {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "Name")?;
        for sketch in &self.sample_names {
            writeln!(f, "{sketch}")?;
        }
        Ok(())
    }
}
