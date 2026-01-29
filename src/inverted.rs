//! The class to support .ski creation, reading and writing, containing an inverted
//! index of multiple sketches
#[cfg(target_arch = "wasm32")]
use std::fmt;
#[cfg(not(target_arch = "wasm32"))]
use std::sync::mpsc;
#[cfg(not(target_arch = "wasm32"))]
use std::{cmp, fmt};

#[cfg(not(target_arch = "wasm32"))]
extern crate needletail;

#[cfg(target_arch = "wasm32")]
use hashbrown::HashMap;
#[cfg(not(target_arch = "wasm32"))]
use hashbrown::{HashMap, HashSet};
#[cfg(not(target_arch = "wasm32"))]
use indicatif::ParallelProgressIterator;
use indicatif::ProgressIterator;
use rayon::prelude::*;
use roaring::{RoaringBitmap, RoaringTreemap};
use serde::{Deserialize, Serialize};

use super::hashing::{
    bloom_filter::KmerFilter, nthash_iterator::NtHashIterator, HashType, RollHash,
};
use crate::distances::distance_matrix::square_to_condensed;
#[cfg(not(target_arch = "wasm32"))]
use crate::io::InputFastx;
#[cfg(target_arch = "wasm32")]
use crate::sketch::Sketch;
#[cfg(not(target_arch = "wasm32"))]
use crate::sketch::{sketch_datafile::SketchArrayWriter, Sketch};
use crate::utils::get_progress_bar;
use anyhow::Error;
use std::fs::File;
use std::io::{BufReader, BufWriter};

#[cfg(target_arch = "wasm32")]
use crate::logw;
#[cfg(target_arch = "wasm32")]
use wasm_bindgen_file_reader::WebSysFile;

type InvSketches = (Vec<Vec<u16>>, Vec<String>);

/// An inverted index and associated metadata
#[derive(Serialize, Deserialize, Default, Clone, PartialEq)]
pub struct Inverted {
    index: Vec<HashMap<u16, RoaringBitmap>>,
    n_samples: usize,
    sample_names: Vec<String>,
    metadata: Option<Vec<String>>,
    labels: Option<Vec<String>>,
    kmer_size: usize,
    sketch_version: String,
    rc: bool,
    hash_type: HashType,
}

impl Inverted {
    #[cfg(not(target_arch = "wasm32"))]
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
        metadata: &Option<Vec<String>>,
        labels: &Option<Vec<String>>,
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
            metadata: metadata.clone(),
            labels: labels.clone(),
            kmer_size: k,
            sketch_version: env!("CARGO_PKG_VERSION").to_string(),
            rc,
            hash_type: seq_type.clone(),
        }
    }

    #[cfg(not(target_arch = "wasm32"))]
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

    #[cfg(target_arch = "wasm32")]
    /// Sketch query files in a compatible manner with the index,
    /// used when querying the index on the fly
    pub fn sketch_queries(
        &self,
        input_files: (&web_sys::File, Option<&web_sys::File>),
        min_count: u16,
        min_qual: u8,
        quiet: bool,
    ) -> InvSketches {
        let file_order: Vec<usize> = if input_files.1.is_some() {
            vec![0, 1]
        } else {
            vec![0]
        };

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

    #[cfg(not(target_arch = "wasm32"))]
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

    #[cfg(target_arch = "wasm32")]
    /// Loads from a web_sys file object. Designed for WebAssembly
    pub fn load(file: &web_sys::File) -> Result<Self, Error> {
        logw("Loading inverted index", Some("info"));

        let ski_file = BufReader::new(WebSysFile::new(file.clone()));
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

    #[cfg(not(target_arch = "wasm32"))]
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
        let mut multientrysamples: HashSet<String> = HashSet::new();
        let mut differentsamples: HashSet<String> = HashSet::new();

        for i in input_files.iter() {
            if differentsamples.contains(&i.0) {
                multientrysamples.insert(i.0.clone());
            } else {
                differentsamples.insert(i.0.clone());
            }
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
                    .map(|((name, fastxvec), genome_idx)| {
                        let mut hash_its: Vec<Box<dyn RollHash>> = match seq_type {
                            HashType::DNA => NtHashIterator::new(fastxvec, rc, min_qual)
                                .into_iter()
                                .map(|it| Box::new(it) as Box<dyn RollHash>)
                                .collect(),
                            _ => unimplemented!("Inverted index only supported for DNA"),
                        };

                        if let Some(hash_it) = hash_its.first_mut() {
                            if hash_it.seq_len() == 0 {
                                panic!("Genome {genome_idx} has no valid sequence");
                            }

                            let mut read_filter = if hash_it.reads() {
                                let mut filter = KmerFilter::new(min_count);
                                filter.init();
                                Some(filter)
                            } else {
                                None
                            };

                            let signs = Sketch::get_signs_no_densify(
                                &mut **hash_it,
                                k,
                                &mut read_filter,
                                sketch_size,
                            );
                            // if densified {
                            //     log::trace!("{name} was densified");
                            // }
                            (*genome_idx, signs, name)
                        } else {
                            panic!("Empty hash iterator for {name}");
                        }
                    })
                    .for_each_with(tx, |tx, result| {
                        let _ = tx.send(result);
                    });
            });
        });

        let mut sketch_results: Vec<Vec<u16>> =
            vec![Vec::with_capacity(sketch_size as usize); differentsamples.len()];
        let mut indexes: HashSet<usize> = HashSet::with_capacity(multientrysamples.len());
        while let Ok((genome_idx, mut sketch, name)) = rx.recv() {
            if differentsamples.contains(name) {
                // Not yet written!
                if !multientrysamples.contains(name) {
                    // Densifying now!
                    Sketch::densify_bin(&mut sketch);
                } else {
                    // We'll need to densify afterwards, let's save the index
                    indexes.insert(genome_idx);
                }
                sketch_results[genome_idx] = sketch.iter().map(|h| *h as u16).collect();
                differentsamples.remove(name);
            } else {
                // already written! We have to merge
                for bin in 0..sketch_size {
                    let saved_sketch = &mut sketch_results[genome_idx][bin as usize];
                    *saved_sketch = cmp::min(*saved_sketch, sketch[bin as usize] as u16);
                }
            }
        }

        for pos in indexes.iter() {
            // TODO: this is clearly improvable, but I'm going to do it temporarily this way...
            let mut tmpvec: Vec<u64> = sketch_results[*pos].iter().map(|h| *h as u64).collect();
            Sketch::densify_bin(&mut tmpvec[..]);
            sketch_results[*pos] = tmpvec.iter().map(|h| *h as u16).collect();
        }

        // Sample names in the correct order
        // (clones names, but reference would be annoying here)
        let mut sample_names: Vec<String> = vec!["".to_string(); sketch_results.len()];
        file_order
            .iter()
            .zip(input_files)
            .for_each(|(idx, (name, _))| sample_names[*idx] = name.to_string());

        (sketch_results, sample_names)
    }

    #[cfg(target_arch = "wasm32")]
    fn sketch_files_inverted(
        input_files: (&web_sys::File, Option<&web_sys::File>),
        file_order: &[usize],
        k: usize,
        sketch_size: u64,
        seq_type: &HashType,
        rc: bool,
        min_count: u16,
        min_qual: u8,
        _quiet: bool,
    ) -> InvSketches {
        let mut hash_its: Vec<Box<dyn RollHash>> = match seq_type {
            HashType::DNA => NtHashIterator::new(input_files, rc, min_qual)
                .into_iter()
                .map(|it| Box::new(it) as Box<dyn RollHash>)
                .collect(),
            _ => unimplemented!("Inverted index only supported for DNA"),
        };

        if let Some(hash_it) = hash_its.first_mut() {
            if hash_it.seq_len() == 0 {
                panic!("Genome 0 has no valid sequence");
            }

            let mut read_filter = if hash_it.reads() {
                let mut filter = KmerFilter::new(min_count);
                filter.init();
                Some(filter)
            } else {
                None
            };

            let (signs, densified) =
                Sketch::get_signs(&mut **hash_it, k, &mut read_filter, sketch_size);
            if densified {
                logw("The query was densified", Some("trace"));
            }

            let mut sketch_results = vec![Vec::new(); 1];

            sketch_results[0] = signs.iter().map(|h| *h as u16).collect();

            (sketch_results, vec!["".to_string(); file_order.len()])
        } else {
            panic!("Empty hash iterator for the query");
        }
    }

    #[cfg(not(target_arch = "wasm32"))]
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

    /// Get the sample names, with the same order as in the index
    #[cfg(target_arch = "wasm32")]
    pub fn get_sample_names(&self) -> &Vec<String> {
        &self.sample_names
    }

    /// Get the metadata, with the same order as in the index
    #[cfg(target_arch = "wasm32")]
    pub fn get_metadata(&self) -> &Option<Vec<String>> {
        &self.metadata
    }

    /// Get the sample labels, with the same order as in the index
    #[cfg(target_arch = "wasm32")]
    pub fn get_sample_labels(&self) -> &Option<Vec<String>> {
        &self.labels
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
