//! Methods to sketch samples, save/load sketches
use std::cmp::Ordering;
use std::fmt;
#[cfg(not(target_arch = "wasm32"))]
use std::sync::mpsc;

#[cfg(not(target_arch = "wasm32"))]
use indicatif::ParallelProgressIterator;
#[cfg(not(target_arch = "wasm32"))]
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use super::hashing::{bloom_filter::KmerFilter, RollHash};

#[cfg(not(target_arch = "wasm32"))]
use super::hashing::{HashType, simdsketch_wrapper::sketch_with_simd};

#[cfg(not(target_arch = "wasm32"))]
use crate::hashing::aahash_iterator::AaHashIterator;
#[cfg(not(target_arch = "wasm32"))]
use crate::io::InputFastx;
#[cfg(feature = "3di")]
use crate::structures::pdb_to_3di;
#[cfg(not(target_arch = "wasm32"))]
use crate::utils::get_progress_bar;

pub mod multisketch;

pub mod sketch_datafile;
#[cfg(not(target_arch = "wasm32"))]
use self::sketch_datafile::SketchArrayWriter;

#[cfg(not(target_arch = "wasm32"))]
use simd_sketch::{SketchParams, Sketch as Sketch_simd, SketchAlg, BitSketch};

/// Total width of all bins (used as sign % sign_mod)
pub const SIGN_MOD: u64 = (1 << 61) - 1;


/// A single sample's sketch
#[derive(Serialize, Deserialize, Debug, Default, Clone, PartialEq)]
pub struct Sketch {
    #[serde(skip)]
    usigs: Vec<u16>,
    name: String,
    index: Option<usize>,
    rc: bool,
    reads: bool,
    seq_length: usize,
    densified: bool,
    acgt: [usize; 4],
    non_acgt: usize,
}

// TODO: should this take hash_it and filter as input?
impl Sketch {
    /// Sketch a sample from a hash generator over its k-mers, transposing
    pub fn new<H: RollHash + ?Sized>(
        seq_hashes: &mut H,
        name: &str,
        kmer_lengths: &[usize],
        sketch_size: u64,
        rc: bool,
        min_count: u16,
    ) -> Self {
        let mut usigs: Vec<u16> = Vec::with_capacity(sketch_size as usize * kmer_lengths.len());

        let mut read_filter = if seq_hashes.reads() {
            let mut filter = KmerFilter::new(min_count);
            filter.init();
            Some(filter)
        } else {
            None
        };

        // Build the sketches across k-mer lengths
        let mut minhash_sum = 0.0;
        let mut densified = false;
        for k in kmer_lengths {
            log::debug!("Running sketching at k={k}");
            let (signs, k_densified) = Self::get_signs(seq_hashes, *k, &mut read_filter, sketch_size);
            densified |= k_densified;
            minhash_sum += (signs[0] as f64) / (SIGN_MOD as f64);

            for h in signs {
                usigs.push(h as u16);
            }
        }
        let (reads, acgt, non_acgt) = seq_hashes.sketch_data();

        // Estimate of sequence length from read data
        let seq_length = if reads {
            ((kmer_lengths.len() as f64) / minhash_sum) as usize
        } else {
            seq_hashes.seq_len()
        };

        Self {
            usigs,
            name: name.to_string(),
            index: None,
            rc,
            reads,
            seq_length,
            densified,
            acgt,
            non_acgt,
        }
    }

    /// Create a Sketch object from a simd_sketch Sketch object
    #[cfg(not(target_arch = "wasm32"))]
    pub fn from_sketch_simd(
        sketches: &mut Vec<Sketch_simd>,
        name: &str,
        reads: bool,
    ) -> Self {
        let testinputsketchs;
        let rc;
        match &sketches[0] {
            Sketch_simd::BottomSketch(_sketch) => panic!("We only support BucketSketch."),
            Sketch_simd::BucketSketch(sketch) => {
                testinputsketchs = &sketch.buckets;
                rc = sketch.rc;
            },
        }

        let testwithwhichtoiter;
        match testinputsketchs {
            BitSketch::B16(thevec) => testwithwhichtoiter = thevec,
            _ => panic!("Only supporting 16 bits as bin size.")
        }
        let sketch_size = testwithwhichtoiter.len() as usize;

        let mut usigs = Vec::with_capacity(sketch_size * sketches.len());
        let mut densified = false;

        for is in sketches {
            let inputsketchs;
            let empty;
            match is {
                Sketch_simd::BottomSketch(_sketch) => panic!("We only support BucketSketch."),
                Sketch_simd::BucketSketch(sketch) => {
                    inputsketchs = &mut sketch.buckets;
                    empty = &sketch.empty;
                }
            }

            match inputsketchs {
                BitSketch::B16(thevec) => {
                    densified |= Self::densify_bin_u16(thevec, empty);
                    usigs.append(thevec);
                }
                _ => panic!("Only supporting 16 bits as bin size.")
            }

            // TODO: get the other stats/values listed below
        }

        Self {
            usigs,
            name: name.to_string(),
            index: None,
            rc,
            reads : reads,
            seq_length : 0,                 // TODO
            densified,
            acgt : [0,0,0,0],               // TODO
            non_acgt : 0,                   // TODO
        }
    }

    /// Get the sketch bins for a sample, but do not transpose
    pub fn get_signs<H: RollHash + ?Sized>(
        seq_hashes: &mut H,
        kmer_size: usize,
        filter: &mut Option<KmerFilter>,
        num_bins: u64,
    ) -> (Vec<u64>, bool) {
        // Setup storage for each k
        let mut signs = vec![u64::MAX; num_bins as usize];
        if let Some(read_filter) = filter {
            read_filter.clear();
        }
        seq_hashes.set_k(kmer_size);

        // Calculate bin minima across all sequence
        let bin_size: u64 = SIGN_MOD.div_ceil(num_bins);
        for hash in seq_hashes.iter() {
            Self::bin_sign(&mut signs, hash % SIGN_MOD, bin_size, filter);
        }
        // Densify
        let densified = Self::densify_bin(&mut signs);
        (signs, densified)
    }

    /// Get the sketch bins for a sample, but do not transpose
    pub fn get_signs_no_densify<H: RollHash + ?Sized>(
        seq_hashes: &mut H,
        kmer_size: usize,
        filter: &mut Option<KmerFilter>,
        num_bins: u64,
    ) -> Vec<u64> {
        // Setup storage for each k
        let mut signs = vec![u64::MAX; num_bins as usize];
        if let Some(read_filter) = filter {
            read_filter.clear();
        }
        seq_hashes.set_k(kmer_size);

        // Calculate bin minima across all sequence
        let bin_size: u64 = SIGN_MOD.div_ceil(num_bins);
        for hash in seq_hashes.iter() {
            Self::bin_sign(&mut signs, hash % SIGN_MOD, bin_size, filter);
        }

        signs
    }

    /// The name of the sample
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Set a position to be saved in a [`multisketch::MultiSketch`]
    pub fn set_index(&mut self, index: usize) {
        self.index = Some(index);
    }

    /// Get the position that has been saved in an .skd
    pub fn get_index(&self) -> usize {
        self.index.unwrap()
    }

    /// Take the (transposed) sketch, emptying it from the [`Sketch`]
    pub fn get_usigs(&mut self) -> Vec<u16> {
        std::mem::take(&mut self.usigs)
    }

    fn bin_sign(signs: &mut [u64], sign: u64, binsize: u64, read_filter: &mut Option<KmerFilter>) {
        let binidx = (sign / binsize) as usize;
        // log::trace!("sign:{sign} idx:{binidx} curr_sign:{}", signs[binidx]);
        if let Some(filter) = read_filter {
            if sign < signs[binidx] && filter.filter(sign) == Ordering::Equal {
                signs[binidx] = sign;
            }
        } else {
            signs[binidx] = signs[binidx].min(sign);
        }
    }

    #[cfg(not(target_arch = "wasm32"))]
    #[inline(always)]
    fn bit_at_pos(x: u64, pos: u64) -> u64 {
        (x & (1_u64 << pos)) >> pos
    }

    #[inline(always)]
    fn universal_hash(s: u64, t: u64) -> u64 {
        let x = s
            .wrapping_mul(1009)
            .wrapping_add(t.wrapping_mul(1000 * 1000 + 3));
        (x.wrapping_mul(48271).wrapping_add(11)) % ((1 << 31) - 1)
    }

    // TODO could use newer method
    // http://proceedings.mlr.press/v115/mai20a.html
    // https://github.com/zhaoxiaofei/bindash/blob/eb4f81e50b3c42a1fdc00901290b35d0fa9a1e8d/src/hashutils.hpp#L109
    /// Densifies an array of bins
    pub fn densify_bin(signs: &mut [u64]) -> bool {
        let mut minval = u64::MAX;
        let mut maxval = 0;
        for sign in &mut *signs {
            minval = minval.min(*sign);
            maxval = maxval.max(*sign);
        }
        if maxval != u64::MAX {
            false
        } else {
            for i in 0..signs.len() {
                let mut j = i;
                let mut n_attempts = 0;
                while signs[j] == u64::MAX {
                    j = (Self::universal_hash(i as u64, n_attempts as u64) as usize) % signs.len();
                    n_attempts += 1;
                }
                signs[i] = signs[j];
            }
            true
        }
    }

    #[cfg(not(target_arch = "wasm32"))]
    #[inline(always)]
    pub(crate) fn empty_mask_to_vec(empty: &[u64], len: usize) -> Vec<bool> {
        let mut is_empty = vec![false; len];
        for (chunk_idx, bits) in empty.iter().enumerate() {
            let base_idx = chunk_idx * u64::BITS as usize;
            for bit_idx in 0..u64::BITS as usize {
                let idx = base_idx + bit_idx;
                if idx >= len {
                    break;
                }
                is_empty[idx] = ((bits >> bit_idx) & 1) == 1;
            }
        }
        is_empty
    }

    #[cfg(not(target_arch = "wasm32"))]
    pub fn densify_bin_u16(signs: &mut [u16], empty: &[u64]) -> bool {
        let mut is_empty = Self::empty_mask_to_vec(empty, signs.len());
        if !is_empty.iter().any(|&empty_bin| empty_bin) {
            return false;
        }
        if is_empty.iter().all(|&empty_bin| empty_bin) {
            panic!("Cannot densify sketch with all bins empty");
        }
        for i in 0..signs.len() {
            let mut j = i;
            let mut n_attempts = 0;
            while is_empty[j] {
                j = (Self::universal_hash(i as u64, n_attempts as u64) as usize) % signs.len();
                n_attempts += 1;
            }
            if is_empty[i] {
                signs[i] = signs[j];
                is_empty[i] = false;
            }
        }
        true
    }
}

impl fmt::Display for Sketch {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(
            f,
            "{}\t{}\t[{}, {}, {}, {}]\t{}\t{}\t{}\t{}",
            self.name,
            self.seq_length,
            self.acgt[0],
            self.acgt[1],
            self.acgt[3],
            self.acgt[2],
            self.non_acgt,
            self.reads,
            !self.rc,
            self.densified
        )
    }
}

#[cfg(not(target_arch = "wasm32"))]
/// Main function to create sketches from a set of input files, which is parallelised
/// over the input files
pub fn sketch_files(
    output_prefix: &str,
    input_files: &[InputFastx],
    concat_fasta: bool,
    #[cfg(feature = "3di")] convert_pdb: bool,
    k: &[usize],
    sketch_size: u64,
    seq_type: &HashType,
    rc: bool,
    min_count: u16,
    min_qual: u8,
    est_coverage: usize,
    quiet: bool,
) -> Vec<Sketch> {
    // TODO: perhaps, although not 100% necessary, some IO logics could be removed now that we are not supporting
    // arbitrary BBITS values, though this might not imply a very large speed improvement (it's just a couple of variables, I think...)
    let bin_stride = 1;
    let kmer_stride = (sketch_size * 16) as usize;
    let sample_stride = kmer_stride * k.len();

    #[cfg(feature = "3di")]
    let struct_strings = if convert_pdb {
        log::info!("Converting PDB files into 3Di representations");
        Some(pdb_to_3di(input_files).expect("Error converting to 3Di"))
    } else {
        None
    };
    #[cfg(not(feature = "3di"))]
    let struct_strings: Option<Vec<String>> = None;

    log::trace!("{struct_strings:?}");

    // Open output file
    let data_filename = format!("{output_prefix}.skd");
    let mut serial_writer =
        SketchArrayWriter::new(&data_filename, bin_stride, kmer_stride, sample_stride);

    // Set up sender (sketching) and receiver (writing)
    let (tx, rx) = mpsc::channel();
    let mut sketches: Vec<Sketch> = Vec::with_capacity(input_files.len());

    let sketchers = if *seq_type == HashType::DNA {
        Some(k.iter().map(|ik| SketchParams {
                alg: SketchAlg::Bucket,
                rc: rc,
                k : *ik,
                s : sketch_size as usize,
                b : 16,
                seed : 0,
                count : min_count as usize,
                coverage: est_coverage,
                filter_empty: true,
                filter_out_n: true,
            }).collect::<Vec<_>>())
    } else {
        None
    };

    let percent = false;
    let progress_bar = get_progress_bar(input_files.len(), percent, quiet);
    // With thanks to https://stackoverflow.com/a/76963325
    rayon::scope(|s| {
        s.spawn(move |_| {
            input_files
                .par_iter()
                .progress_with(progress_bar)
                .enumerate()
                .map(|(idx, (name, fastxvec))| {
                    match seq_type {
                        HashType::AA(_) | HashType::PDB => {
                            // Read in sequence and set up rolling hash by alphabet type
                            let mut hash_its: Vec<Box<dyn RollHash>> = match seq_type {
                                HashType::DNA => panic!("Not supposed to get here!"),
                                HashType::AA(level) => {
                                    AaHashIterator::new(fastxvec, level.clone(), concat_fasta)
                                        .into_iter()
                                        .map(|it| Box::new(it) as Box<dyn RollHash>)
                                        .collect()
                                }
                                HashType::PDB => {
                                    if let Some(di) = &struct_strings {
                                        log::trace!("Length of string: {}", di.len());
                                        AaHashIterator::from_3di_string(di[idx].clone()) // TODO: clone is not ideal
                                            .into_iter()
                                            .map(|it| Box::new(it) as Box<dyn RollHash>)
                                            .collect()
                                    } else {
                                        AaHashIterator::from_3di_file(fastxvec)
                                            .into_iter()
                                            .map(|it| Box::new(it) as Box<dyn RollHash>)
                                            .collect()
                                    }
                                }
                            };

                            hash_its
                                .iter_mut()
                                .enumerate()
                                .map(|(idx, hash_it)| {
                                    let sample_name = if concat_fasta {
                                        format!("{name}_{}", idx + 1)
                                    } else {
                                        name.to_string()
                                    };
                                    if hash_it.seq_len() == 0 {
                                        panic!("{sample_name} has no valid sequence");
                                    }
                                    // Run the sketching
                                    Sketch::new(&mut **hash_it, &sample_name, k, sketch_size, rc, min_count)
                                })
                                .collect::<Vec<Sketch>>()
                        }
                        HashType::DNA => {
                            // Note that we must not pass as reference the sketchers, as those might be used simultaneously by several parallel threads,
                            // and we might be processing reads and assemblies mixed!
                            let (reads, mut sketches_simd) = sketch_with_simd(fastxvec, min_qual, est_coverage, sketchers.clone().unwrap());
                            vec![Sketch::from_sketch_simd(
                                &mut sketches_simd, 
                                &name.to_string(), 
                                reads
                            )]
                        }
                    }
                })
                .for_each_with(tx, |tx, sketch| {
                    // Emit the sketch results to the writer thread
                    let _ = tx.send(sketch);
                });
        });
        // Write each sketch to the .skd file as it comes in
        for sketch_file in rx {
            // Note double loop as single file may contain multiple samples with concat_fasta
            for mut sketch in sketch_file {
                let index = serial_writer.write_sketch(&sketch.get_usigs());
                sketch.set_index(index);
                // Also append (without usigs) to the metadata, which is Vec<Sketch>
                sketches.push(sketch);
            }
        }
    });

    sketches
}
