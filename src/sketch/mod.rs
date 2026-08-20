//! Methods to sketch samples, save/load sketches
use std::cmp::Ordering;
use std::fmt;
#[cfg(not(target_arch = "wasm32"))]
use std::sync::mpsc;

use aligned_vec::{AVec, ConstAlign};
#[cfg(not(target_arch = "wasm32"))]
use indicatif::ParallelProgressIterator;
#[cfg(not(target_arch = "wasm32"))]
use needletail::parse_fastx_file;
#[cfg(not(target_arch = "wasm32"))]
use needletail::parser::Format;
#[cfg(not(target_arch = "wasm32"))]
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use super::hashing::{bloom_filter::KmerFilter, RollHash};

#[cfg(not(target_arch = "wasm32"))]
use super::hashing::nthash_iterator::NtHashIterator;
use super::hashing::HashType;

#[cfg(not(target_arch = "wasm32"))]
use crate::hashing::aahash_iterator::AaHashIterator;
#[cfg(not(target_arch = "wasm32"))]
use crate::io::InputFastx;
#[cfg(not(target_arch = "wasm32"))]
use crate::io::NeedletailIterator;
#[cfg(feature = "3di")]
use crate::structures::pdb_to_3di;
#[cfg(not(target_arch = "wasm32"))]
use crate::utils::get_progress_bar;

pub mod multisketch;

pub mod sketch_datafile;
#[cfg(not(target_arch = "wasm32"))]
use self::sketch_datafile::SketchArrayWriter;

/// Bin bits (lowest of 64-bits to keep)
pub const BIN_BITS: usize = 16;
/// Legacy bin-packing width used by `.skd`/`.skm` databases written before
/// v0.4.0 (see [`multisketch::MIN_SKETCH_VERSION`]), prior to the `BIN_BITS`
/// 14->16 change. Distance-calculation use only — there is no code path
/// that creates new sketches with this width; sketch creation always uses
/// [`BIN_BITS`].
pub const LEGACY_BIN_BITS: usize = 14;
/// Byte alignment used for bit-packed sketch buffers.
pub const SKETCH_ALIGNMENT: usize = 64;

/// Aligned storage for bit-packed sketch buffers.
pub type SketchVec = AVec<u64, ConstAlign<SKETCH_ALIGNMENT>>;

pub(crate) fn aligned_sketch_vec() -> SketchVec {
    SketchVec::new(SKETCH_ALIGNMENT)
}

pub(crate) fn aligned_sketch_vec_with_capacity(capacity: usize) -> SketchVec {
    SketchVec::with_capacity(SKETCH_ALIGNMENT, capacity)
}

pub(crate) fn aligned_sketch_vec_filled(len: usize, value: u64) -> SketchVec {
    let mut vec = aligned_sketch_vec_with_capacity(len);
    vec.resize(len, value);
    vec
}

/// Get the number of elements in the sketch vectors for a given sketch size
///
/// Returns a tuple:
/// - First element is sketch size divided by 64 (used in Jaccard fn)
/// - Second element is the number of bins (rounded up to the
///   nearest 64)
/// - Third element is the number of transposed bins
///
/// # Arguments
///
/// - `sketch_size` -- number of bins wanted.
pub fn num_bins(sketch_size: u64) -> (u64, u64, u64) {
    let sketchsize64 = sketch_size.div_ceil(u64::BITS as u64);
    let signs_size = sketchsize64 * (u64::BITS as u64);
    let usigs_size = sketchsize64 * (BIN_BITS as u64);
    (sketchsize64, signs_size, usigs_size)
}

#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
/// Options and parameters used in sketching
pub struct SketchingOpts {
    /// Sample name
    pub name: String,
    /// Concatenate records in a multifasta file?
    pub concat_fasta: bool,
    /// k-mer sizes to use
    pub k_vals: Vec<usize>,
    /// Sketch size
    pub sketch_size: u64,
    /// Sequence type (DNA, AA)
    pub seq_type: HashType,
    /// Add reverse complements to sketch?
    pub add_rc: bool,
    /// Minimum k-mer count to use for fastq input
    pub min_count: u16,
    /// Minimum quality score to use for fastq input
    pub min_qual: u8,
    /// Is the input reads?
    pub is_reads: bool,
}

impl Default for SketchingOpts {
    fn default() -> SketchingOpts {
        SketchingOpts {
            name: String::new(),
            concat_fasta: false,
            k_vals: Vec::new(),
            sketch_size: crate::cli::DEFAULT_SKETCHSIZE,
            seq_type: HashType::DNA,
            add_rc: crate::cli::DEFAULT_STRAND,
            min_count: crate::cli::DEFAULT_MINCOUNT,
            min_qual: crate::cli::DEFAULT_MINQUAL,
            is_reads: false,
        }
    }
}

/// A single sample's sketch
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct Sketch {
    #[serde(skip, default = "aligned_sketch_vec")]
    usigs: SketchVec,
    name: String,
    index: Option<usize>,
    /// Sequence and reverse complement used
    pub rc: bool,
    /// Built from reads (fastq)
    pub reads: bool,
    /// Assembled sequence length (fasta: exact; fastq: estimate)
    pub seq_length: usize,
    /// Sketch required densification
    pub densified: bool,
    /// Number of [A, C, G, T] bases
    pub acgt: [usize; 4],
    /// Number of non-ACGT bases
    pub non_acgt: usize,
}

impl Default for Sketch {
    fn default() -> Self {
        Self {
            usigs: aligned_sketch_vec(),
            name: String::default(),
            index: None,
            rc: false,
            reads: false,
            seq_length: 0,
            densified: false,
            acgt: [0; 4],
            non_acgt: 0,
        }
    }
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
    ) -> anyhow::Result<Self> {
        let (_sketchsize64, num_bins, usigs_size) = num_bins(sketch_size);
        let flattened_size_u64 = usigs_size as usize * kmer_lengths.len();
        let mut usigs = aligned_sketch_vec_with_capacity(flattened_size_u64);

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
            let (signs, k_densified) = Self::get_signs(seq_hashes, *k, &mut read_filter, num_bins)?;
            densified |= k_densified;
            minhash_sum += (signs[0] as f64) / (u64::MAX as f64);

            // Transpose the bins and save to the sketch map
            log::debug!("Transposing bins");
            let mut kmer_usigs = aligned_sketch_vec_filled(usigs_size as usize, 0);
            Self::fill_usigs(&mut kmer_usigs, &signs);
            usigs.append(&mut kmer_usigs);
        }
        let (reads, acgt, non_acgt) = seq_hashes.sketch_data();

        // Estimate of sequence length from read data
        let seq_length = if reads {
            ((kmer_lengths.len() as f64) / minhash_sum) as usize
        } else {
            seq_hashes.seq_len()
        };

        Ok(Self {
            usigs,
            name: name.to_string(),
            index: None,
            rc,
            reads,
            seq_length,
            densified,
            acgt,
            non_acgt,
        })
    }

    /// Get the sketch bins for a sample, but do not transpose
    pub fn get_signs<H: RollHash + ?Sized>(
        seq_hashes: &mut H,
        kmer_size: usize,
        filter: &mut Option<KmerFilter>,
        num_bins: u64,
    ) -> anyhow::Result<(Vec<u64>, bool)> {
        // Setup storage for each k
        let mut signs = vec![u64::MAX; num_bins as usize];
        if let Some(read_filter) = filter {
            read_filter.clear();
        }
        seq_hashes.set_k(kmer_size)?;

        // Calculate bin minima across all sequence
        for hash in seq_hashes.iter() {
            Self::bin_sign(&mut signs, hash, num_bins, filter);
        }
        // Densify
        let densified = Self::densify_bin(&mut signs);
        Ok((signs, densified))
    }

    /// Get the sketch bins for a sample, but do not transpose
    pub fn get_signs_no_densify<H: RollHash + ?Sized>(
        seq_hashes: &mut H,
        kmer_size: usize,
        filter: &mut Option<KmerFilter>,
        num_bins: u64,
    ) -> anyhow::Result<Vec<u64>> {
        // Setup storage for each k
        let mut signs = vec![u64::MAX; num_bins as usize];
        if let Some(read_filter) = filter {
            read_filter.clear();
        }
        seq_hashes.set_k(kmer_size)?;

        // Calculate bin minima across all sequence
        for hash in seq_hashes.iter() {
            Self::bin_sign(&mut signs, hash, num_bins, filter);
        }

        Ok(signs)
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
    pub fn get_usigs(&mut self) -> SketchVec {
        // Move the old aligned buffer out and install an empty one; elements are not copied.
        std::mem::replace(&mut self.usigs, aligned_sketch_vec())
    }

    fn bin_sign(signs: &mut [u64], sign: u64, num_bins: u64, read_filter: &mut Option<KmerFilter>) {
        let binidx = Self::bin_index(sign, num_bins);
        // log::trace!("sign:{sign} idx:{binidx} curr_sign:{}", signs[binidx]);
        if let Some(filter) = read_filter {
            if sign < signs[binidx] && filter.filter(sign) == Ordering::Equal {
                signs[binidx] = sign;
            }
        } else {
            signs[binidx] = signs[binidx].min(sign);
        }
    }

    #[inline(always)]
    fn bin_index(sign: u64, num_bins: u64) -> usize {
        (((sign as u128) * (num_bins as u128)) >> u64::BITS) as usize
    }

    #[inline(always)]
    fn bit_at_pos(x: u64, pos: u64) -> u64 {
        (x & (1_u64 << pos)) >> pos
    }

    fn fill_usigs(usigs: &mut [u64], signs: &[u64]) {
        debug_assert_eq!(signs.len() % (u64::BITS as usize), 0);
        debug_assert_eq!(usigs.len(), signs.len() / (u64::BITS as usize) * BIN_BITS);

        for (usig_chunk, sign_chunk) in usigs
            .as_chunks_mut::<BIN_BITS>()
            .0
            .iter_mut()
            .zip(signs.as_chunks::<{ u64::BITS as usize }>().0.iter())
        {
            for (bit_pos, usig) in usig_chunk.iter_mut().enumerate() {
                *usig = sign_chunk
                    .iter()
                    .enumerate()
                    .fold(0, |bits, (sign_index, &sign)| {
                        bits | (Self::bit_at_pos(sign, bit_pos as u64) << sign_index)
                    });
            }
        }
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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sketch_vec_is_64_byte_aligned_after_allocation() {
        let mut sketch_vec = aligned_sketch_vec_with_capacity(16);
        sketch_vec.resize(16, 0);
        assert_eq!(sketch_vec.as_ptr() as usize % SKETCH_ALIGNMENT, 0);
    }

    #[test]
    fn bin_index_uses_full_u64_hash_space() {
        assert_eq!(Sketch::bin_index(0, 64), 0);
        assert_eq!(Sketch::bin_index(1_u64 << 63, 64), 32);
        assert_eq!(Sketch::bin_index(u64::MAX, 64), 63);
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
/// Create sketches from an iterator over sequence data
///
/// # Examples
///
/// ## Filter an assembly and sketch only some contigs
/// ```rust
/// use sketchlib::sketch::{Sketch, SketchingOpts};
/// use sketchlib::sketch::sketch_data;
///
/// use std::collections::HashSet;
/// use std::path::{Path, PathBuf};
///
/// // Iterator for needletail records
/// pub struct NeedletailFilterIterator {
///     reader: Box<dyn needletail::FastxReader>,
///     want_ids: HashSet<u32>,
///     current_idx: u32,
/// }
///
/// impl NeedletailFilterIterator {
///     // Construct from needletail readers
///     pub fn new(
///         reader: Box<dyn needletail::FastxReader>,
///         want_ids: HashSet<u32>,
///     ) -> Self {
///         Self {
///             reader,
///             want_ids,
///             current_idx: 0_u32,
///         }
///     }
/// }
///
/// impl Iterator for NeedletailFilterIterator {
///     type Item = (Vec<u8>, Option<Vec<u8>>);
///
///     fn next(
///         &mut self,
///     ) -> Option<(Vec<u8>, Option<Vec<u8>>)> {
///         while let Some(try_record) = self.reader.next() {
///             if self.want_ids.contains(&self.current_idx) {
///                 let record = try_record.expect("Invalid fastX record");
///                 let seq = record.seq();
///                 let qual = record.qual().map(|qual| qual.to_vec());
///                 self.current_idx += 1;
///                 return Some((seq.to_vec(), qual))
///             }
///             self.current_idx += 1;
///         }
///         None
///     }
/// }
///
/// // Sketch a fastX file filtered by the index of reads to include in the sketch
/// pub fn sketch_reads_with_filter(
///     fastx_path: &Path,
///     want_ids: HashSet<u32>,
///     opts: SketchingOpts,
/// ) -> Vec<Sketch> {
///     let reader = needletail::parse_fastx_file(fastx_path).unwrap();
///     let mut filtered_iters = vec![NeedletailFilterIterator::new(reader, want_ids)];
///
///     sketch_data(&mut filtered_iters, opts).unwrap()
/// }
///
/// let fastq_path_str = "tests/test_files_in/14412_3#82.contigs_velvet.fa.gz";
/// let mut fastq_path = PathBuf::from(fastq_path_str);
///
/// let mut opts = SketchingOpts::default();
/// opts.k_vals = vec![21_usize, 31_usize, 51_usize];
/// opts.name = fastq_path_str.to_string();
/// # opts.sketch_size = 1;
///
/// let want_ids: HashSet<u32> = HashSet::from_iter(vec![0_32, 5_u32, 3_u32].into_iter());
/// let sketch = sketch_reads_with_filter(&fastq_path, want_ids, opts);
///
/// # let mut sketch = sketch;
/// # assert_eq!(sketch.len(), 1);
/// # assert_eq!(sketch[0].name(), fastq_path_str);
/// ```
pub fn sketch_data<I: Iterator<Item = (Vec<u8>, Option<Vec<u8>>)>>(
    records_readers: &mut [I],
    opts: SketchingOpts,
    #[cfg(feature = "3di")] convert_pdb: bool,
    #[cfg(feature = "3di")] struct_string: Option<String>,
) -> anyhow::Result<Vec<Sketch>> {
    // Read in sequence and set up rolling hash by alphabet type
    let mut hash_its: Vec<Box<dyn RollHash>> = match opts.seq_type {
        HashType::DNA => NtHashIterator::new(
            records_readers,
            opts.k_vals[0],
            opts.add_rc,
            opts.min_qual,
            opts.is_reads,
        )?
        .into_iter()
        .map(|it| Box::new(it) as Box<dyn RollHash>)
        .collect(),
        HashType::AA(level) => AaHashIterator::new(
            records_readers,
            &opts.name,
            level.clone(),
            opts.concat_fasta,
        )
        .into_iter()
        .map(|it| Box::new(it) as Box<dyn RollHash>)
        .collect(),
        HashType::PDB => {
            #[cfg(feature = "3di")]
            if let Some(di) = &struct_string {
                AaHashIterator::from_3di_string(di.clone()) // TODO: clone is not ideal
                    .into_iter()
                    .map(|it| Box::new(it) as Box<dyn RollHash>)
                    .collect()
            } else {
                AaHashIterator::from_3di_file(records_readers, &opts.name)
                    .into_iter()
                    .map(|it| Box::new(it) as Box<dyn RollHash>)
                    .collect()
            }
            #[cfg(not(feature = "3di"))]
            AaHashIterator::from_3di_file(records_readers, &opts.name)
                .into_iter()
                .map(|it| Box::new(it) as Box<dyn RollHash>)
                .collect()
        }
    };

    let mut sketches: Vec<Sketch> = Vec::with_capacity(hash_its.len());
    for(idx, hash_it) in hash_its.iter_mut().enumerate() {
        let sample_name = if opts.concat_fasta {
            format!("{}_{}", opts.name, idx + 1)
        } else {
            opts.name.to_string()
        };
        if hash_it.seq_len() == 0 {
            return Err(anyhow::anyhow!("{sample_name} has no valid sequence"))
        }
        // Run the sketching
        let sketch = Sketch::new(
            &mut **hash_it,
            &sample_name,
            &opts.k_vals,
            opts.sketch_size,
            opts.add_rc,
            opts.min_count,
        )?;
        sketches.push(sketch);
    }
    Ok(sketches)
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
    quiet: bool,
) -> Vec<Sketch> {
    let bin_stride = 1;
    let kmer_stride = sketch_size as usize * BIN_BITS;
    let sample_stride = kmer_stride * k.len();

    #[cfg(feature = "3di")]
    let struct_strings = if convert_pdb {
        log::info!("Converting PDB files into 3Di representations");
        Some(pdb_to_3di(input_files).expect("Error converting to 3Di"))
    } else {
        None
    };
    #[cfg(feature = "3di")]
    log::trace!("{struct_strings:?}");

    // Open output file
    let data_filename = format!("{output_prefix}.skd");
    let mut serial_writer =
        SketchArrayWriter::new(&data_filename, bin_stride, kmer_stride, sample_stride);

    // Set up sender (sketching) and receiver (writing)
    let (tx, rx) = mpsc::channel();
    let mut sketches: Vec<Sketch> = Vec::with_capacity(input_files.len());

    let percent = false;
    let progress_bar = get_progress_bar(input_files.len(), percent, quiet);
    // With thanks to https://stackoverflow.com/a/76963325
    rayon::scope(|s| {
        s.spawn(move |_| {
            input_files
                .par_iter()
                .progress_with(progress_bar)
                .enumerate()
                .map(|(_idx, (name, fastxvec))| {
                    // Read in sequence and set up rolling hash by alphabet type

                    let reads = if seq_type == &HashType::DNA {
                        // Check if we're working with reads, and initalise the filter if so
                        let mut reader_peek = parse_fastx_file(fastxvec[0].clone())
                            .unwrap_or_else(|_| panic!("Invalid path/file: {}", fastxvec[0]));
                        let seq_peek = reader_peek
                            .next()
                            .expect("Invalid FASTA/Q record")
                            .expect("Invalid FASTA/Q record");
                        let mut reads = false;
                        if seq_peek.format() == Format::Fastq {
                            reads = true;
                            if fastxvec.len() > 2 {
                                panic!("Input files are reads, but there are more than two input files");
                            }
                        }
                        reads
                    } else {
                        false
                    };

                    let opts = SketchingOpts {
                        name: name.clone(),
                        k_vals: k.to_vec(),
                        seq_type: seq_type.clone(),
                        is_reads: reads,
                        concat_fasta,
                        sketch_size,
                        add_rc: rc,
                        min_count,
                        min_qual,
                    };

                    let mut records_readers = fastxvec.iter().map(|file| {
                        let reader = parse_fastx_file(file).unwrap_or_else(|_| panic!("Invalid path/file: {file}"));
                        NeedletailIterator::new(reader)
                    }).collect::<Vec<NeedletailIterator>>();

                    #[cfg(feature = "3di")]
                    let di = struct_strings.as_ref().map(|structs| structs[_idx].clone());

                    sketch_data(
                        &mut records_readers,
                        opts,
                        #[cfg(feature = "3di")]
                        convert_pdb,
                        #[cfg(feature = "3di")]
                        di,
                    ).unwrap()
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
