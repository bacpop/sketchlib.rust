use std::fmt;
use std::sync::Arc;
use std::cmp::Ordering;

extern crate needletail;
use needletail::{parse_fastx_file, parser::Format};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use indicatif::ParallelProgressIterator;
use serde::{Deserialize, Serialize};

use super::hashing::{encode_base, NtHashIterator};
use crate::bloom_filter::KmerFilter;
use crate::hashing::valid_base;
use crate::io::InputFastx;
use crate::sketch_datafile::SketchArrayFile;

/// Character to use for invalid nucleotides
pub const SEQSEP: u8 = 5;

/// Bin bits (lowest of 64-bits to keep)
pub const BBITS: u64 = 14;
/// Total width of all bins (used as sign % sign_mod)
pub const SIGN_MOD: u64 = (1 << 61) - 1;

#[derive(Serialize, Deserialize, Debug, Default, Clone)]
pub struct Sketch {
    #[serde(skip)]
    usigs: Vec<u64>,
    name: String,
    index: Option<usize>,
    rc: bool,
    reads: bool,
    seq_length: usize,
    densified: bool,
    acgt: [usize; 4],
    non_acgt: usize,
}

impl Sketch {
    pub fn new(
        name: &str,
        files: (&str, Option<&String>),
        kmer_lengths: &[usize],
        sketch_size: u64,
        min_qual: u8,
        min_count: u16,
        rc: bool,
    ) -> Self {
        let size_u64 = (sketch_size * BBITS) as usize * kmer_lengths.len();
        let mut sketch = Self {
            usigs: Vec::with_capacity(size_u64),
            name: name.to_string(),
            index: None,
            rc,
            reads: false,
            seq_length: 0,
            densified: false,
            acgt: [0; 4],
            non_acgt: 0,
        };

        // Check if we're working with reads, and initalise the filter if so
        let mut reader_peek =
            parse_fastx_file(files.0).unwrap_or_else(|_| panic!("Invalid path/file: {}", files.0));
        let seq_peek = reader_peek
            .next()
            .expect("Invalid FASTA/Q record")
            .expect("Invalid FASTA/Q record");
        if seq_peek.format() == Format::Fastq {
            sketch.reads = true;
        }
        let mut filter = match sketch.reads {
            true => {
                let mut kmer_filter = KmerFilter::new(min_count);
                kmer_filter.init();
                Some(kmer_filter)
            }
            false => None,
        };

        // Read sequence into memory (as we go through multiple times)
        log::debug!("Preprocessing sequence");
        let mut sequence: Vec<u8> = Vec::new();
        sketch.add_seq(files.0, &mut sequence, min_qual);
        if let Some(filename) = files.1 {
            sketch.add_seq(filename, &mut sequence, min_qual);
        }

        sketch.seq_length = sketch.acgt.iter().sum();
        if sketch.seq_length == 0 {
            panic!("{} has no valid sequence", files.0);
        }

        // Build the sketches across k-mer lengths
        let mut minhash_sum = 0.0;
        let num_bins: u64 = sketch_size * (u64::BITS as u64);
        let bin_size: u64 = (SIGN_MOD + num_bins - 1) / num_bins;
        for k in kmer_lengths {
            // Calculate bin minima across all sequence
            let mut signs = vec![u64::MAX; num_bins as usize];
            let hash_it = NtHashIterator::new(&sequence, *k, rc);
            // let bar = ProgressBar::new(sequence.len() as u64);
            for hash in hash_it {
                // bar.inc(1);
                Self::bin_sign(&mut signs, hash % SIGN_MOD, bin_size, &mut filter);
            }
            // bar.finish();

            // Densify
            sketch.densified |= Self::densify_bin(&mut signs);
            minhash_sum += (signs[0] as f64) / (SIGN_MOD as f64);

            // Transpose the bins and save to the sketch map
            log::debug!("Transposing bins");
            let mut usigs = vec![0; (sketch_size * BBITS) as usize];
            Self::fill_usigs(&mut usigs, &signs);
            sketch.usigs.append(&mut usigs);
        }

        // Estimate of sequence length from read data
        if sketch.reads {
            sketch.seq_length = ((kmer_lengths.len() as f64) / minhash_sum) as usize;
        }

        sketch
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn set_index(&mut self, index: usize) {
        self.index = Some(index);
    }

    pub fn get_index(&self) -> usize {
        self.index.unwrap()
    }

    // Take the (transposed) sketch, emptying it from the [`Sketch`]
    pub fn get_usigs(&mut self) -> Vec<u64> {
        std::mem::take(&mut self.usigs)
    }

    fn add_seq(&mut self, filename: &str, sequence: &mut Vec<u8>, min_qual: u8) {
        let mut reader =
            parse_fastx_file(filename).unwrap_or_else(|_| panic!("Invalid path/file: {filename}"));
        while let Some(record) = reader.next() {
            let seqrec = record.expect("Invalid FASTA/Q record");
            if let Some(quals) = seqrec.qual() {
                for (base, qual) in seqrec.seq().iter().zip(quals) {
                    if *qual >= min_qual {
                        if valid_base(*base) {
                            let encoded_base = encode_base(*base);
                            self.acgt[encoded_base as usize] += 1;
                            sequence.push(encoded_base)
                        } else {
                            self.non_acgt += 1;
                            sequence.push(SEQSEP);
                        }
                    } else {
                        sequence.push(SEQSEP);
                    }
                }
            } else {
                for base in seqrec.seq().iter() {
                    if valid_base(*base) {
                        let encoded_base = encode_base(*base);
                        self.acgt[encoded_base as usize] += 1;
                        sequence.push(encoded_base)
                    } else {
                        self.non_acgt += 1;
                        sequence.push(SEQSEP);
                    }
                }
            }

            sequence.push(SEQSEP);
        }
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

    #[inline(always)]
    fn bit_at_pos(x: u64, pos: u64) -> u64 {
        (x & (1_u64 << pos)) >> pos
    }

    fn fill_usigs(usigs: &mut [u64], signs: &[u64]) {
        eprintln!("{signs:?}");
        for (sign_index, sign) in signs.iter().enumerate() {
            let leftshift = sign_index % (u64::BITS as usize);
            for i in 0..BBITS {
                let orval = Self::bit_at_pos(*sign, i) << leftshift;
                usigs[sign_index / (u64::BITS as usize) * (BBITS as usize) + (i as usize)] |= orval;
            }
        }
        eprintln!("{usigs:?}");
    }

    #[inline(always)]
    fn universal_hash(s: u64, t: u64) -> u64 {
        let x = (1009) * s + (1000 * 1000 + 3) * t;
        (48271 * x + 11) % ((1 << 31) - 1)
    }

    // TODO could use newer method
    //  http://proceedings.mlr.press/v115/mai20a.html
    // https://github.com/zhaoxiaofei/bindash/blob/eb4f81e50b3c42a1fdc00901290b35d0fa9a1e8d/src/hashutils.hpp#L109
    fn densify_bin(signs: &mut [u64]) -> bool {
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

pub fn sketch_files(
    output_prefix: &str,
    input_files: &[InputFastx],
    k: &[usize],
    sketch_size: u64,
    rc: bool,
    min_count: u16,
    min_qual: u8,
) -> Vec<Sketch> {
    let bin_stride = 1;
    let kmer_stride = (sketch_size * BBITS) as usize;
    let sample_stride = kmer_stride * k.len();

    let data_filename = format!("{output_prefix}.skd");
    let serial_writer = Arc::new(SketchArrayFile::new(
        &data_filename,
        bin_stride,
        kmer_stride,
        sample_stride,
    ));

    let mut sketches: Vec<Sketch> = input_files
        .par_iter()
        .progress_count(input_files.len() as u64)
        .map(|(name, fastx1, fastx2)| {
            let mut sketch = Sketch::new(
                name,
                (&fastx1, fastx2.as_ref()),
                k,
                sketch_size,
                min_qual,
                min_count,
                rc,
            );
            let writer = Arc::clone(&serial_writer);
            let index = writer.write_sketch(&sketch.get_usigs());
            sketch.set_index(index);
            sketch
        })
        .collect();
    // Sort to be in the same order as they were saved to the file
    sketches.sort_unstable_by_key(|k| k.get_index());
    sketches
}
