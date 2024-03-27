use std::ops::Mul;
use std::{cmp::Ordering, collections::HashMap};

use serde::{Deserialize, Serialize};
extern crate needletail;
use indicatif::{ParallelProgressIterator, ProgressBar, ProgressIterator};
use needletail::{parse_fastx_file, parser::Format};
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};

use super::hashing::{encode_base, NtHashIterator};
use crate::bloom_filter::KmerFilter;
use crate::hashing::valid_base;
use crate::io::InputFastx;
use crate::multisketch::MultiSketch;

pub const SEQSEP: u8 = 5;

/// Bin bits
pub const BBITS: u64 = 14;
pub const SIGN_MOD: u64 = (1 << 61) - 1;

#[derive(Serialize, Deserialize, Debug, Default)]
pub struct Sketch {
    #[serde(skip)] // In Multisketch
    sketch_size: u64,
    #[serde(skip)]
    kmer_sketches: HashMap<usize, Vec<u64>>,
    index: usize,
    rc: bool,
    reads: bool,
    seq_length: usize,
    missing_bases: usize,
    densified: bool,
    acgt: [usize; 4],
    non_acgt: usize,
}

impl Sketch {
    pub fn new(
        files: (&str, Option<&String>),
        kmer_lengths: &[usize],
        index: usize,
        sketch_size: u64,
        min_qual: u8,
        min_count: u16,
        rc: bool,
    ) -> Self {
        let mut sketch = Self {
            sketch_size: sketch_size,
            kmer_sketches: HashMap::with_capacity(kmer_lengths.len()),
            index,
            rc: rc,
            reads: false,
            seq_length: 0,
            missing_bases: 0,
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
        // TODO more efficient to pass the filter and init if reads, but maybe
        // when multithreaded this is better?
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
            sketch.kmer_sketches.insert(*k, usigs);
        }

        // Estimate of sequence length from read data
        if sketch.reads {
            sketch.seq_length = ((kmer_lengths.len() as f64) / minhash_sum) as usize;
        }

        sketch
    }

    // TODO: would be fun to try putting signs in a roaringbitmap then just
    // doing intersection_len() as this function
    pub fn jaccard_dist(&self, other: &Sketch, k: usize) -> f64 {
        let unionsize = (u64::BITS as u64 * self.sketch_size) as f64;
        let mut samebits: u32 = 0;
        for i in 0..self.sketch_size {
            let mut bits: u64 = !0;
            for j in 0..BBITS {
                bits &= !(self.kmer_sketches[&k][(i * BBITS + j) as usize]
                    ^ other.kmer_sketches[&k][(i * BBITS + j) as usize]);
            }
            samebits += bits.count_ones();
        }
        let maxnbits = self.sketch_size as u32 * u64::BITS;
        let expected_samebits = maxnbits >> BBITS;
        if expected_samebits != 0 {
            samebits as f64
        } else {
            let diff = samebits.saturating_sub(maxnbits);
            let intersize = (diff * maxnbits / (maxnbits - maxnbits)) as f64;
            intersize / unionsize
        }
    }

    pub fn core_acc_dist(&self, other: &Sketch) -> (f64, f64) {
        if self.kmer_sketches.len() < 2 {
            panic!("Need at least two k-mer lengths to calculate core/accessory distances");
        }
        let mut xsum = 0.0;
        let mut ysum = 0.0;
        let mut xysum = 0.0;
        let mut xsquaresum = 0.0;
        let mut ysquaresum = 0.0;
        let mut n = 0.0;
        let tolerance = 5.0 / (self.sketch_size as f64);
        for k in self.kmer_sketches.keys() {
            let y = self.jaccard_dist(other, *k).exp();
            if y < tolerance {
                break;
            }
            let k_fl = *k as f64;
            xsum += k_fl;
            ysum += y;
            xysum += k_fl * y;
            xsquaresum += k_fl * k_fl;
            ysquaresum += y * y;
            n += 1.0;
        }
        Self::simple_linear_regression(xsum, ysum, xysum, xsquaresum, ysquaresum, n)
    }

    pub fn index(&self) -> usize {
        self.index
    }

    pub fn get_usigs(&mut self) -> &HashMap<usize, Vec<u64>> {
        &self.kmer_sketches
    }

    fn simple_linear_regression(
        xsum: f64,
        ysum: f64,
        xysum: f64,
        xsquaresum: f64,
        ysquaresum: f64,
        n: f64,
    ) -> (f64, f64) {
        let xbar = xsum / n;
        let ybar = ysum / n;
        let x_diff = xsquaresum - xsum * xsum / n;
        let y_diff = ysquaresum - ysum * ysum / n;
        let xstddev = ((xsquaresum - xsum * xsum / n) / n).sqrt();
        let ystddev = ((ysquaresum - ysum * ysum / n) / n).sqrt();
        let r = (xysum - xsum * ysum / n) / (x_diff * y_diff).sqrt();
        let beta = r * ystddev / xstddev;
        let alpha = -beta * xbar + ybar;

        let (mut core, mut acc) = (0.0, 0.0);
        if beta < 0.0 {
            core = 1.0 - beta.exp();
        }
        if alpha < 0.0 {
            acc = 1.0 - alpha.exp();
        }
        (core, acc)
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
        x & (1 << pos) >> pos
    }

    fn fill_usigs(usigs: &mut [u64], signs: &[u64]) {
        for (sign_index, sign) in signs.iter().enumerate() {
            let leftshift = sign_index % (u64::BITS as usize);
            for i in 0..BBITS {
                let orval = Self::bit_at_pos(*sign, i) << leftshift;
                usigs[sign_index / (u64::BITS as usize) * (BBITS as usize) + (i as usize)] |= orval;
            }
        }
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

// TODO may want to write these one by one rather than storing them all
pub fn sketch_files(
    sketch_container: &mut MultiSketch,
    input_files: &[InputFastx],
    k: &[usize],
    sketch_size: u64,
    rc: bool,
    min_count: u16,
    min_qual: u8,
) {
    input_files
        .par_iter()
        .enumerate()
        .progress_count(input_files.len() as u64)
        .map(|(idx, fastx)| {
            let mut sketch = Sketch::new(
                (&fastx.1, fastx.2.as_ref()),
                k,
                idx,
                sketch_size,
                min_qual,
                min_count,
                rc,
            );
            // TODO use name &fastx.0 to add to multisketch
            sketch_container.add_sketch(&fastx.0, &mut sketch);
        });
}
