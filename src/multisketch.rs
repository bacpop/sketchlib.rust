use std::error::Error;
use std::fmt;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::mem;

use hashbrown::HashMap;
use serde::{Deserialize, Serialize};

use crate::hashing::HashType;
use crate::sketch::{Sketch, BBITS};
use crate::sketch_datafile::SketchArrayFile;

pub fn jaccard_dist(sketch1: &[u64], sketch2: &[u64], sketch_size: u64) -> f32 {
    let unionsize = (u64::BITS as u64 * sketch_size) as f32;
    let mut samebits: u32 = 0;
    for i in 0..sketch_size {
        let mut bits: u64 = !0;
        for j in 0..BBITS {
            bits &= !(sketch1[(i * BBITS + j) as usize] ^ sketch2[(i * BBITS + j) as usize]);
        }
        samebits += bits.count_ones();
    }
    let maxnbits = sketch_size as u32 * u64::BITS;
    let expected_samebits = maxnbits >> BBITS;

    log::trace!(
        "samebits:{samebits} expected_samebits:{expected_samebits} maxnbits:{maxnbits}"
    );
    if expected_samebits != 0 {
        samebits as f32
    } else {
        let diff = samebits.saturating_sub(expected_samebits);
        let intersize = (diff * maxnbits) as f32 / (maxnbits - expected_samebits) as f32;
        log::trace!("intersize:{intersize} unionsize:{unionsize}");
        intersize / unionsize
    }
}

pub fn core_acc_dist(ref_sketches: &MultiSketch, query_sketches: &MultiSketch, ref_sketch_idx: usize, query_sketch_idx: usize) -> (f32, f32) {
    if ref_sketches.kmer_lengths.len() < 2 {
        panic!("Need at least two k-mer lengths to calculate core/accessory distances");
    }
    let (mut xsum, mut ysum, mut xysum, mut xsquaresum, mut ysquaresum, mut n) =
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    let tolerance = (5.0_f32 / (ref_sketches.sketch_size as f32)).ln();
    for (k_idx, k) in ref_sketches.kmer_lengths.iter().enumerate() {
        let y = jaccard_dist(ref_sketches.get_sketch_slice(ref_sketch_idx, k_idx), query_sketches.get_sketch_slice(query_sketch_idx, k_idx), ref_sketches.sketch_size).ln();
        if y < tolerance {
            break;
        }
        let k_fl = *k as f32;
        xsum += k_fl;
        ysum += y;
        xysum += k_fl * y;
        xsquaresum += k_fl * k_fl;
        ysquaresum += y * y;
        n += 1.0;
    }
    simple_linear_regression(xsum, ysum, xysum, xsquaresum, ysquaresum, n)
}

fn simple_linear_regression(
    xsum: f32,
    ysum: f32,
    xysum: f32,
    xsquaresum: f32,
    ysquaresum: f32,
    n: f32,
) -> (f32, f32) {
    log::trace!(
        "xsum:{xsum} ysum:{ysum} xysum:{xysum} xsquaresum:{xsquaresum} ysquaresum:{ysquaresum}"
    );
    let xbar = xsum / n;
    let ybar = ysum / n;
    let x_diff = xsquaresum - xsum * xsum / n;
    let y_diff = ysquaresum - ysum * ysum / n;
    let xstddev = ((xsquaresum - xsum * xsum / n) / n).sqrt();
    let ystddev = ((ysquaresum - ysum * ysum / n) / n).sqrt();
    let r = (xysum - xsum * ysum / n) / (x_diff * y_diff).sqrt();
    let beta = r * ystddev / xstddev;
    let alpha = -beta * xbar + ybar;
    log::trace!("r:{r} alpha:{alpha} beta:{beta}");

    let (mut core, mut acc) = (0.0, 0.0);
    if beta < 0.0 {
        core = 1.0 - beta.exp();
    }
    if alpha < 0.0 {
        acc = 1.0 - alpha.exp();
    }
    (core, acc)
}

#[derive(Serialize, Deserialize)]
pub struct MultiSketch {
    pub sketch_size: u64,
    kmer_lengths: Vec<usize>,
    sketch_metadata: Vec<Sketch>,
    name_map: HashMap<String, usize>,
    #[serde(skip)]
    // NB: another way to do this is with the ouroboros crate, which allows this to reference self
    // But this requires manual impl for ser and deser, and does the same indirection as an index anyway so not worth it
    block_reindex: Option<Vec<usize>>,
    #[serde(skip)]
    sketch_bins: Vec<u64>,
    bin_stride: usize,
    kmer_stride: usize,
    sample_stride: usize,
    sketch_version: String,
    hash_type: HashType,
}

impl MultiSketch {
    pub fn new(
        names: &mut Vec<String>,
        sketches: &mut Vec<Sketch>,
        sketch_size: u64,
        kmer_lengths: &[usize],
    ) -> Self {
        let mut name_map = HashMap::with_capacity(names.len());
        for (name, sketch) in names.iter_mut().zip(sketches.iter()) {
            name_map.insert(mem::take(name), sketch.get_index());
        }

        let kmer_stride = (sketch_size * BBITS) as usize;
        Self {
            sketch_size,
            kmer_lengths: kmer_lengths.to_vec(),
            sketch_metadata: mem::take(sketches),
            name_map,
            block_reindex: None,
            sketch_bins: Vec::new(),
            bin_stride: 1,
            kmer_stride,
            sample_stride: kmer_stride * kmer_lengths.len(),
            sketch_version: env!("CARGO_PKG_VERSION").to_string(),
            hash_type: HashType::DNA,
        }
    }

    /// Saves the metadata
    pub fn save_metadata(&self, file_prefix: &str) -> Result<(), Box<dyn Error>> {
        let filename = format!("{}.skm", file_prefix);
        log::info!("Saving sketch metadata to {filename}");
        let serial_file = BufWriter::new(File::create(filename)?);
        let mut compress_writer = snap::write::FrameEncoder::new(serial_file);
        ciborium::ser::into_writer(self, &mut compress_writer)?;
        Ok(())
    }

    pub fn load(file_prefix: &str) -> Result<Self, Box<dyn Error>> {
        let filename = format!("{}.skm", file_prefix);
        log::info!("Loading sketch metadata from {filename}");
        let skm_file = BufReader::new(File::open(filename)?);
        let decompress_reader = snap::read::FrameDecoder::new(skm_file);
        let skm_obj: Self = ciborium::de::from_reader(decompress_reader)?;
        Ok(skm_obj)
    }

    pub fn number_samples_loaded(&self) -> usize {
        match &self.block_reindex {
            Some(block_map) => block_map.len(),
            None => self.sketch_metadata.len(),
        }
    }

    pub fn get_k_idx(&self, k: usize) -> Option<usize> {
        self.kmer_lengths
            .iter()
            .enumerate()
            .find_map(|(idx, val)| if *val == k { Some(idx) } else { None })
    }

    pub fn sketch_name(&self, index: usize) -> &str {
        match &self.block_reindex {
            Some(block_map) => self.sketch_metadata[block_map[index]].name(),
            None => self.sketch_metadata[index].name(),
        }
    }

    pub fn read_sketch_data(&mut self, file_prefix: &str) {
        let filename = format!("{}.skd", file_prefix);
        log::debug!(
            "bin_stride:{} kmer_stride:{} sample_stride:{}",
            self.bin_stride,
            self.kmer_stride,
            self.sample_stride
        );
        self.sketch_bins =
            SketchArrayFile::read_all(&filename, self.sample_stride * self.sketch_metadata.len());
    }

    pub fn read_sketch_data_block(&mut self, file_prefix: &str, names: &[String]) {
        // Find the given names in the sketch metadata
        let mut block_reindex = Vec::with_capacity(names.len());
        let mut read_indices = Vec::with_capacity(names.len());
        for name in names {
            if let Some(sketch_idx) = self.name_map.get(name) {
                read_indices.push(self.sketch_metadata[*sketch_idx].get_index());
                block_reindex.push(*sketch_idx);
            } else {
                panic!("Could not find requested sample {name} in sketch metadata");
            }
        }
        self.block_reindex = Some(block_reindex);

        let filename = format!("{}.skd", file_prefix);
        self.sketch_bins =
            SketchArrayFile::read_batch(&filename, &read_indices, self.sample_stride);
    }

    pub fn get_sketch_slice(&self, sketch_idx: usize, k_idx: usize) -> &[u64] {
        debug_assert!(sketch_idx < self.sketch_metadata.len());
        let s1_offset = sketch_idx * self.sample_stride + k_idx * self.kmer_stride;
        let s1_slice =
            &self.sketch_bins[s1_offset..(s1_offset + (self.sketch_size * BBITS) as usize)];
        log::trace!(
            "s1_start:{s1_offset} s1_end:{}",
            s1_offset + s1_slice.len(),
        );
        s1_slice
    }
}

impl fmt::Debug for MultiSketch {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "sketch_version={}\nsequence_type={:?}\nsketch_size={}\nn_samples={}\nkmers={:?}",
            self.sketch_version,
            self.hash_type,
            self.sketch_size * u64::BITS as u64,
            self.sketch_metadata.len(),
            self.kmer_lengths,
        )
    }
}

impl fmt::Display for MultiSketch {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "Name\tSequence length\tBase frequencies\tMissing/ambig bases\tFrom reads\tSingle strand\tDensified")?;
        for sketch in &self.sketch_metadata {
            write!(f, "{sketch}")?;
        }
        Ok(())
    }
}
