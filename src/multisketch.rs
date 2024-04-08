use std::error::Error;
use std::fmt;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::mem;
use std::process::Output;
use std::rc::Rc;

use hashbrown::{HashMap, HashSet};
use serde::{Deserialize, Serialize};

use crate::sketch::{Sketch, BBITS};
use crate::sketch_datafile::SketchArrayFile;

// TODO might be better to have HashMap<String, usize> of name maping
// sketch_metadata is then just Vec<Sketch>
#[derive(Serialize, Deserialize)]
pub struct MultiSketch<'a> {
    sketch_size: u64,
    kmer_lengths: Vec<usize>,
    sketch_metadata: Vec<Sketch>,
    name_map: HashMap<String, usize>, // TODO may not actually need this
    #[serde(skip)]
    block_metadata: Option<Vec<&'a Sketch>>, // TODO: may need to become Arc for multithreading? Also could just list the indices and use indirection into sketch_metadata
    #[serde(skip)]
    sketch_bins: Vec<u64>,
    bin_stride: usize,
    kmer_stride: usize,
    sample_stride: usize,
    sketch_version: String,
}

impl<'a> MultiSketch<'a> {
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
            block_metadata: None,
            sketch_bins: Vec::new(),
            bin_stride: 1,
            kmer_stride,
            sample_stride: kmer_stride * kmer_lengths.len(),
            sketch_version: env!("CARGO_PKG_VERSION").to_string(),
        }
    }

    /// Saves the metadata
    pub fn save_metadata(&self, file_prefix: &str) -> Result<(), Box<dyn Error>> {
        let filename = format!("{}.skm", file_prefix);
        let serial_file = BufWriter::new(File::create(filename)?);
        let mut compress_writer = snap::write::FrameEncoder::new(serial_file);
        ciborium::ser::into_writer(self, &mut compress_writer)?;
        Ok(())
    }

    pub fn load(file_prefix: &str) -> Result<Self, Box<dyn Error>> {
        let filename = format!("{}.skm", file_prefix);
        let skm_file = BufReader::new(File::open(filename)?);
        let decompress_reader = snap::read::FrameDecoder::new(skm_file);
        let skm_obj: Self = ciborium::de::from_reader(decompress_reader)?;
        Ok(skm_obj)
    }

    pub fn number_samples_loaded(&self) -> usize {
        match &self.block_metadata {
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
        match &self.block_metadata {
            Some(block_map) => block_map[index].name(),
            None => self.sketch_metadata[index].name(),
        }
    }

    pub fn read_sketch_data(&mut self, file_prefix: &str) {
        let filename = format!("{}.skd", file_prefix);
        self.sketch_bins =
            SketchArrayFile::read_all(&filename, self.sample_stride * self.sketch_metadata.len());
    }

    pub fn read_sketch_data_block(&'a mut self, file_prefix: &str, names: &[String]) {
        // Find the given names in the sketch metadata
        let mut block_metadata = Vec::with_capacity(names.len());
        let mut read_indices = Vec::with_capacity(names.len());
        for name in names {
            if let Some(sketch_idx) = self.name_map.get(name) {
                let block_sketch = &self.sketch_metadata[*sketch_idx];
                read_indices.push(block_sketch.get_index());
                block_metadata.push(block_sketch);
            } else {
                panic!("Could not find {name} in sketch metadata");
            }
        }
        let filename = format!("{}.skd", file_prefix);
        self.sketch_bins =
            SketchArrayFile::read_batch(&filename, &read_indices, self.sample_stride);
        self.block_metadata = Some(block_metadata);
    }

    // TODO: would be fun to try putting signs in a roaringbitmap then just
    // doing intersection_len() as this function
    pub fn jaccard_dist(&self, sketch1_idx: usize, sketch2_idx: usize, k_idx: usize) -> f64 {
        let s1_offset = sketch1_idx * self.sample_stride + k_idx * self.kmer_stride;
        let s2_offset = sketch2_idx * self.sample_stride + k_idx * self.kmer_stride;
        let s1_slice =
            &self.sketch_bins[s1_offset..(s1_offset + (self.sketch_size * BBITS) as usize)];
        let s2_slice =
            &self.sketch_bins[s2_offset..(s2_offset + (self.sketch_size * BBITS) as usize)];
        let unionsize = (u64::BITS as u64 * self.sketch_size) as f64;
        let mut samebits: u32 = 0;
        for i in 0..self.sketch_size {
            let mut bits: u64 = !0;
            for j in 0..BBITS {
                bits &= !(s1_slice[(i * BBITS + j) as usize] ^ s2_slice[(i * BBITS + j) as usize]);
            }
            samebits += bits.count_ones();
        }
        let maxnbits = self.sketch_size as u32 * u64::BITS;
        let expected_samebits = maxnbits >> BBITS;
        if expected_samebits != 0 {
            samebits as f64
        } else {
            let diff = samebits.saturating_sub(expected_samebits);
            let intersize = (diff * maxnbits) as f64 / (maxnbits - expected_samebits) as f64;
            intersize / unionsize
        }
    }

    pub fn core_acc_dist(&self, sketch1_idx: usize, sketch2_idx: usize) -> (f64, f64) {
        if self.kmer_lengths.len() < 2 {
            panic!("Need at least two k-mer lengths to calculate core/accessory distances");
        }
        let (mut xsum, mut ysum, mut xysum, mut xsquaresum, mut ysquaresum, mut n) =
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        let tolerance = 5.0 / (self.sketch_size as f64);
        for (k_idx, k) in self.kmer_lengths.iter().enumerate() {
            let y = self.jaccard_dist(sketch1_idx, sketch2_idx, k_idx).exp();
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
}

impl<'a> fmt::Debug for MultiSketch<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "sketch_version={}\nsketch_size={}\nn_samples={}\nkmers={:?}",
            self.sketch_version,
            self.sketch_size * u64::BITS as u64,
            self.sketch_metadata.len(),
            self.kmer_lengths,
        )
    }
}

impl<'a> fmt::Display for MultiSketch<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        todo!()
    }
}
