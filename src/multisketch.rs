use std::sync::{Arc};
use std::{collections::HashMap, process::Output};
use std::mem;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::error::Error;

use serde::{Deserialize, Serialize};

use crate::sketch::{Sketch, BBITS};
use crate::sketch_datafile::SketchArrayFile;

#[derive(Serialize, Deserialize, Debug)]
pub struct MultiSketch {
    sketch_size: u64,
    kmer_lengths: Vec<usize>,
    sketch_metadata: HashMap<String, Sketch>,
    #[serde(skip)]
    sketch_data_buffer: Vec<u8>, // reading only
    #[serde(skip)]
    flat_sketch_array: Vec<u64>,
    flat_sketch_array_size: usize,
    #[serde(skip)]
    file_prefix: String,
}

impl MultiSketch {
    pub fn new(names: &[String], sketches: &mut [Sketch], sketch_size: u64, kmer_lengths: &[usize], file_prefix: &str) -> Self {
        let mut sketch_metadata: HashMap<String, Sketch> = HashMap::with_capacity(sketches.len());
        for (name, sketch) in names.iter().zip(sketches) {
            sketch_metadata.insert(name.to_string(), mem::take(sketch));
        }

        Self {
            sketch_size,
            kmer_lengths: kmer_lengths.to_vec(),
            sketch_metadata,
            sketch_data_buffer: Vec::new(),
            flat_sketch_array: Vec::new(),
            flat_sketch_array_size: (sketch_size * BBITS) as usize * kmer_lengths.len() * names.len(),
            file_prefix: file_prefix.to_string(),
        }
    }

    /// Saves the metadata
    pub fn save_metadata(&self) -> Result<(), Box<dyn Error>> {
        let filename = format!("{}.skm", self.file_prefix);
        let serial_file = BufWriter::new(File::create(filename)?);
        let mut compress_writer = snap::write::FrameEncoder::new(serial_file);
        ciborium::ser::into_writer(self, &mut compress_writer)?;
        Ok(())
    }

    pub fn load(file_prefix: &str) -> Result<Self, Box<dyn Error>> {
        // Read the serde part
        let filename = format!("{}.skm", file_prefix);
        let skm_file = BufReader::new(File::open(filename)?);
        let decompress_reader = snap::read::FrameDecoder::new(skm_file);
        let mut skm_obj: Self = ciborium::de::from_reader(decompress_reader)?;
        skm_obj.file_prefix = file_prefix.to_string();
        Ok(skm_obj)
    }

    pub fn read_all(&mut self) {
        // Just stream the whole file and convert to u64 vec
        let filename = format!("{}.skd", self.file_prefix);
        let mut reader = BufReader::new(File::open(filename).expect("Could not read from {filename}"));
        let mut buffer = [0u8; mem::size_of::<u64>()];
        self.flat_sketch_array.reserve_exact(self.flat_sketch_array_size);
        while let Ok(_read) = reader.read_exact(&mut buffer) {
            self.flat_sketch_array.push(u64::from_le_bytes(buffer));
        }
    }

    pub fn read_block(&mut self, names: &[String]) {
        // Find the given names in the sketch metadata
        // Discard any not in the name list (i.e. delete from the hashmap)
        // Open sketch data via memmap2
        // Copy to a u64 with new indices (use Self::load_sketch)
        // update index in hashmap entry
    }

    pub fn to_map(&mut self) -> &HashMap<String, Sketch> {
        &self.sketch_metadata
    }

    fn memmap_file(&mut self, file: &str) {
        todo!()
    }

    fn load_sketch(&mut self, sketch_index: usize) {
        todo!()
    }

    // TODO: would be fun to try putting signs in a roaringbitmap then just
    // doing intersection_len() as this function
    pub fn jaccard_dist(&self, sketch1_idx: usize, sketch2_idx: usize, k: usize) -> f64 {
        todo!();
        // TODO Need to load appropriate blocks
        /*
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
        */
    }

    pub fn core_acc_dist(&self, sketch1_idx: usize, sketch2_idx: usize) -> (f64, f64) {
        if self.kmer_lengths.len() < 2 {
            panic!("Need at least two k-mer lengths to calculate core/accessory distances");
        }
        let mut xsum = 0.0;
        let mut ysum = 0.0;
        let mut xysum = 0.0;
        let mut xsquaresum = 0.0;
        let mut ysquaresum = 0.0;
        let mut n = 0.0;
        let tolerance = 5.0 / (self.sketch_size as f64);
        for k in &self.kmer_lengths {
            let y = self.jaccard_dist(sketch1_idx, sketch2_idx, *k).exp();
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
