use std::{collections::HashMap, process::Output};
use std::mem;
use std::fs::File;
use std::io::{BufReader, BufWriter, Seek, SeekFrom, Write};
use std::error::Error;

use serde::{Deserialize, Serialize};

use crate::sketch::Sketch;

#[derive(Serialize, Deserialize, Debug)]
pub struct MultiSketch {
    sketch_size: u64,
    kmer_lengths: Vec<usize>,
    sketch_metadata: HashMap<String, Sketch>,
    #[serde(skip)]
    sketch_data: Vec<u8>, // reading only
    #[serde(skip)]
    serial_writer: Option<BufWriter<File>>,
    #[serde(skip)]
    metadata_filename: String,
    bin_stride: usize,
    kmer_stride: usize,
    sample_stride: usize,
}

impl MultiSketch {
    pub fn new(n_samples: usize, sketch_size: u64, kmer_lengths: &[usize], output_prefix: &str) -> Self {
        let bin_stride = 1;
        let kmer_stride = sketch_size as usize;
        let sample_stride = kmer_stride * kmer_lengths.len();

        // TODO need to compare speed vs memmapping
        let data_filename = format!("{output_prefix}.skd");
        let serial_file = File::create(data_filename).expect("Couldn't write to {data_filename}");
        serial_file.set_len(n_samples as u64 * kmer_lengths.len() as u64 * sketch_size * (u64::BITS / u8::BITS) as u64).expect("Couldn't set size of {data_filename}");
        let serial_writer = Some(BufWriter::new(serial_file));
        let metadata_filename = format!("{output_prefix}.skm");

        Self {
            sketch_size,
            kmer_lengths: kmer_lengths.to_vec(),
            sketch_metadata: HashMap::with_capacity(n_samples),
            sketch_data: Vec::new(),
            serial_writer,
            metadata_filename,
            bin_stride,
            kmer_stride,
            sample_stride,
        }
    }

    /// Writes the sketch to file
    pub fn add_sketch(&mut self, name: &str, sketch: &mut Sketch) {
        let sketch_idx = sketch.index();
        let sketch_map = sketch.get_usigs();
        let file_writer = self.serial_writer.as_mut().unwrap();
        for kmer_idx in 0..self.kmer_lengths.len() {
            let kmer = self.kmer_lengths[kmer_idx];
            let offset = sketch_idx * self.sample_stride + kmer_idx * self.kmer_stride;
            match sketch_map.get(&kmer) {
                Some(usigs) => Self::write_sketch_data(file_writer, usigs, offset).unwrap(),
                None => panic!("Could not find k-mer {kmer} in sketch {name}")
            }
        }
        self.sketch_metadata.insert(name.to_string(), mem::take(sketch));
    }

    fn write_sketch_data<W: Write + Seek>(writer: &mut W, usigs: &[u64], offset: usize) -> Result<(), Box<dyn Error>>{
        writer.seek(SeekFrom::Start(offset as u64))?;
        for bin_val in usigs {
            writer.write(&bin_val.to_le_bytes())?;
        }
        Ok(())
    }

    /// Saves the metadata
    pub fn save(&self) -> Result<(), Box<dyn Error>> {
        let serial_file = BufWriter::new(File::create(self.metadata_filename.as_str())?);
        let mut compress_writer = snap::write::FrameEncoder::new(serial_file);
        ciborium::ser::into_writer(self, &mut compress_writer)?;
        Ok(())
    }

    pub fn load(file_prefix: &str, names: &[&str]) -> Self {
        // Read the serde part
        // Find the given names in the sketch metadata
        // Discard any not in the name list
        // Open sketch data via memmap2
        // Call private function to fill each sketche's usigs from this file
        todo!()
    }

    pub fn to_map(&mut self, names: &[&str]) -> &HashMap<String, Sketch> {
        &self.sketch_metadata
    }

    fn memmap_file(&mut self, file: &str) {
        todo!()
    }

    fn load_sketch(&mut self, sketch_index: usize) {
        todo!()
    }
}
