use core::panic;
use std::error::Error;
use std::fmt;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::mem;

use hashbrown::HashMap;
use serde::{Deserialize, Serialize};

use crate::hashing::HashType;
use crate::sketch::{self, Sketch, BBITS};
use crate::sketch_datafile::SketchArrayFile;

use std::collections::HashSet;

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
        sketches: &mut Vec<Sketch>,
        sketch_size: u64,
        kmer_lengths: &[usize],
        hash_type: HashType,
    ) -> Self {
        let mut name_map = HashMap::with_capacity(sketches.len());
        for sketch in sketches.iter() {
            name_map.insert(sketch.name().to_string(), sketch.get_index());
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
            hash_type,
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

    pub fn kmer_lengths(&self) -> &[usize] {
        &self.kmer_lengths
    }

    pub fn get_hash_type(&self) -> &HashType {
        &self.hash_type
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
        log::trace!("s1_start:{s1_offset} s1_end:{}", s1_offset + s1_slice.len(),);
        s1_slice
    }

    pub fn remove_sketches(&self, ids: &Vec<String>) {
        // TODO: remove sketch bins which belong to the duplicate ids
        todo!();
    }

    pub fn is_compatible_with(&self, sketch2: &Self) -> bool {
        self.kmer_lengths() == sketch2.kmer_lengths()
            && self.sketch_size == sketch2.sketch_size
            && self.get_hash_type() == sketch2.get_hash_type()
    }

    pub fn merge_sketches(&self, sketch2: &Self) -> Self {
        let mut sketch1_names: HashSet<String> = HashSet::new();
        for sketch in self.sketch_metadata.iter() {
            let name = sketch.name().to_string();
            println!("{}", name);
            sketch1_names.insert(name);
        }

        // First metadata
        let mut merged_metadata = self.sketch_metadata.clone();
        for sketch in sketch2.sketch_metadata.iter() {
            if sketch1_names.contains(sketch.name()) {
                let sample_name = sketch.name();
                panic!("{sample_name} seems to appear in both databases. Cannot merge sketches.");
            } else {
                merged_metadata.push(sketch.clone());
            }
        }
        // then merge sketches, create new multisketch instance
        let mut merged_sketch = Self::new(
            &mut merged_metadata,
            self.sketch_size,
            &self.kmer_lengths,
            self.hash_type.clone(),
        );
        // Merge actual sketche infos
        merged_sketch.sketch_bins = self.sketch_bins.clone();
        for bin in &sketch2.sketch_bins {
            merged_sketch.sketch_bins.push(*bin);
        }
        // Update infos
        merged_sketch.bin_stride = self.bin_stride;
        merged_sketch.kmer_stride = self.kmer_stride;
        merged_sketch.sample_stride = self.sample_stride;

        let mut merged_name_map = self.name_map.clone();
        merged_name_map.reserve(sketch2.sketch_metadata.len());

        // update hashmap
        for sketch in sketch2.sketch_metadata.iter() {
            if !merged_name_map.contains_key(sketch.name()) {
                let new_index = sketch.get_index() + self.sketch_metadata.len();
                merged_name_map.insert(sketch.name().to_string(), new_index);
            }
        }

        // update index
        for sketch in merged_sketch.sketch_metadata.iter_mut() {
            let new_index = merged_name_map[sketch.name()];
            sketch.set_index(new_index);
        }
        merged_sketch.name_map = merged_name_map;

        merged_sketch
    }

    pub fn save_multi_sketches(&self, file_prefix: &str) -> Result<(), Box<dyn Error>> {
        let data_filename = format!("{file_prefix}.skd");
        let mut file = std::fs::File::create(&data_filename)?;
        SketchArrayFile::write_sketch_data(&mut file, &self.sketch_bins)?;
        Ok(())
    }

    pub fn strip_sketch_extension(file_name: &str) -> &str {
        if file_name.ends_with(".skm") || file_name.ends_with(".skd") {
            &file_name[..file_name.len() - 4]
        } else {
            file_name
        }
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

impl PartialEq for MultiSketch {
    fn eq(&self, other: &Self) -> bool {
        if self.sketch_metadata.len() != other.sketch_metadata.len() {
            
            panic!(
                "sketch_metadata length mismatch: self {} != other {}",
                self.sketch_metadata.len(),
                other.sketch_metadata.len()
            );
        }
        for (i, (self_sketch, other_sketch)) in self
            .sketch_metadata
            .iter()
            .zip(other.sketch_metadata.iter())
            .enumerate()
        {
            if self_sketch != other_sketch {
                panic!(
                    "Mismatch in sketch_metadata at index {}: {:?} != {:?}",
                    i, self_sketch, other_sketch
                );
            }
        }
        macro_rules! check_field {
            ($field:ident) => {
                if self.$field != other.$field {
                    panic!(
                        "{} mismatch: {:?} != {:?}",
                        stringify!($field),
                        self.$field,
                        other.$field
                    );
                }
            };
        }
        check_field!(sketch_size);
        check_field!(kmer_lengths);
        check_field!(name_map);
        check_field!(block_reindex);
        check_field!(sketch_bins);
        check_field!(bin_stride);
        check_field!(kmer_stride);
        check_field!(sample_stride);
        check_field!(sketch_version);
        check_field!(hash_type);
        true
    }
}
