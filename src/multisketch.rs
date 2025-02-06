//! The class to support .skm/.skd reading and writing, containing multiple [`Sketch`] objects
use anyhow::anyhow;
use anyhow::bail;
use anyhow::Error;
// use thiserror::Error;
use core::panic;
use std::fmt;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::mem;

use hashbrown::HashMap;
use serde::{Deserialize, Serialize};

use crate::hashing::HashType;
use crate::sketch::{Sketch, BBITS};
use crate::sketch_datafile::SketchArrayFile;

use std::collections::HashSet;

use rusqlite::{params, Connection, Result};

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
    pub fn save_metadata(&self, file_prefix: &str) -> Result<(), Error> {
        let filename = format!("{}.db", file_prefix);
        log::info!("Saving sketch metadata to {filename}");
        // Creates database
        let conn = Connection::open(filename)?;
        // Initialise table in database
        conn.execute(
            "CREATE TABLE sketch_metadata (
                            id INTEGER PRIMARY KEY,
                            name TEXT NOT NULL,
                            length INTEGER
                           )",
            (),
        )?;
        // Iterate over metadata and add to database
        for (index, metadata) in self.sketch_metadata.iter().enumerate() {
            conn.execute(
                "INSERT INTO sketch_metadata (id, name, length) VALUES (?1, ?2, ?3)",
                (index, metadata.name(), metadata.seq_length()),
            )?;
        }
        Ok(())
    }

    pub fn query_metadata<T>(path: String, args: T) -> Result<Option<Vec<usize>>, Error> {
        let conn = Connection::open(path)?;
        conn.execute("SELECT id FROM sketch_metadata WHERE", ())?;
        todo!()
    }

    pub fn load(file_prefix: &str) -> Result<Self, Error> {
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

    pub fn is_compatible_with(&self, sketch2: &Self) -> bool {
        self.kmer_lengths() == sketch2.kmer_lengths()
            && self.sketch_size == sketch2.sketch_size
            && self.get_hash_type() == sketch2.get_hash_type()
    }

    pub fn append_compatibility(&self, name_vec: &[(String, String, Option<String>)]) -> bool {
        let mut compatibility = true;
        let mut duplicate_list = Vec::new();
        for (id, _, _) in name_vec.iter() {
            if self.name_map.contains_key(id) {
                duplicate_list.push(id);
                compatibility = false;
            }
        }

        if !duplicate_list.is_empty() {
            println!("Duplicates found: {:?}", duplicate_list);
        }

        compatibility
    }

    pub fn merge_sketches(&mut self, sketch2: &Self) -> &mut Self {
        // First metadata
        let offset = self.sketch_metadata.len();
        for sketch in sketch2.sketch_metadata.iter() {
            if self.name_map.contains_key(sketch.name()) {
                panic!(
                    "{} appears in both databases. Cannot merge sketches.",
                    sketch.name()
                );
            } else {
                let mut temp_sketch: Sketch = sketch.clone();
                let new_index = temp_sketch.get_index() + offset;
                temp_sketch.set_index(new_index);
                self.name_map
                    .insert(temp_sketch.name().to_string(), new_index);
                self.sketch_metadata.push(temp_sketch);
            }
        }

        self
    }
    pub fn remove_metadata(
        &mut self,
        output_file_name: &str,
        genome_ids_to_remove: &[String],
    ) -> anyhow::Result<()> {
        let mut new_sketch_metadata: Vec<Sketch> =
            Vec::with_capacity(self.sketch_metadata.len() - genome_ids_to_remove.len());
        let mut removed_samples = Vec::new();

        for sketch in &self.sketch_metadata {
            if !genome_ids_to_remove.contains(&(*sketch.name()).to_string()) {
                new_sketch_metadata.push(sketch.clone());
            } else {
                removed_samples.push(sketch.name());
            }
        }

        let set1: HashSet<&str> = removed_samples.iter().map(AsRef::as_ref).collect();
        let set2: HashSet<&str> = genome_ids_to_remove.iter().map(AsRef::as_ref).collect();
        let missing: Vec<&&str> = set2.difference(&set1).collect();
        if !missing.is_empty() {
            bail!(
                "The following samples have not been found in the database: {:?}",
                missing
            );
        }

        self.sketch_metadata = new_sketch_metadata;
        self.save_metadata(output_file_name)?;
        Ok(())
    }

    pub fn remove_genomes(
        &mut self,
        input_prefix: &str,
        output_file: &str,
        genome_ids_to_remove: &[String],
    ) -> anyhow::Result<()> {
        let mut positions_to_remove = Vec::new();
        let mut missing_ids = Vec::new();

        for id in genome_ids_to_remove {
            if let Some(&position) = self.name_map.get(id) {
                positions_to_remove.push(position);
            } else {
                missing_ids.push(id);
            }
        }

        if !missing_ids.is_empty() {
            bail!("The following genome IDs were not found: {:?}", missing_ids);
        }

        // Create a list of indices to keep
        let indices_to_keep: Vec<usize> = (0..self.sketch_metadata.len())
            .filter(|&idx| !positions_to_remove.contains(&idx))
            .collect();

        let input_filename = format!("{}.skd", input_prefix);
        let output_filename = format!("{}.skd", output_file);
        if let Err(e) = SketchArrayFile::write_batch(
            &input_filename,
            &output_filename,
            &indices_to_keep,
            self.sample_stride,
        ) {
            return Err(anyhow!("Error during batch write: {}", e));
        }

        Ok(())
    }

    // This function is called when sketches are merged, not when they are
    // first sketched (this is handled by sketch::sketch_files())
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

// This is only used in the tests
// Ignores name_map
impl PartialEq for MultiSketch {
    fn eq(&self, other: &Self) -> bool {
        let mut metadata_match = true;
        if self.sketch_metadata.len() != other.sketch_metadata.len() {
            metadata_match = false;
            eprintln!(
                "Sketch metadata lengths is mismatching. Self: {}, Other: {}",
                self.sketch_metadata.len(),
                other.sketch_metadata.len()
            );
        }
        for (self_sketch, other_sketch) in self
            .sketch_metadata
            .iter()
            .zip(other.sketch_metadata.iter())
        {
            if self_sketch != other_sketch {
                metadata_match = false;
                eprintln!(
                    "Sketches mismatching. Self: {}, Other: {}",
                    self_sketch, other_sketch
                );
                break;
            }
        }
        if self.sketch_size != other.sketch_size {
            metadata_match = false;
            eprintln!(
                "Sketch sizes are mismatching. Self: {}, Other: {}",
                self.sketch_size, other.sketch_size
            );
        }

        if self.kmer_lengths != other.kmer_lengths {
            metadata_match = false;
            eprintln!(
                "Kmer lengths are mismatching. Self: {:?}, Other: {:?}",
                self.kmer_lengths, other.kmer_lengths
            );
        }

        if self.sketch_bins != other.sketch_bins {
            metadata_match = false;
            eprintln!(
                "Sketch bins are mismatching. Self: {:?}, Other: {:?}",
                self.sketch_bins, other.sketch_bins
            );
        }

        if self.bin_stride != other.bin_stride {
            metadata_match = false;
            eprintln!(
                "Bin strides are mismatching. Self: {}, Other: {}",
                self.bin_stride, other.bin_stride
            );
        }

        if self.kmer_stride != other.kmer_stride {
            metadata_match = false;
            eprintln!(
                "Kmer strides are mismatching. Self: {}, Other: {}",
                self.kmer_stride, other.kmer_stride
            );
        }

        if self.sample_stride != other.sample_stride {
            metadata_match = false;
            eprintln!(
                "Sample strides are mismatching. Self: {}, Other: {}",
                self.sample_stride, other.sample_stride
            );
        }

        if self.sketch_version != other.sketch_version {
            metadata_match = false;
            eprintln!(
                "Sketch versions are mismatching. Self: {:?}, Other: {:?}",
                self.sketch_version, other.sketch_version
            );
        }

        if self.hash_type != other.hash_type {
            metadata_match = false;
            eprintln!(
                "Hash types are mismatching. Self: {:?}, Other: {:?}",
                self.hash_type, other.hash_type
            );
        }
        metadata_match
    }
}
