//! The class to support .ski creation, reading and writing, containing an inverted
//! index of multiple sketches
#[cfg(target_arch = "wasm32")]
use std::fmt;
#[cfg(not(target_arch = "wasm32"))]
use std::sync::mpsc;
#[cfg(not(target_arch = "wasm32"))]
use std::{cmp, fmt};
#[cfg(not(target_arch = "wasm32"))]
use std::{fmt::Write as FmtWrite, io::Write as IoWrite};

#[cfg(not(target_arch = "wasm32"))]
extern crate needletail;

#[cfg(target_arch = "wasm32")]
use hashbrown::HashMap;
#[cfg(not(target_arch = "wasm32"))]
use hashbrown::{HashMap, HashSet};
#[cfg(not(target_arch = "wasm32"))]
use indicatif::ParallelProgressIterator;
use indicatif::ProgressIterator;
use rayon::prelude::*;
use roaring::{RoaringBitmap, RoaringTreemap};
use serde::{Deserialize, Serialize};

use super::hashing::{
    bloom_filter::KmerFilter, nthash_iterator::NtHashIterator, HashType, RollHash,
};
use crate::distances::distance_matrix::square_to_condensed;
use crate::sketch::multisketch::{parse_version, MIN_SKETCH_VERSION};
#[cfg(target_arch = "wasm32")]
use crate::sketch::Sketch;
#[cfg(not(target_arch = "wasm32"))]
use crate::sketch::{sketch_datafile::SketchArrayWriter, Sketch};
use crate::utils::get_progress_bar;
#[cfg(not(target_arch = "wasm32"))]
use crate::{
    cli::InvertedQueryType,
    io::{normalise_metadata, warn_missing_metadata, InputFastx, SampleMetadata},
};
#[cfg(not(target_arch = "wasm32"))]
use anyhow::Context;
use anyhow::{bail, Error};
use std::fs::File;
use std::io::{BufReader, BufWriter};

#[cfg(target_arch = "wasm32")]
use crate::logw;
#[cfg(target_arch = "wasm32")]
use wasm_bindgen_file_reader::WebSysFile;

type InvSketches = (Vec<Vec<u16>>, Vec<String>);

#[cfg(not(target_arch = "wasm32"))]
fn append_sample(output: &mut String, name: &str, metadata: Option<&str>) {
    output.push_str(name);
    if let Some(metadata) = metadata {
        output.push(';');
        output.push_str(metadata);
    }
}

/// Fails if a loaded .ski index was built with a sketchlib version older
/// than the current minimum sketch format version.
fn check_version(sketch_version: &str) -> Result<(), Error> {
    let version_ok =
        parse_version(sketch_version).is_some_and(|version| version >= MIN_SKETCH_VERSION);
    if !version_ok {
        let found_version = if sketch_version.is_empty() {
            "<unknown>"
        } else {
            sketch_version
        };
        log::error!(
            "Inverted index was created with sketchlib v{found_version}, which is older than the minimum supported v{}.{}.{}. Please re-sketch with the current version.",
            MIN_SKETCH_VERSION.0,
            MIN_SKETCH_VERSION.1,
            MIN_SKETCH_VERSION.2
        );
        bail!("Incompatible sketch file version");
    }
    Ok(())
}

/// An inverted index and associated metadata
#[derive(Serialize, Deserialize, Default, Clone, PartialEq)]
pub struct Inverted {
    index: Vec<HashMap<u16, RoaringBitmap>>,
    n_samples: usize,
    sample_names: Vec<String>,
    metadata: Option<Vec<String>>,
    labels: Option<Vec<String>>,
    kmer_size: usize,
    sketch_version: String,
    rc: bool,
    hash_type: HashType,
}

impl Inverted {
    #[cfg(not(target_arch = "wasm32"))]
    /// Sketch files without transposing bins, and invert the index. Files
    /// may be reordered via [`crate::io::reorder_input_files`]
    ///
    /// Optionally write untransposed samples to an skq file with `write_skq`.
    pub fn new(
        input_files: &[InputFastx],
        write_skq: Option<String>,
        file_order: &[usize],
        k: usize,
        sketch_size: u64,
        seq_type: &HashType,
        rc: bool,
        min_count: u16,
        min_qual: u8,
        quiet: bool,
        metadata: &Option<Vec<String>>,
        labels: &Option<Vec<String>>,
    ) -> Self {
        log::info!("Creating sketches");
        let (sketches, names) = Self::sketch_files_inverted(
            input_files,
            file_order,
            k,
            sketch_size,
            seq_type,
            rc,
            min_count,
            min_qual,
            quiet,
        );
        if let Some(skq_file) = write_skq {
            log::info!("Writing bins for use with precluster as {skq_file}");
            let (bin_stride, kmer_stride, sample_stride) = (1, 1, sketch_size as usize);
            let mut skq_writer =
                SketchArrayWriter::new(&skq_file, bin_stride, kmer_stride, sample_stride);
            for sketch in &sketches {
                skq_writer.write_sketch(sketch);
            }
        }
        log::info!("Inverting sketch order");
        Self {
            index: Self::build_inverted_index(&sketches, sketch_size),
            n_samples: names.len(),
            sample_names: names,
            metadata: metadata.clone(),
            labels: labels.clone(),
            kmer_size: k,
            sketch_version: env!("CARGO_PKG_VERSION").to_string(),
            rc,
            hash_type: seq_type.clone(),
        }
    }

    #[cfg(not(target_arch = "wasm32"))]
    /// Sketch query files in a compatible manner with the index,
    /// used when querying the index on the fly
    pub fn sketch_queries(
        &self,
        input_files: &[InputFastx],
        min_count: u16,
        min_qual: u8,
        quiet: bool,
    ) -> InvSketches {
        let file_order: Vec<usize> = (0..input_files.len()).collect();
        Self::sketch_files_inverted(
            input_files,
            &file_order,
            self.kmer_size,
            self.index.len() as u64,
            &self.hash_type,
            self.rc,
            min_count,
            min_qual,
            quiet,
        )
    }

    #[cfg(target_arch = "wasm32")]
    /// Sketch query files in a compatible manner with the index,
    /// used when querying the index on the fly
    pub fn sketch_queries(
        &self,
        input_files: (&web_sys::File, Option<&web_sys::File>),
        min_count: u16,
        min_qual: u8,
        quiet: bool,
    ) -> InvSketches {
        let file_order: Vec<usize> = if input_files.1.is_some() {
            vec![0, 1]
        } else {
            vec![0]
        };

        Self::sketch_files_inverted(
            input_files,
            &file_order,
            self.kmer_size,
            self.index.len() as u64,
            &self.hash_type,
            self.rc,
            min_count,
            min_qual,
            quiet,
        )
    }

    /// Sample names in the index
    pub fn sample_names(&self) -> &Vec<String> {
        &self.sample_names
    }

    /// Number of samples in the index
    pub fn n_samples(&self) -> usize {
        self.sample_names.len()
    }

    /// Sample name at the given index
    pub fn sample_at(&self, idx: usize) -> &str {
        &self.sample_names[idx]
    }

    #[cfg(not(target_arch = "wasm32"))]
    pub(crate) fn metadata_at(&self, idx: usize) -> Option<&str> {
        self.metadata.as_ref().map(|values| values[idx].as_str())
    }

    // Only changes the loaded object, preserving the serialized index format.
    #[cfg(not(target_arch = "wasm32"))]
    pub(crate) fn prepare_metadata_for_output(&mut self, suppress: bool) -> Result<(), Error> {
        if suppress {
            self.metadata = None;
        }
        if let Some(values) = &mut self.metadata {
            if values.len() != self.sample_names.len() {
                bail!("Index metadata length does not match the number of samples");
            }
            let mut missing = 0;
            for (name, value) in self.sample_names.iter().zip(values) {
                missing += usize::from(
                    normalise_metadata(value)
                        .with_context(|| format!("Invalid index metadata for sample {name}"))?,
                );
            }
            warn_missing_metadata(missing, 0);
        }
        Ok(())
    }

    /// Compute and format query rows in parallel while the caller writes them.
    #[cfg(not(target_arch = "wasm32"))]
    pub(crate) fn query_stream<W: IoWrite + Send>(
        &self,
        writer: &mut W,
        queries: &[Vec<u16>],
        query_names: &[String],
        query_metadata: &SampleMetadata,
        query_type: &InvertedQueryType,
        threads: usize,
        quiet: bool,
    ) -> Result<(), Error> {
        let mut header = String::from("Query");
        if *query_type == InvertedQueryType::MatchCount {
            for (idx, name) in self.sample_names.iter().enumerate() {
                header.push('\t');
                append_sample(&mut header, name, self.metadata_at(idx));
            }
        } else {
            header.push_str("\tMatches");
        }
        header.push('\n');
        writer.write_all(header.as_bytes())?;

        let (tx, rx) = mpsc::sync_channel::<String>(threads.max(1) * 2);
        let progress_bar = get_progress_bar(queries.len(), false, quiet);
        rayon::scope(|s| -> Result<(), Error> {
            s.spawn(move |_| {
                queries
                    .par_iter()
                    .zip(query_names)
                    .progress_with(progress_bar)
                    .for_each_with(tx, |tx, (query, name)| {
                        let results = match query_type {
                            InvertedQueryType::MatchCount => {
                                self.query_against_inverted_index(query)
                            }
                            InvertedQueryType::AllBins => self.all_shared_bins(query),
                            InvertedQueryType::AnyBins => self.any_shared_bins(query),
                        };
                        let mut row = String::new();
                        append_sample(&mut row, name, query_metadata.get(name));
                        if *query_type == InvertedQueryType::MatchCount {
                            for count in results {
                                let _ = write!(row, "\t{count}");
                            }
                        } else {
                            for (match_idx, sample_idx) in results.into_iter().enumerate() {
                                row.push(if match_idx == 0 { '\t' } else { ',' });
                                let sample_idx = sample_idx as usize;
                                append_sample(
                                    &mut row,
                                    self.sample_at(sample_idx),
                                    self.metadata_at(sample_idx),
                                );
                            }
                        }
                        row.push('\n');
                        let _ = tx.send(row);
                    });
            });

            // Consuming rx here also disconnects blocked senders on a write error.
            for row in rx {
                writer
                    .write_all(row.as_bytes())
                    .context("Error writing streamed inverted query output")?;
            }
            writer
                .flush()
                .context("Error flushing inverted query output")?;
            Ok(())
        })
    }

    /// K-mer length used for the index
    pub fn kmer(&self) -> usize {
        self.kmer_size
    }

    /// Sketch size used for the index
    pub fn sketch_size(&self) -> usize {
        self.index.len()
    }

    /// Saves to `file_prefix.ski`, using MessagePack as the serialisation format
    pub fn save(&self, file_prefix: &str) -> Result<(), Error> {
        let filename = format!("{file_prefix}.ski");
        log::info!("Saving inverted index to {filename}");
        let serial_file = BufWriter::new(File::create(filename)?);
        let mut compress_writer = snap::write::FrameEncoder::new(serial_file);
        rmp_serde::encode::write(&mut compress_writer, self)?;
        Ok(())
    }

    #[cfg(not(target_arch = "wasm32"))]
    /// Loads from `file_prefix.ski`, using MessagePack as the serialisation format
    // NB MessagePack rather the CBOR uses here because of
    // https://github.com/enarx/ciborium/issues/96
    pub fn load(file_prefix: &str) -> Result<Self, Error> {
        let filename = format!("{file_prefix}.ski");
        log::info!("Loading inverted index from {filename}");
        let ski_file = BufReader::new(File::open(filename)?);
        let decompress_reader = snap::read::FrameDecoder::new(ski_file);
        let ski_obj: Self = rmp_serde::decode::from_read(decompress_reader)?;
        check_version(&ski_obj.sketch_version)?;
        Ok(ski_obj)
    }

    #[cfg(target_arch = "wasm32")]
    /// Loads from a web_sys file object. Designed for WebAssembly
    pub fn load(file: &web_sys::File) -> Result<Self, Error> {
        logw("Loading inverted index", Some("info"));

        let ski_file = BufReader::new(WebSysFile::new(file.clone()));
        let decompress_reader = snap::read::FrameDecoder::new(ski_file);
        let ski_obj: Self = rmp_serde::decode::from_read(decompress_reader)?;
        check_version(&ski_obj.sketch_version)?;
        Ok(ski_obj)
    }

    /// Query a single sample against the index, returning the number of matches
    /// of bins against all samples in the index
    pub fn query_against_inverted_index(&self, query_sigs: &[u16]) -> Vec<u32> {
        let mut match_counts = vec![0; self.sample_names.len()];

        for (bin_idx, query_bin_hash) in query_sigs.iter().enumerate() {
            if let Some(matching_samples) = self.index[bin_idx].get(query_bin_hash) {
                for sample_idx in matching_samples {
                    match_counts[sample_idx as usize] += 1;
                }
            }
        }
        match_counts
    }

    /// Return indexes of samples where all sketch bins are the same as a query
    pub fn all_shared_bins(&self, query_sigs: &[u16]) -> Vec<u32> {
        let mut matching_bits = RoaringBitmap::new();
        matching_bits.insert_range(0..self.sample_names.len() as u32);
        for (bin_idx, query_bin_hash) in query_sigs.iter().enumerate() {
            if let Some(matching_samples) = self.index[bin_idx].get(query_bin_hash) {
                matching_bits &= matching_samples;
            } else {
                matching_bits.clear();
                break;
            }
        }

        matching_bits.iter().collect()
    }

    /// Return indexes of samples where at least one sketch bins is the same as a query
    pub fn any_shared_bins(&self, query_sigs: &[u16]) -> Vec<u32> {
        let mut matching_bits = RoaringBitmap::new();
        for (bin_idx, query_bin_hash) in query_sigs.iter().enumerate() {
            if let Some(matching_samples) = self.index[bin_idx].get(query_bin_hash) {
                matching_bits |= matching_samples;
            }
        }

        matching_bits.iter().collect()
    }

    /// Get the pairwise comparison list used by prefiler
    pub fn any_shared_bin_list(&self, quiet: bool) -> RoaringTreemap {
        let percent = false;
        let progress_bar = get_progress_bar(self.index.len(), percent, quiet);
        self.index
            .iter()
            .progress_with(progress_bar)
            .map(|bin| {
                bin.par_values()
                    .map(|hash_pres| {
                        let mut pair_map_hash = RoaringTreemap::new();
                        let samples_together: Vec<u32> = hash_pres.iter().collect();
                        for (i, sample1_idx) in samples_together.iter().enumerate() {
                            for sample2_idx in samples_together.iter().skip(i + 1) {
                                pair_map_hash.insert(square_to_condensed(
                                    *sample1_idx as usize,
                                    *sample2_idx as usize,
                                    self.n_samples,
                                ) as u64);
                            }
                        }
                        pair_map_hash
                    })
                    .reduce(RoaringTreemap::new, |pair_map_hash_a, pair_map_hash_b| {
                        pair_map_hash_a | pair_map_hash_b
                    }) // Reduction from pair_map_hash produces pair_map_bin
            })
            .reduce(|pair_map_all, pair_map_bin| pair_map_all | pair_map_bin)
            .unwrap()
        // Reduction from pair_map_bin produces pair_map_all
    }

    #[cfg(not(target_arch = "wasm32"))]
    fn sketch_files_inverted(
        input_files: &[InputFastx],
        file_order: &[usize],
        k: usize,
        sketch_size: u64,
        seq_type: &HashType,
        rc: bool,
        min_count: u16,
        min_qual: u8,
        quiet: bool,
    ) -> InvSketches {
        let mut multientrysamples: HashSet<String> = HashSet::new();
        let mut differentsamples: HashSet<String> = HashSet::new();

        for i in input_files.iter() {
            if differentsamples.contains(&i.0) {
                multientrysamples.insert(i.0.clone());
            } else {
                differentsamples.insert(i.0.clone());
            }
        }
        let (tx, rx) = mpsc::channel();

        let percent = false;
        let progress_bar = get_progress_bar(input_files.len(), percent, quiet);
        rayon::scope(|s| {
            s.spawn(move |_| {
                input_files
                    .par_iter()
                    .zip(file_order)
                    .progress_with(progress_bar)
                    .map(|((name, fastxvec), genome_idx)| {
                        let mut hash_its: Vec<Box<dyn RollHash>> = match seq_type {
                            HashType::DNA => {
                            // Check if we're working with reads, and initalise the filter if so
                            let mut reader_peek = needletail::parse_fastx_file(fastxvec[0].clone())
                                .unwrap_or_else(|_| panic!("Invalid path/file: {}", fastxvec[0]));
                            let seq_peek = reader_peek
                                .next()
                                .expect("Invalid FASTA/Q record")
                                .expect("Invalid FASTA/Q record");
                            let mut reads = false;
                            if seq_peek.format() == needletail::parser::Format::Fastq {
                                reads = true;
                                if fastxvec.len() > 2 {
                                    panic!("Input files are reads, but there are more than two input files");
                                }
                            }

                            let mut records_readers = fastxvec.iter().map(|file| {
                                let reader = needletail::parse_fastx_file(file).unwrap_or_else(|_| panic!("Invalid path/file: {file}"));
                                crate::io::NeedletailIterator::new(reader)
                            }).collect::<Vec<crate::io::NeedletailIterator>>();

                                NtHashIterator::new(&mut records_readers, k, rc, min_qual, reads).unwrap()
                                .into_iter()
                                .map(|it| Box::new(it) as Box<dyn RollHash>)
                                .collect()
                            },
                            _ => unimplemented!("Inverted index only supported for DNA"),
                        };

                        if let Some(hash_it) = hash_its.first_mut() {
                            if hash_it.seq_len() == 0 {
                                panic!("Genome {genome_idx} has no valid sequence");
                            }

                            let mut read_filter = if hash_it.reads() {
                                let mut filter = KmerFilter::new(min_count);
                                filter.init();
                                Some(filter)
                            } else {
                                None
                            };

                            let signs = Sketch::get_signs_no_densify(
                                &mut **hash_it,
                                k,
                                &mut read_filter,
                                sketch_size,
                            );
                            // if densified {
                            //     log::trace!("{name} was densified");
                            // }
                            (*genome_idx, signs, name)
                        } else {
                            panic!("Empty hash iterator for {name}");
                        }
                    })
                    .for_each_with(tx, |tx, result| {
                        let _ = tx.send(result);
                    });
            });
        });

        let mut sketch_results: Vec<Vec<u16>> =
            vec![Vec::with_capacity(sketch_size as usize); differentsamples.len()];
        let mut indexes: HashSet<usize> = HashSet::with_capacity(multientrysamples.len());
        while let Ok((genome_idx, mut sketch, name)) = rx.recv() {
            if differentsamples.contains(name) {
                // Not yet written!
                if !multientrysamples.contains(name) {
                    // Densifying now!
                    Sketch::densify_bin(sketch.as_mut().unwrap());
                } else {
                    // We'll need to densify afterwards, let's save the index
                    indexes.insert(genome_idx);
                }
                sketch_results[genome_idx] =
                    sketch.as_ref().unwrap().iter().map(|h| *h as u16).collect();
                differentsamples.remove(name);
            } else {
                // already written! We have to merge
                for bin in 0..sketch_size {
                    let saved_sketch = &mut sketch_results[genome_idx][bin as usize];
                    *saved_sketch =
                        cmp::min(*saved_sketch, sketch.as_ref().unwrap()[bin as usize] as u16);
                }
            }
        }

        for pos in indexes.iter() {
            // TODO: this is clearly improvable, but I'm going to do it temporarily this way...
            let mut tmpvec: Vec<u64> = sketch_results[*pos].iter().map(|h| *h as u64).collect();
            Sketch::densify_bin(&mut tmpvec[..]);
            sketch_results[*pos] = tmpvec.iter().map(|h| *h as u16).collect();
        }

        // Sample names in the correct order
        // (clones names, but reference would be annoying here)
        let mut sample_names: Vec<String> = vec!["".to_string(); sketch_results.len()];
        file_order
            .iter()
            .zip(input_files)
            .for_each(|(idx, (name, _))| sample_names[*idx] = name.to_string());

        (sketch_results, sample_names)
    }

    #[cfg(target_arch = "wasm32")]
    fn sketch_files_inverted(
        input_files: (&web_sys::File, Option<&web_sys::File>),
        file_order: &[usize],
        k: usize,
        sketch_size: u64,
        seq_type: &HashType,
        rc: bool,
        min_count: u16,
        min_qual: u8,
        _quiet: bool,
    ) -> InvSketches {
        let mut hash_its: Vec<Box<dyn RollHash>> = match seq_type {
            HashType::DNA => NtHashIterator::new(input_files, k, rc, min_qual)
                .into_iter()
                .map(|it| Box::new(it) as Box<dyn RollHash>)
                .collect(),
            _ => unimplemented!("Inverted index only supported for DNA"),
        };

        if let Some(hash_it) = hash_its.first_mut() {
            if hash_it.seq_len() == 0 {
                panic!("Genome 0 has no valid sequence");
            }

            let mut read_filter = if hash_it.reads() {
                let mut filter = KmerFilter::new(min_count);
                filter.init();
                Some(filter)
            } else {
                None
            };

            let (signs, densified) =
                Sketch::get_signs(&mut **hash_it, k, &mut read_filter, sketch_size).unwrap();
            if densified {
                logw("The query was densified", Some("trace"));
            }

            let mut sketch_results = vec![Vec::new(); 1];

            sketch_results[0] = signs.iter().map(|h| *h as u16).collect();

            (sketch_results, vec!["".to_string(); file_order.len()])
        } else {
            panic!("Empty hash iterator for the query");
        }
    }

    #[cfg(not(target_arch = "wasm32"))]
    fn build_inverted_index(
        genome_sketches: &[Vec<u16>],
        sketch_size: u64,
    ) -> Vec<HashMap<u16, RoaringBitmap>> {
        // initialize inverted index structure
        let mut inverted_index: Vec<HashMap<u16, RoaringBitmap>> =
            vec![HashMap::new(); sketch_size as usize];

        // process each sketch
        // this could be parallelised over sketch bin, but probably not worth it
        // (NB doing par_iter on the below doesn't quite work as multiple mut borrows
        // of inverted_index needed)
        for (genome_idx, genome_signs) in genome_sketches.iter().enumerate() {
            for (i, hash) in genome_signs.iter().enumerate() {
                // add current genome to the inverted index at the current position
                inverted_index[i]
                    .entry(*hash)
                    .and_modify(|genome_list| {
                        genome_list.insert(genome_idx as u32);
                    })
                    .or_insert_with(|| {
                        let mut rb = RoaringBitmap::new();
                        rb.insert(genome_idx as u32);
                        rb
                    });
            }
        }
        inverted_index
    }

    /// Get the sample names, with the same order as in the index
    #[cfg(target_arch = "wasm32")]
    pub fn get_sample_names(&self) -> &Vec<String> {
        &self.sample_names
    }

    /// Get the metadata, with the same order as in the index
    pub fn get_metadata(&self) -> &Option<Vec<String>> {
        &self.metadata
    }

    /// Get the sample labels, with the same order as in the index
    #[cfg(target_arch = "wasm32")]
    pub fn get_sample_labels(&self) -> &Option<Vec<String>> {
        &self.labels
    }
}

impl fmt::Debug for Inverted {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "sketch_version={}\nsequence_type={:?}\nsketch_size={}\nn_samples={}\nkmer={}\nrc={}\ninverted=true\n",
            self.sketch_version,
            self.hash_type,
            self.index.len(),
            self.sample_names.len(),
            self.kmer_size,
            self.rc,
        )?;
        let mut sizes = Vec::new();
        for bin in self.index.iter() {
            sizes.push(bin.len());
        }
        write!(
            f,
            "max_hashes_per_bin={}\nmin_hashes_per_bin={}\navg_hashes_per_bin={}",
            sizes.iter().max().unwrap(),
            sizes.iter().min().unwrap(),
            sizes.iter().sum::<usize>() as f64 / sizes.len() as f64
        )
    }
}

impl fmt::Display for Inverted {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "Name")?;
        for sketch in &self.sample_names {
            writeln!(f, "{sketch}")?;
        }
        Ok(())
    }
}

#[cfg(all(test, not(target_arch = "wasm32")))]
mod tests {
    use super::*;
    use crate::io::parse_metadata_info;
    use std::io;

    fn index() -> Inverted {
        Inverted {
            index: Inverted::build_inverted_index(&[vec![1, 2], vec![1, 3]], 2),
            n_samples: 2,
            sample_names: vec!["refA".to_owned(), "refB".to_owned()],
            metadata: Some(vec!["description A".to_owned(), "description B".to_owned()]),
            sketch_version: env!("CARGO_PKG_VERSION").to_owned(),
            ..Inverted::default()
        }
    }

    #[test]
    fn loaded_metadata_is_validated_and_normalised_before_output() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("legacy");
        let prefix = prefix.to_str().unwrap();
        for separator in ['\t', ',', ';'] {
            let mut original = index();
            original.metadata.as_mut().unwrap()[0] = format!("old{separator}metadata");
            original.save(prefix).unwrap();
            let mut loaded = Inverted::load(prefix).unwrap();
            assert!(loaded.prepare_metadata_for_output(false).is_err());
            loaded.prepare_metadata_for_output(true).unwrap();
            assert_eq!(loaded.metadata_at(0), None);
            // Suppression affects the loaded object, not the stored metadata.
            assert_eq!(Inverted::load(prefix).unwrap(), original);
        }
        let mut loaded = index();
        loaded.metadata = Some(vec![String::new(), "  ".to_owned()]);
        loaded.prepare_metadata_for_output(false).unwrap();
        assert_eq!(loaded.metadata_at(0), Some("NA"));
        assert_eq!(loaded.metadata_at(1), Some("NA"));
    }

    #[test]
    fn streaming_preserves_rows_and_no_match_output() {
        let mut metadata_file = tempfile::NamedTempFile::new().unwrap();
        writeln!(metadata_file, "query\tquery description").unwrap();
        let metadata = parse_metadata_info(metadata_file.path().to_str(), ["query"]).unwrap();
        let index = index();
        for (query_type, expected) in [
            (
                InvertedQueryType::MatchCount,
                "Query\trefA;description A\trefB;description B\nquery;query description\t2\t1\n",
            ),
            (
                InvertedQueryType::AllBins,
                "Query\tMatches\nquery;query description\trefA;description A\n",
            ),
            (
                InvertedQueryType::AnyBins,
                "Query\tMatches\nquery;query description\trefA;description A,refB;description B\n",
            ),
        ] {
            let mut output = Vec::new();
            index
                .query_stream(
                    &mut output,
                    &[vec![1, 2]],
                    &["query".to_owned()],
                    &metadata,
                    &query_type,
                    1,
                    true,
                )
                .unwrap();
            assert_eq!(String::from_utf8(output).unwrap(), expected);
        }
        for query_type in [InvertedQueryType::AllBins, InvertedQueryType::AnyBins] {
            let mut output = Vec::new();
            index
                .query_stream(
                    &mut output,
                    &[vec![9, 9]],
                    &["query".to_owned()],
                    &metadata,
                    &query_type,
                    1,
                    true,
                )
                .unwrap();
            assert_eq!(
                String::from_utf8(output).unwrap(),
                "Query\tMatches\nquery;query description\n"
            );
        }
    }

    struct RowWriter {
        rows: usize,
        fail_rows: bool,
        fail_flush: bool,
    }

    impl IoWrite for RowWriter {
        fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
            // Header and query rows arrive as complete, already formatted strings.
            assert!(bytes.ends_with(b"\n"));
            assert_eq!(bytes.iter().filter(|&&b| b == b'\n').count(), 1);
            if self.fail_rows && self.rows > 0 {
                return Err(io::Error::new(
                    io::ErrorKind::BrokenPipe,
                    "row write failed",
                ));
            }
            self.rows += 1;
            Ok(bytes.len())
        }

        fn flush(&mut self) -> io::Result<()> {
            if self.fail_flush {
                Err(io::Error::new(io::ErrorKind::BrokenPipe, "flush failed"))
            } else {
                Ok(())
            }
        }
    }

    #[test]
    fn bounded_streaming_completes_and_propagates_write_errors() {
        let index = index();
        let queries = vec![vec![1, 2]; 100];
        let names: Vec<String> = (0..queries.len()).map(|i| format!("q{i}")).collect();
        for threads in [1, 2] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads + 1)
                .build()
                .unwrap();
            for (fail_rows, fail_flush) in [(false, false), (true, false), (false, true)] {
                let mut writer = RowWriter {
                    rows: 0,
                    fail_rows,
                    fail_flush,
                };
                let result = pool.install(|| {
                    index.query_stream(
                        &mut writer,
                        &queries,
                        &names,
                        &SampleMetadata::default(),
                        &InvertedQueryType::AnyBins,
                        threads,
                        true,
                    )
                });
                if fail_rows || fail_flush {
                    assert_eq!(
                        result
                            .unwrap_err()
                            .downcast_ref::<io::Error>()
                            .unwrap()
                            .kind(),
                        io::ErrorKind::BrokenPipe
                    );
                } else {
                    result.unwrap();
                    assert_eq!(writer.rows, queries.len() + 1);
                }
            }
        }
    }
}
