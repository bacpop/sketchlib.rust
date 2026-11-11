//! The class to support .ski creation, reading and writing, containing an inverted
//! index of multiple sketches
#[cfg(target_arch = "wasm32")]
use std::fmt;
use std::sync::mpsc;
#[cfg(not(target_arch = "wasm32"))]
use std::{cmp, fmt};
use std::{fmt::Write as FmtWrite, io::Write as IoWrite};

#[cfg(not(target_arch = "wasm32"))]
extern crate needletail;

#[cfg(target_arch = "wasm32")]
use hashbrown::HashMap;
#[cfg(not(target_arch = "wasm32"))]
use hashbrown::{HashMap, HashSet};
use indicatif::ParallelProgressIterator;
use rayon::prelude::*;
use roaring::RoaringBitmap;
use serde::{Deserialize, Serialize};

use super::hashing::{
    bloom_filter::KmerFilter, nthash_iterator::NtHashIterator, HashType, RollHash,
};
#[cfg(not(target_arch = "wasm32"))]
use crate::io::InputFastx;
use crate::sketch::multisketch::{parse_version, MIN_SKETCH_VERSION};
#[cfg(target_arch = "wasm32")]
use crate::sketch::Sketch;
#[cfg(not(target_arch = "wasm32"))]
use crate::sketch::{sketch_datafile::SketchArrayWriter, Sketch};
use crate::utils::get_progress_bar;
use crate::{cli::InvertedQueryType, io::normalise_metadata};
use anyhow::{bail, Context, Error};
use std::fs::File;
use std::io::{BufReader, BufWriter};

#[cfg(target_arch = "wasm32")]
use crate::logw;
#[cfg(target_arch = "wasm32")]
use wasm_bindgen_file_reader::WebSysFile;

type InvSketches = (Vec<Vec<u16>>, Vec<String>);

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
    #[serde(skip)]
    query_metadata: Option<HashMap<String, String>>,
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
            query_metadata: None,
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

    pub(crate) fn metadata_at(&self, idx: usize) -> Option<&str> {
        self.metadata.as_ref().map(|values| values[idx].as_str())
    }

    /// Metadata for a named query, or `None` when unavailable or suppressed.
    pub fn get_query_metadata(&self, name: &str) -> Option<&str> {
        self.query_metadata.as_ref()?.get(name).map(String::as_str)
    }

    /// Disable query and reference metadata output for this loaded index.
    pub fn suppress_metadata(&mut self) {
        self.metadata = None;
        self.query_metadata = None;
    }

    /// Compute and format query rows in parallel while the caller writes them.
    pub(crate) fn query_stream<W: IoWrite + Send>(
        &self,
        writer: &mut W,
        queries: &[Vec<u16>],
        query_names: &[String],
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
                        append_sample(&mut row, name, self.get_query_metadata(name));
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

    /// Load a MessagePack index and initialize query metadata and output suppression.
    /// Accepts a file prefix on native targets or a browser file on WebAssembly.
    // NB MessagePack rather the CBOR uses here because of
    // https://github.com/enarx/ciborium/issues/96
    pub fn load(
        #[cfg(not(target_arch = "wasm32"))] file_prefix: &str,
        #[cfg(target_arch = "wasm32")] file: &web_sys::File,
        query_metadata: Option<HashMap<String, String>>,
        suppress_metadata: bool,
    ) -> Result<Self, Error> {
        #[cfg(not(target_arch = "wasm32"))]
        let ski_file = {
            let filename = format!("{file_prefix}.ski");
            log::info!("Loading inverted index from {filename}");
            BufReader::new(File::open(filename)?)
        };
        #[cfg(target_arch = "wasm32")]
        let ski_file = {
            logw("Loading inverted index", Some("info"));
            BufReader::new(WebSysFile::new(file.clone()))
        };
        let decompress_reader = snap::read::FrameDecoder::new(ski_file);
        let mut ski_obj: Self = rmp_serde::decode::from_read(decompress_reader)?;
        check_version(&ski_obj.sketch_version)?;
        ski_obj.query_metadata = query_metadata;
        if suppress_metadata {
            ski_obj.suppress_metadata();
        }
        if let Some(values) = &mut ski_obj.metadata {
            if values.len() != ski_obj.sample_names.len() {
                bail!("Index metadata length does not match the number of samples");
            }
            let mut missing = 0;
            for (name, value) in ski_obj.sample_names.iter().zip(values) {
                missing += usize::from(
                    normalise_metadata(value)
                        .with_context(|| format!("Invalid index metadata for sample {name}"))?,
                );
            }
            if missing > 0 {
                log::warn!(
                    "Metadata: {missing} missing values set to NA; 0 records without sample names ignored"
                );
            }
        }
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

    /// Count unique unordered pairs of distinct samples sharing at least one bin.
    /// Reconstructs hashes using two bytes per sample per bin, then retains only
    /// the matching-sample bitmap for each parallel query.
    pub fn count_prefilter_pairs(&self, quiet: bool) -> u64 {
        let n = self.n_samples();
        let n_bins = self.index.len();
        if n < 2 || n_bins == 0 {
            return 0;
        }

        let mut sample_hashes = vec![0_u16; n * n_bins];
        for (bin_idx, bin) in self.index.iter().enumerate() {
            for (&hash, samples) in bin {
                for sample_idx in samples {
                    sample_hashes[sample_idx as usize * n_bins + bin_idx] = hash;
                }
            }
        }

        let progress_bar = get_progress_bar(n, false, quiet);
        sample_hashes
            .par_chunks(n_bins)
            .enumerate()
            .progress_with(progress_bar)
            .map(|(sample_idx, query)| {
                let mut matches = RoaringBitmap::new();
                for (bin, hash) in self.index.iter().zip(query) {
                    if let Some(samples) = bin.get(hash) {
                        matches |= samples;
                    }
                }
                // rank includes self and all lower indexes; count each pair once.
                matches.len() - matches.rank(sample_idx as u32)
            })
            .sum()
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
    fn prefilter_count_matches_unique_pair_oracle() {
        let mut cases = vec![
            Vec::<Vec<u16>>::new(),
            vec![vec![]; 3],
            vec![vec![0, u16::MAX]],
            vec![vec![0, 0], vec![1, 1], vec![2, 2]],
            vec![vec![3, 7]; 6],
            vec![
                vec![1, 11, 21],
                vec![1, 12, 22],
                vec![2, 12, 23],
                vec![3, 13, 21],
            ],
        ];
        for n in [2, 7, 40] {
            for n_bins in [1, 4, 11] {
                cases.push(
                    (0..n)
                        .map(|sample| {
                            (0..n_bins)
                                .map(|bin| ((sample * 17 + bin * 7 + sample * bin * 3) % 11) as u16)
                                .collect()
                        })
                        .collect(),
                );
            }
        }
        for threads in [1, 2] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            for sketches in &cases {
                let n_bins = sketches.first().map_or(0, Vec::len);
                let index = Inverted {
                    index: Inverted::build_inverted_index(sketches, n_bins as u64),
                    sample_names: (0..sketches.len()).map(|i| format!("sample{i}")).collect(),
                    n_samples: sketches.len(),
                    ..Inverted::default()
                };
                // Independent oracle: inspect each distinct pair of input sketches once.
                let expected = sketches
                    .iter()
                    .enumerate()
                    .map(|(i, a)| {
                        sketches
                            .iter()
                            .skip(i + 1)
                            .filter(|b| a.iter().zip(b.iter()).any(|(x, y)| x == y))
                            .count() as u64
                    })
                    .sum::<u64>();
                assert_eq!(
                    pool.install(|| index.count_prefilter_pairs(true)),
                    expected,
                    "{sketches:?}"
                );
            }
        }
    }

    #[test]
    fn streaming_preserves_rows_and_no_match_output() {
        let mut metadata_file = tempfile::NamedTempFile::new().unwrap();
        writeln!(metadata_file, "query\tquery description").unwrap();
        let mut index = index();
        index.query_metadata =
            parse_metadata_info(metadata_file.path().to_str(), ["query"]).unwrap();
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
