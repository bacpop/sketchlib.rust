//! I/O support and memory mapping used by for lower level read/write to .skd and .skq
use memmap2::Mmap;
use num_traits::ToBytes;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::mem;
use std::sync::atomic::{AtomicUsize, Ordering};

use anyhow::{Error, Result};

/// Write to an .skd or .skq file
/// This writer is not parallel safe, serial only!
#[derive(Debug)]
pub struct SketchArrayWriter {
    writer: BufWriter<File>,
    current_index: AtomicUsize,
    // Fields here in case needed, but not currently used
    _bin_stride: usize,
    _kmer_stride: usize,
    _sample_stride: usize,
}

/// Append selected samples to an .skd file providing the selected samples
/// with `sample_indices`
pub fn append_batch(
    input_reader: &SketchArrayReader,
    output_writer: &mut SketchArrayWriter,
    sample_indices: &[usize],
) -> Result<(), Error> {
    if let Some(mmap) = input_reader.mmap_file() {
        for &sample_idx in sample_indices {
            let mut sample_data = Vec::with_capacity(input_reader.sample_stride);
            for bin_idx in 0..input_reader.sample_stride {
                let start_byte =
                    (sample_idx * input_reader.sample_stride + bin_idx) * mem::size_of::<u64>();
                let bin_val =
                    u64::from_le_bytes(*array_ref![mmap, start_byte, mem::size_of::<u64>()]);
                sample_data.push(bin_val);
            }
            output_writer.write_sketch(&sample_data);
        }
    } else {
        unimplemented!("Batch reading of non memmapped files unsupported");
    }
    output_writer.flush()
}

impl SketchArrayWriter {
    /// Open a new sketch file for writing
    pub fn new(
        filename: &str,
        bin_stride: usize,
        kmer_stride: usize,
        sample_stride: usize,
    ) -> Self {
        log::info!("Saving sketch data to {filename}");
        let current_index = AtomicUsize::new(0);
        let writer = BufWriter::new(
            File::create(filename).unwrap_or_else(|_| panic!("Couldn't write to {filename}")),
        );
        Self {
            writer,
            current_index,
            _bin_stride: bin_stride,
            _kmer_stride: kmer_stride,
            _sample_stride: sample_stride,
        }
    }

    /// Write a single sketch to the currently open file, appending to the end
    /// (at the current write position)
    ///
    /// # Returns:
    ///
    /// - Index in file where sketch was written
    pub fn write_sketch<HashSize: ToBytes>(&mut self, usigs_flat: &[HashSize]) -> usize {
        Self::write_sketch_data(&mut self.writer, usigs_flat).unwrap();
        self.current_index.fetch_add(1, Ordering::SeqCst)
    }

    /// Flush the underlying writer
    pub fn flush(&mut self) -> Result<(), Error> {
        self.writer.flush()?;
        Ok(())
    }

    fn write_sketch_data<W: Write, HashSize: ToBytes>(
        writer: &mut W,
        usigs_flat: &[HashSize],
    ) -> Result<(), Box<dyn std::error::Error>> {
        for bin_val in usigs_flat {
            writer.write_all(bin_val.to_le_bytes().as_ref())?;
        }
        Ok(())
    }
}

/// Reads from an .skd or .skq file
#[derive(Debug)]
pub struct SketchArrayReader {
    sketch_reader: BufReader<File>,
    mmap_file: Option<Mmap>,
    /// Number of bins between each sample in the underlying file
    pub sample_stride: usize,
    // Fields here in case needed, but not currently used
    _bin_stride: usize,
    _kmer_stride: usize,
}

impl SketchArrayReader {
    /// Open an .skd or .skq file for reading
    ///
    /// If `mmap` is true, will memorymap the underlying file, which is more
    /// efficient when reading subsets of data from the file
    pub fn open(
        filename: &str,
        mmap: bool,
        bin_stride: usize,
        kmer_stride: usize,
        sample_stride: usize,
    ) -> Self {
        let sketch_reader = BufReader::new(
            File::open(filename).unwrap_or_else(|_| panic!("Could not read from {filename}")),
        );
        let mmap_file = if mmap {
            Some(
                Self::memmap_file(filename)
                    .unwrap_or_else(|_| panic!("Could not memory map {filename}")),
            )
        } else {
            None
        };
        Self {
            sketch_reader,
            mmap_file,
            sample_stride,
            _bin_stride: bin_stride,
            _kmer_stride: kmer_stride,
        }
    }

    fn mmap_file(&self) -> Option<&Mmap> {
        self.mmap_file.as_ref()
    }

    /// Reads a single sketch from an .skq into memory
    pub fn read_all_from_skq(&mut self, total_number_bins_hint: Option<usize>) -> Vec<u16> {
        // Fixed-size buffer for u16
        let mut buffer = [0u8; std::mem::size_of::<u16>()];
        // Stream the whole file and convert to u64 vec
        let mut flat_sketch_array: Vec<u16> = Vec::new();
        if let Some(nbins) = total_number_bins_hint {
            flat_sketch_array.reserve_exact(nbins);
        }
        while let Ok(_read) = self.sketch_reader.read_exact(&mut buffer) {
            flat_sketch_array.push(u16::from_le_bytes(buffer));
        }
        flat_sketch_array
    }

    /// Read all sketch data from an .skd into memory
    pub fn read_all_from_skd(&mut self, total_number_bins: usize) -> Vec<u64> {
        // Fixed-size buffer for u64
        let mut buffer = [0u8; std::mem::size_of::<u64>()];
        // Stream the whole file and convert to u64 vec
        let mut flat_sketch_array: Vec<u64> = Vec::with_capacity(total_number_bins);
        while let Ok(_read) = self.sketch_reader.read_exact(&mut buffer) {
            flat_sketch_array.push(u64::from_le_bytes(buffer));
        }
        flat_sketch_array
    }

    /// Read selected samples from an .skd into memory, providing the selected
    /// samples with `sample_indices`
    pub fn read_batch_from_skd(&self, sample_indices: &[usize], sample_stride: usize) -> Vec<u64> {
        let mut flat_sketch_array: Vec<u64> =
            Vec::with_capacity(sample_stride * sample_indices.len());
        // TODO possible improvement would be to combine adjacent indices into ranges
        if let Some(mmap_array) = &self.mmap_file {
            for sample_idx in sample_indices {
                for bin_idx in 0..sample_stride {
                    let start_byte = (sample_idx * sample_stride + bin_idx) * mem::size_of::<u64>();
                    let bin_val = u64::from_le_bytes(*array_ref![
                        mmap_array,
                        start_byte,
                        mem::size_of::<u64>()
                    ]);
                    // let bin_val = u64::from_le_bytes(mmap[start_byte..(start_byte + mem::size_of::<u64>())].try_into().unwrap());
                    flat_sketch_array.push(bin_val);
                }
            }
        } else {
            unimplemented!("Batch reading of non memmapped files unsupported");
        }

        flat_sketch_array
    }

    // Memory maps a file
    fn memmap_file(filename: &str) -> Result<Mmap, Error> {
        let file = File::open(filename)?;
        let mmap = unsafe { Mmap::map(&file)? };
        Ok(mmap)
    }
}
