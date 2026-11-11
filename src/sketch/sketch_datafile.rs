//! I/O support and memory mapping used by for lower level read/write to .skd and .skq
use memmap2::Mmap;
use num_traits::ToBytes;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::mem;
use std::sync::atomic::{AtomicUsize, Ordering};

use anyhow::{ensure, Context, Error, Result};

use super::{aligned_sketch_vec_filled, aligned_sketch_vec_with_capacity, SketchVec};

/// Mutable byte view of initialized words for reading the on-disk representation.
fn sketch_bytes_mut(words: &mut [u64]) -> &mut [u8] {
    // SAFETY: u64 has no padding or invalid bit patterns. The byte slice covers
    // exactly the initialized allocation and exclusively borrows it.
    unsafe { std::slice::from_raw_parts_mut(words.as_mut_ptr().cast(), mem::size_of_val(words)) }
}

fn decode_little_endian(words: &mut [u64]) {
    #[cfg(target_endian = "big")]
    for word in words {
        *word = u64::from_le(*word);
    }
    #[cfg(target_endian = "little")]
    let _ = words;
}

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
            let mut sample_data = aligned_sketch_vec_with_capacity(input_reader.sample_stride);
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
        Self::try_open(filename, mmap, bin_stride, kmer_stride, sample_stride)
            .unwrap_or_else(|err| panic!("Could not read from {filename}: {err}"))
    }

    /// Fallible opening for loaders that propagate I/O errors.
    pub fn try_open(
        filename: &str,
        mmap: bool,
        bin_stride: usize,
        kmer_stride: usize,
        sample_stride: usize,
    ) -> Result<Self> {
        let sketch_reader = BufReader::new(File::open(filename)?);
        let mmap_file = if mmap {
            Some(Self::memmap_file(filename)?)
        } else {
            None
        };
        Ok(Self {
            sketch_reader,
            mmap_file,
            sample_stride,
            _bin_stride: bin_stride,
            _kmer_stride: kmer_stride,
        })
    }

    fn mmap_file(&self) -> Option<&Mmap> {
        self.mmap_file.as_ref()
    }

    /// Reads a single sketch from an .skq into memory
    pub fn read_all_from_skq(&mut self, total_number_bins: usize) -> Vec<u16> {
        // Fixed-size buffer for u16
        let mut buffer = [0u8; std::mem::size_of::<u16>()];
        // Stream the whole file and convert to u64 vec
        let mut flat_sketch_array: Vec<u16> = Vec::with_capacity(total_number_bins);
        while let Ok(_read) = self.sketch_reader.read_exact(&mut buffer) {
            flat_sketch_array.push(u16::from_le_bytes(buffer));
        }
        flat_sketch_array
    }

    /// Read all sketch data from an .skd into memory
    pub fn read_all_from_skd(&mut self, total_number_bins: usize) -> SketchVec {
        let mut flat_sketch_array = aligned_sketch_vec_filled(total_number_bins, 0);
        self.read_all_into(&mut flat_sketch_array)
            .expect("Could not read complete sketch data");
        flat_sketch_array
    }

    /// Read a complete .skd file directly into its final host allocation.
    pub fn read_all_into(&mut self, destination: &mut [u64]) -> Result<()> {
        ensure!(
            self.sketch_reader.get_ref().metadata()?.len() == mem::size_of_val(destination) as u64,
            "Sketch data size does not match metadata"
        );
        self.sketch_reader
            .read_exact(sketch_bytes_mut(destination))?;
        decode_little_endian(destination);
        Ok(())
    }

    /// Read selected samples from an .skd into memory, providing the selected
    /// samples with `sample_indices`
    pub fn read_batch_from_skd(&self, sample_indices: &[usize], sample_stride: usize) -> SketchVec {
        let mut flat_sketch_array =
            aligned_sketch_vec_filled(sample_stride * sample_indices.len(), 0);
        self.read_batch_into(sample_indices, sample_stride, &mut flat_sketch_array)
            .expect("Could not read selected sketch data");
        flat_sketch_array
    }

    /// Copy selected samples from the mmap directly into the destination.
    pub fn read_batch_into(
        &self,
        sample_indices: &[usize],
        sample_stride: usize,
        destination: &mut [u64],
    ) -> Result<()> {
        ensure!(sample_stride != 0, "Sketch sample stride must be nonzero");
        ensure!(
            sample_indices.len().checked_mul(sample_stride) == Some(destination.len()),
            "Sketch destination size does not match selected samples"
        );
        let mmap = self
            .mmap_file
            .as_ref()
            .context("Subset reading requires an mmap")?;
        let sample_bytes = sample_stride
            .checked_mul(mem::size_of::<u64>())
            .context("Sketch sample size overflow")?;
        for (&sample_idx, output) in sample_indices
            .iter()
            .zip(destination.chunks_mut(sample_stride))
        {
            let start = sample_idx
                .checked_mul(sample_bytes)
                .context("Sketch offset overflow")?;
            let end = start
                .checked_add(sample_bytes)
                .context("Sketch offset overflow")?;
            let source = mmap
                .get(start..end)
                .context("Selected sketch exceeds data file size")?;
            sketch_bytes_mut(output).copy_from_slice(source);
        }
        decode_little_endian(destination);
        Ok(())
    }

    // Memory maps a file
    fn memmap_file(filename: &str) -> Result<Mmap, Error> {
        let file = File::open(filename)?;
        let mmap = unsafe { Mmap::map(&file)? };
        Ok(mmap)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sketch::SKETCH_ALIGNMENT;
    use std::fs;

    #[test]
    fn skd_reads_are_64_byte_aligned() {
        let filename = std::env::temp_dir().join(format!(
            "sketchlib_skd_alignment_{}.skd",
            std::process::id()
        ));
        let filename = filename.to_str().expect("temporary path must be UTF-8");

        {
            let mut writer = SketchArrayWriter::new(filename, 1, 4, 4);
            writer.write_sketch(&[1_u64, 2, 3, 4]);
            writer.flush().expect("failed to flush test sketch");
        }

        let mut reader = SketchArrayReader::open(filename, false, 1, 4, 4);
        let all_sketches = reader.read_all_from_skd(4);
        assert_eq!(all_sketches.as_ptr() as usize % SKETCH_ALIGNMENT, 0);

        let reader = SketchArrayReader::open(filename, true, 1, 4, 4);
        let batch_sketches = reader.read_batch_from_skd(&[0], 4);
        assert_eq!(batch_sketches.as_ptr() as usize % SKETCH_ALIGNMENT, 0);

        drop(batch_sketches);
        drop(reader);
        fs::remove_file(filename).expect("failed to remove test sketch");
    }
}
