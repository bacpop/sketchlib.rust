//! I/O support and memory mapping used by for lower level read/write to .skd
use memmap2::Mmap;
use std::error::Error;
use std::fs::File;
use std::io::Read;
use std::io::{BufReader, BufWriter, Write};
use std::mem;
use std::sync::atomic::{AtomicUsize, Ordering};

#[derive(Debug)]
pub struct SketchArrayFile {
    serial_writer: SketchWriter,
    // Fields here in case needed, but not currently used
    _bin_stride: usize,
    _kmer_stride: usize,
    _sample_stride: usize,
}

// Both the file and index are locked together
#[derive(Debug)]
pub struct SketchWriter {
    writer: BufWriter<File>,
    current_index: AtomicUsize,
}

impl SketchWriter {
    pub fn new(filename: &str) -> Self {
        let current_index = AtomicUsize::new(0);
        log::info!("Saving sketch data to {filename}");
        let writer = BufWriter::new(File::create(filename).expect("Couldn't write to {filename}"));
        Self {
            writer,
            current_index,
        }
    }
}

impl SketchArrayFile {
    pub fn new(
        filename: &str,
        bin_stride: usize,
        kmer_stride: usize,
        sample_stride: usize,
    ) -> Self {
        let serial_writer = SketchWriter::new(filename);
        Self {
            serial_writer,
            _bin_stride: bin_stride,
            _kmer_stride: kmer_stride,
            _sample_stride: sample_stride,
        }
    }

    pub fn write_sketch(&mut self, usigs_flat: &[u64]) -> usize {
        Self::write_sketch_data(&mut self.serial_writer.writer, usigs_flat).unwrap();
        self.serial_writer
            .current_index
            .fetch_add(1, Ordering::SeqCst)
    }

    pub fn read_all(filename: &str, number_bins: usize) -> Vec<u64> {
        // Just stream the whole file and convert to u64 vec
        let mut reader = BufReader::new(
            File::open(filename).unwrap_or_else(|_| panic!("Could not read from {filename}")),
        );
        // TODO: longer buffer may be better (and below too)
        let mut buffer = [0u8; mem::size_of::<u64>()];
        let mut flat_sketch_array: Vec<u64> = Vec::with_capacity(number_bins);
        while let Ok(_read) = reader.read_exact(&mut buffer) {
            flat_sketch_array.push(u64::from_le_bytes(buffer));
        }
        flat_sketch_array
    }

    pub fn read_batch(filename: &str, sample_indices: &[usize], sample_stride: usize) -> Vec<u64> {
        // Just stream the whole file and convert to u64 vec
        let mmap = Self::memmap_file(filename)
            .unwrap_or_else(|_| panic!("Could not memory map {filename}"));
        let mut flat_sketch_array: Vec<u64> =
            Vec::with_capacity(sample_stride * sample_indices.len());
        // TODO possible improvement would be to combine adjacent indices into ranges
        for sample_idx in sample_indices {
            for bin_idx in 0..sample_stride {
                let start_byte = (sample_idx * sample_stride + bin_idx) * mem::size_of::<u64>();
                let bin_val =
                    u64::from_le_bytes(*array_ref![mmap, start_byte, mem::size_of::<u64>()]);
                // let bin_val = u64::from_le_bytes(mmap[start_byte..(start_byte + mem::size_of::<u64>())].try_into().unwrap());
                flat_sketch_array.push(bin_val);
            }
        }
        flat_sketch_array
    }

    pub fn write_batch(
        input_filename: &str,
        output_filename: &str,
        sample_indices: &[usize],
        sample_stride: usize,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mmap = Self::memmap_file(input_filename)
            .unwrap_or_else(|_| panic!("Could not memory map {input_filename}"));

        let output_file = File::create(output_filename)?;
        let mut writer = BufWriter::new(output_file);

        for &sample_idx in sample_indices {
            let mut sample_data = Vec::with_capacity(sample_stride);
            for bin_idx in 0..sample_stride {
                let start_byte = (sample_idx * sample_stride + bin_idx) * mem::size_of::<u64>();
                let bin_val =
                    u64::from_le_bytes(*array_ref![mmap, start_byte, mem::size_of::<u64>()]);
                sample_data.push(bin_val);
            }
            Self::write_sketch_data(&mut writer, &sample_data)?;
        }

        writer.flush()?;
        Ok(())
    }

    fn memmap_file(filename: &str) -> Result<Mmap, Box<dyn Error>> {
        let file = File::open(filename)?;
        let mmap = unsafe { Mmap::map(&file)? };
        Ok(mmap)
    }

    pub fn write_sketch_data<W: Write>(
        writer: &mut W,
        usigs_flat: &[u64],
    ) -> Result<(), Box<dyn Error>> {
        for bin_val in usigs_flat {
            writer.write_all(&bin_val.to_le_bytes())?;
        }
        Ok(())
    }
}
