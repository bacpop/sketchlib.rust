use memmap2::Mmap;
use std::error::Error;
use std::fs::File;
use std::io::Read;
use std::io::{BufReader, BufWriter, Seek, SeekFrom, Write};
use std::mem;
use std::sync::{
    atomic::{AtomicUsize, Ordering},
    Arc, RwLock,
};
use std::{collections::HashMap, process::Output};

#[derive(Debug)]
pub struct SketchArrayFile {
    bin_stride: usize,
    kmer_stride: usize,
    sample_stride: usize,
    serial_writer: Arc<RwLock<SketchWriter>>,
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
        let serial_writer = Arc::new(RwLock::new(SketchWriter::new(filename)));
        Self {
            bin_stride,
            kmer_stride,
            sample_stride,
            serial_writer,
        }
    }

    pub fn write_sketch(&self, usigs_flat: &[u64]) -> usize {
        let writer = Arc::clone(&self.serial_writer);
        // Gets the lock. This will be dropped at the end of the function
        let mut writer = writer.write().unwrap();
        Self::write_sketch_data(&mut writer.writer.get_mut(), usigs_flat).unwrap();
        writer.current_index.fetch_add(1, Ordering::SeqCst)
    }

    pub fn read_all(filename: &str, number_bins: usize) -> Vec<u64> {
        // Just stream the whole file and convert to u64 vec
        let mut reader =
            BufReader::new(File::open(filename).expect(&format!("Could not read from {filename}")));
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
        let mut mmap =
            Self::memmap_file(filename).expect(&format!("Could not memory map {filename}"));
        let mut buffer = [0u8; mem::size_of::<u64>()];
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

    fn memmap_file(filename: &str) -> Result<Mmap, Box<dyn Error>> {
        let file = File::open(filename)?;
        let mmap = unsafe { Mmap::map(&file)? };
        Ok(mmap)
    }

    fn write_sketch_data<W: Write>(
        writer: &mut W,
        usigs_flat: &[u64],
    ) -> Result<(), Box<dyn Error>> {
        for bin_val in usigs_flat {
            writer.write_all(&bin_val.to_le_bytes())?;
        }
        Ok(())
    }
}
