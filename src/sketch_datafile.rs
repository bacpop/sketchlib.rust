
use std::sync::{Arc, RwLock, atomic::{AtomicUsize, Ordering}};
use std::{collections::HashMap, process::Output};
use std::mem;
use std::fs::File;
use std::io::{BufReader, BufWriter, Seek, SeekFrom, Write};
use std::error::Error;

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
    current_index: AtomicUsize
}

impl SketchWriter {
    pub fn new(filename: &str) -> Self {
        let current_index = AtomicUsize::new(0);
        let writer = BufWriter::new(File::create(filename).expect("Couldn't write to {filename}"));
        Self {writer, current_index}
    }
}

impl SketchArrayFile {
    pub fn new(filename: &str, bin_stride: usize, kmer_stride: usize, sample_stride: usize) -> Self {
        let serial_writer = Arc::new(RwLock::new(SketchWriter::new(filename)));
        Self { bin_stride, kmer_stride, sample_stride, serial_writer }
    }

    pub fn write_sketch(&self, usigs_flat: &[u64]) -> usize {
        let writer = Arc::clone(&self.serial_writer);
        // Gets the lock. This will be dropped at the end of the function
        let mut writer = writer.write().unwrap();
        Self::write_sketch_data(&mut writer.writer.get_mut(), usigs_flat).unwrap();
        writer.current_index.fetch_add(1, Ordering::SeqCst)
    }

    fn write_sketch_data<W: Write>(writer: &mut W, usigs_flat: &[u64]) -> Result<(), Box<dyn Error>>{
        for bin_val in usigs_flat {
            writer.write_all(&bin_val.to_le_bytes())?;
        }
        Ok(())
    }
}
