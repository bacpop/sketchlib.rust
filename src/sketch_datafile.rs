
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
    current_index: AtomicUsize,
    serial_writer: Arc<RwLock<BufWriter<File>>>,
}

impl SketchArrayFile {
    pub fn new(filename: &str, bin_stride: usize, kmer_stride: usize, sample_stride: usize) -> Self {
        let current_index = AtomicUsize::new(0);
        let serial_file = File::create(filename).expect("Couldn't write to {filename}");
        let serial_writer = Arc::new(RwLock::new(BufWriter::new(serial_file)));
        Self { bin_stride, kmer_stride, sample_stride, current_index, serial_writer }
    }

    pub fn write_sketch(&self, usigs_flat: &[u64]) -> usize {
        let writer = Arc::clone(&self.serial_writer);
        let mut writer = writer.write().unwrap();
        Self::write_sketch_data(&mut writer.get_mut(), usigs_flat).unwrap();
        self.current_index.fetch_add(1, Ordering::SeqCst)
    }

    fn write_sketch_data<W: Write + Seek>(writer: &mut W, usigs_flat: &[u64]) -> Result<(), Box<dyn Error>>{
        for bin_val in usigs_flat {
            writer.write_all(&bin_val.to_le_bytes())?;
        }
        Ok(())
    }
}
