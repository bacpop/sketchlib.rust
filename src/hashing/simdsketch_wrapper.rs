//! Wrapper around the simd-sketch crate to provide sketching for DNA sequences

use needletail::{parse_fastx_file, parser::Format};
use packed_seq::PackedNSeqVec;
use simd_sketch::{SketchParams, Sketch as Sketch_simd};
use std::path::Path;

/// Wrapper around simd_sketch to make DNA sketches
pub fn sketch_with_simd(
    fastxvec: &Vec<String>,
    min_qual: u8,
    est_coverage: usize,
    mut sketchers: Vec<SketchParams>,
) -> (bool, Vec<Sketch_simd>) {
    let mut reader_peek =
        parse_fastx_file(fastxvec[0].clone()).unwrap_or_else(|_| panic!("Invalid path/file: {}", fastxvec[0]));
    let seq_peek = reader_peek
        .next()
        .expect("Invalid FASTA/Q record")
        .expect("Invalid FASTA/Q record");
    let mut reads = false;
    if seq_peek.format() == Format::Fastq {
        reads = true;
        for is in sketchers.iter_mut() {
            is.coverage = est_coverage;
        }
    } else {
        for is in sketchers.iter_mut() {
            is.coverage = 1;
        }
    }

    let seqs: Vec<PackedNSeqVec>;

    if min_qual == 0 || !reads {
        seqs = fastxvec.into_iter().map(|path| PackedNSeqVec::from_fastx(Path::new(path))).collect();
    } else {
        seqs = fastxvec.into_iter().map(|path| PackedNSeqVec::from_fastq_with_quality(Path::new(path), min_qual)).collect();
    }
    
    // // Run the sketching
    let sketches = sketchers
        .iter()
        .map(|is: &SketchParams| {
            is.build().sketch_seqs(&seqs.iter().map(|x| x.as_slice()).collect::<Vec<_>>())
        }).collect::<Vec<_>>();
        (reads, sketches)
}
