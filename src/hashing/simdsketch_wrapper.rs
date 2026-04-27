//! Wrapper around the simd-sketch crate to provide sketching for DNA sequences

use needletail::{parse_fastx_file, parser::Format};
use packed_seq::PackedNSeqVec;
use simd_sketch::{Sketch as Sketch_simd, SketchParams};

/// Wrapper around simd_sketch to make DNA sketches
pub fn sketch_with_simd(
    fastxvec: &[String],
    min_qual: u8,
    est_coverage: usize,
    sketchers: &[SketchParams],
) -> (bool, Vec<Sketch_simd>) {
    let (reads, seqs) = load_packed_nseqs(fastxvec, min_qual);
    let slices = seqs.iter().map(|seq| seq.as_slice()).collect::<Vec<_>>();

    let sketches = sketchers
        .iter()
        .map(|params| {
            let mut params = *params;
            if reads {
                params.coverage = est_coverage;
            } else {
                params.coverage = 1;
                params.count = 1;
            }
            params.build().sketch_seqs(&slices)
        })
        .collect::<Vec<_>>();
    (reads, sketches)
}

fn load_packed_nseqs(fastxvec: &[String], min_qual: u8) -> (bool, Vec<PackedNSeqVec>) {
    let mut reads = None;
    let mut seqs = Vec::with_capacity(fastxvec.len());

    for path in fastxvec {
        let mut reader =
            parse_fastx_file(path).unwrap_or_else(|_| panic!("Invalid path/file: {path}"));
        let mut packed = PackedNSeqVec::default();
        let mut ascii = Vec::with_capacity(16000);
        let mut qual = Vec::with_capacity(16000);

        while let Some(record) = reader.next() {
            let seqrec = record.expect("Invalid FASTA/Q record");
            let sample_reads = *reads.get_or_insert_with(|| seqrec.format() == Format::Fastq);
            let seq = seqrec.seq();
            ascii.extend_from_slice(&seq);
            ascii.push(b'N');

            if sample_reads && min_qual > 0 {
                let record_qual = seqrec
                    .qual()
                    .expect("FASTQ input with quality filtering must contain quality scores");
                qual.extend_from_slice(record_qual);
                qual.push(0);
                if ascii.len() > 16000 {
                    packed.push_from_ascii_and_quality(&ascii, &qual, min_qual);
                    ascii.clear();
                    qual.clear();
                }
            } else if ascii.len() > 16000 {
                packed.push_ascii(&ascii);
                ascii.clear();
            }
        }

        if !ascii.is_empty() {
            if reads.unwrap_or(false) && min_qual > 0 {
                packed.push_from_ascii_and_quality(&ascii, &qual, min_qual);
            } else {
                packed.push_ascii(&ascii);
            }
        }
        seqs.push(packed);
    }

    (reads.expect("Invalid FASTA/Q record"), seqs)
}

#[cfg(test)]
mod test {
    use std::{
        fs,
        path::{Path, PathBuf},
        time::{SystemTime, UNIX_EPOCH},
    };

    use packed_seq::Seq;
    use simd_sketch::{BitSketch, HashMode, Sketch, SketchAlg};

    use super::*;

    fn temp_path(ext: &str) -> PathBuf {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!(
            "sketchlib_simd_wrapper_{}_{}.{}",
            std::process::id(),
            unique,
            ext
        ))
    }

    fn path_string(path: &Path) -> String {
        path.to_str().unwrap().to_string()
    }

    fn packed_bases(seq: &PackedNSeqVec) -> Vec<u8> {
        seq.as_slice().seq.iter_bp().collect()
    }

    fn ambiguity(seq: &PackedNSeqVec) -> Vec<u8> {
        seq.as_slice().ambiguous.iter_bp().collect()
    }

    #[test]
    fn single_pass_fasta_loader_matches_packed_seq_loader() {
        let path = temp_path("fa");
        fs::write(&path, b">r1\nACGTNN\n>r2\nTTAA\n").unwrap();
        let paths = vec![path_string(&path)];

        let (reads, seqs) = load_packed_nseqs(&paths, 0);
        let expected = PackedNSeqVec::from_fastx(&path);
        fs::remove_file(&path).unwrap();

        assert!(!reads);
        assert_eq!(seqs.len(), 1);
        assert_eq!(packed_bases(&seqs[0]), packed_bases(&expected));
        assert_eq!(ambiguity(&seqs[0]), ambiguity(&expected));
    }

    #[test]
    fn single_pass_fastq_loader_detects_reads_without_quality_masking() {
        let path = temp_path("fq");
        fs::write(&path, b"@r1\nACGT\n+\n!!!!\n").unwrap();
        let paths = vec![path_string(&path)];

        let (reads, seqs) = load_packed_nseqs(&paths, 0);
        let expected = PackedNSeqVec::from_fastx(&path);
        fs::remove_file(&path).unwrap();

        assert!(reads);
        assert_eq!(packed_bases(&seqs[0]), packed_bases(&expected));
        assert_eq!(ambiguity(&seqs[0]), ambiguity(&expected));
    }

    #[test]
    fn single_pass_fastq_loader_applies_quality_masking() {
        let path = temp_path("fq");
        fs::write(&path, b"@r1\nACGT\n+\n!!!!\n").unwrap();
        let paths = vec![path_string(&path)];

        let (reads, seqs) = load_packed_nseqs(&paths, 20);
        let expected = PackedNSeqVec::from_fastq_with_quality(&path, 20);
        fs::remove_file(&path).unwrap();

        assert!(reads);
        assert_eq!(packed_bases(&seqs[0]), packed_bases(&expected));
        assert_eq!(ambiguity(&seqs[0]), ambiguity(&expected));
    }

    #[test]
    fn wrapper_reuses_loaded_sequences_for_multiple_k_values() {
        let path = temp_path("fa");
        fs::write(&path, b">r1\nACGTACGTACGTACGT\n").unwrap();
        let paths = vec![path_string(&path)];
        let sketchers = [5, 7].map(|k| SketchParams {
            alg: SketchAlg::Bucket,
            hash_mode: HashMode::NtHash64,
            rc: true,
            k,
            s: 16,
            b: 16,
            seed: 0,
            count: 3,
            coverage: 99,
            filter_empty: true,
            filter_out_n: true,
        });

        let (reads, sketches) = sketch_with_simd(&paths, 0, 99, &sketchers);
        fs::remove_file(&path).unwrap();

        assert!(!reads);
        assert_eq!(sketches.len(), 2);
        for sketch in sketches {
            let Sketch::BucketSketch(sketch) = sketch else {
                panic!("expected bucket sketch")
            };
            assert_eq!(sketch.count, 1);
            assert!(matches!(sketch.buckets, BitSketch::B16(_)));
        }
    }
}
