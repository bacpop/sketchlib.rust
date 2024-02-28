
use std::collections::HashMap;

extern crate needletail;
use needletail::{parse_fastx_file, parser::Format};
use ska::ska_dict::bloom_filter::KmerFilter;

use super::hashing::NtHashIterator;

/// Bin bits
pub const BBITS: u64 = 14;
pub const SIGN_MOD: u64 = (1 << 61) - 1;

struct Sketch {
    name: String,
    sketch_size: u64,
    kmer_sketches: HashMap<usize, Vec<u64>>,
    rc: bool,
    reads: bool,
    seq_length: usize,
    missing_bases: usize,
    densified: bool,
    kmer_filter: KmerFilter,
}

impl Sketch {
    pub fn new(name: &str, files: (&str, Option<&String>), kmer_lengths: &[usize], sketch_size: u64, min_count: u16, rc: bool) -> Self {
        let mut sketch = Self {
            name: name.to_string(),
            sketch_size: sketch_size,
            kmer_sketches: HashMap::with_capacity(kmer_lengths.len()),
            rc: rc,
            reads: false,
            seq_length: 0,
            missing_bases: 0,
            densified: false,
            kmer_filter: KmerFilter::new(min_count)
        };

        // Check if we're working with reads, and initalise the CM filter if so
        let mut reader_peek =
            parse_fastx_file(files.0).unwrap_or_else(|_| panic!("Invalid path/file: {}", files.0));
        let seq_peek = reader_peek
            .next()
            .expect("Invalid FASTA/Q record")
            .expect("Invalid FASTA/Q record");
        let mut is_reads = false;
        if seq_peek.format() == Format::Fastq {
            sketch.kmer_filter.init();
            is_reads = true;
        }

        // Read sequence into memory (as we go through multiple times)
        // TODO allow it to be streamed from disk
        let mut sequence: Vec<Vec<u8>> = Vec::new();
        let mut quals: Vec<Vec<u8>> = Vec::new();

        sketch.add_seq(files.0, &mut sequence, &mut quals);
        if let Some(filename) = files.1 {
            sketch.add_seq(filename, &mut sequence, &mut quals);
        }

        if sequence.len() == 0 {
            panic!("{} has no valid sequence", files.0);
        }

        // Build the sketches across k-mer lengths
        let num_bins: u64 = sketch_size * (u64::BITS as u64);
        let bin_size: u64 = (SIGN_MOD + num_bins - 1) / num_bins;
        for k in kmer_lengths {
            let usigs = vec![0; (sketch_size * BBITS) as usize];
            let signs = vec![u64::MAX, (sketch_size as u32 * u64::BITS) as u64];
            for seq in &sequence {
                let hash_it = NtHashIterator::new(seq, seq.len(), *k, rc);
                for hash in hash_it {
                    todo!()
                }
            }
        }

        /*
        if (sequence.nseqs() == 0) {
            throw std::runtime_error(name + " contains no sequence");
          }
          _bases = sequence.get_composition();
          _missing_bases = sequence.missing_bases();
        
          double minhash_sum = 0.0;
          for (auto kmer_it = kmers.cbegin(); kmer_it != kmers.cend(); ++kmer_it) {
            double minhash = 0;
            bool densified;
            std::tie(usigs[kmer_it->first], minhash, densified) =
                sketch(sequence, sketchsize64, kmer_it->second, _bbits, codon_phased,
                       _use_rc, min_count, exact);
        
            minhash_sum += minhash;
            _densified |= densified; // Densified at any k-mer length
          }
        
          // 1/N =~ 1/E(Y) where Yi = minhash in [0,1] for k-mer i
          // See
          // https://www.cs.princeton.edu/courses/archive/fall13/cos521/lecnotes/lec4final.pdf
          if (sequence.is_reads()) {
            _seq_size = static_cast<size_t>((double)usigs.size() / minhash_sum);
          } else {
            _seq_size = _bases.total;
          }
        */
        sketch
    }

    // TODO this can count Ns and base composition
    // TODO would probably also make sense to encode the whole thing as 0/1/2/3 here
    fn add_seq(&mut self, filename: &str, sequence: &mut Vec<Vec<u8>>, quals: &mut Vec<Vec<u8>>) {
        let mut reader =
            parse_fastx_file(filename).unwrap_or_else(|_| panic!("Invalid path/file: {filename}"));
        while let Some(record) = reader.next() {
            let seqrec = record.expect("Invalid FASTA/Q record");
            self.seq_length += seqrec.num_bases();
            sequence.push(seqrec.seq().to_vec());
            if let Some(qual_scores) = seqrec.qual() {
                quals.push(qual_scores.to_vec());
            }
        }
    }
}

// TODO may want to write these one by one rather than storing them all
pub fn sketch_files() -> Vec<Sketch> {
    todo!()
}