use needletail::{parse_fastx_file, parser::Format};
use std::cmp::Ordering;

use super::*;

/// Stores forward and (optionally) reverse complement hashes of k-mers in a nucleotide sequence
#[derive(Debug)]
pub struct NtHashIterator {
    k: usize,
    rc: bool,
    fh: u64,
    rh: Option<u64>,
    index: usize,
    seq: Vec<u8>,
    seq_len: usize, // NB seq.len() - 1 due to terminating char
    acgt: [usize; 4],
    non_acgt: usize,
    reads: bool,
}

impl RollHash for NtHashIterator {
    fn set_k(&mut self, k: usize) {
        self.k = k;
        if let Some(new_it) = Self::new_iterator(0, &self.seq, k, self.rc) {
            self.fh = new_it.0;
            self.rh = new_it.1;
            self.index = new_it.2;
        } else {
            panic!("K-mer larger than smallest valid sequence");
        }
    }

    /// Retrieve the current hash (minimum of forward and reverse complement hashes)
    fn curr_hash(&self) -> u64 {
        if let Some(rev) = self.rh {
            u64::min(self.fh, rev)
        } else {
            self.fh
        }
    }

    fn hash_type(&self) -> HashType {
        HashType::DNA
    }

    fn seq_len(&self) -> usize {
        self.acgt.iter().sum()
    }

    fn seq(&self) -> &Vec<u8> {
        &self.seq
    }

    fn sketch_data(&self) -> (bool, [usize; 4], usize) {
        (self.reads, self.acgt, self.non_acgt)
    }

    fn reads(&self) -> bool {
        self.reads
    }
}

impl NtHashIterator {
    /// Creates a new ntHash iterator, by loading DNA sequences into memory
    pub fn new(files: (&str, Option<&String>), rc: bool, min_qual: u8) -> Vec<Self> {
        // Check if we're working with reads, and initalise the filter if so
        let mut reader_peek =
            parse_fastx_file(files.0).unwrap_or_else(|_| panic!("Invalid path/file: {}", files.0));
        let seq_peek = reader_peek
            .next()
            .expect("Invalid FASTA/Q record")
            .expect("Invalid FASTA/Q record");
        let mut reads = false;
        if seq_peek.format() == Format::Fastq {
            reads = true;
        }

        let mut hash_it = Self {
            k: 0,
            rc,
            fh: 0,
            rh: None,
            index: 0,
            seq: Vec::new(),
            seq_len: 0,
            acgt: [0, 0, 0, 0],
            non_acgt: 0,
            reads,
        };

        // Read sequence into memory (as we go through multiple times)
        log::debug!("Preprocessing sequence");
        hash_it.add_dna_seq(files.0, min_qual);
        if let Some(filename) = files.1 {
            hash_it.add_dna_seq(filename, min_qual);
        }
        hash_it.seq_len = hash_it.seq.len() - 1;
        vec![hash_it]
    }

    fn add_dna_seq(&mut self, filename: &str, min_qual: u8) {
        let mut reader =
            parse_fastx_file(filename).unwrap_or_else(|_| panic!("Invalid path/file: {filename}"));
        while let Some(record) = reader.next() {
            let seqrec = record.expect("Invalid FASTA/Q record");
            if let Some(quals) = seqrec.qual() {
                for (base, qual) in seqrec.seq().iter().zip(quals) {
                    if *qual >= min_qual {
                        if valid_base(*base) {
                            let encoded_base = encode_base(*base);
                            self.acgt[encoded_base as usize] += 1;
                            self.seq.push(encoded_base)
                        } else {
                            self.non_acgt += 1;
                            self.seq.push(SEQSEP);
                        }
                    } else {
                        self.seq.push(SEQSEP);
                    }
                }
            } else {
                for base in seqrec.seq().iter() {
                    if valid_base(*base) {
                        let encoded_base = encode_base(*base);
                        self.acgt[encoded_base as usize] += 1;
                        self.seq.push(encoded_base)
                    } else {
                        self.non_acgt += 1;
                        self.seq.push(SEQSEP);
                    }
                }
            }

            self.seq.push(SEQSEP);
        }
    }

    fn new_iterator(
        mut start: usize,
        seq: &[u8],
        k: usize,
        rc: bool,
    ) -> Option<(u64, Option<u64>, usize)> {
        let mut fh = 0_u64;
        'outer: while start < (seq.len() - k) {
            '_inner: for (i, v) in seq[start..(start + k)].iter().enumerate() {
                // If invalid seq
                if *v > 3 {
                    start += i + 1;
                    if start >= seq.len() {
                        return None;
                    }
                    fh = 0;
                    continue 'outer; // Try again from new start
                }
                fh = fh.rotate_left(1_u32);
                fh = swapbits033(fh);
                fh ^= nthash_tables::HASH_LOOKUP[*v as usize];
            }
            break 'outer; // success
        }
        if start >= (seq.len() - k) {
            return None;
        }

        let rh = if rc {
            let mut h = 0_u64;
            for v in seq[start..(start + k)].iter().rev() {
                h = h.rotate_left(1_u32);
                h = swapbits033(h);
                h ^= nthash_tables::RC_HASH_LOOKUP[*v as usize];
            }
            Some(h)
        } else {
            None
        };

        Some((fh, rh, start + k))
    }

    /// Move to the next k-mer by adding a new base, removing a base from the end, efficiently updating the hash.
    fn roll_fwd(&mut self, old_base: u8, new_base: u8) {
        self.fh = self.fh.rotate_left(1_u32);
        self.fh = swapbits033(self.fh);
        self.fh ^= nthash_tables::HASH_LOOKUP[new_base as usize];
        self.fh ^= nthash_tables::MS_TAB_31L[(old_base as usize * 31) + (self.k % 31)]
            | nthash_tables::MS_TAB_33R[(old_base as usize) * 33 + (self.k % 33)];

        if let Some(rev) = self.rh {
            let mut h = rev
                ^ (nthash_tables::MS_TAB_31L[(rc_base(new_base) as usize * 31) + (self.k % 31)]
                    | nthash_tables::MS_TAB_33R[(rc_base(new_base) as usize) * 33 + (self.k % 33)]);
            h ^= nthash_tables::RC_HASH_LOOKUP[old_base as usize];
            h = h.rotate_right(1_u32);
            h = swapbits3263(h);
            self.rh = Some(h);
        };
    }
}

impl Iterator for NtHashIterator {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        match self.index.cmp(&self.seq_len) {
            Ordering::Less => {
                let current = self.curr_hash();
                let new_base = self.seq[self.index];
                // Restart hash if invalid base
                if new_base > 3 {
                    if let Some(new_it) =
                        Self::new_iterator(self.index + 1, &self.seq, self.k, self.rc)
                    {
                        self.fh = new_it.0;
                        self.rh = new_it.1;
                        self.index = new_it.2;
                    } else {
                        self.index = self.seq_len; // End of valid sequence
                    }
                } else {
                    self.roll_fwd(self.seq[self.index - self.k], new_base);
                    self.index += 1;
                }
                Some(current)
            }
            Ordering::Equal => {
                // Final hash, do not roll forward further
                self.index += 1;
                Some(self.curr_hash())
            }
            Ordering::Greater => {
                // End of sequence
                None
            }
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.seq_len, Some(self.seq_len))
    }
}

// This lets you use collect etc
impl ExactSizeIterator for NtHashIterator {}
