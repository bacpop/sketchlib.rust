//! Functions to support `aaHash` generation over sequences
#[cfg(not(target_arch = "wasm32"))]
use needletail::parse_fastx_file;

use std::cmp::Ordering;

use super::*;

/// A valid IUPAC amino acid character (upper- and lowercase)
#[inline(always)]
pub fn valid_aa(aa: u8) -> bool {
    matches!(aa, b'a' | b'c'..=b'i' | b'k'..=b'n' | b'p'..=b't' | b'v' | b'w' | b'y' | b'A' | b'C'..=b'I' | b'K'..=b'N' | b'P'..=b'T' | b'V' | b'W' | b'Y' )
}

// Split a 64-bit word into 33 and 31-bit sub-words and left-rotate them separately.
// Increases period from 64 to 1023.
#[inline(always)]
fn srol(x: u64) -> u64 {
    let m: u64 = ((x & 0x8000000000000000) >> 30) | ((x & 0x100000000) >> 32);
    ((x << 1) & 0xFFFFFFFDFFFFFFFF) | m
}

/// Stores forward and (optionally) reverse complement hashes of k-mers in a nucleotide sequence
#[derive(Debug, Clone)]
pub struct AaHashIterator {
    k: usize,
    level: AaLevel,
    fh: u64,
    index: usize,
    seq: Vec<u8>,
    invalid_count: usize,
}

impl RollHash for AaHashIterator {
    fn set_k(&mut self, k: usize) {
        self.k = k;
        if let Some(new_it) = Self::new_iterator(0, &self.level, &self.seq, k) {
            self.fh = new_it.0;
            self.index = new_it.1;
        } else {
            panic!(
                "K-mer larger than smallest valid sequence, which is:\n{}",
                std::str::from_utf8(&self.seq).unwrap()
            );
        }
    }

    /// Retrieve the current hash (minimum of forward and reverse complement hashes)
    fn curr_hash(&self) -> u64 {
        self.fh
    }

    fn hash_type(&self) -> HashType {
        HashType::AA(self.level.clone())
    }

    fn seq_len(&self) -> usize {
        self.seq.len()
    }

    fn seq(&self) -> &Vec<u8> {
        &self.seq
    }

    fn sketch_data(&self) -> (bool, [usize; 4], usize) {
        (false, [0, 0, 0, 0], self.invalid_count)
    }
}

impl AaHashIterator {
    /// Create a new empty aaHash iterator at the set comparison level
    pub fn default(level: AaLevel) -> Self {
        Self {
            k: 0,
            level,
            fh: 0,
            index: 0,
            seq: Vec::new(),
            invalid_count: 0,
        }
    }

    /// Create a new aaHash iterator from a fasta file, at the set comparison level
    pub fn new(files: &[String], level: AaLevel, concat_fasta: bool) -> Vec<Self> {
        let mut hash_vec = Vec::new();

        // Read sequence into memory (as we go through multiple times)
        log::debug!("Preprocessing sequence");
        let mut seq_hash_it = Self::default(level.clone());
        for file in files.iter() {
            let mut reader =
                parse_fastx_file(file).unwrap_or_else(|_| panic!("Invalid path/file: {file}"));
            loop {
                let record_read = reader.next();
                if let Some(record) = record_read {
                    let seqrec = record.expect("Invalid FASTA/Q record");
                    if seqrec.qual().is_some() {
                        panic!("Unexpected quality information with AA sequences in {file}. Correct sequence type set?");
                    } else {
                        for aa in seqrec.seq().iter() {
                            if valid_aa(*aa) {
                                seq_hash_it.seq.push(*aa)
                            } else {
                                seq_hash_it.invalid_count += 1;
                                seq_hash_it.seq.push(SEQSEP);
                            }
                        }
                    }
                    if concat_fasta {
                        hash_vec.push(seq_hash_it);
                        seq_hash_it = Self::default(level.clone());
                    } else {
                        seq_hash_it.seq.push(SEQSEP);
                    }
                } else {
                    break;
                }
            }
        }
        if !concat_fasta {
            hash_vec.push(seq_hash_it);
        }
        hash_vec
    }

    /// Create a new iterator from a 3di embedding file of a structure
    pub fn from_3di_file(files: &[String]) -> Vec<Self> {
        Self::new(files, AaLevel::Level1, false)
    }

    /// Create a new iterator from a 3di embedding string of a structure
    pub fn from_3di_string(sequence: String) -> Vec<Self> {
        let mut hash_it = Self::default(AaLevel::Level1);
        hash_it.seq = sequence.into_bytes();
        vec![hash_it]
    }

    fn new_iterator(
        mut start: usize,
        level: &AaLevel,
        seq: &[u8],
        k: usize,
    ) -> Option<(u64, usize)> {
        let mut fh = 0_u64;
        'outer: while start < (seq.len() - k) {
            '_inner: for (i, v) in seq[start..(start + k)].iter().enumerate() {
                // If invalid seq
                if !valid_aa(*v) {
                    start += i + 1;
                    if start >= seq.len() {
                        return None;
                    }
                    fh = 0;
                    continue 'outer; // Try again from new start
                }
                fh = srol(fh);
                fh ^= level.aa_seed_table(*v);
            }
            break 'outer; // success
        }
        if start >= (seq.len() - k) {
            return None;
        }

        Some((fh, start + k))
    }

    /// Move to the next k-mer by adding a new base, removing a base from the end, efficiently updating the hash.
    fn roll_fwd(&mut self, old_aa: u8, new_aa: u8) {
        self.fh = srol(self.fh);
        self.fh ^= self.level.aa_seed_table(new_aa);
        self.fh ^= self.level.aa_roll_table(old_aa, self.k);
    }
}

impl Iterator for AaHashIterator {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        match self.index.cmp(&self.seq_len()) {
            Ordering::Less => {
                let current = self.curr_hash();
                let new_aa = self.seq[self.index];
                // Restart hash if invalid base
                if !valid_aa(new_aa) {
                    if let Some(new_it) =
                        Self::new_iterator(self.index + 1, &self.level, &self.seq, self.k)
                    {
                        self.fh = new_it.0;
                        self.index = new_it.1;
                    } else {
                        self.index = self.seq_len(); // End of valid sequence
                    }
                } else {
                    self.roll_fwd(self.seq[self.index - self.k], new_aa);
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
        (self.seq_len(), Some(self.seq_len()))
    }
}

// This lets you use collect etc
impl ExactSizeIterator for AaHashIterator {}
