//! Functions to support `ntHash` generation over sequences
#[cfg(not(target_arch = "wasm32"))]
use needletail::{parse_fastx_file, parser::Format};

#[cfg(target_arch = "wasm32")]
use crate::fastx_wasm::{open_fasta, open_fastq};
#[cfg(target_arch = "wasm32")]
use seq_io::fasta::Record as FastaRecord;
#[cfg(target_arch = "wasm32")]
use seq_io::fastq::Record as FastqRecord;
#[cfg(target_arch = "wasm32")]
use wasm_bindgen_file_reader::WebSysFile;

use std::cmp::Ordering;

use super::*;

fn unpack_byte(mut byte_packed: u8) -> [u8; 4] {
    let mut out = [0; 4];
    for i in 0..4 {
        out[i] = (byte_packed & 0b11_00_00_00) >> 6;
        byte_packed <<= 2;
    }
    out
}

/// Stores forward and (optionally) reverse complement hashes of k-mers in a nucleotide sequence
#[derive(Debug)]
pub struct NtHashIterator {
    k: usize,
    rc: bool,
    fh: u64,
    rh: Option<u64>,
    index: usize,
    front_unpacked: [u8; 4],
    back_unpacked: [u8; 4],
    seq: Vec<u8>,
    offsets: std::iter::Peekable<std::vec::IntoIter<usize>>,
    acgt: [usize; 4],
    non_acgt: usize,
    reads: bool,
}

impl RollHash for NtHashIterator {
    fn set_k(&mut self, k: usize) {
        if k != self.k {
            self.k = k;
            if self.next_iterator(0).is_none() {
                panic!("K-mer larger than smallest valid sequence");
            }
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
        self.seq.len()
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
    #[cfg(not(target_arch = "wasm32"))]
    /// Creates a new ntHash iterator, by loading DNA sequences into memory
    pub fn new(files: &[String], k: usize, rc: bool, min_qual: u8) -> Vec<Self> {
        // Check if we're working with reads, and initalise the filter if so
        let mut reader_peek = parse_fastx_file(files[0].clone())
            .unwrap_or_else(|_| panic!("Invalid path/file: {}", files[0]));
        let seq_peek = reader_peek
            .next()
            .expect("Invalid FASTA/Q record")
            .expect("Invalid FASTA/Q record");
        let mut reads = false;
        if seq_peek.format() == Format::Fastq {
            reads = true;
            if files.len() > 2 {
                panic!("Input files are reads, but there are more than two input files");
            }
        }

        let mut seq = Vec::new();
        let mut offsets = Vec::new();
        let mut acgt = [0, 0, 0, 0];
        let mut non_acgt = 0;

        // Read sequence into memory (as we go through multiple times)
        log::debug!("Preprocessing sequence");
        for file in files.iter() {
            Self::add_dna_seq(
                file,
                min_qual,
                &mut seq,
                &mut offsets,
                &mut acgt,
                &mut non_acgt,
            );
        }

        let mut hash_it = Self {
            k: 0, // this should be changed to an Option<usize>
            rc,
            fh: 0,
            rh: None,
            index: 0,
            front_unpacked: [0; 4],
            back_unpacked: [0; 4],
            seq,
            offsets: offsets.into_iter().peekable(),
            acgt,
            non_acgt,
            reads,
        };
        hash_it.set_k(k);
        vec![hash_it]
    }

    #[cfg(target_arch = "wasm32")]
    /// Creates a new ntHash iterator, by loading DNA sequences into memory. WASM version.
    pub fn new(
        files: (&web_sys::File, Option<&web_sys::File>),
        rc: bool,
        min_qual: u8,
    ) -> Vec<Self> {
        // Check if we're working with reads, and initalise the filter if so

        let file_name = files.0.name();
        let mut file_type = file_name
            .split('.')
            .nth(file_name.split('.').count() - 1)
            .unwrap();
        if file_type == "gz" {
            file_type = file_name
                .split('.')
                .nth(file_name.split('.').count() - 2)
                .unwrap();
        }

        let is_reads: bool;
        if ["fasta", "fa"].contains(&file_type) {
            is_reads = false;
        } else if ["fastq", "fq"].contains(&file_type) {
            is_reads = true;
        } else {
            panic!("Unsupported file type.")
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
            reads: is_reads,
        };

        // Read sequence into memory (as we go through multiple times)
        log::debug!("Preprocessing sequence");
        hash_it.add_dna_seq(files.0, min_qual);
        if let Some(fileobj) = files.1 {
            hash_it.add_dna_seq(fileobj, min_qual);
        }
        hash_it.seq_len = hash_it.seq.len() - 1;
        vec![hash_it]
    }

    #[cfg(not(target_arch = "wasm32"))]
    fn add_dna_seq(
        filename: &str,
        min_qual: u8,
        seq: &mut Vec<u8>,
        offsets: &mut Vec<usize>,
        acgt: &mut [usize; 4],
        non_acgt: &mut usize,
    ) {
        let mut reader =
            parse_fastx_file(filename).unwrap_or_else(|_| panic!("Invalid path/file: {filename}"));

        let mut b = 0;
        let mut i = 0;

        while let Some(record) = reader.next() {
            let seqrec = record.expect("Invalid FASTA/Q record");

            for (base_idx, base) in seqrec.seq().iter().enumerate() {
                if valid_base(*base) && seqrec.qual().is_none_or(|q| q[base_idx] >= min_qual) {
                    let encoded_base = encode_base(*base);
                    acgt[encoded_base as usize] += 1;
                    b |= encoded_base;
                    i += 1;

                    if i > 3 {
                        i = 0;
                        seq.push(b);
                        b = 0;
                    } else {
                        b <<= 2;
                    }
                } else {
                    *non_acgt += 1;
                    offsets.push(seq.len() * 4 + i);
                }
            }
            offsets.push(seq.len() * 4 + i);
        }

        if i != 0 {
            seq.push(b);
        }
    }

    #[cfg(target_arch = "wasm32")]
    fn add_dna_seq(&mut self, file: &web_sys::File, min_qual: u8) {
        let file_name = file.name();
        let mut file_type = file_name
            .split('.')
            .nth(file_name.split('.').count() - 1)
            .unwrap();
        if file_type == "gz" {
            file_type = file_name
                .split('.')
                .nth(file_name.split('.').count() - 2)
                .unwrap();
        }

        let is_reads: bool;
        if ["fasta", "fa"].contains(&file_type) {
            is_reads = false;
        } else if ["fastq", "fq"].contains(&file_type) {
            is_reads = true;
        } else {
            panic!("Unsupported file type.")
        }

        let mut filetoparse = WebSysFile::new(file.clone());
        if !is_reads {
            let mut reader = open_fasta(&mut filetoparse);
            while let Some(seqrec) = reader.next() {
                for base in seqrec.expect("Invalid FASTA record").seq().iter() {
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
        } else {
            let mut reader = open_fastq(&mut filetoparse);
            while let Some(seqrec) = reader.next() {
                let rec = seqrec.expect("Invalid FASTQ record");
                for (base, qual) in rec.seq().iter().zip(rec.qual()) {
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
            }
        }
    }

    // Only valid if start is greater or equal to previous call
    fn next_iterator(&mut self, mut start: usize) -> Option<()> {
        self.fh = 0_u64;

        let mut end = start + self.k;

        // Find next offset from start
        if let Some(&current_offset) = self.offsets.peek() {
            // Not sure if >= start is needed - playing it safe for now
            if current_offset >= start && current_offset < end {
                start = current_offset;
                end = start + self.k;

                while let Some(next_offset) = self.offsets.next() {
                    if next_offset >= end {
                        break;
                    }
                    start = next_offset;
                    end = start + self.k;
                }
            }
        } else {
            // Already past end of sequence before call
            return None;
        }

        if start >= (self.seq_len() - self.k) {
            // No valid k-mer left
            return None;
        }

        // calculate hash from start..end using u8 range correctly
        let byterange_start = start / 4;
        let byterange_end = end / 4;
        let base_start = start % 4;
        let base_end = end % 4 + 4 * (byterange_end - byterange_start);

        let mut bases = self.seq[byterange_start..=byterange_end].iter();
        let mut curr_byte = 0;
        let mut rev_bases = Vec::with_capacity(self.k);
        for base_idx in 0..=base_end {
            if base_idx % 4 == 0 {
                curr_byte = *(bases.next().unwrap());
            }
            let base = (curr_byte & 0b11_00_00_00) >> 6;
            // Store parsed bases for rc_hash below
            if self.rc {
                rev_bases.push(base);
            }
            curr_byte <<= 2;
            if base_idx > base_start {
                // hash
                self.fh = self.fh.rotate_left(1_u32);
                self.fh = swapbits033(self.fh);
                self.fh ^= nthash_tables::HASH_LOOKUP[base as usize];
            }
        }

        self.rh = if self.rc {
            let mut h = 0_u64;
            for b in rev_bases.iter().rev() {
                h = h.rotate_left(1_u32);
                h = swapbits033(h);
                h ^= nthash_tables::RC_HASH_LOOKUP[*b as usize];
            }
            Some(h)
        } else {
            None
        };

        // Set the bytes at start and end of roll for use in roll_fwd
        self.front_unpacked = unpack_byte(self.seq[byterange_end]);
        self.back_unpacked = unpack_byte(self.seq[byterange_start]);

        self.index = end;

        Some(())
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
        match self.index.cmp(&self.seq_len()) {
            Ordering::Less => {
                // Get the current hash to return
                let current = self.curr_hash();

                // Calculate the next hash, which will be returned on the next call
                let new_base = self.front_unpacked[self.index % 4];
                //println!(
                //"self.index = {}, self.k = {}, self.index + 1 = {}",
                //self.index,
                //self.k,
                //self.index + 1,
                //);
                let old_base = self.back_unpacked[((self.index + 1) - self.k) % 4];
                self.roll_fwd(old_base, new_base);

                // Move to the next valid base, which prepares for the next call to
                // calculate the next hash, which will be returned in two calls time
                self.index += 1;
                if self.index == *(self.offsets.peek().unwrap()) {
                    if self.next_iterator(self.index).is_none() {
                        self.index = self.seq_len();
                    }
                }
                // Moved to a new byte at the front - update front_unpacked
                if self.index % 4 == 0 {
                    self.front_unpacked = unpack_byte(self.seq[self.index / 4]);
                }
                // Moved to a new byte at the back - update back_unpacked
                if (self.index - self.k + 1) % 4 == 0 {
                    self.back_unpacked = unpack_byte(self.seq[(self.index - self.k + 1) % 4]);
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
}
