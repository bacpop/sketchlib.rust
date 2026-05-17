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
    for byte_val in &mut out {
        *byte_val = (byte_packed & 0b11_00_00_00) >> 6;
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
    // Bases at the front of the k-mer
    front_unpacked: [u8; 4],
    // Bases at the end of the k-mer
    back_unpacked: [u8; 4],
    // Bitpacked seq
    seq: Vec<u8>,
    // Nucleotide positions of N's and record boundaries (monotonically increasing).
    offsets: Vec<usize>,
    // Index into `offsets` for the current k-value pass; reset by set_k().
    offset_idx: usize,
    acgt: [usize; 4],
    non_acgt: usize,
    reads: bool,
}

impl RollHash for NtHashIterator {
    fn set_k(&mut self, k: usize) {
        if k != self.k {
            self.k = k;
            self.offset_idx = 0; // rewind: offsets must be re-traversed for each k
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
            k: 0, // TODO neater if an Option<usize> rather than using 0 as a guard
            rc,
            fh: 0,
            rh: None,
            index: 0,
            front_unpacked: [0; 4],
            back_unpacked: [0; 4],
            seq,
            offsets,
            offset_idx: 0,
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
        k: usize,
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
            front_unpacked: [0; 4],
            back_unpacked: [0; 4],
            seq: Vec::new(),
            offsets: Vec::new(),
            offset_idx: 0,
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
        hash_it.set_k(k);
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

        // Current byte
        let mut b = 0;
        // Nucleotide index within byte 0-3
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
            // Left-align the partial byte so bases are at the MSB positions,
            // matching the layout of full bytes and the unpack_byte function.
            seq.push(b << ((3 - i) * 2));
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
        let mut b = 0;
        let mut i = 0;
        if !is_reads {
            let mut reader = open_fasta(&mut filetoparse);
            while let Some(seqrec) = reader.next() {
                for base in seqrec.expect("Invalid FASTA record").seq().iter() {
                    self.pack_dna_base(*base, true, &mut b, &mut i);
                }
                self.offsets.push(self.seq.len() * 4 + i);
            }
        } else {
            let mut reader = open_fastq(&mut filetoparse);
            while let Some(seqrec) = reader.next() {
                let rec = seqrec.expect("Invalid FASTQ record");
                for (base, qual) in rec.seq().iter().zip(rec.qual()) {
                    self.pack_dna_base(*base, *qual >= min_qual, &mut b, &mut i);
                }
                self.offsets.push(self.seq.len() * 4 + i);
            }
        }

        if i != 0 {
            self.seq.push(b << ((3 - i) * 2));
        }
    }

    #[cfg(target_arch = "wasm32")]
    fn pack_dna_base(&mut self, base: u8, quality_ok: bool, b: &mut u8, i: &mut usize) {
        if quality_ok && valid_base(base) {
            let encoded_base = encode_base(base);
            self.acgt[encoded_base as usize] += 1;
            *b |= encoded_base;
            *i += 1;

            if *i > 3 {
                *i = 0;
                self.seq.push(*b);
                *b = 0;
            } else {
                *b <<= 2;
            }
        } else {
            self.non_acgt += 1;
            self.offsets.push(self.seq.len() * 4 + *i);
        }
    }

    // Only valid if start is greater or equal to previous call
    fn next_iterator(&mut self, mut start: usize) -> Option<()> {
        self.fh = 0_u64;

        let mut end = start + self.k;

        // Find next offset from start
        if self.offset_idx >= self.offsets.len() {
            // Already past end of sequence before call
            return None;
        }
        // Advance past offsets that fall within the current window, moving start
        // past each N/boundary until we find a window with no offset inside it.
        while let Some(&current_offset) = self.offsets.get(self.offset_idx) {
            if current_offset < start || current_offset >= end {
                break;
            }
            self.offset_idx += 1;
            start = current_offset;
            end = start + self.k;
        }

        if start + self.k > self.seq_len() {
            // No valid k-mer left
            return None;
        }

        // calculate hash from start..end using u8 range correctly
        let byterange_start = start / 4;
        let byterange_end = end / 4;
        let base_start = start % 4;

        // Pre-shift first byte so the target start position is at the MSB
        let mut curr_byte = self.seq[byterange_start] << (base_start * 2);
        let mut byte_idx = byterange_start;
        let mut rev_bases = Vec::with_capacity(self.k);

        for hash_idx in 0..self.k {
            let byte_pos = (start + hash_idx) / 4;
            if byte_pos != byte_idx {
                byte_idx = byte_pos;
                curr_byte = self.seq[byte_idx];
            }
            let base = (curr_byte & 0b11_00_00_00) >> 6;
            if self.rc {
                rev_bases.push(base);
            }
            curr_byte <<= 2;
            self.fh = self.fh.rotate_left(1_u32);
            self.fh = swapbits033(self.fh);
            self.fh ^= nthash_tables::HASH_LOOKUP[base as usize];
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

    // Helper function to generate iterator from sequence string in tests
    #[cfg(test)]
    fn from_seq(seq_str: &str, k: usize, rc: bool) -> Self {
        let mut seq: Vec<u8> = Vec::new();
        let mut offsets: Vec<usize> = Vec::new();
        let mut acgt = [0usize; 4];
        let mut non_acgt = 0usize;
        let mut b: u8 = 0;
        let mut i: usize = 0;

        for ch in seq_str.bytes() {
            if valid_base(ch) {
                let e = encode_base(ch);
                acgt[e as usize] += 1;
                b |= e;
                i += 1;
                if i > 3 {
                    i = 0;
                    seq.push(b);
                    b = 0;
                } else {
                    b <<= 2;
                }
            } else {
                non_acgt += 1;
                offsets.push(seq.len() * 4 + i);
            }
        }
        offsets.push(seq.len() * 4 + i); // end-of-record sentinel
        if i != 0 {
            seq.push(b << ((3 - i) * 2));
        }

        let mut hash_it = Self {
            k: 0,
            rc,
            fh: 0,
            rh: None,
            index: 0,
            front_unpacked: [0; 4],
            back_unpacked: [0; 4],
            seq,
            offsets,
            offset_idx: 0,
            acgt,
            non_acgt,
            reads: false,
        };
        hash_it.set_k(k);
        hash_it
    }
}

impl Iterator for NtHashIterator {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        match self.index.cmp(&self.seq_len()) {
            Ordering::Less => {
                // Get the current hash to return
                let current = self.curr_hash();

                // If the next base is an offset (N or record boundary), reinitialize
                // from that position rather than rolling forward (which would corrupt
                // the pre-computed hash for the k-mer we're about to return).
                if self
                    .offsets
                    .get(self.offset_idx)
                    .is_some_and(|&off| self.index == off)
                {
                    if self.next_iterator(self.index).is_none() {
                        // next_iterator reset fh to 0; skip Equal (which would emit a
                        // spurious zero hash) and jump straight to Greater.
                        self.index = self.seq_len() + 1;
                    }
                } else {
                    // Roll forward to compute the hash for the next k-mer
                    let new_base = self.front_unpacked[self.index % 4];
                    let old_base = self.back_unpacked[(self.index - self.k) % 4];
                    self.roll_fwd(old_base, new_base);
                    self.index += 1;
                    // Moved to a new byte at the front - update front_unpacked
                    if self.index.is_multiple_of(4) && self.index < self.seq_len() {
                        self.front_unpacked = unpack_byte(self.seq[self.index / 4]);
                    }
                    // Moved to a new byte at the back - update back_unpacked.
                    // old_base for the next roll is at position index-k; reload when that
                    // position is the first base of a new packed byte.
                    if (self.index - self.k).is_multiple_of(4) {
                        self.back_unpacked = unpack_byte(self.seq[(self.index - self.k) / 4]);
                    }
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

#[cfg(test)]
mod tests {
    use super::super::nthash_tables;
    use super::*;

    // These are reference implementations for how the rolling works in ntHash, and
    // in v0.2.4 of sketchlib. The bitpacking is quite tricky, so this is essentially
    // an internal consistency test which could be removed in future versions
    fn ref_new_iterator(
        mut start: usize,
        seq: &[u8],
        k: usize,
        rc: bool,
    ) -> Option<(u64, Option<u64>, usize)> {
        let mut fh = 0_u64;
        'outer: while start + k <= seq.len() {
            for (i, v) in seq[start..(start + k)].iter().enumerate() {
                if *v > 3 {
                    start += i + 1;
                    fh = 0;
                    continue 'outer;
                }
                fh = fh.rotate_left(1);
                fh = swapbits033(fh);
                fh ^= nthash_tables::HASH_LOOKUP[*v as usize];
            }
            break 'outer;
        }
        if start + k > seq.len() {
            return None;
        }
        let rh = if rc {
            let mut h = 0_u64;
            for v in seq[start..(start + k)].iter().rev() {
                h = h.rotate_left(1);
                h = swapbits033(h);
                h ^= nthash_tables::RC_HASH_LOOKUP[*v as usize];
            }
            Some(h)
        } else {
            None
        };
        Some((fh, rh, start + k))
    }

    fn ref_roll_fwd(fh: &mut u64, rh: &mut Option<u64>, old_base: u8, new_base: u8, k: usize) {
        *fh = fh.rotate_left(1);
        *fh = swapbits033(*fh);
        *fh ^= nthash_tables::HASH_LOOKUP[new_base as usize];
        *fh ^= nthash_tables::MS_TAB_31L[(old_base as usize * 31) + (k % 31)]
            | nthash_tables::MS_TAB_33R[(old_base as usize) * 33 + (k % 33)];
        if let Some(rev) = rh {
            let mut h = *rev
                ^ (nthash_tables::MS_TAB_31L[(rc_base(new_base) as usize * 31) + (k % 31)]
                    | nthash_tables::MS_TAB_33R[(rc_base(new_base) as usize) * 33 + (k % 33)]);
            h ^= nthash_tables::RC_HASH_LOOKUP[old_base as usize];
            h = h.rotate_right(1);
            h = swapbits3263(h);
            *rh = Some(h);
        }
    }

    fn ref_hashes(seq_str: &str, k: usize, rc: bool) -> Vec<u64> {
        let mut seq: Vec<u8> = seq_str
            .bytes()
            .map(|b| {
                if valid_base(b) {
                    encode_base(b)
                } else {
                    SEQSEP
                }
            })
            .collect();
        seq.push(SEQSEP); // end-of-record sentinel
        let seq_len = seq.len() - 1;
        let mut out = Vec::new();
        let Some((mut fh, mut rh, mut index)) = ref_new_iterator(0, &seq, k, rc) else {
            return out;
        };
        let canonical = |fh: u64, rh: Option<u64>| rh.map_or(fh, |r| u64::min(fh, r));
        loop {
            match index.cmp(&seq_len) {
                std::cmp::Ordering::Less => {
                    out.push(canonical(fh, rh));
                    let nb = seq[index];
                    if nb > 3 {
                        match ref_new_iterator(index + 1, &seq, k, rc) {
                            Some((f, r, i)) => {
                                fh = f;
                                rh = r;
                                index = i;
                            }
                            None => break,
                        }
                    } else {
                        ref_roll_fwd(&mut fh, &mut rh, seq[index - k], nb, k);
                        index += 1;
                    }
                }
                std::cmp::Ordering::Equal => {
                    out.push(canonical(fh, rh));
                    break;
                }
                std::cmp::Ordering::Greater => break,
            }
        }
        out
    }

    #[test]
    /// set_k() must rewind the offsets so multiple k values each see all N positions.
    /// This is the key multi-k regression: the sketch command calls set_k() once per k value
    /// using the same iterator, so offsets must not be consumed across k-value passes.
    fn nthash_bitpacked_offsets_rewind_across_k_values() {
        // All three segments are ≥ 7 bases so every k value has valid k-mers.
        let seq = "ACGTACGTANACGTACGTNNTACGTACGT";
        for k in [3usize, 5, 7] {
            let expected = ref_hashes(seq, k, true);
            let actual: Vec<u64> = NtHashIterator::from_seq(seq, k, true).collect();
            assert_eq!(actual.len(), expected.len(), "count mismatch at k={k}");
            assert_eq!(actual, expected, "hash mismatch at k={k}");
        }
        // Simulate the sketch command: one iterator, multiple set_k calls in sequence.
        let mut it = NtHashIterator::from_seq(seq, 3, true);
        let hashes_k3: Vec<u64> = it.by_ref().collect();
        it.set_k(5);
        let hashes_k5: Vec<u64> = it.by_ref().collect();
        it.set_k(7);
        let hashes_k7: Vec<u64> = it.by_ref().collect();
        assert_eq!(hashes_k3, ref_hashes(seq, 3, true), "k=3 hashes wrong");
        assert_eq!(
            hashes_k5,
            ref_hashes(seq, 5, true),
            "k=5 hashes wrong after set_k(5)"
        );
        assert_eq!(
            hashes_k7,
            ref_hashes(seq, 7, true),
            "k=7 hashes wrong after set_k(7)"
        );
    }

    #[test]
    fn nthash_bitpacked_matches_reference_50bp() {
        // Use actual test sequences from distance test fixtures
        let seq = "CTAGGGCCCTTTCCCGGATATAAACGCCAGGTTGAATCCGCATTTGGAGG";
        for k in [3usize, 17, 31] {
            let expected = ref_hashes(seq, k, true);
            let actual: Vec<u64> = NtHashIterator::from_seq(seq, k, true).collect();
            assert_eq!(
                actual.len(),
                expected.len(),
                "Hash count differs for k={k}: got {}, expected {}",
                actual.len(),
                expected.len()
            );
            assert_eq!(actual, expected, "Hashes differ for k={k}");
        }
    }

    #[test]
    fn nthash_bitpacked_matches_reference_no_n() {
        let seq = "ACGTACGTACGT";
        let k = 4;
        let expected = ref_hashes(seq, k, false);
        assert!(!expected.is_empty());
        let actual: Vec<u64> = NtHashIterator::from_seq(seq, k, false).collect();
        assert_eq!(actual, expected, "k={k} forward hashes differ");
    }

    #[test]
    /// N within the last k-1 positions of a segment must not emit a spurious zero hash.
    fn nthash_bitpacked_no_spurious_hash_after_terminal_n() {
        // k=4: N at position 6 is within the last 3 positions of "ACGTACG" segment.
        // A naive impl resets fh=0 then emits it via the Equal branch.
        let seq = "ACGTACGNACGT";
        for k in [4usize, 5] {
            let expected = ref_hashes(seq, k, true);
            let actual: Vec<u64> = NtHashIterator::from_seq(seq, k, true).collect();
            assert_eq!(
                actual.len(),
                expected.len(),
                "spurious hash count for k={k}"
            );
            assert_eq!(actual, expected, "spurious hash value for k={k}");
        }
    }

    #[test]
    fn nthash_bitpacked_matches_reference_with_n() {
        let seq = "ACGTANACGT";
        let k = 4;
        let expected = ref_hashes(seq, k, false);
        assert!(!expected.is_empty());
        let actual: Vec<u64> = NtHashIterator::from_seq(seq, k, false).collect();
        assert_eq!(actual, expected, "k={k} hashes with N differ");
    }

    #[test]
    fn nthash_bitpacked_matches_reference_rc() {
        let seq = "ACGTACGTACGT";
        let k = 4;
        let expected = ref_hashes(seq, k, true);
        let actual: Vec<u64> = NtHashIterator::from_seq(seq, k, true).collect();
        assert_eq!(actual, expected, "k={k} RC hashes differ");
    }
}
