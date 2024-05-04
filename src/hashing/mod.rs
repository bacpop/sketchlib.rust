use std::{borrow::Cow, cmp::Ordering};

use serde::{Deserialize, Serialize};
use clap::ValueEnum;

mod nthash_tables;

#[derive(Clone, Debug, Serialize, Deserialize, ValueEnum)]
pub enum HashType {
    DNA,
    AA,
    Structure,
}

/// Table from bits 0-3 to ASCII (use [`decode_base()`] not this table).
const LETTER_CODE: [u8; 4] = [b'A', b'C', b'T', b'G'];

/// Encode an ASCII char to bits 0-3.
#[inline(always)]
pub fn encode_base(base: u8) -> u8 {
    (base >> 1) & 0x3
}

/// Decode bits 0-3 to ASCII.
#[inline(always)]
pub fn decode_base(bitbase: u8) -> u8 {
    LETTER_CODE[bitbase as usize]
}

/// Reverse complement an encoded base.
#[inline(always)]
pub fn rc_base(base: u8) -> u8 {
    base ^ 2
}

#[inline(always)]
pub fn valid_base(mut base: u8) -> bool {
    base |= 0x20; // to lower
    matches!(base, b'a' | b'c' | b'g' | b't')
}

#[inline(always)]
pub fn swapbits033(v: u64) -> u64 {
    let x = (v ^ (v >> 33)) & 1;
    v ^ (x | (x << 33))
}

#[inline(always)]
pub fn swapbits3263(v: u64) -> u64 {
    let x = ((v >> 32) ^ (v >> 63)) & 1;
    v ^ ((x << 32) | (x << 63))
}

// TODO aaHash for proteins

// TODO generic hash for structure alphabet
// https://github.com/eldruin/wyhash-rs

pub trait RollHash: Iterator<Item = u64> {
    // fn new(seq: &'a [u8], k: usize, rc: bool) -> Self;
    fn curr_hash(&self) -> u64;
    fn hash_type(&self) -> HashType;
}

/// Stores forward and (optionally) reverse complement hashes of k-mers in a nucleotide sequence
#[derive(Debug)]
pub struct NtHashIterator<'a> {
    k: usize,
    rc: bool,
    fh: u64,
    rh: Option<u64>,
    index: usize,
    seq: Cow<'a, [u8]>,
    seq_len: usize, // NB seq.len() - 1 due to terminating char
}

impl<'a> RollHash for NtHashIterator<'a> {
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
}

impl<'a> NtHashIterator<'a> {
    /// Creates a new iterator over a sequence with a given k-mer size
    pub fn new(seq: &'a [u8], k: usize, rc: bool) -> Self {
        if let Some(new_it) = Self::new_iterator(0, seq, k, rc) {
            Self {
                k,
                rc,
                fh: new_it.0,
                rh: new_it.1,
                index: new_it.2,
                seq: Cow::Borrowed(seq),
                seq_len: seq.len() - 1,
            }
        } else {
            panic!("K-mer larger than smallest valid sequence");
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

impl<'a> Iterator for NtHashIterator<'a> {
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
impl<'a> ExactSizeIterator for NtHashIterator<'a> {}
