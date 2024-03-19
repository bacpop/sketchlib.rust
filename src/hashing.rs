use std::borrow::Cow;

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

const HASH_LOOKUP: [u64; 4] = [
    0x3c8b_fbb3_95c6_0474,
    0x3193_c185_62a0_2b4c,
    0x2955_49f5_4be2_4456,
    0x2032_3ed0_8257_2324,
];
const RC_HASH_LOOKUP: [u64; 4] = [
    0x2955_49f5_4be2_4456,
    0x2032_3ed0_8257_2324,
    0x3c8b_fbb3_95c6_0474,
    0x3193_c185_62a0_2b4c,
];

// TODO aaHash for proteins

// TODO generic hash for structure alphabet

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

impl<'a> NtHashIterator<'a> {
    /// Creates a new iterator over a sequence with a given k-mer size
    pub fn new(seq: &'a [u8], k: usize, rc: bool) -> NtHashIterator {
        if let Some(new_it) = Self::new_iterator(0, seq, k, rc) {
            Self {
                k,
                rc,
                fh: new_it.0,
                rh: new_it.1,
                index: new_it.2 + k,
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
        let mut fh = 0;
        while start < (seq.len() - k) {
            for (i, v) in seq[start..(start + k)].iter().enumerate() {
                // If invalid seq
                if *v > 3 {
                    start += i + 1;
                    if start >= seq.len() {
                        return None;
                    }
                    fh = 0;
                    break;
                }
                fh ^= HASH_LOOKUP[*v as usize].rotate_left((k - i - 1) as u32);
            }
            break; // success
        }
        if start >= (seq.len() - k) {
            return None;
        }

        let rh = if rc {
            let mut h = 0;
            for (i, v) in seq[start..(start + k)].iter().rev().enumerate() {
                h ^= RC_HASH_LOOKUP[*v as usize].rotate_left((k - i - 1) as u32);
            }
            Some(h)
        } else {
            None
        };

        Some((fh, rh, start + k + 1))
    }

    /// Move to the next k-mer by adding a new base, removing a base from the end, efficiently updating the hash.
    fn roll_fwd(&mut self, old_base: u8, new_base: u8) {
        self.fh = self.fh.rotate_left(1)
            ^ HASH_LOOKUP[old_base as usize].rotate_left(self.k as u32)
            ^ HASH_LOOKUP[new_base as usize];

        if let Some(rev) = self.rh {
            self.rh = Some(
                rev.rotate_right(1)
                    ^ RC_HASH_LOOKUP[old_base as usize].rotate_right(1)
                    ^ RC_HASH_LOOKUP[new_base as usize].rotate_left(self.k as u32 - 1),
            )
        };
    }

    /// Retrieve the current hash (minimum of forward and reverse complement hashes)
    pub fn curr_hash(&self) -> u64 {
        if let Some(rev) = self.rh {
            u64::min(self.fh, rev)
        } else {
            self.fh
        }
    }
}

impl<'a> Iterator for NtHashIterator<'a> {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index < self.seq_len {
            let current = self.curr_hash();
            let new_base = self.seq[self.index];
            // Restart hash if invalid base
            if new_base > 3 {
                if let Some(new_it) = Self::new_iterator(self.index + 1, &self.seq, self.k, self.rc) {
                    self.fh = new_it.0;
                    self.rh = new_it.1;
                    self.index = new_it.2;
                } else {
                    self.index = self.seq_len; // End of valid sequence
                }
            } else {
                self.roll_fwd(self.seq[self.index - self.k - 1], new_base);
                self.index += 1;
            }
            Some(current)
        } else {
            None // End of sequence
        }
    }
}
