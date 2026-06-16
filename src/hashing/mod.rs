//! [nthash](https://github.com/bcgsc/ntHash) and aahash iterators
use clap::ValueEnum;
use serde::{Deserialize, Serialize};

#[cfg(not(target_arch = "wasm32"))]
pub mod aahash_iterator;
#[cfg(not(target_arch = "wasm32"))]
mod aahash_tables;
pub mod bloom_filter;
pub mod nthash_iterator;
mod nthash_tables;

/// Character to use for invalid nucleotides
pub const SEQSEP: u8 = 5;
/// Default aaHash 'level'
pub const DEFAULT_LEVEL: AaLevel = AaLevel::Level1;

/// aaHash levels
#[derive(Clone, Debug, PartialEq, PartialOrd, Serialize, Deserialize, ValueEnum)]
pub enum AaLevel {
    /// Level1: All amino acids are different
    Level1,
    /// Level2: Groups T,S; D,E; Q,K,R; V,I,L,M; W,F,Y
    Level2,
    /// Level3: Additionally groups A with T,S; N with D,E
    Level3,
}

/// Type of sequence hashed
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, PartialOrd, Default)]
pub enum HashType {
    #[default]
    /// Nucleotides
    DNA,
    /// Amino acids at set [`AaLevel`]
    AA(AaLevel),
    /// Structures using 3di embedding
    PDB,
}

// TODO: for PDB need to use 3Di sequences. These are from an embedding
// With foldseek, can run:
// foldseek createdb 5uak.pdb 5Uak_DB
// foldseek lndb 5uak_DB 5uak_DB_ss_h
// foldseek convert2fasta 5uak_DB_ss_h 5uak.fasta
// The Cpp code looks like a bit of a pain to build in:
// https://github.com/steineggerlab/foldseek/blob/master/lib/3di/structureto3di.cpp
// This python port is an alternative:
// https://github.com/althonos/mini3di
// (this seems to give different results to the above)
// Seems pretty easy to run python from within rust:
// https://pyo3.rs/v0.21.2/python-from-rust/calling-existing-code
// Or a language model:
// https://github.com/mheinzinger/ProstT5

// NB: this is needed because ValueEnum (for clap) only works with non-unit types
// So here set a default for the level and set it properly later (in lib.rs)
impl clap::ValueEnum for HashType {
    fn value_variants<'a>() -> &'a [Self] {
        &[
            HashType::DNA,
            HashType::AA(DEFAULT_LEVEL),
            // HashType::AA(AaLevel::Level2),
            // HashType::AA(AaLevel::Level3),
            HashType::PDB,
        ]
    }
    fn to_possible_value<'a>(&self) -> ::std::option::Option<clap::builder::PossibleValue> {
        match self {
            Self::DNA => Some(clap::builder::PossibleValue::new("dna")),
            Self::AA(_) => Some(clap::builder::PossibleValue::new("aa")),
            // Self::AA(AaLevel::Level1) => Some(clap::builder::PossibleValue::new("aa_1")),
            // Self::AA(AaLevel::Level2) => Some(clap::builder::PossibleValue::new("aa_2")),
            // Self::AA(AaLevel::Level3) => Some(clap::builder::PossibleValue::new("aa_3")),
            Self::PDB => Some(clap::builder::PossibleValue::new("pdb")),
        }
    }
}

/// Encode an ASCII char to bits 0-3.
#[inline(always)]
fn encode_base(base: u8) -> u8 {
    (base >> 1) & 0x3
}

/// Reverse complement an encoded base.
#[inline(always)]
fn rc_base(base: u8) -> u8 {
    base ^ 2
}

/// A valid base a,c,g,t,u/A,C,G,T,U
#[inline(always)]
fn valid_base(mut base: u8) -> bool {
    base |= 0x20; // to lower
    matches!(base, b'a' | b'c' | b'g' | b't' | b'u')
}

#[inline(always)]
fn swapbits033(v: u64) -> u64 {
    let x = (v ^ (v >> 33)) & 1;
    v ^ (x | (x << 33))
}

#[inline(always)]
fn swapbits3263(v: u64) -> u64 {
    let x = ((v >> 32) ^ (v >> 63)) & 1;
    v ^ ((x << 32) | (x << 63))
}

// TODO generic hash for structure alphabet
// https://github.com/eldruin/wyhash-rs

/// Rolling functions supported by both ntHash and aaHash
pub trait RollHash: Iterator<Item = u64> {
    /// Set the k-mer size
    fn set_k(&mut self, k: usize);
    /// Get the current hash
    fn curr_hash(&self) -> u64;
    /// The type of sequence being hashed
    fn hash_type(&self) -> HashType;
    /// Total length of the sequence
    fn seq_len(&self) -> usize;
    /// The underlying sequence as bytes
    fn seq(&self) -> &Vec<u8>;
    /// Sequence metadata (reads, \[a,c,g,t\] counts, non-acgt bases)
    fn sketch_data(&self) -> (bool, [usize; 4], usize);

    /// An iterator over the hashes
    fn iter(&mut self) -> Box<dyn Iterator<Item = u64> + '_> {
        Box::new(self)
    }

    /// Whether the underlying sequences is fastq/read data
    fn reads(&self) -> bool {
        false
    }
}
