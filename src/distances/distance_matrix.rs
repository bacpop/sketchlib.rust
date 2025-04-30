//! Functions and traits for calculating and storing distances
use std::cmp::Ordering;
use std::fmt;

// use ordered_float::NotNan;

use crate::sketch::multisketch::MultiSketch;

/// Index k in long form, given ref vs ref sample pair (i, j) and total samples n in square form
#[inline(always)]
pub fn square_to_condensed(i: usize, j: usize, n: usize) -> usize {
    debug_assert!(j > i);
    n * i - ((i * (i + 1)) >> 1) + j - 1 - i
}

/// Index k in long form, given ref vs query sample pair (i, j) and total samples n in rectangular form
#[inline(always)]
pub fn ref_query_index(i: usize, j: usize, n: usize) -> usize {
    debug_assert!(j > i);
    i * n + j
}

/// Indexes (i, j) in rectangular form, given k and total samples n in long form
#[inline(always)]
pub fn calc_query_indices(k: usize, n: usize) -> (usize, usize) {
    let i = k / n;
    let j = k % n;
    debug_assert!(i < n);
    debug_assert!(j < n);
    (i, j)
}

/// Indexes j in square form, given k, index i (use [`calc_row_idx`]) and total samples n in long form
#[inline(always)]
pub fn calc_col_idx(k: usize, i: usize, n: usize) -> usize {
    debug_assert!(i < n);
    let k_i64 = k as i64;
    let i_i64 = i as i64;
    let n_i64 = n as i64;
    (k_i64 + i_i64 + 1 - n_i64 * (n_i64 - 1) / 2 + (n_i64 - i_i64) * ((n_i64 - i_i64) - 1) / 2)
        as usize
}

/// Indexes i in square form, given k and total samples n in long form
#[inline(always)]
pub fn calc_row_idx(k: usize, n: usize) -> usize {
    let k_i64 = k as i64;
    let n_i64 = n as i64;
    n - 2
        - (((-8 * k_i64 + 4 * n_i64 * (n_i64 - 1) - 7) as f64).sqrt() / 2.0 - 0.5).floor() as usize
}

/// Types of distance, single k-mer (Jaccard) or multi-k (core/accessory)
#[derive(PartialEq, PartialOrd)]
pub enum DistType {
    /// Jaccard distance (k-mer index, k-mer size, ANI on/off)
    Jaccard(usize, f32, bool),
    /// Core and accessory distances
    CoreAcc,
}

impl fmt::Display for DistType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            DistType::CoreAcc => write!(f, "Distances: core/accessory regression"),
            DistType::Jaccard(_, k, ani) => {
                let k = k as usize;
                if ani {
                    write!(f, "Distances: ANI at k={k}")
                } else {
                    write!(f, "Distances: Jaccard distances at k={k}")
                }
            }
        }
    }
}

/// Metadata shared by all distance matrices
pub trait Distances<'a> {
    /// Distance type
    fn jaccard(&self) -> &DistType;

    /// Whether calculating ANI
    fn ani(&self) -> bool {
        match self.jaccard() {
            DistType::Jaccard(_, _, ani_on) => *ani_on,
            _ => false,
        }
    }

    /// If calcualting Jaccard distances, the index of the k-mer, and the k-mer as a float
    fn k_vals(&self) -> Option<(usize, f32)> {
        match self.jaccard() {
            DistType::Jaccard(k_idx, k_val, _) => Some((*k_idx, *k_val)),
            _ => None,
        }
    }

    /// Number of distance columns in the output
    fn n_dist_cols(&self) -> usize {
        match self.jaccard() {
            DistType::CoreAcc => 2,
            DistType::Jaccard(_, _, _) => 1,
        }
    }

    /// Names of the sketch files in the distance matrix
    fn sketch_names(sketches: &'a MultiSketch) -> Vec<&'a str> {
        let n_samples = sketches.number_samples_loaded();
        let mut names = Vec::with_capacity(n_samples);
        for idx in 0..n_samples {
            names.push(sketches.sketch_name(idx));
        }
        names
    }
}

/// A dense distance matrix in long form, which can represent ref vs ref
/// or ref vs query depending on whether `query_names` is set
pub struct DistanceMatrix<'a> {
    /// The total number of distances (upper triangle only)
    pub n_distances: usize,
    jaccard: DistType,
    distances: Vec<f32>,
    ref_names: Vec<&'a str>,
    query_names: Option<Vec<&'a str>>,
}

impl<'a> DistanceMatrix<'a> {
    /// Create a new distance matrix for the given [`MultiSketch`] objects
    /// with the parameters set by the [`DistType`]
    pub fn new(
        ref_sketches: &'a MultiSketch,
        query_sketches: Option<&'a MultiSketch>,
        jaccard: DistType,
    ) -> Self {
        let n_distances;
        let query_names = if let Some(query) = query_sketches {
            n_distances = ref_sketches.number_samples_loaded() * query.number_samples_loaded();
            Some(Self::sketch_names(query))
        } else {
            n_distances = ref_sketches.number_samples_loaded()
                * (ref_sketches.number_samples_loaded() - 1)
                / 2;
            None
        };

        // Pre-allocate distances
        let mut distances = vec![0.0; n_distances];
        if jaccard == DistType::CoreAcc {
            distances.append(&mut vec![0.0; n_distances]);
        }

        Self {
            n_distances,
            distances,
            ref_names: Self::sketch_names(ref_sketches),
            query_names,
            jaccard,
        }
    }

    /// Reference to the underlying distances which can be written to
    pub fn dists_mut(&mut self) -> &mut Vec<f32> {
        &mut self.distances
    }
}

impl<'a> Distances<'a> for DistanceMatrix<'a> {
    fn jaccard(&self) -> &DistType {
        &self.jaccard
    }
}

impl fmt::Display for DistanceMatrix<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut dist_idx = 0;
        if let Some(queries) = &self.query_names {
            for ref_name in &self.ref_names {
                for query_name in queries {
                    write!(f, "{ref_name}\t{query_name}\t{}", self.distances[dist_idx])?;
                    if self.jaccard == DistType::CoreAcc {
                        write!(f, "\t{}", self.distances[dist_idx + 1])?;
                        dist_idx += 1;
                    }
                    writeln!(f)?;
                    dist_idx += 1;
                }
            }
        } else {
            for (i, ref_name) in self.ref_names.iter().enumerate() {
                for j in (i + 1)..self.ref_names.len() {
                    write!(
                        f,
                        "{ref_name}\t{}\t{}",
                        self.ref_names[j], self.distances[dist_idx]
                    )?;
                    if self.jaccard == DistType::CoreAcc {
                        write!(f, "\t{}", self.distances[dist_idx + 1])?;
                        dist_idx += 1;
                    }
                    writeln!(f)?;
                    dist_idx += 1;
                }
            }
        }
        Ok(())
    }
}

/// Sparse distance struct which contains the index of the query sample
/// and the Jaccard/ANI distance
#[derive(Debug, Clone)]
pub struct SparseJaccard(pub usize, pub f32);
impl Ord for SparseJaccard {
    fn cmp(&self, other: &Self) -> Ordering {
        self.1.partial_cmp(&other.1).unwrap()
                                              // Could also use
                                              /*
                                              NotNan::new(other.1)
                                                  .unwrap()
                                                  .cmp(&NotNan::new(self.1).unwrap())
                                              */
    }
}
impl PartialOrd for SparseJaccard {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl PartialEq for SparseJaccard {
    fn eq(&self, other: &Self) -> bool {
        self.1 == other.1
    }
}
impl Eq for SparseJaccard {}

/// Sparse distance struct for a single distance entry which contains the
/// index of the query sample and the core and accessory distances.

// TODO: could either change the field to compare on, or add Euclidean dists
#[derive(Debug, Clone)]
pub struct SparseCoreAcc(pub usize, pub f32, pub f32);

impl Ord for SparseCoreAcc {
    fn cmp(&self, other: &Self) -> Ordering {
        self.1.partial_cmp(&other.1).unwrap()
        // Could also use
        /*
        NotNan::new(self.1)
            .unwrap()
            .cmp(&NotNan::new(other.1).unwrap())
        */
    }
}
impl PartialOrd for SparseCoreAcc {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl PartialEq for SparseCoreAcc {
    fn eq(&self, other: &Self) -> bool {
        self.1 == other.1
    }
}
impl Eq for SparseCoreAcc {}

/// Underlying distance objects for sparse distances which wrap distance and
/// query index for both types of distance
pub enum DistVec {
    /// Jaccard distances (and index)
    Jaccard(Vec<SparseJaccard>),
    /// Core-accessory distances (and index)
    CoreAcc(Vec<SparseCoreAcc>),
}

/// A sparse distance matrix with a maximum of `knn` distances for each sample
pub struct SparseDistanceMatrix<'a> {
    /// Total number of distances
    pub n_distances: usize,
    /// Maximum number of distances kept per sample: k smallest distances
    pub knn: usize,
    jaccard: DistType,
    distances: DistVec,
    ref_names: Vec<&'a str>,
}

impl<'a> SparseDistanceMatrix<'a> {
    /// Create a new sparse distance matrix for a [`MultiSketch`] keeping the
    /// minimum `knn` distances with distance options set by `jaccard`
    pub fn new(ref_sketches: &'a MultiSketch, knn: usize, jaccard: DistType) -> Self {
        let n_distances = ref_sketches.number_samples_loaded() * knn;

        // Pre-allocate distances
        let distances = match jaccard {
            DistType::CoreAcc => DistVec::CoreAcc(vec![SparseCoreAcc(0, 0.0, 0.0); n_distances]),
            DistType::Jaccard(_, _, _) => {
                DistVec::Jaccard(vec![SparseJaccard(0, 0.0); n_distances])
            }
        };

        Self {
            n_distances,
            knn,
            jaccard,
            distances,
            ref_names: Self::sketch_names(ref_sketches),
        }
    }

    /// Mutable reference to the underlying distance storage
    pub fn dists_mut(&mut self) -> &mut DistVec {
        &mut self.distances
    }
}

impl<'a> Distances<'a> for SparseDistanceMatrix<'a> {
    fn jaccard(&self) -> &DistType {
        &self.jaccard
    }
}

impl fmt::Display for SparseDistanceMatrix<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut ref_name_iter = self.ref_names.iter();
        let mut ref_name = ref_name_iter.next().unwrap();
        let mut k = 0;
        match &self.distances {
            DistVec::Jaccard(dists) => {
                for dist_item in dists {
                    k += 1;
                    if k > self.knn {
                        ref_name = ref_name_iter.next().unwrap();
                        k = 1;
                    }
                    let query_name = self.ref_names[dist_item.0];
                    // If fewer items than knn, padded with ref=query and dist=1.0 values
                    // which are skipped here
                    // TODO: more rust-like way of doing this would be to have
                    // SparseJaccard as an enum with an empty value
                    if dist_item.1 < 1.0_f32 || query_name != *ref_name {
                        writeln!(f, "{ref_name}\t{query_name}\t{}", dist_item.1)?;
                    }
                }
            }
            DistVec::CoreAcc(dists) => {
                for dist_item in dists {
                    k += 1;
                    if k > self.knn {
                        ref_name = ref_name_iter.next().unwrap();
                        k = 1;
                    }
                    writeln!(
                        f,
                        "{ref_name}\t{}\t{}\t{}",
                        self.ref_names[dist_item.0], dist_item.1, dist_item.2,
                    )?;
                }
            }
        }
        Ok(())
    }
}
