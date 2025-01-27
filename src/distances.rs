//! Functions and traits for calculating and storing distances
use std::cmp::Ordering;
use std::fmt;

// use ordered_float::NotNan;

use crate::multisketch::MultiSketch;

#[inline(always)]
pub fn square_to_condensed(i: usize, j: usize, n: usize) -> usize {
    debug_assert!(j > i);
    n * i - ((i * (i + 1)) >> 1) + j - 1 - i
}

#[inline(always)]
pub fn ref_query_index(i: usize, j: usize, n: usize) -> usize {
    debug_assert!(j > i);
    i * n + j
}

#[inline(always)]
pub fn calc_query_indices(k: usize, n: usize) -> (usize, usize) {
    let i = k / n;
    let j = k % n;
    debug_assert!(i < n);
    debug_assert!(j < n);
    (i, j)
}

#[inline(always)]
pub fn calc_col_idx(k: usize, i: usize, n: usize) -> usize {
    debug_assert!(i < n);
    let k_i64 = k as i64;
    let i_i64 = i as i64;
    let n_i64 = n as i64;
    (k_i64 + i_i64 + 1 - n_i64 * (n_i64 - 1) / 2 + (n_i64 - i_i64) * ((n_i64 - i_i64) - 1) / 2)
        as usize
}

#[inline(always)]
pub fn calc_row_idx(k: usize, n: usize) -> usize {
    let k_i64 = k as i64;
    let n_i64 = n as i64;
    n - 2
        - (((-8 * k_i64 + 4 * n_i64 * (n_i64 - 1) - 7) as f64).sqrt() / 2.0 - 0.5).floor() as usize
}

#[derive(PartialEq, PartialOrd)]
pub enum DistType {
    Jaccard(usize, bool),
    CoreAcc,
}

impl fmt::Display for DistType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            DistType::CoreAcc => write!(f, "Distances: core/accessory regression"),
            DistType::Jaccard(k, ani) => {
                if ani {
                    write!(f, "Distances: ANI at k={k}")
                } else {
                    write!(f, "Distances: Jaccard distances at k={k}")
                }
            }
        }
    }
}

pub trait Distances<'a> {
    fn jaccard(&self) -> &DistType;

    fn n_dist_cols(&self) -> usize {
        match self.jaccard() {
            DistType::CoreAcc => 2,
            DistType::Jaccard(_, _) => 1,
        }
    }

    fn sketch_names(sketches: &'a MultiSketch) -> Vec<&'a str> {
        let n_samples = sketches.number_samples_loaded();
        let mut names = Vec::with_capacity(n_samples);
        for idx in 0..n_samples {
            names.push(sketches.sketch_name(idx));
        }
        names
    }
}

pub struct DistanceMatrix<'a> {
    pub n_distances: usize,
    jaccard: DistType,
    distances: Vec<f32>,
    ref_names: Vec<&'a str>,
    query_names: Option<Vec<&'a str>>,
}

impl<'a> DistanceMatrix<'a> {
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

    pub fn dists_mut(&mut self) -> &mut Vec<f32> {
        &mut self.distances
    }
}

impl<'a> Distances<'a> for DistanceMatrix<'a> {
    fn jaccard(&self) -> &DistType {
        &self.jaccard
    }
}

impl<'a> fmt::Display for DistanceMatrix<'a> {
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

#[derive(Debug, Clone)]
pub struct SparseJaccard(pub usize, pub f32);
impl Ord for SparseJaccard {
    fn cmp(&self, other: &Self) -> Ordering {
        other.1.partial_cmp(&self.1).unwrap() // NB: backwards
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
        other.1.partial_cmp(&self.1) // NB: backwards
    }
}
impl PartialEq for SparseJaccard {
    fn eq(&self, other: &Self) -> bool {
        self.1 == other.1
    }
}
impl Eq for SparseJaccard {}

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
        self.1.partial_cmp(&other.1)
    }
}
impl PartialEq for SparseCoreAcc {
    fn eq(&self, other: &Self) -> bool {
        self.1 == other.1
    }
}
impl Eq for SparseCoreAcc {}

pub enum DistVec {
    Jaccard(Vec<SparseJaccard>),
    CoreAcc(Vec<SparseCoreAcc>),
}

pub struct SparseDistanceMatrix<'a> {
    pub n_distances: usize,
    pub knn: usize,
    jaccard: DistType,
    distances: DistVec,
    ref_names: Vec<&'a str>,
}

impl<'a> SparseDistanceMatrix<'a> {
    pub fn new(ref_sketches: &'a MultiSketch, knn: usize, jaccard: DistType) -> Self {
        let n_distances = ref_sketches.number_samples_loaded() * knn;

        // Pre-allocate distances
        let distances = match jaccard {
            DistType::CoreAcc => DistVec::CoreAcc(vec![SparseCoreAcc(0, 0.0, 0.0); n_distances]),
            DistType::Jaccard(_, _) => DistVec::Jaccard(vec![SparseJaccard(0, 0.0); n_distances]),
        };

        Self {
            n_distances,
            knn,
            jaccard,
            distances,
            ref_names: Self::sketch_names(ref_sketches),
        }
    }

    pub fn dists_mut(&mut self) -> &mut DistVec {
        &mut self.distances
    }
}

impl<'a> Distances<'a> for SparseDistanceMatrix<'a> {
    fn jaccard(&self) -> &DistType {
        &self.jaccard
    }
}

impl<'a> fmt::Display for SparseDistanceMatrix<'a> {
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
                    writeln!(
                        f,
                        "{ref_name}\t{}\t{}",
                        self.ref_names[dist_item.0], dist_item.1,
                    )?;
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
