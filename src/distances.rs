use std::fmt;

use crate::multisketch::MultiSketch;

#[inline(always)]
pub fn square_to_condensed(i: usize, j: usize, n: usize) -> usize {
    debug_assert!(j > i);
    return n * i - ((i * (i + 1)) >> 1) + j - 1 - i;
}

#[inline(always)]
pub fn ref_query_index(i: usize, j: usize, n: usize) -> usize {
    debug_assert!(j > i);
    i * n + j
}

#[inline(always)]
pub fn calc_col_idx(k: usize, i: usize, n: usize) -> usize {
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
    Jaccard(usize),
    CoreAcc,
}

impl fmt::Display for DistType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            &DistType::CoreAcc => write!(f, "Distances: core/accessory regression"),
            &DistType::Jaccard(k) => write!(f, "Distances: Jaccard distances at k={k}"),
        }
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

    pub fn n_dist_cols(&self) -> usize {
        match self.jaccard {
            DistType::CoreAcc => 2,
            DistType::Jaccard(_) => 1,
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
                    write!(f, "\n")?;
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
                    write!(f, "\n")?;
                    dist_idx += 1;
                }
            }
        }
        Ok(())
    }
}
