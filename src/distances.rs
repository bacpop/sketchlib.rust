use std::fmt;

use crate::multisketch::MultiSketch;


pub struct DistanceMatrix<'a> {
    pub n_distances: usize,
    distances: Vec<f64>,
    other_distances: Vec<f64>,
    ref_names: Vec<&'a str>,
    query_names: Option<Vec<&'a str>>,
}

impl<'a> DistanceMatrix<'a> {
    pub fn new(ref_sketches: &'a MultiSketch, query_sketches: Option<&'a MultiSketch>, jaccard: bool) -> Self {
        if let Some(query) = query_sketches {
            let n_distances = ref_sketches.number_samples_loaded() * query.number_samples_loaded();
            let mut other_distances = Vec::new();
            if jaccard {
                other_distances.reserve(n_distances);
            }

            Self { n_distances, distances: Vec::with_capacity(n_distances), other_distances, ref_names: Self::sketch_names(ref_sketches), query_names: Some(Self::sketch_names(query))}
        } else {
            let n_distances = ref_sketches.number_samples_loaded() * (ref_sketches.number_samples_loaded() - 1) / 2;
            let mut other_distances = Vec::new();
            if jaccard {
                other_distances.reserve(n_distances);
            }

            Self { n_distances, distances: Vec::with_capacity(n_distances), other_distances, ref_names: Self::sketch_names(ref_sketches), query_names: None}
        }
    }

    pub fn add_jaccard_dist(&mut self, dist: f64) {
        self.distances.push(dist);
    }

    pub fn add_core_acc_dist(&mut self, core_dist: f64, acc_dist: f64) {
        self.distances.push(core_dist);
        self.other_distances.push(acc_dist);
    }

    fn sketch_names(sketches: &MultiSketch) -> Vec<&str> {
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
                    if self.other_distances.len() > 0 {
                        write!(f, "\t{}", self.other_distances[dist_idx])?;
                    }
                    write!(f, "\n")?;
                    dist_idx += 1;
                }
            }
        } else {
            for (i, ref_name) in self.ref_names.iter().enumerate() {
                for j in i..self.ref_names.len() {
                    write!(f, "{ref_name}\t{}\t{}", self.ref_names[j], self.distances[dist_idx])?;
                    if self.other_distances.len() > 0 {
                        write!(f, "\t{}", self.other_distances[dist_idx])?;
                    }
                    write!(f, "\n")?;
                    dist_idx += 1;
                }
            }
        }
        Ok(())
    }
}