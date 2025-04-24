//! Functions to calculate distances between sample sets
use std::collections::BinaryHeap;

use indicatif::ParallelProgressIterator;
use rayon::prelude::*;

use crate::get_progress_bar;
use crate::multisketch::MultiSketch;

pub mod distance_matrix;
use self::distance_matrix::*;
mod jaccard;
use self::jaccard::*;

/// Chunk size in parallel distance calculations
const CHUNK_SIZE: usize = 1000;
// Distance progress bars use percent rather than number of comparisons
const BAR_PERCENT: bool = true;

type DistK = (DistType, Option<usize>, f32);

/// Set type of distances to use and set up k-mer index
pub fn set_k(sketches: &MultiSketch, kmer: Option<usize>, ani: bool) -> DistK {
    let k_idx;
    let mut k_f32 = 0.0;
    let dist_type = if let Some(k) = kmer {
        k_idx = sketches.get_k_idx(k);
        k_f32 = k as f32;
        DistType::Jaccard(k, ani)
    } else {
        k_idx = None;
        DistType::CoreAcc
    };
    log::info!("{dist_type}");
    (dist_type, k_idx, k_f32)
}

// Add to distances, only keep the best knn or fewer
#[inline(always)]
fn push_heap<T: PartialOrd + Ord>(heap: &mut BinaryHeap<T>, dist_item: T, knn: usize) {
    if heap.len() < knn || dist_item < *heap.peek().unwrap() {
        heap.push(dist_item);
        if heap.len() > knn {
            heap.pop();
        }
    }
}

// Notes and ideas
//      Possible improvement would be to load sketch slices when i, j change
//      This would require a change to core_acc where multiple k-mer lengths are loaded at once
//      Overall this would be nicer I think (not sure about speed)
//
//      Streaming out of distances, when a sample is 'ready'

/// Self query mode (dense, all distances)
pub fn self_dists_all(
    sketches: &MultiSketch,
    n: usize,
    k_idx: Option<usize>,
    k_f32: f32,
    dist_type: DistType,
    ani: bool,
    quiet: bool,
) -> DistanceMatrix {
    let mut distances = DistanceMatrix::new(&sketches, None, dist_type);
    let par_chunk = CHUNK_SIZE * distances.n_dist_cols();
    let progress_bar = get_progress_bar(par_chunk, BAR_PERCENT, quiet);
    distances
        .dists_mut()
        .par_chunks_mut(par_chunk)
        .progress_with(progress_bar)
        .enumerate()
        .for_each(|(chunk_idx, dist_slice)| {
            // Get first i, j index for the chunk
            let start_dist_idx = chunk_idx * CHUNK_SIZE;
            let mut i = calc_row_idx(start_dist_idx, n);
            let mut j = calc_col_idx(start_dist_idx, i, n);
            for dist_idx in 0..CHUNK_SIZE {
                if let Some(k) = k_idx {
                    let mut dist = jaccard_dist(
                        sketches.get_sketch_slice(i, k),
                        sketches.get_sketch_slice(j, k),
                        sketches.sketchsize64,
                    );
                    dist = if ani {
                        ani_pois(dist, k_f32)
                    } else {
                        1.0_f32 - dist
                    };
                    dist_slice[dist_idx] = dist;
                } else {
                    let dist = core_acc_dist(&sketches, &sketches, i, j);
                    dist_slice[dist_idx * 2] = dist.0;
                    dist_slice[dist_idx * 2 + 1] = dist.1;
                }

                // Move to next index in upper triangle
                j += 1;
                if j >= n {
                    i += 1;
                    j = i + 1;
                    // End of all dists reached (final chunk)
                    if i >= (n - 1) {
                        break;
                    }
                }
            }
        });
    distances
}

/// Self query mode (dense, all distances)
pub fn self_dists_knn(
    sketches: &MultiSketch,
    n: usize,
    knn: usize,
    k_idx: Option<usize>,
    k_f32: f32,
    dist_type: DistType,
    ani: bool,
    quiet: bool,
) -> SparseDistanceMatrix {
    let mut sp_distances = SparseDistanceMatrix::new(&sketches, knn, dist_type);
    let progress_bar = get_progress_bar(n * knn, BAR_PERCENT, quiet);
    match sp_distances.dists_mut() {
        DistVec::Jaccard(distances) => {
            let k = k_idx.unwrap();
            distances
                .par_chunks_mut(knn)
                .progress_with(progress_bar)
                .enumerate()
                .for_each(|(i, row_dist_slice)| {
                    let mut heap = BinaryHeap::with_capacity(knn + 1);
                    let i_sketch = sketches.get_sketch_slice(i, k);
                    for j in 0..n {
                        if i == j {
                            continue;
                        }
                        let mut dist = jaccard_dist(
                            i_sketch,
                            sketches.get_sketch_slice(j, k),
                            sketches.sketchsize64,
                        );
                        dist = if ani {
                            ani_pois(dist, k_f32)
                        } else {
                            1.0_f32 - dist
                        };
                        let dist_item = SparseJaccard(j, dist);
                        push_heap(&mut heap, dist_item, knn);
                    }
                    debug_assert_eq!(row_dist_slice.len(), heap.len());
                    row_dist_slice.clone_from_slice(&heap.into_sorted_vec());
                });
        }
        DistVec::CoreAcc(distances) => {
            distances
                .par_chunks_mut(knn)
                .progress_with(progress_bar)
                .enumerate()
                .for_each(|(i, row_dist_slice)| {
                    let mut heap = BinaryHeap::with_capacity(knn + 1);
                    for j in 0..n {
                        if i == j {
                            continue;
                        }
                        let dists = core_acc_dist(&sketches, &sketches, i, j);
                        let dist_item = SparseCoreAcc(j, dists.0, dists.1);
                        push_heap(&mut heap, dist_item, knn);
                    }
                    debug_assert_eq!(row_dist_slice.len(), heap.len());
                    row_dist_slice.clone_from_slice(&heap.into_sorted_vec());
                });
        }
    }
    sp_distances
}

/// Self query mode (dense, all distances)
pub fn self_query_dists_all<'a>(
    ref_sketches: &'a MultiSketch,
    query_sketches: &'a MultiSketch,
    n: usize,
    nq: usize,
    k_idx: Option<usize>,
    k_f32: f32,
    dist_type: DistType,
    ani: bool,
    quiet: bool,
) -> DistanceMatrix<'a> {
    let mut distances =
        DistanceMatrix::new(&ref_sketches, Some(&query_sketches), dist_type);
    let par_chunk = CHUNK_SIZE * distances.n_dist_cols();
    let progress_bar = get_progress_bar(par_chunk, BAR_PERCENT, quiet);
    distances
        .dists_mut()
        .par_chunks_mut(par_chunk)
        .progress_with(progress_bar)
        .enumerate()
        .for_each(|(chunk_idx, dist_slice)| {
            // Get first i, j index for the chunk
            let start_dist_idx = chunk_idx * CHUNK_SIZE;
            let (mut i, mut j) = calc_query_indices(start_dist_idx, nq);
            for dist_idx in 0..CHUNK_SIZE {
                if let Some(k) = k_idx {
                    let mut dist = jaccard_dist(
                        ref_sketches.get_sketch_slice(i, k),
                        query_sketches.get_sketch_slice(j, k),
                        ref_sketches.sketchsize64,
                    );
                    dist = if ani {
                        ani_pois(dist, k_f32)
                    } else {
                        1.0_f32 - dist
                    };
                    dist_slice[dist_idx] = dist;
                } else {
                    let dist = core_acc_dist(&ref_sketches, &query_sketches, i, j);
                    dist_slice[dist_idx * 2] = dist.0;
                    dist_slice[dist_idx * 2 + 1] = dist.1;
                }

                // Move to next index
                j += 1;
                if j >= nq {
                    i += 1;
                    j = 0;
                    // End of all dists reached (final chunk)
                    if i >= n {
                        break;
                    }
                }
            }
        });
        distances
    }
