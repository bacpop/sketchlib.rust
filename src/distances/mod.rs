//! Functions to calculate distances between sample sets
use std::collections::BinaryHeap;

use anyhow::{Context, Error};
use hashbrown::HashMap;
use indicatif::ParallelProgressIterator;
use rayon::prelude::*;

use crate::get_progress_bar;
use crate::inverted::Inverted;
use crate::sketch::multisketch::MultiSketch;

pub mod distance_matrix;
use self::distance_matrix::*;
mod jaccard;
use self::jaccard::*;

/// Chunk size in parallel distance calculations
const CHUNK_SIZE: usize = 1000;
// Distance progress bars use percent rather than number of comparisons
const BAR_PERCENT: bool = true;

/// Set type of distances to use and set up k-mer index
pub fn set_k(sketches: &MultiSketch, kmer: Option<usize>, ani: bool) -> Result<DistType, Error> {
    let k_idx;
    let dist_type = if let Some(k) = kmer {
        k_idx = sketches
            .get_k_idx(k)
            .with_context(|| format!("K-mer size {k} not found in file"))?;
        DistType::Jaccard(k_idx, k as f32, ani)
    } else {
        DistType::CoreAcc
    };
    log::info!("{dist_type}");
    Ok(dist_type)
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
    dist_type: DistType,
    quiet: bool,
) -> DistanceMatrix {
    let mut distances = DistanceMatrix::new(sketches, None, dist_type);
    let k_vals = distances.k_vals();
    let ani = distances.ani();
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
                if let Some((k_idx, k_f32)) = k_vals {
                    let j_index = jaccard_index(
                        sketches.get_sketch_slice(i, k_idx),
                        sketches.get_sketch_slice(j, k_idx),
                        sketches.sketchsize64,
                    );
                    let dist = if ani {
                        ani_pois(j_index, k_f32)
                    } else {
                        1.0_f32 - j_index
                    };
                    dist_slice[dist_idx] = dist;
                } else {
                    let dist = core_acc_dist(sketches, sketches, i, j);
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
    dist_type: DistType,
    quiet: bool,
) -> SparseDistanceMatrix {
    let mut sp_distances = SparseDistanceMatrix::new(sketches, knn, dist_type);
    let k_vals = sp_distances.k_vals();
    let ani = sp_distances.ani();
    let progress_bar = get_progress_bar(n, BAR_PERCENT, quiet);
    match sp_distances.dists_mut() {
        DistVec::Jaccard(distances) => {
            let (k_idx, k_f32) = k_vals.unwrap();
            distances
                .par_chunks_mut(knn)
                .progress_with(progress_bar)
                .enumerate()
                .for_each(|(i, row_dist_slice)| {
                    let mut heap = BinaryHeap::with_capacity(knn + 1);
                    let i_sketch = sketches.get_sketch_slice(i, k_idx);
                    for j in 0..n {
                        if i == j {
                            continue;
                        }
                        let mut dist = jaccard_index(
                            i_sketch,
                            sketches.get_sketch_slice(j, k_idx),
                            sketches.sketchsize64,
                        );
                        dist = if ani {
                            // This is just done so the heap sorts correctly (as want to keep higher ANI)
                            1.0_f32 - ani_pois(dist, k_f32)
                        } else {
                            1.0_f32 - dist
                        };
                        let dist_item = SparseJaccard(j, dist);
                        push_heap(&mut heap, dist_item, knn);
                    }
                    debug_assert_eq!(row_dist_slice.len(), heap.len());
                    if ani {
                        // Undo the above transform
                        heap.into_sorted_vec().iter().zip(row_dist_slice).for_each(|(inverse_ani, output_ani)| {
                            *output_ani = SparseJaccard(inverse_ani.0, 1.0_f32 - inverse_ani.1);
                        });
                    } else {
                        row_dist_slice.clone_from_slice(&heap.into_sorted_vec());
                    }
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
                        let dists = core_acc_dist(sketches, sketches, i, j);
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
    dist_type: DistType,
    quiet: bool,
) -> DistanceMatrix<'a> {
    let mut distances = DistanceMatrix::new(ref_sketches, Some(query_sketches), dist_type);
    let k_vals = distances.k_vals();
    let ani = distances.ani();
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
                if let Some((k_idx, k_f32)) = k_vals {
                    let j_index = jaccard_index(
                        ref_sketches.get_sketch_slice(i, k_idx),
                        query_sketches.get_sketch_slice(j, k_idx),
                        ref_sketches.sketchsize64,
                    );
                    let dist = if ani {
                        ani_pois(j_index, k_f32)
                    } else {
                        1.0_f32 - j_index
                    };
                    dist_slice[dist_idx] = dist;
                } else {
                    let dist = core_acc_dist(ref_sketches, query_sketches, i, j);
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

/// Same as [`self_dists_knn`], but also using an inverted_index to precluster
/// to reduce the number of comparisons
pub fn self_dists_knn_precluster<'a>(
    sketches: &'a MultiSketch,
    inverted_index: &Inverted,
    skq_bins: &[u16],
    skq_stride: usize,
    n: usize,
    knn: usize,
    dist_type: DistType,
    quiet: bool,
) -> SparseDistanceMatrix<'a> {
    // Check that sample sets in ski and skm are the same and create i,j lookup
    let mut skq_lookup = HashMap::with_capacity(n);
    for (skq_index, skq_sample) in inverted_index.sample_names().iter().enumerate() {
        skq_lookup.insert(skq_sample.as_str(), skq_index);
    }
    let mut not_found = Vec::new();
    let mut skq_index_lookup = Vec::with_capacity(n);
    for skd_sample_idx in 0..sketches.number_samples_loaded() {
        let sample_name = sketches.sketch_name(skd_sample_idx);
        match skq_lookup.get(sample_name) {
            Some(skq_index) => skq_index_lookup.push(skq_index),
            None => not_found.push(sample_name),
        };
    }
    if !not_found.is_empty() {
        panic!("The following samples in the .skd could not be found in the .ski:\n{not_found:?}");
    }

    let mut sp_distances = SparseDistanceMatrix::new(sketches, knn, dist_type);
    let k_vals = sp_distances.k_vals();
    let ani = sp_distances.ani();
    let progress_bar = get_progress_bar(n, BAR_PERCENT, quiet);
    match sp_distances.dists_mut() {
        DistVec::Jaccard(distances) => {
            let (k_idx, k_f32) = k_vals.unwrap();
            distances
                .par_chunks_mut(knn)
                .progress_with(progress_bar)
                .enumerate()
                .for_each(|(i, row_dist_slice)| {
                    // Prefilter step here
                    let skq_offset = i * skq_stride;
                    let flat_i_sketch = &skq_bins[skq_offset..(skq_offset + skq_stride)];
                    let prefiltered_samples = inverted_index.any_shared_bins(flat_i_sketch);
                    // Standard search
                    let mut heap = BinaryHeap::with_capacity(knn + 1);
                    let i_sketch = sketches.get_sketch_slice(i, k_idx);
                    for j in prefiltered_samples {
                        let j = j as usize;
                        if i == j {
                            continue;
                        }
                        let mut dist = jaccard_index(
                            i_sketch,
                            sketches.get_sketch_slice(j, k_idx),
                            sketches.sketchsize64,
                        );
                        dist = if ani {
                            1.0_f32 - ani_pois(dist, k_f32)
                        } else {
                            1.0_f32 - dist
                        };
                        let dist_item = SparseJaccard(j, dist);
                        push_heap(&mut heap, dist_item, knn);
                    }
                    let mut dist_vec = heap.into_sorted_vec();
                    if ani {
                        // Undo the above transform
                        dist_vec.iter_mut().for_each(|inverse_ani| {
                            *inverse_ani = SparseJaccard(inverse_ani.0, 1.0_f32 - inverse_ani.1);
                        });
                    }
                    // If there are fewer prefiltered dists than knn, add null values at the end
                    if dist_vec.len() < row_dist_slice.len() {
                        // TODO: more rust-like way of doing this would be to have
                        // SparseJaccard as an enum with an empty value
                        dist_vec.append(&mut vec![
                            SparseJaccard(i, 1.0);
                            row_dist_slice.len() - dist_vec.len()
                        ]);
                    }
                    row_dist_slice.clone_from_slice(&dist_vec);
                });
        }
        DistVec::CoreAcc(_) => {
            unimplemented!("Prefilter only available for single k-mer distances");
        }
    }
    sp_distances
}
