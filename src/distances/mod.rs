//! Functions to calculate distances between sample sets
use std::collections::BinaryHeap;
use std::fmt::Write as FmtWrite;
use std::io::Write as IoWrite;
use std::sync::mpsc;

use anyhow::{bail, Context, Error};
use hashbrown::HashMap;
use indicatif::ParallelProgressIterator;
use rayon::prelude::*;

use crate::cli::RetainUnmatched;
use crate::get_progress_bar;
use crate::inverted::Inverted;
use crate::sketch::multisketch::MultiSketch;
use crate::sketch::{BIN_BITS, LEGACY_BIN_BITS};

/// Error message used whenever a reference and query database with
/// mismatched sketch generations (legacy pre-v0.4 vs. current) are compared.
/// Legacy and new-format databases use incompatible bin-packing schemes, so
/// this is an explicit rejection rather than a silent/lossy conversion.
fn mismatched_generation_message(ref_legacy: bool, query_legacy: bool) -> String {
    format!(
        "Cannot compare reference and query databases with different sketch generations (reference is_legacy={ref_legacy}, query is_legacy={query_legacy}): legacy (pre-v0.4) and new-format databases use incompatible bin-packing schemes and cannot be directly compared. Please re-sketch both databases with the current version."
    )
}

pub mod distance_matrix;
use self::distance_matrix::*;
pub mod jaccard;
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
        DistType::Jaccard(k_idx, k as f64, ani)
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
///
/// Computes all pairwise distances within a single set of `sketches`, iterating
/// the upper triangle: `i` and `j` both index into `sketches`, with `i` the outer
/// (row) index and `j` the inner (column) index, `j` always `> i`.
pub fn self_dists_all<'a>(
    sketches: &'a MultiSketch,
    n: usize,
    dist_type: DistType,
    quiet: bool,
    completeness_vec: Option<&Vec<f64>>,
    completeness_cutoff: f64,
) -> DistanceMatrix<'a> {
    if sketches.is_legacy_format() {
        self_dists_all_generic::<LEGACY_BIN_BITS>(
            sketches,
            n,
            dist_type,
            quiet,
            completeness_vec,
            completeness_cutoff,
        )
    } else {
        self_dists_all_generic::<BIN_BITS>(
            sketches,
            n,
            dist_type,
            quiet,
            completeness_vec,
            completeness_cutoff,
        )
    }
}

fn self_dists_all_generic<'a, const BB: usize>(
    sketches: &'a MultiSketch,
    n: usize,
    dist_type: DistType,
    quiet: bool,
    completeness_vec: Option<&Vec<f64>>,
    completeness_cutoff: f64,
) -> DistanceMatrix<'a> {
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
                if let Some((k_idx, k_f64)) = k_vals {
                    // If completeness_vec is Some, extract the value at index i (or j) from the inner vector.
                    // If completeness_vec is None, the result will also be None.
                    // This uses Option::map to safely access the completeness value for each sample.
                    let c1 = completeness_vec.map(|cv| cv[i]);
                    let c2 = completeness_vec.map(|cv| cv[j]);
                    let j_index = jaccard_index_generic::<BB>(
                        sketches.get_sketch_slice(i, k_idx),
                        sketches.get_sketch_slice(j, k_idx),
                        sketches.sketchsize64,
                        c1,
                        c2,
                        completeness_cutoff,
                    );
                    let dist = if ani {
                        ani_pois(j_index, k_f64) as f32
                    } else {
                        (1.0_f64 - j_index) as f32
                    };
                    dist_slice[dist_idx] = dist;
                } else {
                    let dist = core_acc_dist_generic::<BB>(
                        sketches,
                        sketches,
                        i,
                        j,
                        completeness_vec,
                        completeness_vec,
                        completeness_cutoff,
                    );
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

/// Self query mode (dense, all distances), streaming text output directly to `writer`
/// instead of materializing a [`DistanceMatrix`]. This keeps peak memory bounded
/// (independent of `n`) for large all-vs-all comparisons.
///
/// Line format matches [`DistanceMatrix`]'s `Display` impl exactly, but pair order is
/// not guaranteed to match it (chunks may be written as they complete, from whichever
/// thread finishes first) — every pair is still written exactly once.
pub fn self_dists_all_stream<W: IoWrite + Send>(
    writer: &mut W,
    sketches: &MultiSketch,
    n: usize,
    dist_type: DistType,
    quiet: bool,
    completeness_vec: Option<&Vec<f64>>,
    completeness_cutoff: f64,
    threads: usize,
) -> Result<(), Error> {
    if sketches.is_legacy_format() {
        self_dists_all_stream_generic::<LEGACY_BIN_BITS, W>(
            writer,
            sketches,
            n,
            dist_type,
            quiet,
            completeness_vec,
            completeness_cutoff,
            threads,
        )
    } else {
        self_dists_all_stream_generic::<BIN_BITS, W>(
            writer,
            sketches,
            n,
            dist_type,
            quiet,
            completeness_vec,
            completeness_cutoff,
            threads,
        )
    }
}

fn self_dists_all_stream_generic<const BB: usize, W: IoWrite + Send>(
    writer: &mut W,
    sketches: &MultiSketch,
    n: usize,
    dist_type: DistType,
    quiet: bool,
    completeness_vec: Option<&Vec<f64>>,
    completeness_cutoff: f64,
    threads: usize,
) -> Result<(), Error> {
    let ani = matches!(dist_type, DistType::Jaccard(_, _, true));
    let k_vals = match dist_type {
        DistType::Jaccard(k_idx, k_val, _) => Some((k_idx, k_val)),
        DistType::CoreAcc => None,
    };
    let ref_names = <DistanceMatrix as Distances>::sketch_names(sketches);

    let n_distances = n * (n - 1) / 2;
    let n_chunks = n_distances.div_ceil(CHUNK_SIZE);
    let progress_bar = get_progress_bar(n_chunks, BAR_PERCENT, quiet);

    // Bounded channel: caps in-flight chunk memory at a small constant regardless of
    // n, rather than relying on producer/writer speeds happening to balance out.
    let channel_bound = threads.max(1) * 2;
    let (tx, rx) = mpsc::sync_channel::<String>(channel_bound);

    rayon::scope(|s| -> Result<(), Error> {
        s.spawn(move |_| {
            (0..n_chunks)
                .into_par_iter()
                .progress_with(progress_bar)
                .for_each_with(tx, |tx, chunk_idx| {
                    let start = chunk_idx * CHUNK_SIZE;
                    let end = (start + CHUNK_SIZE).min(n_distances);
                    let mut i = calc_row_idx(start, n);
                    let mut j = calc_col_idx(start, i, n);
                    let mut buf = String::with_capacity((end - start) * 24);

                    for _ in start..end {
                        if let Some((k_idx, k_f64)) = k_vals {
                            let c1 = completeness_vec.map(|cv| cv[i]);
                            let c2 = completeness_vec.map(|cv| cv[j]);
                            let j_index = jaccard_index_generic::<BB>(
                                sketches.get_sketch_slice(i, k_idx),
                                sketches.get_sketch_slice(j, k_idx),
                                sketches.sketchsize64,
                                c1,
                                c2,
                                completeness_cutoff,
                            );
                            let dist = if ani {
                                ani_pois(j_index, k_f64) as f32
                            } else {
                                (1.0_f64 - j_index) as f32
                            };
                            let _ = writeln!(buf, "{}\t{}\t{dist}", ref_names[i], ref_names[j]);
                        } else {
                            let d = core_acc_dist_generic::<BB>(
                                sketches,
                                sketches,
                                i,
                                j,
                                completeness_vec,
                                completeness_vec,
                                completeness_cutoff,
                            );
                            let _ = writeln!(
                                buf,
                                "{}\t{}\t{}\t{}",
                                ref_names[i], ref_names[j], d.0, d.1
                            );
                        }

                        // Move to next index in upper triangle
                        j += 1;
                        if j >= n {
                            i += 1;
                            j = i + 1;
                        }
                    }
                    let _ = tx.send(buf);
                });
        });

        for chunk_text in rx {
            writer
                .write_all(chunk_text.as_bytes())
                .context("Error writing streamed distance output")?;
        }
        Ok(())
    })
}

/// Self query mode (dense, all distances)
pub fn self_dists_knn<'a>(
    sketches: &'a MultiSketch,
    n: usize,
    knn: usize,
    dist_type: DistType,
    quiet: bool,
    completeness_vec: Option<&Vec<f64>>,
    completeness_cutoff: f64,
) -> SparseDistanceMatrix<'a> {
    if sketches.is_legacy_format() {
        self_dists_knn_generic::<LEGACY_BIN_BITS>(
            sketches,
            n,
            knn,
            dist_type,
            quiet,
            completeness_vec,
            completeness_cutoff,
        )
    } else {
        self_dists_knn_generic::<BIN_BITS>(
            sketches,
            n,
            knn,
            dist_type,
            quiet,
            completeness_vec,
            completeness_cutoff,
        )
    }
}

fn self_dists_knn_generic<'a, const BB: usize>(
    sketches: &'a MultiSketch,
    n: usize,
    knn: usize,
    dist_type: DistType,
    quiet: bool,
    completeness_vec: Option<&Vec<f64>>,
    completeness_cutoff: f64,
) -> SparseDistanceMatrix<'a> {
    let mut sp_distances = SparseDistanceMatrix::new(sketches, knn, dist_type);
    let k_vals = sp_distances.k_vals();
    let ani = sp_distances.ani();
    let progress_bar = get_progress_bar(n, BAR_PERCENT, quiet);
    match sp_distances.dists_mut() {
        DistVec::Jaccard(distances) => {
            let (k_idx, k_f64) = k_vals.unwrap();
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
                        // If completeness_vec is Some, extract the value at index i (or j) from the inner vector.
                        // If completeness_vec is None, the result will also be None.
                        // This uses Option::map to safely access the completeness value for each sample.
                        let c1 = completeness_vec.map(|cv| cv[i]);
                        let c2 = completeness_vec.map(|cv| cv[j]);
                        let dist = jaccard_index_generic::<BB>(
                            i_sketch,
                            sketches.get_sketch_slice(j, k_idx),
                            sketches.sketchsize64,
                            c1,
                            c2,
                            completeness_cutoff,
                        );
                        let dist_f32 = if ani {
                            // This is just done so the heap sorts correctly (as want to keep higher ANI)
                            (1.0_f64 - ani_pois(dist, k_f64)) as f32
                        } else {
                            (1.0_f64 - dist) as f32
                        };
                        let dist_item = SparseJaccard(j, dist_f32);
                        push_heap(&mut heap, dist_item, knn);
                    }
                    debug_assert_eq!(row_dist_slice.len(), heap.len());
                    if ani {
                        // Undo the above transform
                        heap.into_sorted_vec().iter().zip(row_dist_slice).for_each(
                            |(inverse_ani, output_ani)| {
                                *output_ani = SparseJaccard(inverse_ani.0, 1.0_f32 - inverse_ani.1);
                            },
                        );
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
                        let dists = core_acc_dist_generic::<BB>(
                            sketches,
                            sketches,
                            i,
                            j,
                            completeness_vec,
                            completeness_vec,
                            completeness_cutoff,
                        );
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

/// Cross-query mode (dense, all distances)
///
/// Computes all pairwise distances between `ref_sketches` and `query_sketches`:
/// `i` indexes `ref_sketches` (outer/row), `j` indexes `query_sketches`
/// (inner/column), chunked by `n_query` — i.e. ref-outer, query-inner ordering,
/// matching the row-major layout consumed by [`DistanceMatrix`].
pub fn cross_dists_all<'a>(
    ref_sketches: &'a MultiSketch,
    query_sketches: &'a MultiSketch,
    n: usize,
    n_query: usize,
    dist_type: DistType,
    quiet: bool,
    ref_completeness_vec: Option<&Vec<f64>>,
    query_completeness_vec: Option<&Vec<f64>>,
    completeness_cutoff: f64,
) -> DistanceMatrix<'a> {
    let (ref_legacy, query_legacy) = (
        ref_sketches.is_legacy_format(),
        query_sketches.is_legacy_format(),
    );
    if ref_legacy != query_legacy {
        panic!(
            "{}",
            mismatched_generation_message(ref_legacy, query_legacy)
        );
    }
    if ref_legacy {
        cross_dists_all_generic::<LEGACY_BIN_BITS>(
            ref_sketches,
            query_sketches,
            n,
            n_query,
            dist_type,
            quiet,
            ref_completeness_vec,
            query_completeness_vec,
            completeness_cutoff,
        )
    } else {
        cross_dists_all_generic::<BIN_BITS>(
            ref_sketches,
            query_sketches,
            n,
            n_query,
            dist_type,
            quiet,
            ref_completeness_vec,
            query_completeness_vec,
            completeness_cutoff,
        )
    }
}

fn cross_dists_all_generic<'a, const BB: usize>(
    ref_sketches: &'a MultiSketch,
    query_sketches: &'a MultiSketch,
    n: usize,
    n_query: usize,
    dist_type: DistType,
    quiet: bool,
    ref_completeness_vec: Option<&Vec<f64>>,
    query_completeness_vec: Option<&Vec<f64>>,
    completeness_cutoff: f64,
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
            let (mut i, mut j) = calc_query_indices(start_dist_idx, n_query);
            for dist_idx in 0..CHUNK_SIZE {
                if let Some((k_idx, k_f64)) = k_vals {
                    let c1 = ref_completeness_vec.map(|cv| cv[i]);
                    let c2 = query_completeness_vec.map(|cv| cv[j]);
                    let j_index = jaccard_index_generic::<BB>(
                        ref_sketches.get_sketch_slice(i, k_idx),
                        query_sketches.get_sketch_slice(j, k_idx),
                        ref_sketches.sketchsize64,
                        c1,
                        c2,
                        completeness_cutoff,
                    );
                    let dist = if ani {
                        ani_pois(j_index, k_f64) as f32
                    } else {
                        (1.0_f64 - j_index) as f32
                    };
                    dist_slice[dist_idx] = dist;
                } else {
                    let dist = core_acc_dist_generic::<BB>(
                        ref_sketches,
                        query_sketches,
                        i,
                        j,
                        ref_completeness_vec,
                        query_completeness_vec,
                        completeness_cutoff,
                    );
                    dist_slice[dist_idx * 2] = dist.0;
                    dist_slice[dist_idx * 2 + 1] = dist.1;
                }

                // Move to next index
                j += 1;
                if j >= n_query {
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

/// Cross-query mode (dense, all distances), streaming variant of [`cross_dists_all`].
/// Same ordering/format contract as [`self_dists_all_stream`].
pub fn cross_dists_all_stream<W: IoWrite + Send>(
    writer: &mut W,
    ref_sketches: &MultiSketch,
    query_sketches: &MultiSketch,
    n: usize,
    n_query: usize,
    dist_type: DistType,
    quiet: bool,
    ref_completeness_vec: Option<&Vec<f64>>,
    query_completeness_vec: Option<&Vec<f64>>,
    completeness_cutoff: f64,
    threads: usize,
) -> Result<(), Error> {
    let (ref_legacy, query_legacy) = (
        ref_sketches.is_legacy_format(),
        query_sketches.is_legacy_format(),
    );
    if ref_legacy != query_legacy {
        bail!(mismatched_generation_message(ref_legacy, query_legacy));
    }
    if ref_legacy {
        cross_dists_all_stream_generic::<LEGACY_BIN_BITS, W>(
            writer,
            ref_sketches,
            query_sketches,
            n,
            n_query,
            dist_type,
            quiet,
            ref_completeness_vec,
            query_completeness_vec,
            completeness_cutoff,
            threads,
        )
    } else {
        cross_dists_all_stream_generic::<BIN_BITS, W>(
            writer,
            ref_sketches,
            query_sketches,
            n,
            n_query,
            dist_type,
            quiet,
            ref_completeness_vec,
            query_completeness_vec,
            completeness_cutoff,
            threads,
        )
    }
}

fn cross_dists_all_stream_generic<const BB: usize, W: IoWrite + Send>(
    writer: &mut W,
    ref_sketches: &MultiSketch,
    query_sketches: &MultiSketch,
    n: usize,
    n_query: usize,
    dist_type: DistType,
    quiet: bool,
    ref_completeness_vec: Option<&Vec<f64>>,
    query_completeness_vec: Option<&Vec<f64>>,
    completeness_cutoff: f64,
    threads: usize,
) -> Result<(), Error> {
    let ani = matches!(dist_type, DistType::Jaccard(_, _, true));
    let k_vals = match dist_type {
        DistType::Jaccard(k_idx, k_val, _) => Some((k_idx, k_val)),
        DistType::CoreAcc => None,
    };
    let ref_names = <DistanceMatrix as Distances>::sketch_names(ref_sketches);
    let query_names = <DistanceMatrix as Distances>::sketch_names(query_sketches);

    let n_distances = n * n_query;
    let n_chunks = n_distances.div_ceil(CHUNK_SIZE);
    let progress_bar = get_progress_bar(n_chunks, BAR_PERCENT, quiet);

    let channel_bound = threads.max(1) * 2;
    let (tx, rx) = mpsc::sync_channel::<String>(channel_bound);

    rayon::scope(|s| -> Result<(), Error> {
        s.spawn(move |_| {
            (0..n_chunks)
                .into_par_iter()
                .progress_with(progress_bar)
                .for_each_with(tx, |tx, chunk_idx| {
                    let start = chunk_idx * CHUNK_SIZE;
                    let end = (start + CHUNK_SIZE).min(n_distances);
                    let (mut i, mut j) = calc_query_indices(start, n_query);
                    let mut buf = String::with_capacity((end - start) * 24);

                    for _ in start..end {
                        if let Some((k_idx, k_f64)) = k_vals {
                            let c1 = ref_completeness_vec.map(|cv| cv[i]);
                            let c2 = query_completeness_vec.map(|cv| cv[j]);
                            let j_index = jaccard_index_generic::<BB>(
                                ref_sketches.get_sketch_slice(i, k_idx),
                                query_sketches.get_sketch_slice(j, k_idx),
                                ref_sketches.sketchsize64,
                                c1,
                                c2,
                                completeness_cutoff,
                            );
                            let dist = if ani {
                                ani_pois(j_index, k_f64) as f32
                            } else {
                                (1.0_f64 - j_index) as f32
                            };
                            let _ = writeln!(buf, "{}\t{}\t{dist}", ref_names[i], query_names[j]);
                        } else {
                            let d = core_acc_dist_generic::<BB>(
                                ref_sketches,
                                query_sketches,
                                i,
                                j,
                                ref_completeness_vec,
                                query_completeness_vec,
                                completeness_cutoff,
                            );
                            let _ = writeln!(
                                buf,
                                "{}\t{}\t{}\t{}",
                                ref_names[i], query_names[j], d.0, d.1
                            );
                        }

                        // Move to next index
                        j += 1;
                        if j >= n_query {
                            i += 1;
                            j = 0;
                        }
                    }
                    let _ = tx.send(buf);
                });
        });

        for chunk_text in rx {
            writer
                .write_all(chunk_text.as_bytes())
                .context("Error writing streamed distance output")?;
        }
        Ok(())
    })
}

/// Cross-query mode with kNN filtering.
///
/// For each query genome, computes distances to all `n` reference genomes and
/// retains only the `knn` nearest neighbours using a priority queue (max-heap
/// capped at `knn`). Output has `n_query × knn` entries — one row per query genome.
///
/// This is the cross-database analogue of [`self_dists_knn`].
pub fn cross_dists_knn<'a>(
    ref_sketches: &'a MultiSketch,
    query_sketches: &'a MultiSketch,
    n: usize,
    n_query: usize,
    knn: usize,
    dist_type: DistType,
    quiet: bool,
    ref_completeness_vec: Option<&Vec<f64>>,
    query_completeness_vec: Option<&Vec<f64>>,
    completeness_cutoff: f64,
) -> SparseDistanceMatrix<'a> {
    let (ref_legacy, query_legacy) = (
        ref_sketches.is_legacy_format(),
        query_sketches.is_legacy_format(),
    );
    if ref_legacy != query_legacy {
        panic!(
            "{}",
            mismatched_generation_message(ref_legacy, query_legacy)
        );
    }
    if ref_legacy {
        cross_dists_knn_generic::<LEGACY_BIN_BITS>(
            ref_sketches,
            query_sketches,
            n,
            n_query,
            knn,
            dist_type,
            quiet,
            ref_completeness_vec,
            query_completeness_vec,
            completeness_cutoff,
        )
    } else {
        cross_dists_knn_generic::<BIN_BITS>(
            ref_sketches,
            query_sketches,
            n,
            n_query,
            knn,
            dist_type,
            quiet,
            ref_completeness_vec,
            query_completeness_vec,
            completeness_cutoff,
        )
    }
}

fn cross_dists_knn_generic<'a, const BB: usize>(
    ref_sketches: &'a MultiSketch,
    query_sketches: &'a MultiSketch,
    n: usize,
    n_query: usize,
    knn: usize,
    dist_type: DistType,
    quiet: bool,
    ref_completeness_vec: Option<&Vec<f64>>,
    query_completeness_vec: Option<&Vec<f64>>,
    completeness_cutoff: f64,
) -> SparseDistanceMatrix<'a> {
    if n == 0 {
        panic!("Reference database has no loaded samples");
    }
    if n_query == 0 {
        panic!("Query database has no loaded samples");
    }
    // Can't have more neighbours than there are reference genomes
    let knn = knn.min(n);
    let mut sp_distances =
        SparseDistanceMatrix::new_cross_query(ref_sketches, query_sketches, knn, dist_type);
    let k_vals = sp_distances.k_vals();
    let ani = sp_distances.ani();
    let progress_bar = get_progress_bar(n_query, BAR_PERCENT, quiet);
    match sp_distances.dists_mut() {
        DistVec::Jaccard(distances) => {
            let (k_idx, k_f64) = k_vals.unwrap();
            distances
                .par_chunks_mut(knn)
                .progress_with(progress_bar)
                .enumerate()
                .for_each(|(qi, row_dist_slice)| {
                    let mut heap = BinaryHeap::with_capacity(knn + 1);
                    let qi_sketch = query_sketches.get_sketch_slice(qi, k_idx);
                    for ri in 0..n {
                        let c1 = query_completeness_vec.map(|cv| cv[qi]);
                        let c2 = ref_completeness_vec.map(|cv| cv[ri]);
                        let dist = jaccard_index_generic::<BB>(
                            qi_sketch,
                            ref_sketches.get_sketch_slice(ri, k_idx),
                            ref_sketches.sketchsize64,
                            c1,
                            c2,
                            completeness_cutoff,
                        );
                        let dist_f32 = if ani {
                            (1.0_f64 - ani_pois(dist, k_f64)) as f32
                        } else {
                            (1.0_f64 - dist) as f32
                        };
                        push_heap(&mut heap, SparseJaccard(ri, dist_f32), knn);
                    }
                    if ani {
                        heap.into_sorted_vec().iter().zip(row_dist_slice).for_each(
                            |(inverse_ani, output_ani)| {
                                *output_ani = SparseJaccard(inverse_ani.0, 1.0_f32 - inverse_ani.1);
                            },
                        );
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
                .for_each(|(qi, row_dist_slice)| {
                    let mut heap = BinaryHeap::with_capacity(knn + 1);
                    for ri in 0..n {
                        let dists = core_acc_dist_generic::<BB>(
                            ref_sketches,
                            query_sketches,
                            ri,
                            qi,
                            ref_completeness_vec,
                            query_completeness_vec,
                            completeness_cutoff,
                        );
                        push_heap(&mut heap, SparseCoreAcc(ri, dists.0, dists.1), knn);
                    }
                    row_dist_slice.clone_from_slice(&heap.into_sorted_vec());
                });
        }
    }
    sp_distances
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
    completeness_vec: Option<&Vec<f64>>,
    completeness_cutoff: f64,
    retain_unmatched: &Option<RetainUnmatched>,
) -> SparseDistanceMatrix<'a> {
    // Check that sample sets in ski and skm are the same and create i,j lookup
    let mut skq_lookup = HashMap::with_capacity(n);
    for (skq_index, skq_sample) in inverted_index.sample_names().iter().enumerate() {
        skq_lookup.insert(skq_sample.as_str(), skq_index);
    }

    let mut not_found = Vec::new();
    // Map: skd index to ski index
    // Needed to convert i (from skd ordering) to ski ordering
    let mut skq_index_lookup = Vec::with_capacity(n);
    for skd_sample_idx in 0..sketches.number_samples_loaded() {
        let sample_name = sketches.sketch_name(skd_sample_idx);
        match skq_lookup.get(sample_name) {
            Some(skq_index) => skq_index_lookup.push(*skq_index),
            None => not_found.push(sample_name),
        };
    }
    if !not_found.is_empty() {
        panic!("The following samples in the .skd could not be found in the .ski:\n{not_found:?}");
    }

    // Reverse map: ski index to skd index
    // Needed to convert j (from inverted index, ski ordering) back to skd ordering
    let mut skd_index_from_ski: Vec<usize> = vec![0; n];
    for (skd_idx, &ski_idx) in skq_index_lookup.iter().enumerate() {
        skd_index_from_ski[ski_idx] = skd_idx;
    }

    let mut sp_distances = SparseDistanceMatrix::new(sketches, knn, dist_type);
    let k_vals = sp_distances.k_vals();
    let ani = sp_distances.ani();
    let progress_bar = get_progress_bar(n, BAR_PERCENT, quiet);
    match sp_distances.dists_mut() {
        DistVec::Jaccard(distances) => {
            let (k_idx, k_f64) = k_vals.unwrap();
            distances
                .par_chunks_mut(knn)
                .progress_with(progress_bar)
                .enumerate()
                .for_each(|(i, row_dist_slice)| {
                    // Prefilter step here
                    let skq_offset = skq_index_lookup[i] * skq_stride;
                    let flat_i_sketch = &skq_bins[skq_offset..(skq_offset + skq_stride)];
                    let prefiltered_samples = inverted_index.any_shared_bins(flat_i_sketch);
                    // Standard search
                    let mut heap = BinaryHeap::with_capacity(knn + 1);
                    let i_sketch = sketches.get_sketch_slice(i, k_idx);
                    for j in prefiltered_samples {
                        let j = j as usize;
                        if skq_index_lookup[i] == j {
                            continue;
                        }
                        let skd_j = skd_index_from_ski[j];
                        // If completeness_vec is Some, extract the value at index i (or j) from the inner vector.
                        // If completeness_vec is None, the result will also be None.
                        // This uses Option::map to safely access the completeness value for each sample.
                        let c1 = completeness_vec.map(|cv| cv[i]);
                        let c2 = completeness_vec.map(|cv| cv[skd_j]);
                        let dist = jaccard_index(
                            i_sketch,
                            sketches.get_sketch_slice(skd_j, k_idx),
                            sketches.sketchsize64,
                            c1,
                            c2,
                            completeness_cutoff,
                        );
                        let dist_f32 = if ani {
                            (1.0_f64 - ani_pois(dist, k_f64)) as f32
                        } else {
                            (1.0_f64 - dist) as f32
                        };
                        let dist_item = SparseJaccard(skd_j, dist_f32);
                        push_heap(&mut heap, dist_item, knn);
                    }
                    let mut dist_vec = heap.into_sorted_vec();

                    // Handle unmatched genomes (no prefiltered matches)
                    if dist_vec.is_empty() {
                        match retain_unmatched {
                            Some(RetainUnmatched::Singleton) => {
                                let mut singleton_vec = vec![SparseJaccard(i, 0.0)];
                                singleton_vec.append(&mut vec![
                                    SparseJaccard(i, 1.0);
                                    row_dist_slice.len() - 1
                                ]);
                                row_dist_slice.clone_from_slice(&singleton_vec);
                                return;
                            }
                            Some(RetainUnmatched::Bruteforce) => {
                                let mut bf_heap = BinaryHeap::with_capacity(knn + 1);
                                for j in 0..n {
                                    if i == j {
                                        continue;
                                    }
                                    let c1 = completeness_vec.map(|cv| cv[i]);
                                    let c2 = completeness_vec.map(|cv| cv[j]);
                                    let dist = jaccard_index(
                                        i_sketch,
                                        sketches.get_sketch_slice(j, k_idx),
                                        sketches.sketchsize64,
                                        c1,
                                        c2,
                                        completeness_cutoff,
                                    );
                                    let dist_f32 = if ani {
                                        (1.0_f64 - ani_pois(dist, k_f64)) as f32
                                    } else {
                                        (1.0_f64 - dist) as f32
                                    };
                                    let dist_item = SparseJaccard(j, dist_f32);
                                    push_heap(&mut bf_heap, dist_item, knn);
                                }
                                dist_vec = bf_heap.into_sorted_vec();
                            }
                            None => {}
                        }
                    }

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

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    /// Replays the chunk generation/stepping logic used by `self_dists_all_stream`
    /// (self/upper-triangle mode), without any sketch data, to verify the index
    /// arithmetic visits every pair exactly once.
    fn self_mode_chunk_pairs(n: usize) -> Vec<(usize, usize)> {
        let n_distances = n * n.saturating_sub(1) / 2;
        let n_chunks = n_distances.div_ceil(CHUNK_SIZE);
        let mut pairs = Vec::with_capacity(n_distances);
        for chunk_idx in 0..n_chunks {
            let start = chunk_idx * CHUNK_SIZE;
            let end = (start + CHUNK_SIZE).min(n_distances);
            let mut i = calc_row_idx(start, n);
            let mut j = calc_col_idx(start, i, n);
            for _ in start..end {
                pairs.push((i, j));
                j += 1;
                if j >= n {
                    i += 1;
                    j = i + 1;
                }
            }
        }
        pairs
    }

    /// Replays the chunk generation/stepping logic used by `cross_dists_all_stream`
    /// (rectangular ref-vs-query mode).
    fn cross_mode_chunk_pairs(n: usize, n_query: usize) -> Vec<(usize, usize)> {
        let n_distances = n * n_query;
        let n_chunks = n_distances.div_ceil(CHUNK_SIZE);
        let mut pairs = Vec::with_capacity(n_distances);
        for chunk_idx in 0..n_chunks {
            let start = chunk_idx * CHUNK_SIZE;
            let end = (start + CHUNK_SIZE).min(n_distances);
            let (mut i, mut j) = calc_query_indices(start, n_query);
            for _ in start..end {
                pairs.push((i, j));
                j += 1;
                if j >= n_query {
                    i += 1;
                    j = 0;
                }
            }
        }
        pairs
    }

    #[test]
    fn self_mode_chunk_boundaries_cover_every_pair_once() {
        // n values chosen to exercise: no pairs, a single pair, small n, and
        // n_distances landing exactly on / either side of the CHUNK_SIZE (1000)
        // boundary (n=46 -> 1035 pairs: one full chunk of 1000 + a 35-pair tail).
        for n in [0usize, 1, 2, 3, 4, 45, 46, 47, 63, 64, 100] {
            let pairs = self_mode_chunk_pairs(n);
            let expected_count = n * n.saturating_sub(1) / 2;
            assert_eq!(pairs.len(), expected_count, "wrong pair count for n={n}");

            let unique: HashSet<_> = pairs.iter().copied().collect();
            assert_eq!(
                unique.len(),
                pairs.len(),
                "duplicate pair detected for n={n}"
            );

            for &(i, j) in &pairs {
                assert!(i < j && j < n, "pair ({i},{j}) out of range for n={n}");
            }
            for i in 0..n {
                for j in (i + 1)..n {
                    assert!(unique.contains(&(i, j)), "missing pair ({i},{j}) for n={n}");
                }
            }
        }
    }

    #[test]
    fn cross_mode_chunk_boundaries_cover_every_pair_once() {
        for (n, n_query) in [
            (0usize, 0usize),
            (1, 1),
            (1, 5),
            (5, 1),
            (4, 4),
            (32, 32),
            (31, 33),
            (10, 100),
        ] {
            let pairs = cross_mode_chunk_pairs(n, n_query);
            let expected_count = n * n_query;
            assert_eq!(
                pairs.len(),
                expected_count,
                "wrong pair count for n={n}, n_query={n_query}"
            );

            let unique: HashSet<_> = pairs.iter().copied().collect();
            assert_eq!(
                unique.len(),
                pairs.len(),
                "duplicate pair detected for n={n}, n_query={n_query}"
            );

            for &(i, j) in &pairs {
                assert!(
                    i < n && j < n_query,
                    "pair ({i},{j}) out of range for n={n}, n_query={n_query}"
                );
            }
            for i in 0..n {
                for j in 0..n_query {
                    assert!(
                        unique.contains(&(i, j)),
                        "missing pair ({i},{j}) for n={n}, n_query={n_query}"
                    );
                }
            }
        }
    }
}
