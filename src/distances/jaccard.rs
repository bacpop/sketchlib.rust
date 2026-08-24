//! Implementation of Jaccard, core and accessory distance calculations
use crate::sketch::multisketch::MultiSketch;
use crate::sketch::{BIN_BITS, LEGACY_BIN_BITS};

pub(super) fn jaccard_index_generic<const BB: usize>(
    sketch1: &[u64],
    sketch2: &[u64],
    sketchsize64: u64,
    c1: Option<f64>,
    c2: Option<f64>,
    completeness_cutoff: f64,
) -> f64 {
    let samebits = same_bits_dispatch::<BB>(sketch1, sketch2) as f64;
    let unionsize = (u64::BITS as u64 * sketchsize64) as f64;
    // Correction for random matches
    let expected_random = unionsize / (1u64 << BB) as f64;
    let mut jaccard_index =
        ((samebits - expected_random) / (unionsize - expected_random)).clamp(0.0, 1.0);

    log::trace!("samebits:{samebits} expected_random:{expected_random} jaccard:{jaccard_index}");

    // Apply completeness correction if both completeness values are provided
    if let (Some(c1_val), Some(c2_val)) = (c1, c2) {
        if c1_val * c2_val >= completeness_cutoff {
            // Cap the corrected Jaccard index at 1.0 to prevent negative distances
            jaccard_index = completeness_correction(jaccard_index, c1_val, c2_val).min(1.0);
        }
    }

    jaccard_index
}

/// Returns the Jaccard index between two samples
pub fn jaccard_index(
    sketch1: &[u64],
    sketch2: &[u64],
    sketchsize64: u64,
    c1: Option<f64>,
    c2: Option<f64>,
    completeness_cutoff: f64,
) -> f64 {
    jaccard_index_generic::<BIN_BITS>(sketch1, sketch2, sketchsize64, c1, c2, completeness_cutoff)
}

/// Legacy analogue of [`jaccard_index`], for pre-v0.4 (`BIN_BITS` = 14)
/// databases — see [`LEGACY_BIN_BITS`] and
/// [`crate::sketch::multisketch::MultiSketch::is_legacy_format`].
/// Random-match correction is applied identically to `jaccard_index`, just
/// using the 14-bit random-match denominator. Dispatched automatically by
/// `distances::self_dists_*`/`distances::cross_dists_*` based on the loaded
/// database's format.
// `distances::mod` dispatches via `jaccard_index_generic::<BB>` directly
// (no extra monomorphization hop in the per-pair hot loop), so this named
// entry point is only exercised by tests — kept for API discoverability.
#[allow(dead_code)]
pub(crate) fn jaccard_index_legacy(
    sketch1: &[u64],
    sketch2: &[u64],
    sketchsize64: u64,
    c1: Option<f64>,
    c2: Option<f64>,
    completeness_cutoff: f64,
) -> f64 {
    jaccard_index_generic::<LEGACY_BIN_BITS>(
        sketch1,
        sketch2,
        sketchsize64,
        c1,
        c2,
        completeness_cutoff,
    )
}

/// Computes the "same bits" popcount for a given bin-packing width `BB`.
/// Dispatches to the NEON-accelerated path only when `BB == BIN_BITS` (the
/// current/new format); every other width — in practice only
/// [`LEGACY_BIN_BITS`] — always uses the portable scalar implementation, on
/// every architecture including aarch64. `BB == BIN_BITS` is resolved per
/// monomorphization, so for the legacy instantiation this is unconditionally
/// `false` and `jaccard_same_bits`/the NEON kernel is never invoked.
#[inline(always)]
fn same_bits_dispatch<const BB: usize>(sketch1: &[u64], sketch2: &[u64]) -> u32 {
    if BB == BIN_BITS {
        jaccard_same_bits(sketch1, sketch2)
    } else {
        jaccard_same_bits_general::<BB>(sketch1, sketch2)
    }
}

#[inline(always)]
fn jaccard_same_bits(sketch1: &[u64], sketch2: &[u64]) -> u32 {
    #[cfg(target_arch = "aarch64")]
    unsafe {
        jaccard_neon_unroll2_inner(sketch1, sketch2)
    }

    #[cfg(not(target_arch = "aarch64"))]
    {
        jaccard_same_bits_general::<BIN_BITS>(sketch1, sketch2)
    }
}

#[inline(always)]
pub(crate) fn jaccard_same_bits_general<const BB: usize>(sketch1: &[u64], sketch2: &[u64]) -> u32 {
    debug_assert_eq!(sketch1.len(), sketch2.len());
    debug_assert_eq!(sketch1.len() % BB, 0);
    sketch1
        .as_chunks::<BB>()
        .0
        .iter()
        .zip(sketch2.as_chunks::<BB>().0.iter())
        .map(|(chunk1, chunk2)| {
            let mut bits: u64 = !0;
            chunk1.iter().zip(chunk2.iter()).for_each(|(&s1, &s2)| {
                bits &= !(s1 ^ s2);
            });
            bits.count_ones()
        })
        .sum()
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn jaccard_neon_unroll2_inner(a: &[u64], b: &[u64]) -> u32 {
    use std::arch::aarch64::*;

    debug_assert_eq!(a.len(), b.len());
    debug_assert_eq!(a.len() % BIN_BITS, 0);

    let chunk = BIN_BITS;
    let n_chunks = a.len() / chunk;
    let n_pairs = n_chunks / 2;
    let mut total = 0u32;

    #[inline(always)]
    unsafe fn one_chunk(ap: *const u8, bp: *const u8) -> u8 {
        let a0 = vld1q_u8(ap);
        let a1 = vld1q_u8(ap.add(16));
        let a2 = vld1q_u8(ap.add(32));
        let a3 = vld1q_u8(ap.add(48));
        let a4 = vld1q_u8(ap.add(64));
        let a5 = vld1q_u8(ap.add(80));
        let a6 = vld1q_u8(ap.add(96));
        let a7 = vld1q_u8(ap.add(112));
        let b0 = vld1q_u8(bp);
        let b1 = vld1q_u8(bp.add(16));
        let b2 = vld1q_u8(bp.add(32));
        let b3 = vld1q_u8(bp.add(48));
        let b4 = vld1q_u8(bp.add(64));
        let b5 = vld1q_u8(bp.add(80));
        let b6 = vld1q_u8(bp.add(96));
        let b7 = vld1q_u8(bp.add(112));
        let x0 = veorq_u8(a0, b0);
        let x1 = veorq_u8(a1, b1);
        let x2 = veorq_u8(a2, b2);
        let x3 = veorq_u8(a3, b3);
        let x4 = veorq_u8(a4, b4);
        let x5 = veorq_u8(a5, b5);
        let x6 = veorq_u8(a6, b6);
        let x7 = veorq_u8(a7, b7);
        let or01 = vorrq_u8(x0, x1);
        let or23 = vorrq_u8(x2, x3);
        let or45 = vorrq_u8(x4, x5);
        let or67 = vorrq_u8(x6, x7);
        let or_all = vorrq_u8(vorrq_u8(or01, or23), vorrq_u8(or45, or67));
        let not_or = vmvnq_u8(or_all);
        vaddv_u8(vcnt_u8(vand_u8(vget_low_u8(not_or), vget_high_u8(not_or))))
    }

    for i in 0..n_pairs {
        let base = i * 2 * chunk;
        total += one_chunk(a.as_ptr().add(base) as _, b.as_ptr().add(base) as _) as u32;
        total += one_chunk(
            a.as_ptr().add(base + chunk) as _,
            b.as_ptr().add(base + chunk) as _,
        ) as u32;
    }
    if !n_chunks.is_multiple_of(2) {
        let base = (n_chunks - 1) * chunk;
        total += one_chunk(a.as_ptr().add(base) as _, b.as_ptr().add(base) as _) as u32;
    }
    total
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scalar_same_bits_counts_matching_bins() {
        let sketch1 = [u64::MAX; BIN_BITS * 2];
        let mut sketch2 = [u64::MAX; BIN_BITS * 2];
        sketch2[BIN_BITS] = 0;

        assert_eq!(
            jaccard_same_bits_general::<BIN_BITS>(&sketch1, &sketch1),
            128
        );
        assert_eq!(
            jaccard_same_bits_general::<BIN_BITS>(&sketch1, &sketch2),
            64
        );
    }

    #[cfg(target_arch = "aarch64")]
    #[test]
    fn neon_same_bits_matches_scalar() {
        let mut sketch1 = [0u64; BIN_BITS * 3];
        let mut sketch2 = [0u64; BIN_BITS * 3];
        for i in 0..sketch1.len() {
            sketch1[i] = (i as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
            sketch2[i] = sketch1[i] ^ ((i as u64) << (i % 17));
        }

        let scalar = jaccard_same_bits_general::<BIN_BITS>(&sketch1, &sketch2);
        let neon = unsafe { jaccard_neon_unroll2_inner(&sketch1, &sketch2) };
        assert_eq!(neon, scalar);
    }

    #[test]
    fn random_match_correction_applied_for_large_sketch_size() {
        // sketchsize64 chosen so that unionsize = sketchsize64 * 64 exceeds
        // 2^BIN_BITS, the regime the random-match correction is meant for.
        let sketchsize64: u64 = 2048;
        let unionsize = (u64::BITS as u64 * sketchsize64) as f64;
        assert!(unionsize > (1u64 << BIN_BITS) as f64);

        let n_chunks = sketchsize64 as usize;
        let sketch1 = vec![u64::MAX; BIN_BITS * n_chunks];
        let mut sketch2 = sketch1.clone();
        // Mismatch exactly one chunk (one 64-bit "column") out of n_chunks.
        for word in sketch2[..BIN_BITS].iter_mut() {
            *word = 0;
        }

        let samebits = jaccard_same_bits_general::<BIN_BITS>(&sketch1, &sketch2) as f64;
        let expected_random = unionsize / (1u64 << BIN_BITS) as f64;
        let expected_jaccard =
            ((samebits - expected_random) / (unionsize - expected_random)).clamp(0.0, 1.0);
        // The correction should be non-trivial (i.e. not equal to the raw,
        // uncorrected ratio) for this sketch size.
        assert!((expected_jaccard - samebits / unionsize).abs() > f64::EPSILON);

        let jaccard = jaccard_index(&sketch1, &sketch2, sketchsize64, None, None, 0.0);
        assert_eq!(jaccard, expected_jaccard);

        // Identical sketches should give a jaccard index of exactly 1.0.
        let identical = jaccard_index(&sketch1, &sketch1, sketchsize64, None, None, 0.0);
        assert_eq!(identical, 1.0);
    }

    #[test]
    fn legacy_scalar_same_bits_counts_matching_bins() {
        // Mirrors `scalar_same_bits_counts_matching_bins`, but for the legacy
        // 14-bit bin-packing width. Chunk width affects how many words feed
        // into each chunk, not the per-chunk `count_ones()` range, so the
        // expected popcounts are numerically identical to the BIN_BITS case.
        let sketch1 = [u64::MAX; LEGACY_BIN_BITS * 2];
        let mut sketch2 = [u64::MAX; LEGACY_BIN_BITS * 2];
        sketch2[LEGACY_BIN_BITS] = 0;

        assert_eq!(
            jaccard_same_bits_general::<LEGACY_BIN_BITS>(&sketch1, &sketch1),
            128
        );
        assert_eq!(
            jaccard_same_bits_general::<LEGACY_BIN_BITS>(&sketch1, &sketch2),
            64
        );
    }

    #[test]
    fn legacy_random_match_correction_applied_for_large_sketch_size() {
        // Mirrors `random_match_correction_applied_for_large_sketch_size`,
        // but exercising `jaccard_index_legacy` (LEGACY_BIN_BITS = 14).
        let sketchsize64: u64 = 2048;
        let unionsize = (u64::BITS as u64 * sketchsize64) as f64;
        assert!(unionsize > (1u64 << LEGACY_BIN_BITS) as f64);

        let n_chunks = sketchsize64 as usize;
        let sketch1 = vec![u64::MAX; LEGACY_BIN_BITS * n_chunks];
        let mut sketch2 = sketch1.clone();
        for word in sketch2[..LEGACY_BIN_BITS].iter_mut() {
            *word = 0;
        }

        let samebits = jaccard_same_bits_general::<LEGACY_BIN_BITS>(&sketch1, &sketch2) as f64;
        let expected_random = unionsize / (1u64 << LEGACY_BIN_BITS) as f64;
        let expected_jaccard =
            ((samebits - expected_random) / (unionsize - expected_random)).clamp(0.0, 1.0);
        assert!((expected_jaccard - samebits / unionsize).abs() > f64::EPSILON);

        let jaccard = jaccard_index_legacy(&sketch1, &sketch2, sketchsize64, None, None, 0.0);
        assert_eq!(jaccard, expected_jaccard);

        let identical = jaccard_index_legacy(&sketch1, &sketch1, sketchsize64, None, None, 0.0);
        assert_eq!(identical, 1.0);
    }
}

/// Converts between Jaccard distance and ANI, using a Poisson model of mutations
#[inline(always)]
pub fn ani_pois(jaccard: f64, k: f64) -> f64 {
    0.0_f64.max(1.0 + 1.0 / k * (((2.0 * jaccard) / (1.0 + jaccard)).ln()))
}

/// Completeness correction for MAGs
#[inline(always)]
pub fn completeness_correction(jaccard: f64, c1: f64, c2: f64) -> f64 {
    jaccard / (c1 * c2 / (c1 + c2 - c1 * c2))
}

pub(super) fn core_acc_dist_generic<const BB: usize>(
    ref_sketches: &MultiSketch,
    query_sketches: &MultiSketch,
    ref_sketch_idx: usize,
    query_sketch_idx: usize,
    ref_completeness_vec: Option<&Vec<f64>>,
    query_completeness_vec: Option<&Vec<f64>>,
    completeness_cutoff: f64,
) -> (f32, f32) {
    if ref_sketches.kmer_lengths().len() < 2 {
        panic!("Need at least two k-mer lengths to calculate core/accessory distances");
    }
    let (mut xsum, mut ysum, mut xysum, mut xsquaresum, mut ysquaresum, mut n) =
        (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
    let tolerance = (2.0_f64 / ((ref_sketches.sketch_size * u64::BITS as u64) as f64)).ln();
    //let tolerance = -100.0_f32;
    for (k_idx, k) in ref_sketches.kmer_lengths().iter().enumerate() {
        let c1 = ref_completeness_vec.map(|cv| cv[ref_sketch_idx]);
        let c2 = query_completeness_vec.map(|cv| cv[query_sketch_idx]);
        let y = jaccard_index_generic::<BB>(
            ref_sketches.get_sketch_slice(ref_sketch_idx, k_idx),
            query_sketches.get_sketch_slice(query_sketch_idx, k_idx),
            ref_sketches.sketchsize64,
            c1,
            c2,
            completeness_cutoff,
        )
        .ln();
        if y < tolerance {
            break;
        }
        let k_fl = *k as f64;
        xsum += k_fl;
        ysum += y;
        xysum += k_fl * y;
        xsquaresum += k_fl * k_fl;
        ysquaresum += y * y;
        n += 1.0;
    }
    simple_linear_regression(xsum, ysum, xysum, xsquaresum, ysquaresum, n)
}

/// Core and accessory distances between two sketches, using the PopPUNK regression
/// model
// `distances::mod` dispatches via `core_acc_dist_generic::<BB>` directly;
// kept as a named, behaviorally-unchanged entry point.
#[allow(dead_code)]
pub fn core_acc_dist(
    ref_sketches: &MultiSketch,
    query_sketches: &MultiSketch,
    ref_sketch_idx: usize,
    query_sketch_idx: usize,
    ref_completeness_vec: Option<&Vec<f64>>,
    query_completeness_vec: Option<&Vec<f64>>,
    completeness_cutoff: f64,
) -> (f32, f32) {
    core_acc_dist_generic::<BIN_BITS>(
        ref_sketches,
        query_sketches,
        ref_sketch_idx,
        query_sketch_idx,
        ref_completeness_vec,
        query_completeness_vec,
        completeness_cutoff,
    )
}

/// Legacy analogue of [`core_acc_dist`] — see [`jaccard_index_legacy`].
#[allow(dead_code)]
pub(crate) fn core_acc_dist_legacy(
    ref_sketches: &MultiSketch,
    query_sketches: &MultiSketch,
    ref_sketch_idx: usize,
    query_sketch_idx: usize,
    ref_completeness_vec: Option<&Vec<f64>>,
    query_completeness_vec: Option<&Vec<f64>>,
    completeness_cutoff: f64,
) -> (f32, f32) {
    core_acc_dist_generic::<LEGACY_BIN_BITS>(
        ref_sketches,
        query_sketches,
        ref_sketch_idx,
        query_sketch_idx,
        ref_completeness_vec,
        query_completeness_vec,
        completeness_cutoff,
    )
}

// Linear regression for calculating core/accessory distances from matches, with some
// sensible bounds for bad fits
fn simple_linear_regression(
    xsum: f64,
    ysum: f64,
    xysum: f64,
    xsquaresum: f64,
    ysquaresum: f64,
    n: f64,
) -> (f32, f32) {
    log::trace!(
        "xsum:{xsum} ysum:{ysum} xysum:{xysum} xsquaresum:{xsquaresum} ysquaresum:{ysquaresum}"
    );
    // No matches
    if ysum.is_nan() || ysum == f64::NEG_INFINITY || n < 3.0 {
        return (1.0, 1.0);
    }

    let xbar = xsum / n;
    let ybar = ysum / n;
    let x_diff = xsquaresum - xsum * xsum / n;
    let y_diff = ysquaresum - ysum * ysum / n;
    let xstddev = ((xsquaresum - xsum * xsum / n) / n).sqrt();
    let ystddev = ((ysquaresum - ysum * ysum / n) / n).sqrt();
    let r = (xysum - xsum * ysum / n) / (x_diff * y_diff).sqrt();
    let beta = r * ystddev / xstddev;
    let alpha = -beta * xbar + ybar;
    log::trace!("r:{r} alpha:{alpha} beta:{beta}");

    let (mut core, mut acc) = (0.0_f64, 0.0_f64);
    if beta < 0.0 {
        core = 1.0 - beta.exp();
    }
    if alpha < 0.0 {
        acc = 1.0 - alpha.exp();
    }
    (core as f32, acc as f32)
}
