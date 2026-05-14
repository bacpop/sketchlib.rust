//! Implementation of Jaccard, core and accessory distance calculations
use crate::sketch::multisketch::MultiSketch;
use crate::sketch::BIN_BITS;

/// Returns the Jaccard index between two samples
pub fn jaccard_index(
    sketch1: &[u64],
    sketch2: &[u64],
    sketchsize64: u64,
    c1: Option<f64>,
    c2: Option<f64>,
    completeness_cutoff: f64,
) -> f64 {
    let unionsize = (u64::BITS as u64 * sketchsize64) as f64;
    let samebits = jaccard_same_bits(sketch1, sketch2);
    let raw_jaccard = samebits as f64 / unionsize;
    let mut jaccard_index = random_match_correction(raw_jaccard);

    log::trace!("samebits:{samebits} raw_jaccard:{raw_jaccard} jaccard:{jaccard_index}");

    // Apply completeness correction if both completeness values are provided
    if let (Some(c1_val), Some(c2_val)) = (c1, c2) {
        if c1_val * c2_val >= completeness_cutoff {
            jaccard_index = completeness_correction(jaccard_index, c1_val, c2_val);
            // Cap the corrected Jaccard index at 1.0 to prevent negative distances
            jaccard_index = jaccard_index.min(1.0);
        }
    }

    jaccard_index
}

#[inline(always)]
fn jaccard_same_bits(sketch1: &[u64], sketch2: &[u64]) -> u32 {
    #[cfg(target_arch = "aarch64")]
    unsafe {
        jaccard_neon_unroll2_inner(sketch1, sketch2)
    }

    #[cfg(not(target_arch = "aarch64"))]
    {
        jaccard_same_bits_general(sketch1, sketch2)
    }
}

#[cfg_attr(target_arch = "aarch64", allow(dead_code))]
pub(crate) fn jaccard_same_bits_general(sketch1: &[u64], sketch2: &[u64]) -> u32 {
    debug_assert_eq!(sketch1.len(), sketch2.len());
    debug_assert_eq!(sketch1.len() % BIN_BITS, 0);
    sketch1
        .chunks_exact(BIN_BITS)
        .zip(sketch2.chunks_exact(BIN_BITS))
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
    if n_chunks % 2 != 0 {
        let base = (n_chunks - 1) * chunk;
        total += one_chunk(a.as_ptr().add(base) as _, b.as_ptr().add(base) as _) as u32;
    }
    total
}

#[inline(always)]
fn random_match_correction(jaccard: f64) -> f64 {
    let bb = (1u32 << (BIN_BITS as u32)) as f64;
    ((bb * jaccard - 1.0).max(0.0)) / (bb - 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn random_match_correction_maps_expected_random_to_zero() {
        let expected_random = 1.0 / ((1u32 << (BIN_BITS as u32)) as f64);
        assert_eq!(random_match_correction(expected_random), 0.0);
        assert_eq!(random_match_correction(1.0), 1.0);
    }

    #[test]
    fn scalar_same_bits_counts_matching_bins() {
        let sketch1 = [u64::MAX; BIN_BITS * 2];
        let mut sketch2 = [u64::MAX; BIN_BITS * 2];
        sketch2[BIN_BITS] = 0;

        assert_eq!(jaccard_same_bits_general(&sketch1, &sketch1), 128);
        assert_eq!(jaccard_same_bits_general(&sketch1, &sketch2), 64);
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

        let scalar = jaccard_same_bits_general(&sketch1, &sketch2);
        let neon = unsafe { jaccard_neon_unroll2_inner(&sketch1, &sketch2) };
        assert_eq!(neon, scalar);
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

/// Core and accessory distances between two sketches, using the PopPUNK regression
/// model
pub fn core_acc_dist(
    ref_sketches: &MultiSketch,
    query_sketches: &MultiSketch,
    ref_sketch_idx: usize,
    query_sketch_idx: usize,
    completeness_vec: Option<&Vec<f64>>,
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
        let c1 = completeness_vec.map(|cv| cv[ref_sketch_idx]);
        let c2 = completeness_vec.map(|cv| cv[query_sketch_idx]);
        let y = jaccard_index(
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
    } else if r > 0.0 {
        core = 1.0;
    }
    if alpha < 0.0 {
        acc = 1.0 - alpha.exp();
    }
    (core as f32, acc as f32)
}
