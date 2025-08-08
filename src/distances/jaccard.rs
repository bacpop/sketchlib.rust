//! Implementation of Jaccard, core and accessory distance calculations
use crate::sketch::multisketch::MultiSketch;
use crate::sketch::BBITS;

/// Returns the Jaccard index between two samples
pub fn jaccard_index(sketch1: &[u64], sketch2: &[u64], sketchsize64: u64) -> f64 {
    let unionsize = (u64::BITS as u64 * sketchsize64) as f64;
    let samebits: u32 = sketch1
        .chunks_exact(BBITS as usize)
        .zip(sketch2.chunks_exact(BBITS as usize))
        .map(|(chunk1, chunk2)| {
            let mut bits: u64 = !0;
            chunk1.iter().zip(chunk2.iter()).for_each(|(&s1, &s2)| {
                bits &= !(s1 ^ s2);
            });
            bits.count_ones()
        })
        .sum();
    let maxnbits = sketchsize64 as u32 * u64::BITS;
    let expected_samebits = maxnbits >> BBITS;

    log::trace!("samebits:{samebits} expected_samebits:{expected_samebits} maxnbits:{maxnbits}");
    let diff = samebits.saturating_sub(expected_samebits);
    let intersize = (diff as f64 * maxnbits as f64) / (maxnbits - expected_samebits) as f64;
    log::trace!("intersize:{intersize} unionsize:{unionsize}");
    intersize / unionsize
}

/// Converts between Jaccard distance and ANI, using a Poisson model of mutations
#[inline(always)]
pub fn ani_pois(jaccard: f64, k: f64) -> f64 {
    0.0_f64.max(1.0 + 1.0 / k * (((2.0 * jaccard) / (1.0 + jaccard)).ln()))
}

/// Completeness correction for MAGs
#[inline(always)]
pub fn completeness_correction(jaccard: f64, c1: &f64, c2: &f64) -> f64 {
    jaccard / (c1 * c2 / (c1 + c2 - c1 * c2))
}

/// Core and accessory distances between two sketches, using the PopPUNK regression
/// model
pub fn core_acc_dist(
    ref_sketches: &MultiSketch,
    query_sketches: &MultiSketch,
    ref_sketch_idx: usize,
    query_sketch_idx: usize,
) -> (f32, f32) {
    if ref_sketches.kmer_lengths().len() < 2 {
        panic!("Need at least two k-mer lengths to calculate core/accessory distances");
    }
    let (mut xsum, mut ysum, mut xysum, mut xsquaresum, mut ysquaresum, mut n) =
        (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
    let tolerance = (2.0_f64 / ((ref_sketches.sketch_size * u64::BITS as u64) as f64)).ln();
    //let tolerance = -100.0_f32;
    for (k_idx, k) in ref_sketches.kmer_lengths().iter().enumerate() {
        let y = jaccard_index(
            ref_sketches.get_sketch_slice(ref_sketch_idx, k_idx),
            query_sketches.get_sketch_slice(query_sketch_idx, k_idx),
            ref_sketches.sketchsize64,
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
