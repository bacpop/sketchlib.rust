//! Implementation of Jaccard, core and accessory distance calculations
use crate::sketch::multisketch::MultiSketch;

/// Returns the Jaccard index between two samples
pub fn jaccard_index(
    sketch1: &[u16],
    sketch2: &[u16],
    sketchsize: u64,
    c1: Option<f64>,
    c2: Option<f64>,
    completeness_cutoff: f64,
) -> f64 {
    // Ragnar distance, could be optimised
    let both_empty = std::iter::zip(sketch1, sketch2).map(|(a, b)| ((a | b) == 0) as u32).sum::<u32>() as f64;
    let mut jaccard_index = 1.0 - std::iter::zip(sketch1, sketch2)
        .map(|(a, b)| (a != b) as u32)
        .sum::<u32>() as f64
        / (sketchsize as f64 - both_empty);
    
    let bb = (1usize << 16) as f64;

    // Correction for accidental matches.
    // Take a max with 0 to avoid correcting into a negative jaccard similarity
    // for uncorrelated sketches.
    jaccard_index = (bb * jaccard_index - 1.0).max(0.0) / (bb - 1.0);

    // John's sketchlib-distance
    
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
    let tolerance = (2.0_f64 / ((ref_sketches.sketch_size * u64::BITS as u64) as f64)).ln(); // ASK JOHN
    //let tolerance = -100.0_f32;
    for (k_idx, k) in ref_sketches.kmer_lengths().iter().enumerate() {
        let c1 = completeness_vec.map(|cv| cv[ref_sketch_idx]);
        let c2 = completeness_vec.map(|cv| cv[query_sketch_idx]);
        let y = jaccard_index(
            ref_sketches.get_sketch_slice(ref_sketch_idx, k_idx),
            query_sketches.get_sketch_slice(query_sketch_idx, k_idx),
            ref_sketches.sketch_size,
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
