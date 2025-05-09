//! Implementation of Jaccard, core and accessory distance calculations
use crate::sketch::multisketch::MultiSketch;
use crate::sketch::BBITS;

/// Returns the Jaccard index between two samples
pub fn jaccard_index(sketch1: &[u64], sketch2: &[u64], sketchsize64: u64) -> f32 {
    let unionsize = (u64::BITS as u64 * sketchsize64) as f32;
    let mut samebits: u32 = 0;
    unsafe {
        for i in 0..sketchsize64 {
            let mut bits: u64 = !0;
            for j in 0..BBITS {
                bits &= !(sketch1.get_unchecked((i * BBITS + j) as usize) ^ sketch2.get_unchecked((i * BBITS + j) as usize));
            }
            samebits += bits.count_ones();
        }
    }
    let maxnbits = sketchsize64 as u32 * u64::BITS;
    let expected_samebits = maxnbits >> BBITS;

    log::trace!("samebits:{samebits} expected_samebits:{expected_samebits} maxnbits:{maxnbits}");
    let diff = samebits.saturating_sub(expected_samebits);
    // Do float multiplication here. diff * maxnbits for large s will overflow u32
    // f32 is sufficient for s=10M (I think f64 would get to one bit diff precision larger than this)
    let intersize = (diff as f64 * maxnbits as f64) / (maxnbits - expected_samebits) as f64;
    log::trace!("intersize:{intersize} unionsize:{unionsize}");
    intersize as f32 / unionsize
}

/// Returns the Jaccard index between two samples
pub fn jaccard_index2(sketch1: &[u16], sketch2: &[u16], sketchsize: usize) -> f32 {
    let unionsize = sketchsize as f32;
    let mut samebits: u32 = 0;
    for i in 0..sketchsize {
        if sketch1[i] == sketch2[i] {
            samebits += 1;
        }
    }
    let maxnbits = sketchsize as u32;
    let expected_samebits = maxnbits >> 16;

    log::trace!("samebits:{samebits} expected_samebits:{expected_samebits} maxnbits:{maxnbits}");
    let diff = samebits.saturating_sub(expected_samebits);
    // Do float multiplication here. diff * maxnbits for large s will overflow u32
    // f32 is sufficient for s=10M (I think f64 would get to one bit diff precision larger than this)
    let intersize = (diff as f64 * maxnbits as f64) / (maxnbits - expected_samebits) as f64;
    log::trace!("intersize:{intersize} unionsize:{unionsize}");
    intersize as f32 / unionsize
}

/// Converts between Jaccard distance and ANI, using a Poisson model of mutations
#[inline(always)]
pub fn ani_pois(jaccard: f32, k: f32) -> f32 {
    0.0_f32.max(1.0 + 1.0 / k * (((2.0 * jaccard) / (1.0 + jaccard)).ln()))
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
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    let tolerance = (2.0_f32 / ((ref_sketches.sketch_size * u64::BITS as u64) as f32)).ln();
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
        let k_fl = *k as f32;
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
    xsum: f32,
    ysum: f32,
    xysum: f32,
    xsquaresum: f32,
    ysquaresum: f32,
    n: f32,
) -> (f32, f32) {
    log::trace!(
        "xsum:{xsum} ysum:{ysum} xysum:{xysum} xsquaresum:{xsquaresum} ysquaresum:{ysquaresum}"
    );
    // No matches
    if ysum.is_nan() || ysum == f32::NEG_INFINITY || n < 3.0 {
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

    let (mut core, mut acc) = (0.0, 0.0);
    if beta < 0.0 {
        core = 1.0 - beta.exp();
    } else if r > 0.0 {
        core = 1.0;
    }
    if alpha < 0.0 {
        acc = 1.0 - alpha.exp();
    }
    (core, acc)
}
