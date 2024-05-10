use crate::multisketch::MultiSketch;
use crate::sketch::BBITS;

pub fn jaccard_dist(sketch1: &[u64], sketch2: &[u64], sketch_size: u64) -> f32 {
    let unionsize = (u64::BITS as u64 * sketch_size) as f32;
    let mut samebits: u32 = 0;
    for i in 0..sketch_size {
        let mut bits: u64 = !0;
        for j in 0..BBITS {
            bits &= !(sketch1[(i * BBITS + j) as usize] ^ sketch2[(i * BBITS + j) as usize]);
        }
        samebits += bits.count_ones();
    }
    let maxnbits = sketch_size as u32 * u64::BITS;
    let expected_samebits = maxnbits >> BBITS;

    log::trace!("samebits:{samebits} expected_samebits:{expected_samebits} maxnbits:{maxnbits}");
    if expected_samebits != 0 {
        samebits as f32
    } else {
        let diff = samebits.saturating_sub(expected_samebits);
        let intersize = (diff * maxnbits) as f32 / (maxnbits - expected_samebits) as f32;
        log::trace!("intersize:{intersize} unionsize:{unionsize}");
        intersize / unionsize
    }
}

#[inline(always)]
pub fn ani_pois(jaccard: f32, k: f32) -> f32 {
    0.0_f32.max(1.0 + 1.0 / k * (((2.0 * jaccard) / (1.0 + jaccard)).ln()))
}

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
    let tolerance = (5.0_f32 / (ref_sketches.sketch_size as f32)).ln();
    for (k_idx, k) in ref_sketches.kmer_lengths().iter().enumerate() {
        let y = jaccard_dist(
            ref_sketches.get_sketch_slice(ref_sketch_idx, k_idx),
            query_sketches.get_sketch_slice(query_sketch_idx, k_idx),
            ref_sketches.sketch_size,
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

    let (mut core, mut acc) = (1.0, 1.0);
    if beta < 0.0 {
        core = 1.0 - beta.exp();
    }
    if alpha < 0.0 {
        acc = 1.0 - alpha.exp();
    }
    (core, acc)
}
