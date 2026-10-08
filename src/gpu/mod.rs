use std::fmt::Write as FmtWrite;
use std::io::Write as IoWrite;

use anyhow::Error;

use cuda_core::{CudaContext, DeviceBuffer, LaunchConfig1D, simt::LaunchConfig};
use cuda_device::{DisjointSlice, kernel, launch_bounds, launch_contract, thread, device, gpu_printf};
use cuda_host::cuda_module;

use crate::sketch::multisketch::MultiSketch;
use crate::sketch::BIN_BITS;
use crate::distances::distance_matrix::*;
use crate::distances::jaccard::ani_pois;

const BBITS: u64 = BIN_BITS as u64;

#[device]
fn jaccard_index(
    sketch1: &[u64],
    sketch2: &[u64],
    sketch1_stride: u64,
    sketch2_stride: u64,
    sketchsize64: u64,
) -> f64 {
    let unionsize = (u64::BITS as u64 * sketchsize64) as f64;
    let mut samebits: u32 = 0;
    for i in 0..sketchsize64 {
        let mut bits: u64 = !0;
        for j in 0..BBITS {
            bits &= !(sketch1[((i * BBITS + j) * sketch1_stride) as usize] ^ sketch2[((i * BBITS + j) * sketch2_stride) as usize]);
        }
        samebits += bits.count_ones();
    }
    // Correction for random matches
    let expected_random = unionsize / (1u64 << BBITS) as f64;
    let jaccard_index =
        ((samebits as f64 - expected_random) / (unionsize - expected_random)).clamp(0.0, 1.0);

    jaccard_index
}

#[cuda_module]
mod kernels {
    use super::*;

    #[kernel]
    //#[launch_bounds(256)]
    //#[launch_contract(domain = 1, block = (256, 1, 1))]
    pub fn dists(
        sketches: &[u64],
        sketch_strides: (usize, usize, usize), // (sample, k-mer, bin)
        sketchsize64: u64,
        kmers: &[usize],
        dist_type: DistType,
        k_vals: Option<(usize, f64)>,
        ani: bool,
        mut dist: DisjointSlice<f32>
    ) {
        let idx = thread::index_1d();
        let idx_raw = idx.get();

        // let i = calc_row_idx(idx_raw, dist.len());

        let n = dist.len();

        let k_i64 = idx_raw as i64;
        let n_i64 = n as i64;

        let square = (-8 * k_i64 + 4 * n_i64 * (n_i64 - 1) - 7) as f64;
        let sqrt_square = (square).sqrt();
        let sqrt_square_floor = sqrt_square / 2.0 - 0.5;

        let i = n - 2 - (sqrt_square_floor).floor() as usize;


        let j = calc_col_idx(idx_raw, i, dist.len());

        gpu_printf!("{}, {}, {}, {}\n", square, sqrt_square, sqrt_square_floor, i);

        //if idx_raw == 0 {
            //gpu_printf!("{}\n", dist.len());
        //}

        if let Some((k_idx, k_f64)) = k_vals {
            let i_idx = i * sketch_strides.0 + k_idx * sketch_strides.1;
            let j_idx = j * sketch_strides.0 + k_idx * sketch_strides.1;

            let j_index = jaccard_index(&sketches[i_idx..], &sketches[j_idx..], sketch_strides.2 as u64, sketch_strides.2 as u64, sketchsize64);

            if let Some(dist_elem) = dist.get_mut(idx) {
                *dist_elem = if ani {
                    ani_pois(j_index, k_f64) as f32
                } else {
                    (1.0_f64 - j_index) as f32
                };
            }
        } else {
            unimplemented!("path iterating over k-mers and doing regression")
        }

    }
}

pub fn self_dists_all_stream<'a, W: IoWrite + Send>(
    writer: &mut W,
    sketches: &'a MultiSketch,
    n: usize,
    dist_type: DistType,
    quiet: bool,
    completeness_vec: Option<&Vec<f64>>,
    completeness_cutoff: f64,
    gpu_device_id: usize,
) -> Result<DistanceMatrix<'a>, Error> {
    let ctx = CudaContext::new(gpu_device_id)?;
    let stream = ctx.default_stream();

    let mut distances = DistanceMatrix::new(sketches, None, dist_type);
    let k_vals = distances.k_vals();
    let ani = distances.ani();
    let n_distances = distances.n_distances;

    let sketch_dev = DeviceBuffer::from_host(&stream, sketches.get_sketch_all())?;
    let kmers_dev = DeviceBuffer::from_host(&stream, sketches.kmer_lengths())?;

    let mut dist_dev = DeviceBuffer::<f32>::zeroed(&stream, n_distances)?;

    // SAFETY: this package owns the embedded device bundle produced for the
    // kernels module above.
    let module = unsafe { kernels::load(&ctx)? };

    //let prepared = module.prepare_dists(
        //LaunchConfig1D::new((n_distances as u32).div_ceil(256), 256, 0)
    //)?;

    let launch_config = LaunchConfig::for_num_elems(n_distances as u32);

    unsafe {
        module.dists(
            &stream,
            //&prepared,
            launch_config,
            &sketch_dev,
            sketches.get_strides(),
            sketches.sketchsize64,
            &kmers_dev,
            dist_type,
            k_vals,
            ani,
            &mut dist_dev,
        )?;
    }

    let dist_mut = distances.dists_mut();
    *dist_mut = dist_dev.to_host_vec(&stream)?;

    Ok(distances)
}
